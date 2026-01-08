#!/usr/bin/env python3
"""
annotate_fastq_reads_with_str_coords.py

Annotates FASTQ reads with STR coordinates using CRAM alignments (primary only).
Handles both forward and reverse strand reads correctly.
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict

try:
    import pysam
except ImportError as e:
    sys.exit("This script requires pysam. Install with: pip install pysam\n" + str(e))


def parse_args():
    p = argparse.ArgumentParser(description="Annotate FASTQ reads with STR coordinates using CRAM and STRaglr TSV")
    p.add_argument("--tsv", required=True, help="STRaglr TSV file")
    p.add_argument("--fastq", required=True, nargs='+', help="Input FASTQ.gz files")
    p.add_argument("--cram", required=True, help="CRAM file (indexed)")
    p.add_argument("--crai", required=True, help="CRAM index (.crai)")
    p.add_argument("--fasta", required=True, help="Reference fasta (indexed)")
    p.add_argument("--outdir", default="annotate_out", help="Output directory")
    p.add_argument("--min-mapq", type=int, default=0, help="Minimum MAPQ for primary alignment")
    return p.parse_args()


# --- TSV utilities ---

def read_tsv(tsv_path):
    rows = []
    header = None
    with open(tsv_path) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.lower().startswith("#chrom"):
                    header = line.lstrip("#").strip().split()
                continue
            cols = line.rstrip("\n").split("\t")
            rows.append(cols)
    return header, rows


def tsv_rows_for_region(rows, chrom, start, end):
    matches = []
    for cols in rows:
        if len(cols) < 3:
            continue
        try:
            rchrom = cols[0]
            rstart = int(cols[1])
            rend = int(cols[2])
        except Exception:
            continue
        if rchrom == chrom and rstart == start and rend == end:
            matches.append(cols)
    return matches


def build_tsv_index(rows, header):
    read_to_rows = defaultdict(list)
    for cols in rows:
        if header:
            d = {header[i]: (cols[i] if i < len(cols) else "") for i in range(len(header))}
            readname = d.get("read_name") or d.get("read") or d.get("readname") or d.get("read_id") or ""
            if not readname and len(cols) > 7:
                readname = cols[7]
            read_to_rows[readname].append(d)
        else:
            readname = cols[7] if len(cols) > 7 else ""
            read_to_rows[readname].append(cols)
    return read_to_rows


# --- FASTQ utilities ---

def load_fastq_read_ids_and_sequences(fastq_path):
    reads = {}
    order = []
    opener = gzip.open if fastq_path.endswith(".gz") else open
    with opener(fastq_path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            fh.readline()  # +
            fh.readline()  # qual
            rid = header.split()[0][1:]
            reads[rid] = seq
            order.append(rid)
    return reads, order


def infer_region_from_fastq(fq_path):
    import re
    fname = os.path.basename(fq_path)
    m = re.match(r"^(chr[\w]+)_(\d+)_(\d+)_(\d+(?:\.\d+)?)\.fastq(?:\.gz)?$", fname)
    if not m:
        raise ValueError(f"FASTQ name does not match expected pattern: {fname}")
    chrom, start, end, allele = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
    return chrom, start, end, allele


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


# --- CRAM utilities ---

def collect_cram_alignments_for_region(cram_path, crai_path, chrom, start, end, fastq_read_ids, min_mapq=0):
    mapping = defaultdict(list)
    cram = pysam.AlignmentFile(cram_path, "rc", index_filename=crai_path)
    ref0_start = start - 1
    ref0_end_excl = end
    found_reads = set()

    for aln in cram.fetch(chrom, ref0_start, ref0_end_excl):
        if aln.is_unmapped or aln.mapping_quality < min_mapq:
            continue
        if aln.query_name in fastq_read_ids:
            mapping[aln.query_name].append(aln)
            found_reads.add(aln.query_name)
        if len(found_reads) == len(fastq_read_ids):
            break
    cram.close()
    return mapping


def choose_primary_alignment(alignments, min_mapq=0):
    prim = [a for a in alignments if not a.is_secondary and not a.is_supplementary]
    if not prim:
        return None
    prim_high = [a for a in prim if a.mapping_quality >= min_mapq]
    if prim_high:
        return sorted(prim_high, key=lambda x: x.mapping_quality, reverse=True)[0]
    return prim[0]


def compute_read_coords_for_ref_region(aln, ref0_start, ref0_end):
    """
    Map reference region to read coordinates.
    Returns (read_start, read_end) in read sequence, inclusive.
    """
    pairs = aln.get_aligned_pairs(matches_only=True)
    rpos_list = [r for r, ref in pairs if ref0_start <= ref <= ref0_end]
    if rpos_list:
        return min(rpos_list), max(rpos_list)
    return None, None


# --- Main script ---

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    header, tsv_rows = read_tsv(args.tsv)
    tsv_read_index = build_tsv_index(tsv_rows, header)

    for fq in args.fastq:
        try:
            chrom, start, end, allele = infer_region_from_fastq(fq)
        except ValueError as e:
            print("ERROR:", e)
            continue

        region_rows = tsv_rows_for_region(tsv_rows, chrom, start, end)
        tsv_info_by_read = defaultdict(list)
        reads_in_region_tsv = set()
        for cols in region_rows:
            if header:
                d = {header[i]: (cols[i] if i < len(cols) else "") for i in range(len(header))}
                rid = d.get("read_name") or d.get("read") or d.get("readname") or d.get("read_id") or ""
            else:
                rid = cols[7] if len(cols) > 7 else ""
            if rid:
                tsv_info_by_read[rid].append(d if header else cols)
                reads_in_region_tsv.add(rid)

        reads_seq, reads_order = load_fastq_read_ids_and_sequences(fq)
        fastq_read_ids = set(reads_order)

        cram_map = collect_cram_alignments_for_region(
            args.cram, args.crai, chrom, start, end, fastq_read_ids, min_mapq=args.min_mapq
        )

        out_base = os.path.basename(fq).replace(".fastq.gz", "").replace(".fastq", "")
        out_tsv = os.path.join(args.outdir, out_base + "_annotation.tsv")

        with open(out_tsv, "wt") as outfh:
            cols_out = [
                "read_id", "in_tsv", "tsv_alleles", "primary_found",
                "primary_ref_start", "primary_ref_end", "read_pos_start", "read_pos_end",
                "seq_between_coords", "other_alignments", "other_alignment_flags",
                "strand", "notes"
            ]
            outfh.write("\t".join(cols_out) + "\n")

            for rid in reads_order:
                seq = reads_seq[rid]
                in_tsv = "yes" if rid in reads_in_region_tsv else "no"
                tsv_alleles = []
                if rid in tsv_info_by_read:
                    for entry in tsv_info_by_read[rid]:
                        if isinstance(entry, dict):
                            tsv_alleles.append(entry.get("allele", ""))
                        else:
                            tsv_alleles.append(entry[13] if len(entry) > 13 else "")
                tsv_alleles_str = ";".join(tsv_alleles)

                aln_list = cram_map.get(rid, [])
                other_flags = []
                for a in aln_list:
                    if a.is_secondary:
                        other_flags.append("sec")
                    if a.is_supplementary:
                        other_flags.append("supp")
                other_flags_str = ",".join(sorted(set(other_flags))) if other_flags else ""
                other_count = len(other_flags)

                primary = choose_primary_alignment(aln_list, min_mapq=args.min_mapq)
                primary_found = "yes" if primary else "no"
                pref_start = str(primary.reference_start) if primary else "NA"
                pref_end = str(primary.reference_end) if primary else "NA"

                read_pos_start = "NA"
                read_pos_end = "NA"
                seq_between = ""
                strand = "NA"
                notes = []

                if primary:
                   rstart, rend = compute_read_coords_for_ref_region(primary, start - 1, end - 1)
                   if rstart is not None and rend is not None:
                       # Determine slice boundaries
                       left = min(rstart, rend)
                       right = max(rstart, rend)

                       # bounds check
                       left = max(left, 0)
                       right = min(right, len(seq) - 1)

                       # slice sequence
                       seq_between = seq[left:right + 1]

                       # assign strand and reverse-complement if needed
                       if primary.is_reverse:
                           seq_oriented = revcomp(seq)
                           strand = "-"
                       else:
                           seq_oriented = seq
                           strand = "+"

                       seq_between = seq_oriented[left:right + 1]

                       # save read positions for TSV
                       read_pos_start = str(left)
                       read_pos_end = str(right)

                       # --- DEBUG block similar to old script ---
                       read_leftmost = 0
                       read_rightmost = len(seq) - 1
                       # Take a slice after the repeat for context (careful with bounds)
                       seq_between2 = seq[right + 1] if right + 1 < len(seq) else ""
                       total_seq = seq[read_leftmost:read_rightmost]
                       print(f"DEBUG: read {rid}, rstart={rstart}, rend={rend}, "
                             f"left={left}, right={right}, seq_between={seq_between}, "
                             f"seq_between2={seq_between2}, total_seq={total_seq}, strand={strand}")


                if in_tsv == "yes" and len(tsv_info_by_read.get(rid, [])) > 1:
                    notes.append("multiple_tsv_rows")
                if primary_found == "no" and aln_list:
                    notes.append("no_primary_but_other_alignments")
                if primary_found == "no" and not aln_list:
                    notes.append("no_alignment_in_cram")
                notes_str = ";".join(notes) if notes else ""

                outfh.write("\t".join([
                    rid, in_tsv, tsv_alleles_str, primary_found, pref_start, pref_end,
                    read_pos_start, read_pos_end, seq_between,
                    str(other_count), other_flags_str, strand, notes_str
                ]) + "\n")

        print(f"Wrote annotation TSV: {out_tsv}")


if __name__ == "__main__":
    main()
