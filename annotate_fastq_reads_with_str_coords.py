#!/usr/bin/env python3
#script created with use of chatGPT 12-12-205 Lennart Johansson
"""
annotate_fastq_reads_with_str_coords.py

For each input FASTQ (name format: chrom_start_end_allele.fastq[.gz]),
find matching TSV rows in the STRaglr TSV, use the CRAM alignments
(ONLY primary alignments) to compute where the reference region maps into
the read, and output a TSV per FASTQ summarizing:

Columns output:
  read_id
  in_tsv (yes/no)
  tsv_alleles (semicolon-separated if multiple TSV rows matched)
  primary_found (yes/no)
  primary_ref_start (0-based, aln.reference_start) or NA
  primary_ref_end   (0-based exclusive, aln.reference_end) or NA
  read_pos_start    (0-based index in read where region begins) or NA
  read_pos_end      (0-based index in read where region ends, inclusive) or NA
  seq_between_coords (the read sequence between read_pos_start..read_pos_end inclusive; empty if NA)
  other_alignments   (number of alignments overlapping region that are secondary/supplementary)
  other_alignment_flags (comma-separated list like "supp,sec" if present)
  notes

This script uses the output of split_fastq_by_straglr_allele.py
Multiple fastq files can be given as input separated by a space

Usage:
  python3 annotate_fastq_reads_with_str_coords.py --tsv straglr.tsv \
      --fastq chr1_57367053_57367108_12.2.fastq.gz \
      --cram sample.cram --fasta ref.fa
"""

import argparse
import gzip
import os
import re
import sys
from collections import defaultdict

try:
    import pysam
except Exception as e:
    sys.exit("This script requires pysam. Install with: pip install pysam\n" + str(e))


def parse_args():
    p = argparse.ArgumentParser(description="Annotate FASTQ reads with STR coordinates using CRAM and STRaglr TSV")
    p.add_argument("--tsv", required=True, help="STRaglr TSV file (can contain many loci)")
    p.add_argument("--fastq", required=True, nargs='+', help="Input FASTQ.gz files (named chrom_start_end_allele.fastq(.gz))")
    p.add_argument("--cram", required=True, help="CRAM file for the sample (indexed .crai required)")
    p.add_argument("--crai", required=True, help="CRAM index file (.crai) explicitly specified")
    p.add_argument("--fasta", required=True, help="Reference fasta that matches the CRAM (indexed .fai required)")
    p.add_argument("--outdir", default="annotate_out", help="Output directory (default: annotate_out)")
    p.add_argument("--min-mapq", type=int, default=0, help="Minimum MAPQ when considering primary alignment (default: 0)")
    return p.parse_args()


def read_tsv(tsv_path):
    """
    Parse TSV. Return:
      header (list or None), and list of rows (list of lists)
    If a header line starting with '#chrom' exists, we use it and strip leading '#'.
    """
    rows = []
    header = None
    with open(tsv_path) as fh:
        for line in fh:
            if line.startswith("#"):
                # potential header like "#chrom  start   end ..."
                if line.lower().startswith("#chrom"):
                    header = line.lstrip("#").strip().split()
                continue
            cols = line.rstrip("\n").split("\t")
            rows.append(cols)
    return header, rows


def tsv_rows_for_region(rows, chrom, start, end):
    """
    Given parsed TSV rows, return subset where chrom,start,end match.
    Assumes TSV start/end columns correspond to columns 0/1/2 as in examples.
    Input 'start' and 'end' are integers (1-based in TSV, but we compare raw strings).
    We'll compare numerically to be robust.
    """
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
    """
    Build maps for quick access: read_name -> list of TSV rows (as dict if header available)
    If header is present, return dicts with keys from header; else fallback indices.
    """
    read_to_rows = defaultdict(list)
    for cols in rows:
        if header:
            # create dict mapping header->value (safe for different column counts)
            d = {}
            for i, key in enumerate(header):
                d[key] = cols[i] if i < len(cols) else ""
            readname = d.get("read_name") or d.get("read") or d.get("readname") or d.get("read_id") or d.get("read.name") or ""
            if not readname:
                # fallback to col 8 (index 7)
                if len(cols) > 7:
                    readname = cols[7]
            read_to_rows[readname].append(d)
        else:
            # fallback: col 8 -> index 7, allele in index 13 as used previously
            readname = cols[7] if len(cols) > 7 else ""
            read_to_rows[readname].append(cols)
    return read_to_rows


import re
import os

def infer_region_from_fastq(fq_path):
    import re, os
    fname = os.path.basename(fq_path)

    # Expecting: chrXX_start_end_allele.fastq.gz
    # Example: chr16_66490400_66490453_30.4.fastq.gz
    m = re.match(r"^(chr[\w]+)_(\d+)_(\d+)_(\d+(?:\.\d+)?)\.fastq\.gz$", fname)
    if not m:
        raise ValueError(f"FASTQ name does not match expected pattern: {fname}")

    chrom = m.group(1)
    start = int(m.group(2))
    end   = int(m.group(3))
    allele = m.group(4)

    return chrom, start, end, allele


def load_fastq_read_ids_and_sequences(fastq_path):
    """
    Return dict read_id -> sequence (string). We'll iterate reads sequentially in main loop,
    but for convenience also return an ordered list of read ids (to preserve FASTQ order).
    """
    reads = {}
    order = []
    opener = gzip.open if fastq_path.endswith(".gz") else open
    with opener(fastq_path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            plus = fh.readline()
            qual = fh.readline()
            rid = header.split()[0][1:]
            reads[rid] = seq
            order.append(rid)
    return reads, order


import time
from collections import defaultdict
import os
import pysam

def collect_cram_alignments_for_region(cram_path, crai_path, chrom, start, end, fastq_read_ids, min_mapq=0):
    """
    Fetch alignments overlapping [start, end] for only reads in fastq_read_ids.
    Returns dict: read_id -> list of AlignedSegment objects.
    Includes debug prints and early stopping if all reads are found.
    """
    print(f"Opening CRAM {cram_path} with CRAI {crai_path} ...")
    start_time = time.time()
    cram = pysam.AlignmentFile(cram_path, "rc", index_filename=crai_path)
    print(f"CRAM opened in {time.time() - start_time:.2f} seconds")

    ref0_start = start - 1
    ref0_end_excl = end
    mapping = defaultdict(list)
    found_reads = set()
    count = 0

    print(f"Fetching alignments for {chrom}:{start}-{end} ...")
    fetch_start = time.time()
    for aln in cram.fetch(chrom, ref0_start, ref0_end_excl):
        count += 1
        if count % 1000 == 0:
            print(f"  Processed {count} alignments so far...")
        # skip unmapped or low MAPQ
        if aln.is_unmapped or aln.mapping_quality < min_mapq:
            continue
        if aln.query_name in fastq_read_ids:
            mapping[aln.query_name].append(aln)
            found_reads.add(aln.query_name)
        # early exit if all reads are found
        if len(found_reads) == len(fastq_read_ids):
            print(f"  All {len(fastq_read_ids)} FASTQ reads found, stopping fetch early.")
            break

    cram.close()
    print(f"Finished fetching alignments in {time.time() - fetch_start:.2f} seconds; "
          f"Total alignments scanned: {count}, reads found: {len(found_reads)}")

    missing_reads = fastq_read_ids - found_reads
    if missing_reads:
        print(f"Warning: {len(missing_reads)} FASTQ reads not found in CRAM region: {missing_reads}")

    return mapping


def choose_primary_alignment(alignments, min_mapq=0):
    """
    Return the primary alignment among list (or None).
    Primary = not secondary and not supplementary.
    If multiple primary (rare), pick highest mapq >= min_mapq otherwise first.
    """
    prim = [a for a in alignments if (not a.is_secondary and not a.is_supplementary)]
    if not prim:
        return None
    # prefer a primary with mapq >= min_mapq
    prim_high = [a for a in prim if a.mapping_quality >= min_mapq]
    if prim_high:
        # choose highest mapq
        return sorted(prim_high, key=lambda x: x.mapping_quality, reverse=True)[0]
    # otherwise return first
    return prim[0]


def compute_read_coords_for_ref_region(aln, ref0_start, ref0_end):
    """
    Given a pysam AlignedSegment (primary), compute read-based start/end indices (0-based)
    that cover reference coordinates in [ref0_start..ref0_end] inclusive.

    Returns (read_start, read_end) as 0-based indices into read sequence (inclusive).
    If mapping cannot assign (e.g., deletion over region), we attempt a best-guess:
      - use nearest read pos to left or right similarly to previous script.
    If still cannot determine, returns (None, None).
    """
    pairs = aln.get_aligned_pairs(matches_only=False)  # list of (read_pos or None, ref_pos or None)
    # Build mapping ref_pos -> read_pos for non-None
    ref_to_rpos = {}
    for rpos, rref in pairs:
        if rref is not None and rpos is not None:
            ref_to_rpos[rref] = rpos

    # collect read positions for ref positions inside region
    rpos_list = [r for ref, r in ref_to_rpos.items() if ref0_start <= ref <= ref0_end]
    if rpos_list:
        read_start = min(rpos_list)
        read_end = max(rpos_list)
        return read_start, read_end

    # no direct mapping inside region (e.g., region entirely deleted in read). Try nearest flanking positions:
    refs = sorted(ref_to_rpos.keys())
    if not refs:
        return None, None
    left = None
    right = None
    for r in refs:
        if r < ref0_start:
            left = r
        elif r > ref0_end and right is None:
            right = r
            break
    if left is not None:
        # insert after left read pos
        read_start = ref_to_rpos[left] + 1
        # read_end: if right exists, take ref_to_rpos[right] - 1 else set to read_start (best-effort)
        if right is not None:
            read_end = ref_to_rpos[right] - 1
        else:
            read_end = read_start
        return read_start, read_end
    elif right is not None:
        read_end = ref_to_rpos[right] - 1
        read_start = read_end
        return read_start, read_end

    return None, None


def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Reading TSV...")
    header, tsv_rows = read_tsv(args.tsv)
    print(f"Loaded {len(tsv_rows)} TSV rows; header columns: {header if header else 'none detected'}")

    # build tsv index (by read name) for quick lookup if needed later
    tsv_read_index = build_tsv_index(tsv_rows, header)

    # open fasta for region reference retrieval
    fasta = pysam.FastaFile(args.fasta)

    for fq in args.fastq:
        print(f"\nProcessing FASTQ: {fq}")
        try:
            chrom, start, end, allele_from_fname = infer_region_from_fastq(fq)
        except ValueError as e:
            print("ERROR:", e)
            continue
        print(f" Inferred region: {chrom}:{start}-{end} allele={allele_from_fname}")

        # find TSV rows matching region
        region_rows = tsv_rows_for_region(tsv_rows, chrom, start, end)
        print(f" Found {len(region_rows)} TSV rows for this region in the master TSV")

        # build set of read IDs from region rows (use header if available)
        reads_in_region_tsv = set()
        tsv_info_by_read = defaultdict(list)
        for cols in region_rows:
            if header:
                # build dict
                d = {}
                for i, key in enumerate(header):
                    d[key] = cols[i] if i < len(cols) else ""
                readname = d.get("read_name") or d.get("read") or d.get("readname") or d.get("read_id") or ""
                reads_in_region_tsv.add(readname)
                tsv_info_by_read[readname].append(d)
            else:
                # fallback
                readname = cols[7] if len(cols) > 7 else ""
                reads_in_region_tsv.add(readname)
                tsv_info_by_read[readname].append(cols)

        # load fastq reads (sequences)
        print(" Loading FASTQ reads into memory...")
        reads_seq, reads_order = load_fastq_read_ids_and_sequences(fq)
        print(f"  FASTQ contains {len(reads_order)} reads")

        import gzip

        # Load read IDs from the FASTQ
        fastq_read_ids = set()
        with gzip.open(fq, "rt") as f:
            while True:
                header = f.readline()
                if not header:
                    break
                fastq_read_ids.add(header.split()[0][1:])  # Remove the leading '@'
                f.readline()  # skip sequence line
                f.readline()  # skip '+'
                f.readline()  # skip quality
        print(f" FASTQ contains {len(fastq_read_ids)} reads")



        # collect cram alignments overlapping region
        print(" Collecting CRAM alignments overlapping region...")
        cram_map = collect_cram_alignments_for_region(
            args.cram,       # path to CRAM
            args.crai,       # path to CRAI
            chrom,           # chromosome from FASTQ filename
            start,           # start position
            end,             # end position
            fastq_read_ids,  # only fetch these reads
            min_mapq=args.min_mapq
        )
        print(f"  Found alignments for {len(cram_map)} distinct read names in the CRAM that overlap region")

        # prepare output TSV file
        out_base = os.path.basename(fq).replace(".fastq.gz", "").replace(".fastq", "")
        out_tsv = os.path.join(args.outdir, out_base + "_annotation.tsv")
        log_file = os.path.join(args.outdir, out_base + "_annotation.log")

        with open(out_tsv, "wt") as outfh, open(log_file, "wt") as logfh:
            # header
            cols = [
                "read_id", "in_tsv", "tsv_alleles", "primary_found", "primary_ref_start",
                "primary_ref_end", "read_pos_start", "read_pos_end", "seq_between_coords",
                "other_alignments", "other_alignment_flags", "notes"
            ]
            outfh.write("\t".join(cols) + "\n")

            # iterate reads in FASTQ order
            for rid in reads_order:
                seq = reads_seq.get(rid, "")
                in_tsv = "yes" if rid in reads_in_region_tsv else "no"
                tsv_alleles = []
                if rid in tsv_info_by_read:
                    # if header present, try to extract allele key else fallback to index 13
                    for entry in tsv_info_by_read[rid]:
                        if isinstance(entry, dict):
                            a = entry.get("allele", "")
                        else:
                            a = entry[13] if len(entry) > 13 else ""
                        tsv_alleles.append(a)
                tsv_alleles_str = ";".join(tsv_alleles) if tsv_alleles else ""

                # find alignments for this read overlapping region
                aln_list = cram_map.get(rid, [])
                # classify other alignments flags
                other_flags = []
                other_count = 0
                for a in aln_list:
                    if a.is_secondary:
                        other_flags.append("sec")
                        other_count += 1
                    if a.is_supplementary:
                        other_flags.append("supp")
                        other_count += 1
                other_flags_str = ",".join(sorted(set(other_flags))) if other_flags else ""

                # choose primary among alignments
                primary = choose_primary_alignment(aln_list, min_mapq=args.min_mapq)

                primary_found = "yes" if primary is not None else "no"
                pref_start = pref_end = "NA"
                read_pos_start = read_pos_end = "NA"
                seq_between = ""

                if primary is not None:
                    pref_start = str(primary.reference_start)  # 0-based
                    pref_end = str(primary.reference_end)      # 0-based exclusive
                    # compute read coords covering ref region
                    ref0_start = start - 1
                    ref0_end = end - 1
                    rp = compute_read_coords_for_ref_region(primary, ref0_start, ref0_end)
                    if rp[0] is None:
                        read_pos_start = "NA"
                        read_pos_end = "NA"
                        seq_between = ""
                    else:
                        rstart, rend = rp
                        # ensure bounds inside the read
                        # === END-ANCHORED STR EXTRACTION USING STRaglr COPY NUMBER ===

                        # trust the right boundary from alignment
                        read_end = rend

                        repeat_unit = entry.get("target_repeat", "")
                        copy_number = entry.get("copy_number", "")

                        if repeat_unit and copy_number:
                            try:
                                ru_len = len(repeat_unit)
                                n_units = float(copy_number)
                                str_len = int(round(n_units * ru_len))

                                read_start = read_end - str_len + 1
                                if read_start < 0:
                                    read_start = 0

                                read_pos_start = str(read_start)
                                read_pos_end = str(read_end)
                                seq_between = seq[read_start:read_end + 1]

                            except Exception:
                                read_pos_start = "NA"
                                read_pos_end = "NA"
                                seq_between = ""
                        else:
                            read_pos_start = "NA"
                            read_pos_end = "NA"
                            seq_between = ""

                else:
                    # no primary - report presence of other alignments in log
                    if aln_list:
                        logfh.write(f"{rid}\tNO_PRIMARY\t{len(aln_list)}_alignments\tdetails_flags:{other_flags_str}\n")

                notes = []
                if in_tsv == "yes":
                    # check whether read id appears multiple times in TSV rows for region
                    if len(tsv_info_by_read.get(rid, [])) > 1:
                        notes.append("multiple_tsv_rows")
                if primary_found == "no" and aln_list:
                    notes.append("no_primary_but_other_alignments")
                if primary_found == "no" and not aln_list:
                    notes.append("no_alignment_in_cram")
                notes_str = ";".join(notes) if notes else ""

                outfh.write("\t".join([
                    rid, in_tsv, tsv_alleles_str, primary_found, pref_start,
                    pref_end, read_pos_start, read_pos_end, seq_between,
                    str(other_count), other_flags_str, notes_str
                ]) + "\n")

        print(f" Wrote annotation TSV: {out_tsv}")
        print(f" Wrote log: {log_file}")

    fasta.close()
    print("\nAll done.")


if __name__ == "__main__":
    main()
