#!/usr/bin/env python3
#script created with use of chatGPT 15-12-205 Lennart Johansson

"""
Replace a subsequence in FASTQ reads based on coordinates from an annotation TSV.

- Uses read_pos_start / read_pos_end (1-based, inclusive) from TSV
- Verifies that the subsequence in the read is sufficiently similar to the
  reference sequence provided in the TSV (seq_between_coords)
- Replaces that subsequence with --replace-sequence
- Handles insertions (adds '?' qualities) and deletions (removes qualities)
- Writes a modified FASTQ and a log TSV of all changes
"""

import argparse
import gzip
import sys
import os
from difflib import SequenceMatcher


def parse_args():
    p = argparse.ArgumentParser(
        description="Replace STR sequence in FASTQ reads using annotation TSV"
    )
    p.add_argument("--fastq", required=True, help="Input FASTQ(.gz)")
    p.add_argument("--tsv", required=True, help="Annotation TSV (from previous step)")
    p.add_argument("--replace-sequence", required=True,
                   help="Sequence to insert instead of original")
    p.add_argument("--out-fastq", required=True, help="Output modified FASTQ(.gz)")
    p.add_argument("--log-tsv", required=True,
                   help="Log TSV with original and replaced sequences")
    p.add_argument("--min-identity", type=float, default=0.85,
                   help="Minimum sequence identity required (default: 0.85)")
    return p.parse_args()


def open_maybe_gzip(path, mode):
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


def read_annotation_tsv(tsv_path):
    """Return dict: read_id -> dict with positions and reference sequence"""
    data = {}
    with open(tsv_path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}

        required = [
            "read_id",
            "read_pos_start",
            "read_pos_end",
            "seq_between_coords",
        ]
        for r in required:
            if r not in idx:
                sys.exit(f"ERROR: TSV missing required column '{r}'")

        for line in f:
            cols = line.rstrip("\n").split("\t")
            read_id = cols[idx["read_id"]]
            data[read_id] = {
                "start": int(cols[idx["read_pos_start"]]),
                "end": int(cols[idx["read_pos_end"]]),
                "ref_seq": cols[idx["seq_between_coords"]],
            }
    return data


def seq_identity(a, b):
    return SequenceMatcher(None, a, b).ratio()


def main():
    args = parse_args()

    ann = read_annotation_tsv(args.tsv)
    replace_seq = args.replace_sequence

    log_fh = open(args.log_tsv, "w")
    log_fh.write(
        "read_id\tread_pos_start\tread_pos_end\toriginal_sequence\treplaced_sequence\n"
    )

    fin = open_maybe_gzip(args.fastq, "rt")
    fout = open_maybe_gzip(args.out_fastq, "wt")

    while True:
        h = fin.readline()
        if not h:
            break
        seq = fin.readline().rstrip("\n")
        plus = fin.readline()
        qual = fin.readline().rstrip("\n")

        read_id = h.split()[0][1:]

        if read_id not in ann:
            fout.write(h)
            fout.write(seq + "\n")
            fout.write(plus)
            fout.write(qual + "\n")
            continue

        info = ann[read_id]
        s0 = info["start"] - 1
        e0 = info["end"]

        read_subseq = seq[s0:e0]
        ref_subseq = info["ref_seq"]

        identity = seq_identity(read_subseq, ref_subseq)
        if identity < args.min_identity:
            sys.exit(
                f"ERROR: sequence in reads does not match reference sequence "
                f"for read {read_id} (identity={identity:.3f})"
            )

        new_seq = seq[:s0] + replace_seq + seq[e0:]

        old_len = e0 - s0
        new_len = len(replace_seq)

        if new_len > old_len:
            new_qual_mid = qual[s0:e0] + "?" * (new_len - old_len)
        else:
            new_qual_mid = qual[s0:s0 + new_len]

        new_qual = qual[:s0] + new_qual_mid + qual[e0:]

        fout.write(h)
        fout.write(new_seq + "\n")
        fout.write(plus)
        fout.write(new_qual + "\n")

        log_fh.write(
            f"{read_id}\t{info['start']}\t{info['end']}\t"
            f"{read_subseq}\t{replace_seq}\n"
        )

    fin.close()
    fout.close()
    log_fh.close()

    print("Done.")


if __name__ == "__main__":
    main()
