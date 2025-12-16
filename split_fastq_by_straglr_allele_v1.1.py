#!/usr/bin/env python3
#script created with use of chatGPT 10-12-205 Lennart Johansson
import gzip
import argparse
import os
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Split FASTQ by locus/allele from STRaglr TSV")
    parser.add_argument("--tsv", required=True, help="STRaglr TSV file")
    parser.add_argument("--fastq", required=True, help="Input FASTQ.gz file")
    return parser.parse_args()

def load_read_ids(tsv_file):
    """
    Returns a dict:
    {read_id: output_filename}
    """
    read2file = {}
    # keep track of all output files we need to open later
    output_files = set()
    with open(tsv_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip().split("\t")
            if len(cols) < 14:
                continue
            status = cols[14] if len(cols) > 14 else "full"
            if status != "full":
                continue
            chrom, start, end, read_id, allele = cols[0], cols[1], cols[2], cols[7], cols[13]
            outfile = f"{chrom}_{start}_{end}_{allele}.fastq"
            read2file[read_id] = outfile
            output_files.add(outfile)
    return read2file, output_files

def split_fastq(fastq_file, read2file, output_files, outdir):
    os.makedirs(outdir, exist_ok=True)

    out_handles = {
        fname: gzip.open(os.path.join(outdir, fname + ".gz"), "wt")
        for fname in output_files
    }

    with gzip.open(fastq_file, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()

            read_id = header.split()[0][1:]
            if read_id in read2file:
                out_fname = read2file[read_id]
                out_handles[out_fname].write(f"{header}{seq}{plus}{qual}")

    for h in out_handles.values():
        h.close()

def parse_args():
    parser = argparse.ArgumentParser(description="Split FASTQ by locus/allele from STRaglr TSV")
    parser.add_argument("--tsv", required=True, help="STRaglr TSV file")
    parser.add_argument("--fastq", required=True, help="Input FASTQ.gz file")
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for split FASTQ files"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    read2file, output_files = load_read_ids(args.tsv)
    split_fastq(args.fastq, read2file, output_files, args.outdir)
    print(f"Created {len(output_files)} FASTQ.gz files in {args.outdir}.")

if __name__ == "__main__":
    main()
