#!/usr/bin/env python3
# script created with use of chatGPT 16-12-2025 Lennart Johansson

import argparse
import gzip

def parse_args():
    p = argparse.ArgumentParser(
        description="Replace reads in an original FASTQ with adapted reads from a synthetic FASTQ"
    )
    p.add_argument("--original-fastq", required=True)
    p.add_argument("--synthetic-fastq", required=True)
    p.add_argument("--out-fastq", required=True)
    p.add_argument("--log-tsv", required=True)
    p.add_argument("--log-detailed-tsv", required=True)
    return p.parse_args()


def open_fastq(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def read_fastq_to_dict(fastq_path):
    """
    Load FASTQ into dict:
      read_id -> (sequence, quality)
    """
    reads = {}
    with open_fastq(fastq_path, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().rstrip()
            f.readline()  # plus
            qual = f.readline().rstrip()

            read_id = header.split()[0][1:]
            reads[read_id] = (seq, qual)

    return reads


def main():
    args = parse_args()

    print("Loading synthetic FASTQ...")
    synthetic_reads = read_fastq_to_dict(args.synthetic_fastq)
    remaining = set(synthetic_reads.keys())

    print(f" Loaded {len(remaining)} synthetic reads")

    replaced_count = 0

    out_fh = open_fastq(args.out_fastq, "wt")
    log_fh = open(args.log_tsv, "w")
    detailed_log_fh = open(args.log_detailed_tsv, "w")

    log_fh.write("read_id\treplaced\n")
    detailed_log_fh.write(
        "original_fastq\toutput_fastq\tread_id\toriginal_sequence\treplaced_sequence\n"
    )

    print("Processing original FASTQ...")

    with open_fastq(args.original_fastq, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break

            seq = f.readline().rstrip()
            plus = f.readline()
            qual = f.readline().rstrip()

            # FAST PATH: all replacements already done
            if not remaining:
                out_fh.write(header)
                out_fh.write(seq + "\n")
                out_fh.write(plus)
                out_fh.write(qual + "\n")
                continue

            read_id = header.split()[0][1:]

            if read_id in remaining:
                new_seq, new_qual = synthetic_reads[read_id]

                out_fh.write(header)
                out_fh.write(new_seq + "\n")
                out_fh.write(plus)
                out_fh.write(new_qual + "\n")

                log_fh.write(f"{read_id}\tyes\n")
                detailed_log_fh.write(
                    f"{args.original_fastq}\t{args.out_fastq}\t"
                    f"{read_id}\t{seq}\t{new_seq}\n"
                )

                remaining.remove(read_id)
                replaced_count += 1
            else:
                out_fh.write(header)
                out_fh.write(seq + "\n")
                out_fh.write(plus)
                out_fh.write(qual + "\n")

    out_fh.close()
    log_fh.close()
    detailed_log_fh.close()

    print("Finished.")
    print(f" Replaced reads: {replaced_count}")
    print(f" Output FASTQ: {args.out_fastq}")
    print(f" Simple log: {args.log_tsv}")
    print(f" Detailed log: {args.log_detailed_tsv}")


if __name__ == "__main__":
    main()

