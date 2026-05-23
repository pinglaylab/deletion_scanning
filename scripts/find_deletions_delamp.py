"""
Script to identify translocations in a BAM file from nanopore sequencing.
Outputs a TSV file with translocation breakpoints and a BAM file with the reads supporting translocations.

Read name structures supported:

1) 4-field:
   <readid>_<barcode>_<cell>_<UMI>
   Example:
     36edcb02-f99a-4ec7-8ddc-c70716c0fa14_GTTTAAGAAATT_GTCCTAAAGTGAAGTA_AATAATGGCTGCTTTT

2) 2-field (single UMI only):
   <readid>_<UMI>
   Example:
     7dbdf249-64ae-4cb2-b1ec-405c2196197a_CCGTACCCGG

Usage:
    python translocation_finder.py input.bam output.tsv output.bam
"""

import pysam
import argparse
import csv
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Identify translocations from a BAM file.")
    parser.add_argument("bamfile", help="Input BAM file")
    parser.add_argument("output_tsv", help="Output TSV file with translocation breakpoints")
    parser.add_argument("output_bam", help="Output BAM file with reads supporting translocations")
    return parser.parse_args()

def parse_read_name(read_name: str):
    """
    Returns (barcode, cell, UMI) from read_name.

    Supported formats:
      - <readid>_<barcode>_<cell>_<UMI>  (>=4 underscore-delimited parts; uses last 3 fields)
      - <readid>_<UMI>                  (exactly 2 parts; barcode/cell empty)

    If parsing fails, returns (None, None, None).
    """
    parts = read_name.split("_")

    # New format: <readid>_<UMI>
    if len(parts) == 2:
        umi = parts[1]
        return "", "", umi

    # Old format (or anything with extra underscores): take the last 3 fields as barcode/cell/UMI
    if len(parts) >= 4:
        barcode, cell, umi = parts[-3], parts[-2], parts[-1]
        return barcode, cell, umi

    return None, None, None

def main():
    args = parse_args()

    bam = pysam.AlignmentFile(args.bamfile, "rb")
    header = bam.header

    # Group alignments by read name.
    read_groups = defaultdict(list)
    for aln in bam:
        read_groups[aln.query_name].append(aln)
    bam.close()

    with open(args.output_tsv, "w", newline="") as tsvfile, \
         pysam.AlignmentFile(args.output_bam, "wb", header=header) as bam_out:

        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(["read_name", "start", "end", "chr_2", "strand_2", "barcode", "cell", "UMI"])

        for read_name, alignments in read_groups.items():
            ref_aln = None
            new_aln = None

            # Pick one alignment on 'reference' and one alignment on a non-reference chromosome.
            for aln in alignments:
                if aln.reference_name == "reference":
                    ref_aln = aln
                else:
                    new_aln = aln

            if ref_aln is None or new_aln is None:
                continue

            # Parse barcode/cell/UMI
            barcode, cell, umi = parse_read_name(read_name)
            if umi is None:  # malformed name
                continue

            # Write all alignments of the read to the output BAM.
            for aln in alignments:
                bam_out.write(aln)

            # Breakpoints
            start = ref_aln.reference_end  # 0-based exclusive end
            if new_aln.is_reverse:
                end = new_aln.reference_start+1
            else:
                end = new_aln.reference_end

            chr_2 = new_aln.reference_name
            strand_2 = '-' if new_aln.is_reverse else '+'

            writer.writerow([read_name, start, end, chr_2, strand_2, barcode, cell, umi])

if __name__ == "__main__":
    main()
koeppelj@nexus3:/net/shendure/vol10/projects/jonas/nanopore/delamp$ cat sample_fastq_q30.py 
#!/usr/bin/env python3
import argparse, gzip, random, sys

def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def phred_scores(qual: str):
    # Standard FASTQ Phred+33
    return [ord(c) - 33 for c in qual.rstrip("\n")]

def passes(scores, mode, q, frac):
    if not scores:
        return False
    if mode == "mean":
        return (sum(scores) / len(scores)) >= q
    if mode == "min":
        return min(scores) >= q
    if mode == "frac":
        good = sum(s >= q for s in scores) / len(scores)
        return good >= frac
    raise ValueError("unknown mode")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fastq", help="FASTQ or FASTQ.GZ")
    ap.add_argument("k", type=int, help="Number of reads to sample")
    ap.add_argument("--q", type=int, default=30, help="Q threshold (default 30)")
    ap.add_argument("--mode", choices=["mean", "min", "frac"], default="mean",
                    help="How to apply Q threshold across a read")
    ap.add_argument("--frac", type=float, default=0.9,
                    help="For --mode frac: required fraction of bases >= Q (default 0.9)")
    ap.add_argument("--seed", type=int, default=None)
    args = ap.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    reservoir = []
    seen_passing = 0

    with open_maybe_gz(args.fastq) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline()
            p = fh.readline()
            ql = fh.readline()
            if not ql:
                break

            scores = phred_scores(ql)
            if not passes(scores, args.mode, args.q, args.frac):
                continue

            seen_passing += 1
            rec = (h, s, p, ql)

            if len(reservoir) < args.k:
                reservoir.append(rec)
            else:
                j = random.randrange(seen_passing)
                if j < args.k:
                    reservoir[j] = rec

    if len(reservoir) < args.k:
        print(f"WARNING: only {len(reservoir)} reads passed filter; requested {args.k}.",
              file=sys.stderr)

    out = sys.stdout
    for h, s, p, ql in reservoir:
        out.write(h); out.write(s); out.write(p); out.write(ql)

if __name__ == "__main__":
    main()