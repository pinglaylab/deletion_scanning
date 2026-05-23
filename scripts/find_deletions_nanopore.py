"""
Script to identify translocations in a BAM file from nanopore sequencing.
Outputs a TSV file with translocation breakpoints and a BAM file with the reads supporting translocations.

For each read that aligns partly to an extra chromosome called 'reference'
and partly to another chromosome, this script outputs a table with:
  - read_name
  - start: last nucleotide on the 'reference' alignment
  - end: if new alignment is on '+' strand, the first nucleotide; if on '-' strand, the last nucleotide of that alignment
  - chr_2: new chromosome name
  - strand_2: strand of the new chromosome alignment ('+' or '-')
  - barcode: extracted as the last underscore-delimited field from the read name

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

def main():
    args = parse_args()
    
    # Open the input BAM file and store its header for writing the output BAM.
    bam = pysam.AlignmentFile(args.bamfile, "rb")
    header = bam.header
    # Group alignments by read name.
    read_groups = defaultdict(list)
    for aln in bam:
        # Consider both primary and supplementary alignments.
        read_groups[aln.query_name].append(aln)
    bam.close()

    # Open the output TSV and BAM files.
    with open(args.output_tsv, "w", newline="") as tsvfile, \
         pysam.AlignmentFile(args.output_bam, "wb", header=header) as bam_out:
        
        writer = csv.writer(tsvfile, delimiter="\t")
        # Write header row.
        writer.writerow(["read_name", "start", "end", "chr_2", "strand_2", "barcode"])
        
        # Process each read group.
        for read_name, alignments in read_groups.items():
            ref_aln = None
            new_aln = None
            # Identify the alignment on the 'reference' chromosome and the one on the new chromosome.
            for aln in alignments:
                if aln.reference_name == "reference":
                    ref_aln = aln
                else:
                    new_aln = aln
            # Proceed only if both alignments are present.
            if ref_aln is None or new_aln is None:
                continue

            # Write all alignments of the read to the output BAM file.
            for aln in alignments:
                bam_out.write(aln)

            # Compute breakpoint positions.
            start = ref_aln.reference_end  # Last mapped base from the 'reference' alignment.
            # Determine "end" based on strand: first nucleotide if '+'; last nucleotide if '-'
            if new_aln.is_reverse:
                end = new_aln.reference_end  # Last nucleotide (1-indexed, since reference_end is 0-indexed exclusive).
            else:
                end = new_aln.reference_start + 1  # First nucleotide (1-indexed).

            chr_2 = new_aln.reference_name
            strand_2 = '-' if new_aln.is_reverse else '+'
            barcode = read_name.split("_")[-1]

            writer.writerow([read_name, start, end, chr_2, strand_2, barcode])

if __name__ == "__main__":
    main()