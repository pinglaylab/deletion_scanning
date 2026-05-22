#!/usr/bin/env python3
import sys
import re
import pandas as pd
import operator
from collections import Counter

# define standard chromosomes
regular_chr = list(range(1,23)) + ["X", "Y", "M"]
chrom = ["chr" + str(s) for s in regular_chr]

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])

def clean_up(input_file):
    bedfile = open(input_file, "r")
    shortname = input_file.split(".")[0]
    temp = open(shortname + ".txt.temp", 'w')
    
    for line in bedfile:
        fields = line.rstrip("\n").split("\t")
        chr, start, end, name = fields[0:4]
        strand = fields[5]
        cigar = fields[7]
        seq = fields[11]

        # Only keep records where at least 40 bp are mapped.
        if int(end) - int(start) < 40:
            continue

        # extract cell barcode (after underscore) from the UMI_barcode string in name
        # name format: ...:UMI_cellbarcode
        umi_cell = name.rsplit(':', 1)[1]
        barcode = umi_cell.split('_', 1)[1]

        # Adjust position and sequence depending on strand
        if strand == "+":
            pos = end
        else:  # strand == "-"
            pos = start
            seq = reverse_complement(seq)

        # Only keep reads on standard chromosomes
        if chr in chrom:
            # Write fields: barcode, chr, pos, start, end, strand, cigar, seq
            temp.write("\t".join([barcode, chr, pos, start, end, strand, cigar, seq]) + "\n")

    temp.close()
    bedfile.close()


def sort_data(input_file):
    shortname = input_file.split(".")[0]
    filename = shortname + ".txt.temp"
    df = pd.read_table(filename, header=None, sep="\t")
    df.columns = ["barcode", "chr", "pos", "start", "end", "strand", "cigar", "seq"]

    # Sort by chr, pos, strand, barcode
    df = df.sort_values(["chr", "pos", "strand", "barcode"], ascending=True)
    df.to_csv(shortname + ".sort.txt.temp", sep="\t", header=False, index=False)


def duplicate_removal(input_file):
    shortname = input_file.split(".")[0]
    infile = open(shortname + ".sort.txt.temp", 'r')
    output = open(shortname + ".out.txt", 'w')
    
    barcode_dic = {}

    # Initialize with the first line
    first_line = infile.readline().rstrip("\n")
    if not first_line:
        infile.close()
        output.close()
        return

    fields = first_line.split("\t")
    prev_chr, prev_pos, prev_str = fields[1], fields[2], fields[5]
    barcode_dic[fields[0]] = 1
    dup_count = 1

    for line in infile:
        fields = line.rstrip("\n").split("\t")
        barcode, chr, pos, strand = fields[0], fields[1], fields[2], fields[5]

        if chr == prev_chr and pos == prev_pos and strand == prev_str:
            dup_count += 1
            barcode_dic[barcode] = barcode_dic.get(barcode, 0) + 1
        else:
            # write previous group
            sorted_bcs = sorted(barcode_dic.items(), key=operator.itemgetter(1), reverse=True)
            for bc, count in sorted_bcs:
                output.write("\t".join([prev_chr, prev_pos, prev_str, str(dup_count), bc, str(count)]) + "\n")
            # reset for new group
            barcode_dic = {barcode: 1}
            dup_count = 1
            prev_chr, prev_pos, prev_str = chr, pos, strand

    # write last group
    sorted_bcs = sorted(barcode_dic.items(), key=operator.itemgetter(1), reverse=True)
    for bc, count in sorted_bcs:
        output.write("\t".join([prev_chr, prev_pos, prev_str, str(dup_count), bc, str(count)]) + "\n")

    infile.close()
    output.close()

if __name__ == "__main__":
    input_file = sys.argv[1]
    clean_up(input_file)
    sort_data(input_file)
    duplicate_removal(input_file)
