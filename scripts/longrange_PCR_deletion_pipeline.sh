# (1) Extract the amplicon
seqkit amplicon -I --max-mismatch 3 -j 8 -F "GTAAGCCGGTCCACCAACAAAG" -R "CCTTAGCCGCTAATAGGTGAGC" delamp_c4_GFPn2k_reform.fastq | seqkit seq -m 100 > fastq/delamp_c4_GFPn2k.fastq
seqkit amplicon -I --max-mismatch 3 -j 8 -F "GTAAGCCGGTCCACCAACAAAG" -R "CCTTAGCCGCTAATAGGTGAGC" delamp_c4_wt_reform.fastq | seqkit seq -m 100 > fastq/delamp_c4_wt.fastq


# (2) Remove primer
cutadapt \
  -g "GTAAGCCGGTCCACCAACAAAG" \
  --discard-untrimmed \
  --buffer-size=32000000 \
  --length 20000 \
  -j 10 \
  -o fastq/delamp_c4_GFPn2k_trimmed.fastq.gz \
  fastq/delamp_c4_GFPn2k.fastq

cutadapt \
  -g "GTAAGCCGGTCCACCAACAAAG" \
  --discard-untrimmed \
  --buffer-size=32000000 \
  --length 20000 \
  -j 10 \
  -o fastq/delamp_c4_wt_trimmed.fastq.gz \
  fastq/delamp_c4_wt.fastq

# (3) Extract UMI
fastp \
        --umi \
        --umi_loc=read1 \
        --umi_len=10 \
        --umi_delim "_" \
        --disable_adapter_trimming \
        --disable_length_filtering \
        -i fastq/delamp_c4_GFPn2k_trimmed.fastq.gz \
        -o fastq/delamp_c4_GFPn2k_UMI.fastq.gz \
        -Q \
        --thread 8

fastp \
        --umi \
        --umi_loc=read1 \
        --umi_len=10 \
        --umi_delim "_" \
        --disable_adapter_trimming \
        --disable_length_filtering \
        -i fastq/delamp_c4_wt_trimmed.fastq.gz \
        -o fastq/delamp_c4_wt_UMI.fastq.gz \
        -Q \
        --thread 8

# (4) Align to genome
minimap2 -a -x map-hifi -s40 -t 8 /net/shendure/vol8/projects/jonas/nobackup/genome/genome_transposon.fa fastq/delamp_c4_GFPn2k_UMI.fastq.gz > mapped/delamp_c4_GFPn2k.sam
    samtools sort -@ 8 mapped/delamp_c4_GFPn2k.sam -o mapped/delamp_c4_GFPn2k.bam
    samtools index -@ 8 mapped/delamp_c4_GFPn2k.bam
    rm mapped/delamp_c4_GFPn2k.sam

minimap2 -a -x map-hifi -s40 -t 8 /net/shendure/vol8/projects/jonas/nobackup/genome/genome_transposon.fa fastq/delamp_c4_wt_UMI.fastq.gz > mapped/delamp_c4_wt.sam
    samtools sort -@ 8 mapped/delamp_c4_wt.sam -o mapped/delamp_c4_wt.bam
    samtools index -@ 8 mapped/delamp_c4_wt.bam
    rm mapped/delamp_c4_wt.sam

# (5) Call split reads and align filtered bam files for inspection
python find_deletions_delamp.py mapped/delamp_c4_GFPn2k.bam deletions/delamp_c4_GFPn2k.tsv deletions/delamp_c4_GFPn2k.bam
python find_deletions_delamp.py mapped/delamp_c4_wt.bam deletions/delamp_c4_wt.tsv deletions/delamp_c4_wt.bam

samtools sort -o deletions/delamp_c4_GFPn2k_del.bam deletions/delamp_c4_GFPn2k.bam
samtools index deletions/delamp_c4_GFPn2k_del.bam

samtools sort -o deletions/delamp_c4_wt_del.bam deletions/delamp_c4_wt.bam
samtools index deletions/delamp_c4_wt_del.bam
