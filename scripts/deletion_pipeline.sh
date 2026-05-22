rawdata="/path/to/rawdata/"

# ==== Deletion pipeline ====
# (1) Attach UMI to read header
while read -r SAMPLE; do
    R1="$rawdata/${SAMPLE}"_S*_L001_R1_001.fastq.gz
    R2="$rawdata/${SAMPLE}"_S*_L001_R2_001.fastq.gz

    fastp \
        --umi \
        --umi_loc=read2 \
        --umi_len=10 \
        --disable_adapter_trimming \
        --trim_poly_g \
        -i $R1 \
        -I $R2 \
        -o "fastq/${SAMPLE}_fastp_R1.fastq" \
        -O "fastq/${SAMPLE}_fastp_R2.fastq" \
        --average_qual 15 \
        --thread 10
done < $rawdata/samples.txt

# (2) Trim CS2
while read -r SAMPLE; do
  cutadapt \
    -G "^CCTTAGCCGCTAATAGGTGAGC" \
    --discard-untrimmed \
    -j 10 \
    -o "fastq/${SAMPLE}_trimmed_R1.fastq" \
    -p "fastq/${SAMPLE}_trimmed_R2.fastq" \
    "fastq/${SAMPLE}_fastp_R1.fastq" \
    "fastq/${SAMPLE}_fastp_R2.fastq"
done < $rawdata/samples.txt

# (3) Find the beacon barcode and attach to the read header
while read -r SAMPLE; do
  fastp \
    --umi \
    --umi_loc=read2 \
    --umi_len=12 \
    --umi_delim "_" \
    --disable_adapter_trimming \
    --disable_length_filtering \
    --trim_poly_g \
    -i "fastq/${SAMPLE}_trimmed_R1.fastq" \
    -I "fastq/${SAMPLE}_trimmed_R2.fastq" \
    -o "fastq/${SAMPLE}_bc_R1.fastq" \
    -O "fastq/${SAMPLE}_bc_R2.fastq" \
    --average_qual 15 \
    --thread 4 && \
    rm fastq/${SAMPLE}_trimmed_R1.fastq && \
    rm fastq/${SAMPLE}_trimmed_R2.fastq
done < $rawdata/samples.txt

# (4) Check for correct sequence beyond barcode
while read -r SAMPLE; do
  cutadapt \
    -G "^GCTTTAAGGCC" \
    --discard-untrimmed \
    -j 10 \
    -o "fastq/${SAMPLE}_bc_GFP_R1.fastq.gz" \
    -p "fastq/${SAMPLE}_bc_GFP_R2.fastq.gz" \
    "fastq/${SAMPLE}_bc_R1.fastq" \
    "fastq/${SAMPLE}_bc_R2.fastq" && \
    rm fastq/${SAMPLE}_bc_R1.fastq && \
    rm fastq/${SAMPLE}_bc_R2.fastq
done < $rawdata/samples.txt

# (4) Align to the genome
while read SAMPLE; do
    bwa mem -Y -t 12 /path/to/genome/genome_PB_GFP_v3.fa fastq/${SAMPLE}_bc_GFP_R1.fastq.gz > mapped/${SAMPLE}.sam.temp
done < $rawdata/samples.txt

# (5) Sort and index the SAM files
while read -r SAMPLE; do
    samtools sort -@ 8 -o mapped/${SAMPLE}.bam mapped/${SAMPLE}.sam.temp && \
    samtools index mapped/${SAMPLE}.bam && \
    rm mapped/${SAMPLE}.sam.temp
done < $rawdata/samples.txt

# (6) Extract all read locations with samtobed
while read -r SAMPLE; do
    samtools view "mapped/${SAMPLE}.bam" | sam2bed --reduced > "mapped/${SAMPLE}.bed" && \
    gzip mapped/${SAMPLE}.bed
done < "$rawdata/samples.txt"