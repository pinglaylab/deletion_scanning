rawdata="/path/to/rawdata/"

# ==== Mapping pipeline 1 ====
  # (1) UMI extraction and quality filtering with fastp.
  # The UMI is the first 10 bases of read2.
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

# (2) Trim CS1 #TTGCTAGGACCGGCCTTAAAGC 
while read -r SAMPLE; do
  cutadapt \
    -G "^TTGCTAGGACCGGCCTTAAAGC" \
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
    --thread 4
done < $rawdata/samples.txt

# (4) Check for correct sequence beyond barcode (will be different for constructs) # GTGGCCACGCA GCTCAC
while read -r SAMPLE; do
  cutadapt \
    -G "^GTGGCCACGCA" \
    --discard-untrimmed \
    -j 10 \
    -o "fastq/${SAMPLE}_bc_target_R1.fastq" \
    -p "fastq/${SAMPLE}_bc_target_R2.fastq" \
    "fastq/${SAMPLE}_bc_R1.fastq" \
    "fastq/${SAMPLE}_bc_R2.fastq"
done < $rawdata/samples.txt

# (5) Only keep reads that map into the PB/LV/SB ITRs
# Lenti TGCTAGAGATTTTCCACACT
# SB    CAGTTGAAGTCGGAAGTTTA
# PB    CCCTAGAAAGATAGTCTGCG

while read -r SAMPLE; do
  cutadapt \
    --cores 10 \
    --discard-untrimmed \
    -e 0.1 \
    -O 10 \
    -m 30 \
    -a CCCTAGAAAGATAGTCTGCG \
    -o "fastq/${SAMPLE}_ITR.fastq.gz" \
    "fastq/${SAMPLE}_bc_target_R1.fastq"
done < $rawdata/samples.txt

# Alignment
# trimmed reads
while read SAMPLE; do
    bwa mem -Y -t 12 /net/shendure/vol8/projects/jonas/nobackup/genome/genome.fa fastq/${SAMPLE}_ITR.fastq.gz > mapped/${SAMPLE}.sam.temp
done < $rawdata/samples.txt

# Combined processing loop: sort, dedup UMIs, index, and convert to bed
while read SAMPLE; do
    echo "Processing sample: ${SAMPLE}"
    
    echo "  Step 1: Sorting SAM file..."
    samtools sort -@ 8 -o mapped/${SAMPLE}.sam mapped/${SAMPLE}.sam.temp && \
    
    echo "  Step 2: Deduplicating UMIs..."
    umi_tools dedup -I mapped/${SAMPLE}.sam --umi-separator=: -S mapped/${SAMPLE}_dedup.bam && \
    
    echo "  Step 3: Converting BAM to SAM..."
    samtools view mapped/${SAMPLE}_dedup.bam > mapped/${SAMPLE}_dedup.sam && \
    
    echo "  Step 4: Indexing BAM file..."
    samtools index mapped/${SAMPLE}_dedup.bam && \
    
    echo "  Step 5: Converting SAM to BED..."
    sam2bed < mapped/${SAMPLE}_dedup.sam > mapped/${SAMPLE}.bed && \
    
    echo "  Step 6: Cleaning up temporary files..."
    rm mapped/${SAMPLE}.sam.temp mapped/${SAMPLE}.sam && \
    
    echo "  Completed processing for sample: ${SAMPLE}"
    echo "----------------------------------------"
done < $rawdata/samples.txt

# identify insertion sites
while read SAMPLE; do
    python /path/to/collapse_barcodes.py mapped/${SAMPLE}.bed
done < $rawdata/samples.txt

