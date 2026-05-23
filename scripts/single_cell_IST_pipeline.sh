rawdata="/path/to/rawdata/"

# (1) Trim CS2 (22nt)
while read -r SAMPLE; do
    R1=raw/"${SAMPLE}"_S*_L001_R1_001.fastq.gz
    R2=raw/"${SAMPLE}"_S*_L001_R2_001.fastq.gz

    cutadapt \
        -U 22 \
        -j 10 \
        -o "fastq/${SAMPLE}_trimmed_R1.fastq.gz" \
        -p "fastq/${SAMPLE}_trimmed_R2.fastq.gz" \
        $R1 \
        $R2
done < $rawdata/samples_scSad.txt

# (2) Add beacon barcode to read name
while read -r SAMPLE; do
    fastp \
        --umi \
        --umi_loc=read2 \
        --umi_len=12 \
        --umi_delim "_" \
        --disable_adapter_trimming \
        --disable_length_filtering \
        --trim_poly_g \
        -i fastq/${SAMPLE}_trimmed_R1.fastq.gz \
        -I fastq/${SAMPLE}_trimmed_R2.fastq.gz \
        -o fastq/${SAMPLE}_bc_R1.fastq.gz \
        -O fastq/${SAMPLE}_bc_R2.fastq.gz \
        -Q \
        --thread 8
done < $rawdata/samples_scSad.txt


# (3) Extract cell barcode
while read -r SAMPLE; do
    fastp \
        --umi \
        --umi_loc=read1 \
        --umi_len=16 \
        --umi_delim "_" \
        --disable_adapter_trimming \
        --disable_length_filtering \
        -i fastq/${SAMPLE}_bc_R1.fastq.gz \
        -I fastq/${SAMPLE}_bc_R2.fastq.gz \
        -o fastq/${SAMPLE}_CB_R1.fastq \
        -O fastq/${SAMPLE}_CB_R2.fastq \
        -Q \
        --thread 4
done < $rawdata/samples_scSad.txt

# (4) Extract UMI
while read -r SAMPLE; do
    fastp \
        --umi \
        --umi_loc=read1 \
        --umi_len=12 \
        --umi_delim "_" \
        --disable_adapter_trimming \
        --disable_length_filtering \
        -i fastq/${SAMPLE}_CB_R1.fastq \
        -I fastq/${SAMPLE}_CB_R2.fastq \
        -o fastq/${SAMPLE}_CB_UMI_R1.fastq \
        -O fastq/${SAMPLE}_CB_UMI_R2.fastq \
        -Q \
        --thread 4
done < $rawdata/samples_scSad.txt

# (5) Remove poly(A) tail
while read -r SAMPLE; do
    cutadapt --poly-a \
        --minimum-length 20 \
        -o fastq/${SAMPLE}_R2.fastq -p fastq/${SAMPLE}.fastq \
        fastq/${SAMPLE}_CB_UMI_R2.fastq fastq/${SAMPLE}_CB_UMI_R1.fastq
done < $rawdata/samples_scSad.txt

# (6) Align to genome
while read -r SAMPLE; do
    bwa mem -Y -t 12 /net/shendure/vol8/projects/jonas/nobackup/genome/genome_PB_GFP_v3.fa fastq/${SAMPLE}_CB_UMI_R2.fastq > mapped/${SAMPLE}.sam.tmp
done < $rawdata/samples_scSad.txt

# (7) Sort BAM files and convert to bed
while read -r SAMPLE; do
    samtools sort -@ 8 -o mapped/${SAMPLE}.bam mapped/${SAMPLE}.sam.tmp && samtools index mapped/${SAMPLE}.bam && \
    samtools sort -@ 8 -o mapped/${SAMPLE}.sam mapped/${SAMPLE}.sam.tmp && \
    sam2bed < mapped/${SAMPLE}.sam > mapped/${SAMPLE}.bed && \
    gzip mapped/${SAMPLE}.bed 
done < $rawdata/samples_scSad.txt