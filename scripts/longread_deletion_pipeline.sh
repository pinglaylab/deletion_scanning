
# ==== (1) Processing of basecalled fastq ====
while read SAMPLE; do
    # extract the amplicon
    seqkit amplicon -I --max-mismatch 3 -j 8 \
    -R "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG" \
    -F "CCTTAGCCGCTAATAGGTGAGC" \
    ${SAMPLE}.fastq \
    -r 33:-23 > ${SAMPLE}_seqkit.fastq 

    fastp \
    --umi \
    --umi_loc=read1 \
    --umi_len=12 \
    --umi_delim "_" \
    -i ${SAMPLE}_seqkit.fastq \
    -o ${SAMPLE}_bc.fastq \
    --average_qual 15 \
    --thread 4

    rm ${SAMPLE}_seqkit.fastq 
done < samples.txt

# map the reads to the genome
while read SAMPLE; do
    minimap2 -a -x sr genome_PB.mmi ${SAMPLE}_bc.fastq > mapped/${SAMPLE}_bc.sam
    samtools sort mapped/${SAMPLE}_bc.sam -o mapped/${SAMPLE}_bc.bam
    samtools index mapped/${SAMPLE}_bc.bam
    rm mapped/${SAMPLE}_bc.sam
done < samples.txt

# find the deletions
while read SAMPLE; do
    python ./find_deletions_nanopore.py mapped/${SAMPLE}_bc.bam deletions/${SAMPLE}_deletions.tsv deletions/${SAMPLE}_deletions.bam
done < samples_d27_d3.txt

# sort deletion bams
while read SAMPLE; do
    samtools sort -o deletions/${SAMPLE}_del.bam deletions/${SAMPLE}_deletions.bam
    samtools index deletions/${SAMPLE}_del.bam
    rm deletions/${SAMPLE}_deletions.bam
done < samples_d27_d3.txt
