#!/bin/bash
#SBATCH -J pilon
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 200

bwa index ./p_ctg.fa
bwa mem -t 200 ./p_ctg.fa \
./short_reads.fq.gz \
./short_reads.fq.gz \
| samtools sort -@ 200 -O bam -o align.bam
samtools index -@ 200 align.bam
sambamba markdup -t 200 align.bam \
align_markdup.bam
samtools view -q 30 -b -@ 200 align_markdup.bam \
> lign_filter.bam
samtools index -@ 200 align_filter.bam
java -Xmx2000g -jar pilon-1.23.jar \
--genome ./p_ctg.fa \
--frags align_filter.bam \
--fix snps,indels \
--output polished --vcf --threads 200
