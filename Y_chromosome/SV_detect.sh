#!/bin/bash
#SBATCH -J chrY
#SBATCH -N 1
#SBATCH -n 30

breed=
for i in $breed
minimap2 \
-a ref.fasta \
$i.fasta.gz --MD -R '@RG\tID:'$i'\tSM:'$i'' -t 30 | \
samtools view -@ 30 -bS | samtools sort -@ 30 \
-o $i.sorted.bam
samtools index $i.sorted.bam
sniffles --input $i.sorted.bam \
--threads 30 --snf $i.snf
done

name=$(ls *.snf | xargs -I {} printf {}" ")
sniffles --reference ref.fasta \
--input $name \
--vcf multisample.vcf
