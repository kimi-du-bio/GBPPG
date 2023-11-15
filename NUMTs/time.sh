#!/bin/bash
#SBATCH -J time
#SBATCH -N 1
#SABTCH -n 2

cat all.raw.len1 | while read line
do
	chr=$(echo $line | awk '{print $21}')
	start=$(echo $line | awk '{print $22}')
	end=$(echo $line | awk '{print $23}')
	id=$(echo $line | awk '{print $1}')
	breed=$(echo $line | awk '{print $25}')

	samtools faidx $breed.fasta \
	$chr:$start\-$end > fasta/$id.fasta
done
cd fasta/
for i in $(ls *.fasta);do cat Duroc_Warthog.fa $i > fasta_cat/$i;done
cd fasta_cat/
cat all.raw.len1 | awk '{print $1}' | while read line;do clustalo -i $line.fasta  -o $line.aln.fa --outfmt=fa --force;done
cat all.raw.len1 | awk '{print $1}' | while read line;do snp-sites -v -o $line.aln.fa.vcf $line.aln.fa;done
