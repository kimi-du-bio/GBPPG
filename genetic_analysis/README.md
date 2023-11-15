# SV genotyping using graph-based genome
`
sh graphtyper.sh
`

# Heritability analysis and Genome-wide association studies

## Data preparation

### 16s

### Gene expression

To quantify the expression of all genes, all expression results was normalized.

```R
Rscript script/prepare_gene_expression.R
```

### Genetic variants

The SNPs and SVs were performed quality control separately. In brief, the genetic variants were firstly removed as follows: (1) mean depth < 5×, (2) missing rate > 95%, and (3) minor allele frequency < 0.1. Only genetic variants on autosomes and X chromosomes were used.

```bash
vcftools --gzvcf raw.snp.1filtered.vcf.gz --min-meanDP 5 --recode -c | bgzip -@ 40 > snp.vcf.gz
plink --vcf snp.vcf.gz --make-bed --out snp

plink --bfile snp --mind 0.05 --maf 0.1 --hwe 1e-6 --make-bed --out snp_qc
```

```bash
vcftools --gzvcf raw.sv.1filtered.vcf.gz --min-meanDP 5 --recode -c | bgzip -@ 40 > sv.vcf.gz
plink --vcf sv.vcf.gz --make-bed --out sv

plink --bfile sv --mind 0.05 --maf 0.1 --hwe 1e-6 --make-bed --out sv_qc
```

For SNPs, we simultaneously performed LD-based pruning using PLINK.

```bash
plink --bfile snp_qc --indep-pairwise 50 10 0.1 --out snp_purge
```

Then, We merged SNPs and SVs for combining analysis.

```bash
plink --bfile snp_purge --bmerge sv_qc.bed sv_qc.bim sv_qc.fam --make-bed --out merge
```

## Heritability estimation

The LDAK-thin model was applied to estimate the proportion of phenotypic variance explained by genetic variants.

```bash
ldak --bfile merge --cut-weights merge --window-prune 0.98 --section-length 100

ldak --bfile merge --calc-weights-all merge

ldak  --calc-kins-direct LDAK-Thin --bfile merge --weights merge/weights.all --power -.5

ldak --pheno pheno.txt --mpheno -1  --grm LDAK-Thin  --reml 4 --constrain YES
```

## Genome-wide association study

We performed genome-wide association studies based on the framework of the Mixed linear model (MLM) method. Furthermore, we used the leave-one-chromosome-out (LOCO) method and the association studies implemented in GCTA (v1.93.2beta).

```bash
ldak --bfile merge --thin exp --window-prune .5 --window-kb 100

parallel -j18 ldak --calc-kins-direct exp{} --bfile merge --ignore-weights YES --power -0.5 --extract exp.in --chr {} ::: 1..18

parallel '(echo exp{} >> list.All)' ::: 1..18

ldak --add-grm expAll --mgrm list.All

parallel '(echo "expAll exp{}" > list.{}; ldak --sub-grm expN{} --mgrm list.{})' ::: 1..18

parallel gcta --mlma --bfile merge --chr {} --pheno pheno-16s.txt --out chr{} --mpheno -1 --grm G-merge/expN{} ::: 1..18
```

## eQTLs

We calculated the eQTLs using the MatrixEQTL package.

```bash
Rscript eQTL.R
```
