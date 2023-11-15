#!/usr/bin/env Rscript

library("MatrixEQTL")

OUT_PREFIX <- "sv"
TISSUE <- "data_sv"

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

## Location of the package with the data files.
base.dir = "/storage/data/"

# Genotype file name
SNP_file_name = paste0(base.dir, TISSUE, "/SNP.txt");
snps_location_file_name = paste0(base.dir, TISSUE, "/snpsloc.txt");

# Gene expression file name
expression_file_name="/storage/exp/GE.txt"
gene_location_file_name = "/storage/data/geneloc.txt"

# Covariates file name
# Set to character() for no covariates
covariates_file_name=character()

# Output file name
output_file_name_cis = paste0(OUT_PREFIX, ".cis_eqtl.txt");
output_file_name_tra = paste0(OUT_PREFIX, ".trans_eqtl.txt");

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-2;
pvOutputThreshold_tra = 1e-5;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();

# Distance for local gene-SNP pairs
cisDist = 50000;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = " ";
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = " "
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);


# detach('package:MatrixEQTL', unload=TRUE)
# install.packages('/gscmnt/gc2719/halllab/users/cchiang/src/MatrixEQTL', repos = NULL, type="source")
# library('MatrixEQTL')

me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = FALSE);
