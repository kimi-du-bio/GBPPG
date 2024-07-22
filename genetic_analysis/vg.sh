#!/bin/bash

## construct the graph-based pangenome
vg autoindex --workflow giraffe \
-i Presence.fasta \
-t 50 -r reference.fasta \
-v SV.vcf.gz \
-p pig

## map reads to the pangenome using giraffe
vg giraffe -p -t 56 -Z pig.giraffe.gbz \
-d pig.dist -m pig.min \
-f reads_1.fq.gz -f reads_2.fq.gz > sample.gam

vg pack -Q 5 -s 5 -t 56 -x pig.giraffe.gbz -g sample.gam -o sample.pack

##genotype
vg call -t 56 pig.giraffe.gbz -r pig.snarls -k sample.pack --progress -z > sample.vcf
