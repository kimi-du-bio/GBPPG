### Here lists the pipeline of gene clustering and pan-genome construction

run pan_genome.sh to perform Markov-based gene clustering and construction of pan-genome.

The input are amino-acid sequences of protein-coding genes of the 13 tomato genomes.

The output will be a n × m gene matrix, where n is the number of gene families (pipeline produces it) and m is the number of genomes (13 here).

```
sh orthofinder.sh
python3 convert_orthogroup_2_pan_matrix.py Orthogroups.tsv
```
