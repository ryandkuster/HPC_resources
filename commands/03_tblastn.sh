#!/bin/bash

gunzip data/reference/GCF_000013425.1_ASM1342v1_protein.faa.gz

module load blast-plus/2.12.0
tblastn -query data/reference/GCF_000013425.1_ASM1342v1_protein.faa -evalue .1 -db data/derived/e_coli_db -out data/derived/results.txt -outfmt 6 > logs/tblastn.log 2>&1
cut -f 3 data/derived/results.txt | cut -d "." -f 1 | sort | uniq -c | sort -nk1,1 > results/percent_id.txt
