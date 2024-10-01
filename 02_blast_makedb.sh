#!/bin/bash

gunzip data/reference/GCF_000005845.2_ASM584v2_genomic.fna.gz

module load blast-plus/2.12.0
makeblastdb -in data/reference/GCF_000005845.2_ASM584v2_genomic.fna -dbtype nucl -out data/derived/e_coli_db > logs/makeblastdb.log 2>&1
