#!/bin/bash

# Raw Reads Directory = $1
# Transciptome        = $2
# Experiment          = $3
# Results Dir         = $4

# Index the transcriptome
bin/kallisto index -i human_GRCh38.idx $2

# Run Kallisto Quantification for each set of fastq files
bin/kallisto quant --single -l 180 -s 20 -i human_GRCh38.idx -o kallisto/SRR493366 -b 100 tutorial/reads/SRR493366.fastq
bin/kallisto quant --single -l 180 -s 20 -i human_GRCh38.idx -o kallisto/SRR493367 -b 100 tutorial/reads/SRR493367.fastq
bin/kallisto quant --single -l 180 -s 20 -i human_GRCh38.idx -o kallisto/SRR493368 -b 100 tutorial/reads/SRR493368.fastq
bin/kallisto quant --single -l 180 -s 20 -i human_GRCh38.idx -o kallisto/SRR493369 -b 100 tutorial/reads/SRR493369.fastq
bin/kallisto quant --single -l 180 -s 20 -i human_GRCh38.idx -o kallisto/SRR493370 -b 100 tutorial/reads/SRR493370.fastq
bin/kallisto quant --single -l 180 -s 20 -i human_GRCh38.idx -o kallisto/SRR493371 -b 100 tutorial/reads/SRR493371.fastq

# Run Slueth
echo "bin/sleuth.R kallisto $3"
bin/sleuth.R kallisto $3


