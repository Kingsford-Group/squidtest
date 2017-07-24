#!/bin/bash

argv=("$@")
OutputDir=${argv[0]}

mkdir -p ${OutputDir}/Ensemble75
mkdir -p ${OutputDir}/RNAseq

# download genome sequence and annotation
cd ${OutputDir}/Ensemble75
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
echo "Gunzipping reference files"
gunzip Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz

# download RNA-seq
echo "Downloading RNA-seq fastq files"
cd ${OutputDir}/RNAseq
fastq-dump --split-files SRR2532344 # HCC1954
fastq-dump --split-files SRR925710 # HCC1954
fastq-dump --split-files SRR2532336 # HCC1395

cat SRR2532344_1.fastq SRR925710_1.fastq > RNAHCC1954_1.fastq
cat SRR2532344_2.fastq SRR925710_2.fastq > RNAHCC1954_2.fastq

gzip -k RNAHCC1954_1.fastq
gzip -k RNAHCC1954_2.fastq
gzip -k SRR2532336_1.fastq
gzip -k SRR2532336_2.fastq