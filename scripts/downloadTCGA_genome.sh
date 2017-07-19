#!/bin/bash
threads=${argv[0]}
OutputDir=${argv[1]}

mkdir -p ${OutputDir}/
mkdir -p ${OutputDir}/StarIndex

# download genome sequence and annotation
cd ${OutputDir}/
wget ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz
gunzip Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz

# download RNA-seq
STAR --runThreadN ${threads} --runMode genomeGenerate --genomeDir StarIndex/ --genomeFastaFiles Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.87.gtf