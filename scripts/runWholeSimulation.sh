#!/bin/bash

args=("$@")
OutDir=${args[0]}

mkdir -p OutDir
mkdir -p $OutDir/GRCh38
mkdir -p $OutDir/GRCh38/Chromosomes
mkdir -p $OutDir/GRCh38/Annotation

# download genome and reference
cd $OutDir/GRCh38/Chromosomes
for ((i=1; i<=5; i++)); do
	wget ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${i}.fa.gz
	gunzip Homo_sapiens.GRCh38.dna.chromosome.${i}.fa.gz
done
cd ../Annotation
wget ftp://ftp.ensembl.org/pub/release-89/gtf/homo_sapiens/Homo_sapiens.GRCh38.89.gtf.gz
gunzip Homo_sapiens.GRCh38.89.gtf.gz
awk '{if(length($1)==1 && $1>=1 && $1<=5)) print $0}' Homo_sapiens.GRCh38.89.gtf > genes.gtf
awk 'BEGIN{FS="\t";OFS="\t"}{split($NF,a," ");pfx="";s="";for(i=1;i<=length(a);i+=2){if(a[i]=="transcript_id"){pfx=a[i]" "a[i+1]}else{s=s" "a[i]" "a[i+1]}}if(pfx==""){print "[WARN] line "NR" without transcript_id!" > "/dev/stderr"}else{$NF=pfx""s;print$0} }' genes.gtf > genes_clean.gtf

GenomeDir=${OutDir}/GRCh38/Chromosomes
AnnotationFile=${OutDir}/GRCh38/Annotation/genes_clean.gtf

for ((i=200}; i<900; i+=300)); do
	mkdir -p ${OutDir}/SVRNA${i}_1
	mkdir -p ${OutDir}/SVRNA${i}_2
	mkdir -p ${OutDir}/SVRNA${i}_3
	mkdir -p ${OutDir}/SVRNA${i}_4
	./bin/SimSVGenome.sh -g ${GenomeDir} -p ${OutDir}/SVRNA${i}_1 -a ${AnnotationFile} -seq 21 -inv ${i} -ins ${i} -del ${i} -dup ${i} -tra 2 -RNA
	./bin/SimSVGenome.sh -g ${GenomeDir} -p ${OutDir}/SVRNA${i}_2 -a ${AnnotationFile} -seq 25 -inv ${i} -ins ${i} -del ${i} -dup ${i} -tra 2 -RNA
	./bin/SimSVGenome.sh -g ${GenomeDir} -p ${OutDir}/SVRNA${i}_3 -a ${AnnotationFile} -seq 28 -inv ${i} -ins ${i} -del ${i} -dup ${i} -tra 2 -RNA
	./bin/SimSVGenome.sh -g ${GenomeDir} -p ${OutDir}/SVRNA${i}_4 -a ${AnnotationFile} -seq 30 -inv ${i} -ins ${i} -del ${i} -dup ${i} -tra 2 -RNA
done

for ((i=200; i<900; i+=300)); do
	./bin/SVcalling.sh -p ${OutDir}/SVRNA${i}_1
	./bin/SVcalling.sh -p ${OutDir}/SVRNA${i}_2
	./bin/SVcalling.sh -p ${OutDir}/SVRNA${i}_3
	./bin/SVcalling.sh -p ${OutDir}/SVRNA${i}_4
done
