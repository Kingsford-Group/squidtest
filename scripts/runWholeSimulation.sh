#!/bin/bash

OutDir=./
GenomeDir=./genome
AnnotationFile=./annotation/annot.gtf
BamtoolsDir=./bamtools

mkdir -p OutDir


for ((i=200}; i<900; i+=300)); do
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
