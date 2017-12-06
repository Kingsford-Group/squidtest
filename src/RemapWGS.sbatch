#!/bin/bash

#SBATCH -p RM-shared
#SBATCH -t 12:00:00
#SBATCH --ntasks-per-node 12

CaDir=/pghbio/cure/congma/TCGAData/blca

count=0
while read -r line; do
	((count++))

	if [ -e $CaDir/${line:0:12}/${line:0:15}/WGS_1.fastq.gz ] && [ -e $CaDir/${line:0:12}/${line:0:15}/WGS_2.fastq.gz ] && [ ! -e $CaDir/${line:0:12}/${line:0:15}/StarAlign/ValidateSQUID/Aligned.bam ]; then
		echo -e $count"\t"$line

		if [ ! -e $CaDir/${line:0:12}/${line:0:15}/StarAlign/sqfilter_fusedjunc.fa ]; then
			python /pghbio/cure/congma/Code/FusedGenome.py 1 $CaDir/${line:0:12}/${line:0:15}/StarAlign/Aligned.sortedByCoord.out.bam $CaDir/${line:0:12}/${line:0:15}/StarAlign/sqfilter_sv.txt /pghbio/cure/congma/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa $CaDir/${line:0:12}/${line:0:15}/StarAlign/sqfilter_fusedjunc.fa
		fi

		mkdir -p $CaDir/${line:0:12}/${line:0:15}/StarAlign/ValidateSQUID/
		if [ ! -e $CaDir/${line:0:12}/${line:0:15}/StarAlign/sqfilter_fusedjunc.fa ]; then
			ln -s $CaDir/${line:0:12}/${line:0:15}/StarAlign/sqfilter_fusedjunc.fa $CaDir/${line:0:12}/${line:0:15}/StarAlign/ValidateSQUID/
		fi

		speedseq align -t 12 -o $CaDir/${line:0:12}/${line:0:15}/StarAlign/ValidateSQUID/Aligned -R "@RG\tID:id\tSM:samplename\tLB:lib" $CaDir/${line:0:12}/${line:0:15}/StarAlign/ValidateSQUID/sqfilter_fusedjunc.fa $CaDir/${line:0:12}/${line:0:15}/WGS_1.fastq.gz $CaDir/${line:0:12}/${line:0:15}/WGS_2.fastq.gz
	
	fi
done < $CaDir/Barcodes.txt
