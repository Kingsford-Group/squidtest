#!/bin/bash

#SBATCH -p RM-shared
#SBATCH -t 12:00:00
#SBATCH --ntasks-per-node 12

CaDir=/pghbio/cure/congma/TCGAData/prad

count=0
while read -r line; do
	((count++))

	if [ -e $CaDir/${line:0:12}/${line:0:15}/StarAlign/ValidateSQUID/Aligned.bam ] && [ -e $CaDir/${line:0:12}/${line:0:15}/WGS_original.bam ] && [ ! -e $CaDir/${line:0:12}/${line:0:15}/StarAlign/ValidateSQUID/Better.bam ]; then
		echo -e $count"\t"$line

		python3 CompareReMap.py $CaDir/${line:0:12}/${line:0:15}/StarAlign/ValidateSQUID/sqfilter_fusedjunc.fa $CaDir/${line:0:12}/${line:0:15}/StarAlign/ValidateSQUID/Aligned.bam $CaDir/${line:0:12}/${line:0:15}/WGS_original.bam $CaDir/${line:0:12}/${line:0:15}/StarAlign/ValidateSQUID/Better.bam
		x=$(samtools view $CaDir/${line:0:12}/${line:0:15}/StarAlign/ValidateSQUID/Better.bam | cut -f3 | sort | uniq | wc -l)
		y=$(wc -l $CaDir/${line:0:12}/${line:0:15}/StarAlign/sqfilter_sv.txt)
		echo -e $x"\t"$y
	fi

	if [[ ${line:12:3} == "-11" ]]; then
		continue
	fi
	
	if [ -e $CaDir/${line:0:12}/${line:0:15}/StarAlign/ValidateSQUID/Better.bam ]; then
		#echo -e $count"\t"$line
		python3 WGSValidationAccy.py $CaDir/${line:0:12}/${line:0:15}/StarAlign/sqfilter_fusedjunc.fa $CaDir/${line:0:12}/${line:0:15}/StarAlign/ValidateSQUID/Better.bam
	fi
done < $CaDir/Barcodes.txt
