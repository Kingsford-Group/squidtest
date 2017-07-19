#!/bin/bash

argv=("$@")
DataFolder=${argv[0]}
threads=${argv[1]}

echo "DataFolder: "${DataFolder}
echo "threads: "${threads}

for ((i=2; i<9; i+=3)); do
	for ((j=1; j<5; j++)); do
		mkdir -p $DataFolder/SVRNA${i}00_${j}/Alignments

		# Aligning reads with STAR
		if [ ! -d $DataFolder/SVRNA${i}00_${j}/WholeGenome/STAR_genome_rearranged ]; then
			mkdir -p $DataFolder/SVRNA${i}00_${j}/WholeGenome/STAR_genome_rearranged
			STAR --runThreadN ${threads} --runMode genomeGenerate --genomeDir $DataFolder/SVRNA${i}00_${j}/WholeGenome/STAR_genome_rearranged --genomeFastaFiles $DataFolder/SVRNA${i}00_${j}/WholeGenome/genome_rearranged.fa
		fi
		mkdir -p $DataFolder/SVRNA${i}00_${j}/Alignments/Star_rearranged
		STAR --runThreadN ${threads} --genomeDir $DataFolder/SVRNA${i}00_${j}/WholeGenome/STAR_genome_rearranged/ --readFilesIn $DataFolder/SVRNA${i}00_${j}/Reads/RNA1.fq.gz $DataFolder/SVRNA${i}00_${j}/Reads/RNA2.fq.gz --readFilesCommand gunzip -c --outFileNamePrefix $DataFolder/SVRNA${i}00_${j}/Alignments/Star_rearranged/ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --chimSegmentMin 20 --outSAMstrandField intronMotif --limitBAMsortRAM 21943468974
		samtools view -Shb $DataFolder/SVRNA${i}00_${j}/Alignments/Star_rearranged/Chimeric.out.sam -o $DataFolder/SVRNA${i}00_${j}/Alignments/Star_rearranged/Chimeric.out.bam
		samtools merge $DataFolder/SVRNA${i}00_${j}/Alignments/Star_rearranged/Merged_unsort.bam $DataFolder/SVRNA${i}00_${j}/Alignments/Star_rearranged/Aligned.sortedByCoord.out.bam $DataFolder/SVRNA${i}00_${j}/Alignments/Star_rearranged/Chimeric.out.bam
		samtools sort $DataFolder/SVRNA${i}00_${j}/Alignments/Star_rearranged/Merged_unsort.bam -o $DataFolder/SVRNA${i}00_${j}/Alignments/Star_rearranged/Merged.bam
		samtools index $DataFolder/SVRNA${i}00_${j}/Alignments/Star_rearranged/Merged.bam

		# Aligning reads with SpeedSeq
		mkdir -p $DataFolder/SVRNA${i}00_${j}/Alignments/SpeedSeq
		speedseq align -R "@RG\tID:id\tSM:samplename\tLB:lib" -t ${threads} -o $DataFolder/SVRNA${i}00_${j}/Alignments/SpeedSeq/Aligned $DataFolder/SVRNA${i}00_${j}/WholeGenome/genome_rearranged.fa $DataFolder/SVRNA${i}00_${j}/Reads/RNA1.fq.gz $DataFolder/SVRNA${i}00_${j}/Reads/RNA2.fq.gz

		# Call TSVs with different methods
		./bin/SVcalling.sh -p ${DataFolder}/SVRNA${i}00_${j} -threads ${threads}
	done
done