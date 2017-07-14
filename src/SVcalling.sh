#!/bin/bash

args=("$@")
for (( i=0; i<${#args[@]}; )); do
	if [ ${args[$i]} = "-p" ]; then
		if [[ ${args[$(($i+1))]} = /* || ${args[$(($i+1))]} = ~* ]]; then
			ProjectDir=${args[$(($i+1))]}
		else
			ProjectDir=$(pwd)/${args[$(($i+1))]}
		fi
		((i+=2))
	fi
done

# ILP Rearrangement
echo "Running SQUID on STAR alignment"
for (( i=3; i<10; i++)); do
	squid -b $ProjectDir/Alignments/Star_rearranged/Aligned.sortedByCoord.out.bam -c $ProjectDir/Alignments/Star_rearranged/Chimeric.out.bam -w $i -o $ProjectDir/Alignments/Star_rearranged/squid_w${i}
done
echo "Running SQUID on SpeedSeq alignment"
for (( i=3; i<10; i++)); do
	squid --bwa -b $ProjectDir/Alignments/SpeedSeq/Aligned.bam -w $i -o $ProjectDir/Alignments/SpeedSeq/squid_w${i}
done

# SV metrics (ILP, Delly, Lumpy)
cd $ProjectDir
mkdir -p $ProjectDir/SVcall
mkdir -p $ProjectDir/SVcall/DELLY_BWA
mkdir -p $ProjectDir/SVcall/DELLY_STAR
mkdir -p $ProjectDir/SVcall/LUMPY
mkdir -p $ProjectDir/SVcall/TransAbyss
mkdir -p $ProjectDir/SVcall/Trinity

# Delly_STAR
delly call -t DEL -o $ProjectDir/SVcall/DELLY_STAR/delly_DEL.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/Star_rearranged/Merged.bam
delly call -t DUP -o $ProjectDir/SVcall/DELLY_STAR/delly_DUP.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/Star_rearranged/Merged.bam
delly call -t TRA -o $ProjectDir/SVcall/DELLY_STAR/delly_TRA.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/Star_rearranged/Merged.bam
delly call -t INV -o $ProjectDir/SVcall/DELLY_STAR/delly_INV.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/Star_rearranged/Merged.bam
delly call -t INS -o $ProjectDir/SVcall/DELLY_STAR/delly_INS.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/Star_rearranged/Merged.bam
bcftools merge --force-samples -m id -O v -o SVcall/DELLY_STAR/delly_SV.vcf SVcall/DELLY_STAR/delly_*.bcf
python3 bin/VerifySVpred.py 2 WholeGenome/SV_newpos.txt SVcall/DELLY_STAR/delly_SV.vcf SVcall/DELLY_STAR/hit.txt
# Delly_BWA
delly call -t DEL -o $ProjectDir/SVcall/DELLY_BWA/delly_DEL.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/SpeedSeq/Aligned.bam
delly call -t DUP -o $ProjectDir/SVcall/DELLY_BWA/delly_DUP.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/SpeedSeq/Aligned.bam
delly call -t TRA -o $ProjectDir/SVcall/DELLY_BWA/delly_TRA.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/SpeedSeq/Aligned.bam
delly call -t INV -o $ProjectDir/SVcall/DELLY_BWA/delly_INV.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/SpeedSeq/Aligned.bam
delly call -t INS -o $ProjectDir/SVcall/DELLY_BWA/delly_INS.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/SpeedSeq/Aligned.bam
bcftools merge --force-samples -m id -O v -o $ProjectDir/SVcall/DELLY_BWA/delly_SV.vcf $ProjectDir/SVcall/DELLY_BWA/delly_*.bcf
python3 bin/VerifySVpred.py 2 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/SVcall/DELLY_BWA/delly_SV.vcf $ProjectDir/SVcall/DELLY_BWA/hit.txt

# Lumpy
for (( i=3; i<10; i++)); do
	lumpyexpress -B $ProjectDir/Alignments/SpeedSeq/Aligned.bam -D $ProjectDir/Alignments/SpeedSeq/Aligned.discordants.bam -S $ProjectDir/Alignments/SpeedSeq/Aligned.splitters.bam -o $ProjectDir/SVcall/LUMPY/lumpy_thresh$i.vcf -m $i
	python3 bin/VerifySVpred.py 2 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/SVcall/LUMPY/lumpy_thresh$i.vcf $ProjectDir/SVcall/LUMPY/hit_thresh$i.txt
done

# preparing for gmap
mkdir -p $ProjectDir/WholeGenome/Chromosomes/
if [ ! -e $ProjectDir/WholeGenome/Chromosomes/chr1.fa ]; then
	python3 bin/ReadGenome.py $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/WholeGenome/Chromosomes/
fi
if [ ! -d $ProjectDir/WholeGenome/Chromosomes/"genome"${ProjectDir##*/} ]; then
	gmap_build -D $ProjectDir/WholeGenome/Chromosomes/ -d "genome"${ProjectDir##*/} $ProjectDir/WholeGenome/Chromosomes/*.fa
fi

# transabyss + mummer3
transabyss -k 25 --pe $ProjectDir/Reads/RNA1.fq.gz $ProjectDir/Reads/RNA2.fq.gz --outdir $ProjectDir/SVcall/TransAbyss/k25 --name k25 --threads 4 --island 0
transabyss -k 32 --pe $ProjectDir/Reads/RNA1.fq.gz $ProjectDir/Reads/RNA2.fq.gz --outdir $ProjectDir/SVcall/TransAbyss/k32 --name k32 --threads 4 --island 0
transabyss-merge --mink 25 --maxk 32 --prefixes k25. k32. --out $ProjectDir/SVcall/TransAbyss/mergedassembly.fa $ProjectDir/SVcall/TransAbyss/k25/k25-final.fa $ProjectDir/SVcall/TransAbyss/k32/k32-final.fa
nucmer -p $ProjectDir/SVcall/TransAbyss/nucm $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/SVcall/TransAbyss/mergedassembly.fa
show-coords $ProjectDir/SVcall/TransAbyss/nucm.delta > $ProjectDir/SVcall/TransAbyss/nucm_text.txt
./bin/NucmerSV2 $ProjectDir/SVcall/TransAbyss/nucm_text.txt $ProjectDir/SVcall/TransAbyss/nucm_res.bedpe
python3 bin/VerifySVpred.py 5 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/SVcall/TransAbyss/nucm_res.bedpe $ProjectDir/SVcall/TransAbyss/nucm_res_hit.bedpe

# transabyss + gmap
gmap -D $ProjectDir/WholeGenome/Chromosomes/ -d "genome"${ProjectDir##*/} $ProjectDir/SVcall/TransAbyss/mergedassembly.fa > $ProjectDir/SVcall/TransAbyss/gmap.out
./bin/GmapSV $ProjectDir/SVcall/TransAbyss/gmap.out $ProjectDir/SVcall/TransAbyss/gmap_res.bedpe
python3 bin/VerifySVpred.py 5 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/SVcall/TransAbyss/gmap_res.bedpe $ProjectDir/SVcall/TransAbyss/gmap_res_hit.bedpe

# trinity + mummer3
Trinity --seqType fq --max_memory 50G --left $ProjectDir/Reads/RNA1.fq.gz --right $ProjectDir/Reads/RNA2.fq.gz --CPU 6 --output $ProjectDir/SVcall/Trinity
nucmer -p $ProjectDir/SVcall/Trinity/nucm $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/SVcall/Trinity/Trinity.fasta
show-coords $ProjectDir/SVcall/Trinity/nucm.delta > $ProjectDir/SVcall/Trinity/nucm_text.txt
./bin/NucmerSV2 $ProjectDir/SVcall/Trinity/nucm_text.txt $ProjectDir/SVcall/Trinity/nucm_res.bedpe
python3 bin/VerifySVpred.py 5 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/SVcall/Trinity/nucm_res.bedpe $ProjectDir/SVcall/Trinity/nucm_res_hit.bedpe

# trinity + gmap
gmap -D $ProjectDir/WholeGenome/Chromosomes/ -d "genome"${ProjectDir##*/} $ProjectDir/SVcall/Trinity/Trinity.fasta > $ProjectDir/SVcall/Trinity/gmap.out
./bin/GmapSV $ProjectDir/SVcall/Trinity/gmap.out $ProjectDir/SVcall/Trinity/gmap_res.bedpe
python3 bin/VerifySVpred.py 5 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/SVcall/Trinity/gmap_res.bedpe $ProjectDir/SVcall/Trinity/gmap_res_hit.bedpe