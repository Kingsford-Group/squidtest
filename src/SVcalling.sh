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
mkdir -p $ProjectDir/TSVcall
mkdir -p $ProjectDir/TSVcall/SQUID_STAR
mkdir -p $ProjectDir/TSVcall/SQUID_BWA
for (( i=3; i<10; i++)); do
	squid -b $ProjectDir/Alignments/Star_rearranged/Aligned.sortedByCoord.out.bam -c $ProjectDir/Alignments/Star_rearranged/Chimeric.out.bam -w $i -o $ProjectDir/TSVcall/SQUID_STAR/squid_w${i}
	python3 bin/VerifySVpred.py 1 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/TSVcall/SQUID_STAR/squid_w${i}_sv.txt $ProjectDir/TSVcall/SQUID_STAR/squid_w${i}_sv_hit.txt
done
echo "Running SQUID on SpeedSeq alignment"
for (( i=3; i<10; i++)); do
	squid --bwa -b $ProjectDir/Alignments/SpeedSeq/Aligned.bam -w $i -o $ProjectDir/TSVcall/SQUID_BWA/squid_w${i}
	python3 bin/VerifySVpred.py 1 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/TSVcall/SQUID_BWA/squid_w${i}_sv.txt $ProjectDir/TSVcall/SQUID_BWA/squid_w${i}_sv_hit.txt
done

# SV metrics (ILP, Delly, Lumpy)
cd $ProjectDir
mkdir -p $ProjectDir/TSVcall/DELLY_BWA
mkdir -p $ProjectDir/TSVcall/DELLY_STAR
mkdir -p $ProjectDir/TSVcall/LUMPY
mkdir -p $ProjectDir/TSVcall/TransABySS
mkdir -p $ProjectDir/TSVcall/Trinity

# Delly_STAR
delly call -t DEL -o $ProjectDir/TSVcall/DELLY_STAR/delly_DEL.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/Star_rearranged/Merged.bam
delly call -t DUP -o $ProjectDir/TSVcall/DELLY_STAR/delly_DUP.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/Star_rearranged/Merged.bam
delly call -t TRA -o $ProjectDir/TSVcall/DELLY_STAR/delly_TRA.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/Star_rearranged/Merged.bam
delly call -t INV -o $ProjectDir/TSVcall/DELLY_STAR/delly_INV.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/Star_rearranged/Merged.bam
delly call -t INS -o $ProjectDir/TSVcall/DELLY_STAR/delly_INS.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/Star_rearranged/Merged.bam
bcftools merge --force-samples -m id -O v -o TSVcall/DELLY_STAR/delly_SV.vcf TSVcall/DELLY_STAR/delly_*.bcf
python3 bin/VerifySVpred.py 2 WholeGenome/SV_newpos.txt TSVcall/DELLY_STAR/delly_SV.vcf TSVcall/DELLY_STAR/hit.txt
# Delly_BWA
delly call -t DEL -o $ProjectDir/TSVcall/DELLY_BWA/delly_DEL.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/SpeedSeq/Aligned.bam
delly call -t DUP -o $ProjectDir/TSVcall/DELLY_BWA/delly_DUP.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/SpeedSeq/Aligned.bam
delly call -t TRA -o $ProjectDir/TSVcall/DELLY_BWA/delly_TRA.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/SpeedSeq/Aligned.bam
delly call -t INV -o $ProjectDir/TSVcall/DELLY_BWA/delly_INV.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/SpeedSeq/Aligned.bam
delly call -t INS -o $ProjectDir/TSVcall/DELLY_BWA/delly_INS.bcf -g $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/Alignments/SpeedSeq/Aligned.bam
bcftools merge --force-samples -m id -O v -o $ProjectDir/TSVcall/DELLY_BWA/delly_SV.vcf $ProjectDir/TSVcall/DELLY_BWA/delly_*.bcf
python3 bin/VerifySVpred.py 2 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/TSVcall/DELLY_BWA/delly_SV.vcf $ProjectDir/TSVcall/DELLY_BWA/hit.txt

# Lumpy
for (( i=3; i<10; i++)); do
	lumpyexpress -B $ProjectDir/Alignments/SpeedSeq/Aligned.bam -D $ProjectDir/Alignments/SpeedSeq/Aligned.discordants.bam -S $ProjectDir/Alignments/SpeedSeq/Aligned.splitters.bam -o $ProjectDir/TSVcall/LUMPY/lumpy_thresh$i.vcf -m $i
	python3 bin/VerifySVpred.py 2 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/TSVcall/LUMPY/lumpy_thresh$i.vcf $ProjectDir/TSVcall/LUMPY/hit_thresh$i.txt
done

# preparing for gmap
mkdir -p $ProjectDir/WholeGenome/Chromosomes/
if [ ! -e $ProjectDir/WholeGenome/Chromosomes/chr1.fa ]; then
	python3 bin/ReadGenome.py $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/WholeGenome/Chromosomes/
fi
if [ ! -d $ProjectDir/WholeGenome/Chromosomes/"genome"${ProjectDir##*/} ]; then
	gmap_build -D $ProjectDir/WholeGenome/Chromosomes/ -d "genome"${ProjectDir##*/} $ProjectDir/WholeGenome/Chromosomes/*.fa
fi

# TransABySS + mummer3
TransABySS -k 25 --pe $ProjectDir/Reads/RNA1.fq.gz $ProjectDir/Reads/RNA2.fq.gz --outdir $ProjectDir/TSVcall/TransABySS/k25 --name k25 --threads 4 --island 0
TransABySS -k 32 --pe $ProjectDir/Reads/RNA1.fq.gz $ProjectDir/Reads/RNA2.fq.gz --outdir $ProjectDir/TSVcall/TransABySS/k32 --name k32 --threads 4 --island 0
TransABySS-merge --mink 25 --maxk 32 --prefixes k25. k32. --out $ProjectDir/TSVcall/TransABySS/mergedassembly.fa $ProjectDir/TSVcall/TransABySS/k25/k25-final.fa $ProjectDir/TSVcall/TransABySS/k32/k32-final.fa
nucmer -p $ProjectDir/TSVcall/TransABySS/nucm $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/TSVcall/TransABySS/mergedassembly.fa
show-coords $ProjectDir/TSVcall/TransABySS/nucm.delta > $ProjectDir/TSVcall/TransABySS/nucm_text.txt
./bin/NucmerSV2 $ProjectDir/TSVcall/TransABySS/nucm_text.txt $ProjectDir/TSVcall/TransABySS/nucm_res.bedpe
python3 bin/VerifySVpred.py 5 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/TSVcall/TransABySS/nucm_res.bedpe $ProjectDir/TSVcall/TransABySS/nucm_res_hit.bedpe

# TransABySS + gmap
gmap -D $ProjectDir/WholeGenome/Chromosomes/ -d "genome"${ProjectDir##*/} $ProjectDir/TSVcall/TransABySS/mergedassembly.fa > $ProjectDir/TSVcall/TransABySS/gmap.out
./bin/GmapSV $ProjectDir/TSVcall/TransABySS/gmap.out $ProjectDir/TSVcall/TransABySS/gmap_res.bedpe
python3 bin/VerifySVpred.py 5 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/TSVcall/TransABySS/gmap_res.bedpe $ProjectDir/TSVcall/TransABySS/gmap_res_hit.bedpe

# trinity + mummer3
Trinity --seqType fq --max_memory 50G --left $ProjectDir/Reads/RNA1.fq.gz --right $ProjectDir/Reads/RNA2.fq.gz --CPU 6 --output $ProjectDir/TSVcall/Trinity
nucmer -p $ProjectDir/TSVcall/Trinity/nucm $ProjectDir/WholeGenome/genome_rearranged.fa $ProjectDir/TSVcall/Trinity/Trinity.fasta
show-coords $ProjectDir/TSVcall/Trinity/nucm.delta > $ProjectDir/TSVcall/Trinity/nucm_text.txt
./bin/NucmerSV2 $ProjectDir/TSVcall/Trinity/nucm_text.txt $ProjectDir/TSVcall/Trinity/nucm_res.bedpe
python3 bin/VerifySVpred.py 5 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/TSVcall/Trinity/nucm_res.bedpe $ProjectDir/TSVcall/Trinity/nucm_res_hit.bedpe

# trinity + gmap
gmap -D $ProjectDir/WholeGenome/Chromosomes/ -d "genome"${ProjectDir##*/} $ProjectDir/TSVcall/Trinity/Trinity.fasta > $ProjectDir/TSVcall/Trinity/gmap.out
./bin/GmapSV $ProjectDir/TSVcall/Trinity/gmap.out $ProjectDir/TSVcall/Trinity/gmap_res.bedpe
python3 bin/VerifySVpred.py 5 $ProjectDir/WholeGenome/SV_newpos.txt $ProjectDir/TSVcall/Trinity/gmap_res.bedpe $ProjectDir/TSVcall/Trinity/gmap_res_hit.bedpe