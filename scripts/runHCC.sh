#!/bin/bash

argv=("$@")

FusionCatcherExe=""
FusionCatcherDataDir=""
ChimerascanDir="" # chimerascan installation directory
JaffaDir=""
DefuseDir=""
PicardDir=""
OutputDir="" # downloaded reference and RNA-seq folder, the same as output folder
DataHCC=""
threads=2
ExeDir=$(pwd)

for ((i=0; i<${#argv}; i+=2)); do
	if [[ ${argv[$i]} == "--fusioncatcher" ]]; then
		FusionCatcherExe=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--fusioncatcherdata" ]]; then
		FusionCatcherDataDir=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--chimeriscan" ]]; then
		ChimerascanDir=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--jaffa" ]]; then
		JaffaDir=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--defuse" ]]; then
		DefuseDir=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--picard" ]]; then
		PicardDir=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--datadownload" ]]; then
		OutputDir=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--dataHCC" ]]; then
		DataHCC=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--threads" ]]; then
		threads=${argv[(($i+1))]}
	fi
done

echo "INTEGRATE executable are in your path."

GenomeFasta=${OutputDir}/Ensemble75/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa
AnnotationGTF=${OutputDir}/Ensemble75/Homo_sapiens.GRCh37.75.gtf

# preparing STAR index
echo "Generating STAR index"
mkdir -p ${OutputDir}/Ensemble75/STARIndex
STAR --runThreadN ${threads} --runMode genomeGenerate --genomeDir ${OutputDir}/Ensemble75/STARIndex --genomeFastaFiles ${GenomeFasta} --sjdbGTFfile ${AnnotationGTF}
# preparing chimerascan Bowtie index
echo "Generating chimerascan index"
mkdir -p ${OutputDir}/Ensemble75/CSBowtieIndex
gunzip ${DataHCC}/Homo_sapiens.GRCh37.75_chimerascan.txt.gz
python ${ChimerascanDir}/bin/chimerascan_index.py ${GenomeFasta} ${DataHCC}/Homo_sapiens.GRCh37.75_chimerascan.txt ${OutputDir}/Ensemble75/CSBowtieIndex
# preparing INTEGRATE index
echo "Generating INTEGRATE index"
gunzip ${DataHCC}/annot.ensembl.txt.gz
mkdir -p ${OutputDir}/Ensemble75/IntegrateIndex
${IntegrateIntallDir}/bin/Integrate mkbwt -dir ${OutputDir}/Ensemble75/IntegrateIndex ${GenomeFasta}

cd $OutputDir
# align reads using STAR
echo "STAR alignment"
mkdir -p ${OutputDir}/StarAlign/HCC1954
STAR  --runThreadN ${threads} --genomeDir ${OutputDir}/Ensemble75/STARIndex --readFilesIn ${OutputDir}/RNAseq/RNAHCC1954_1.fastq ${OutputDir}/RNAseq/RNAHCC1954_2.fastq --outFileNamePrefix $OutputDir/StarAlign/HCC1954/ --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --chimSegmentMin 15 --outReadsUnmapped Fastx
samtools view -Shb ${OutputDir}/StarAlign/HCC1954/Chimeric.out.sam -o ${OutputDir}/StarAlign/HCC1954/Chimeric.out.bam

mkdir -p ${OutputDir}/StarAlign/HCC1395
STAR  --runThreadN ${threads} --genomeDir ${OutputDir}/Ensemble75/STARIndex --readFilesIn ${OutputDir}/RNAseq/SRR2532336_1.fastq ${OutputDir}/RNAseq/SRR2532336_2.fastq --outFileNamePrefix $OutputDir/StarAlign/HCC1395/ --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --chimSegmentMin 15 --outReadsUnmapped Fastx
samtools view -Shb ${OutputDir}/StarAlign/HCC1395/Chimeric.out.sam -o ${OutputDir}/StarAlign/HCC1395/Chimeric.out.bam

# SQUID
echo "predicting TSV with SQUID"
mkdir -p ${OutputDir}/TSVprediction/squid_1954
squid -b ${OutputDir}/StarAlign/HCC1954/Aligned.sortedByCoord.out.bam -c ${OutputDir}/StarAlign/HCC1954/Chimeric.out.bam -o ${OutputDir}/TSVprediction/squid_1954/squidtsv
awk 'BEGIN{FS="\t";} {if(substr($0,0,1)=="#" || !($1~"M" || $1~"G" || $4~"M" || $4~"G")) print $0;}' ${OutputDir}/TSVprediction/squid_1954/squidtsv_sv.txt > ${OutputDir}/TSVprediction/squid_1954/squidtsv_sv_final.txt

mkdir -p ${OutputDir}/TSVprediction/squid_1395
squid -b ${OutputDir}/StarAlign/HCC1395/Aligned.sortedByCoord.out.bam -c ${OutputDir}/StarAlign/HCC1395/Chimeric.out.bam -o ${OutputDir}/TSVprediction/squid_1395/squidtsv
awk 'BEGIN{FS="\t";} {if(substr($0,0,1)=="#" || !($1~"M" || $1~"G" || $4~"M" || $4~"G")) print $0;}' ${OutputDir}/TSVprediction/squid_1395/squidtsv_sv.txt > ${OutputDir}/TSVprediction/squid_1395/squidtsv_sv_final.txt

# fusioncatcher
mkdir -p ${OutputDir}/TSVprediction/fusioncatcher_1954
mkdir -p ${OutputDir}/TSVprediction/fusioncatcher_1954/data
ln -s ${OutputDir}/RNAseq/RNAHCC1954_1.fastq ${OutputDir}/TSVprediction/fusioncatcher_1954/data/
ln -s ${OutputDir}/RNAseq/RNAHCC1954_2.fastq ${OutputDir}/TSVprediction/fusioncatcher_1954/data/
${FusionCatcherExe} -d ${FusionCatcherDataDir} -i ${OutputDir}/TSVprediction/fusioncatcher_1954/data/ -o ${OutputDir}/TSVprediction/fusioncatcher_1954/

mkdir -p ${OutputDir}/TSVprediction/fusioncatcher_1395
mkdir -p ${OutputDir}/TSVprediction/fusioncatcher_1395/data
ln -s ${OutputDir}/RNAseq/SRR2532336_1.fastq ${OutputDir}/TSVprediction/fusioncatcher_1395/data/
ln -s ${OutputDir}/RNAseq/SRR2532336_2.fastq ${OutputDir}/TSVprediction/fusioncatcher_1395/data/
${FusionCatcherExe} -d ${FusionCatcherDataDir} -i ${OutputDir}/TSVprediction/fusioncatcher_1395/data/ -o ${OutputDir}/TSVprediction/fusioncatcher_1395/

# JAFFA
mkdir -p ${OutputDir}/TSVprediction/jaffa_1954
mkdir -p ${OutputDir}/TSVprediction/jaffa_1954/data
ln -s ${OutputDir}/RNAseq/RNAHCC1954_1.fastq.gz ${OutputDir}/TSVprediction/jaffa_1954/data
ln -s ${OutputDir}/RNAseq/RNAHCC1954_2.fastq.gz ${OutputDir}/TSVprediction/jaffa_1954/data
cd ${OutputDir}/TSVprediction/jaffa_1954
${JaffaDir}/tools/bin/bpipe run ${JaffaDir}/JAFFA_assembly.groovy ./data/RNAHCC1954_*.fastq.gz
cd ${ExeDir}

mkdir -p ${OutputDir}/TSVprediction/jaffa_1395
mkdir -p ${OutputDir}/TSVprediction/jaffa_1395/data
ln -s ${OutputDir}/RNAseq/SRR2532336_1.fastq.gz ${OutputDir}/TSVprediction/jaffa_1395/data/
ln -s ${OutputDir}/RNAseq/SRR2532336_2.fastq.gz ${OutputDir}/TSVprediction/jaffa_1395/data/
cd ${OutputDir}/TSVprediction/jaffa_1395
${JaffaDir}/tools/bin/bpipe run ${JaffaDir}/JAFFA_assembly.groovy ./data/SRR2532336_*.fastq.gz
cd ${ExeDir}

# deFuse
mkdir -p ${OutputDir}/TSVprediction/defuse_1954
${DefuseDir}/scripts/defuse_run.pl -c ${DefuseDir}/scripts/config.txt -d ${DefuseDir}/RefData/ -1 ${OutputDir}/RNAseq/RNAHCC1954_1.fastq -2 ${OutputDir}/RNAseq/RNAHCC1954_2.fastq -o ${OutputDir}/TSVprediction/defuse_1954 -p ${threads}

mkdir -p ${OutputDir}/TSVprediction/defuse_1395
${DefuseDir}/scripts/defuse_run.pl -c ${DefuseDir}/scripts/config.txt -d ${DefuseDir}/RefData/ -1 ${OutputDir}/RNAseq/SRR2532336_1.fastq -2 ${OutputDir}/RNAseq/SRR2532336_2.fastq -o ${OutputDir}/TSVprediction/defuse_1954 -p ${threads}

# Chimerascan
mkdir -p ${OutputDir}/TSVprediction/chimerascan_1954
python ${ChimerascanDir}/bin/chimerascan_run.py -v -p ${threads} ${OutputDir}/Ensemble75/CSBowtieIndex ${OutputDir}/RNAseq/RNAHCC1954_1.fastq ${OutputDir}/RNAseq/RNAHCC1954_2.fastq ${OutputDir}/TSVprediction/chimerascan_1954/

mkdir -p ${OutputDir}/TSVprediction/chimerascan_1395
python ${ChimerascanDir}/bin/chimerascan_run.py -v -p ${threads} ${OutputDir}/Ensemble75/CSBowtieIndex ${OutputDir}/RNAseq/SRR2532336_1.fastq ${OutputDir}/RNAseq/SRR2532336_2.fastq ${OutputDir}/TSVprediction/chimerascan_1395/

# INTEGRATE
samtools merge $OutputDir/StarAlign/HCC1954/MergedAlign.bam $OutputDir/StarAlign/HCC1954/Aligned.sortedByCoord.out.bam $OutputDir/StarAlign/HCC1954/Chimeric.out.bam
samtools sort $OutputDir/StarAlign/HCC1954/MergedAlign.bam -o $OutputDir/StarAlign/HCC1954/MergedAlign_sort.bam
samtools index $OutputDir/StarAlign/HCC1954/MergedAlign_sort.bam
java -jar ${PicardDir}/picard.jar FastqToSam F1=$OutputDir/StarAlign/HCC1954/Unmapped.out.mate1 F2=$OutputDir/StarAlign/HCC1954/Unmapped.out.mate2 O=$OutputDir/StarAlign/HCC1954/Unmapped.bam SM=HCC1954 
mkdir -p ${OutputDir}/TSVprediction/integrate_1954
Integrate fusion ${GenomeFasta} ${DataHCC}/annot.ensembl.txt ${OutputDir}/Ensemble75/IntegrateIndex $OutputDir/StarAlign/HCC1954//MergedAlign_sort.bam $OutputDir/StarAlign/HCC1954/Unmapped.bam

samtools merge $OutputDir/StarAlign/HCC1395/MergedAlign.bam $OutputDir/StarAlign/HCC1395/Aligned.sortedByCoord.out.bam $OutputDir/StarAlign/HCC1395/Chimeric.out.bam
samtools sort $OutputDir/StarAlign/HCC1395/MergedAlign.bam -o $OutputDir/StarAlign/HCC1395/MergedAlign_sort.bam
samtools index $OutputDir/StarAlign/HCC1395/MergedAlign_sort.bam
java -jar ${PicardDir}/picard.jar FastqToSam F1=$OutputDir/StarAlign/HCC1395/Unmapped.out.mate1 F2=$OutputDir/StarAlign/HCC1395/Unmapped.out.mate2 O=$OutputDir/StarAlign/HCC1395/Unmapped.bam SM=HCC1395
mkdir -p ${OutputDir}/TSVprediction/integrate_1395
Integrate fusion ${GenomeFasta} ${DataHCC}/annot.ensembl.txt ${OutputDir}/Ensemble75/IntegrateIndex $OutputDir/StarAlign/HCC1395//MergedAlign_sort.bam $OutputDir/StarAlign/HCC1395/Unmapped.bam
