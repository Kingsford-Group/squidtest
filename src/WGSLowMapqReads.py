#!/bin/python

import sys
import pysam
import subprocess

Nucleotide={'A':'T','C':'G','G':'C','T':'A','R':'Y','Y':'R','S':'W','W':'S','K':'M','M':'K','B':'V','V':'B','D':'H','H':'D', 'N':'N', '.':'.','-':'-'}
def ReverseComp(strs):
	global Nucleotide
	newstrs=""
	for e in strs:
		newstrs+=Nucleotide[e.upper()]
	newstrs=newstrs[::-1]
	return newstrs

def CollectReadName(bamfile, NMthresh):
	ReadNames=[]
	fp=pysam.AlignmentFile(bamfile)
	for read in fp.fetch():
		if read.is_secondary:
			continue
		if read.is_unmapped or read.mate_is_unmapped:
			ReadNames.append(read.query_name)
		elif read.get_tag("NM")>NMthresh:
			ReadNames.append(read.query_name)
	fp.close()
	ReadNames=set(ReadNames)
	ReadNameMap={}
	for e in ReadNames:
		ReadNameMap[e]=[0]
	return ReadNameMap

def WriteBam4Reads(ReadNameMap, inbam, outprefix):
	infp=pysam.AlignmentFile(inbam)
	for read in infp.fetch():
		if len(read.query_sequence)<70:
			continue
		if read.query_name in ReadNameMap:
			v=ReadNameMap[read.query_name][0]
			if v==0:
				if read.is_read1:
					ReadNameMap[read.query_name][0]+=1
					ReadNameMap[read.query_name].append(read)
				else:
					ReadNameMap[read.query_name][0]+=2
					ReadNameMap[read.query_name].append(read)
			elif v==1 and read.is_read2:
				ReadNameMap[read.query_name][0]+=2
				ReadNameMap[read.query_name].append(read)
			elif v==2 and read.is_read1:
				ReadNameMap[read.query_name][0]+=1
				ReadNameMap[read.query_name].append(read)
	infp.close()

	outr1=open(outprefix+"_1.fastq", 'w')
	outr2=open(outprefix+"_2.fastq", "w")
	for key, value in ReadNameMap.iteritems():
		if len(value)==3:
			for i in range(1,3):
				if value[i].is_read1:
					if value[i].is_reverse:
						outr1.write("@{}\n{}\n+\n{}\n".format(read.query_name, ReverseComp(read.query_sequence), "".join([chr(x+33) for x in read.query_qualities])))
					else:
						outr1.write("@{}\n{}\n+\n{}\n".format(read.query_name, read.query_sequence, "".join([chr(x+33) for x in read.query_qualities])))
				else:
					if value[i].is_reverse:
						outr2.write("@{}\n{}\n+\n{}\n".format(read.query_name, ReverseComp(read.query_sequence), "".join([chr(x+33) for x in read.query_qualities])))
					else:
						outr2.write("@{}\n{}\n+\n{}\n".format(read.query_name, read.query_sequence, "".join([chr(x+33) for x in read.query_qualities])))
	outr1.close()
	outr2.close()

def WriteDiscordant(inbam, outprefix):
	global MinAlignLen
	global thresh
	count=0
	Reads={}
	tmpReads={}
	Added=set()
	infp=pysam.AlignmentFile(inbam)
	outfp=pysam.AlignmentFile(outprefix+"_original.bam", "wb", template=infp)
	outr1=open(outprefix+"_1.fastq", 'w')
	outr2=open(outprefix+"_2.fastq", "w")

	for read in infp:
		count+=1
		if count%10000000==0:
			print(count)
		if read.is_secondary or read.is_supplementary or read.is_proper_pair:
			continue
		readdis=False
		if read.is_unmapped or read.mate_is_unmapped:
			readdis=True
		elif read.reference_id!=read.next_reference_id or read.is_reverse==read.mate_is_reverse or abs(read.reference_start-read.next_reference_start)>thresh:
			readdis=True
		elif read.is_reverse and read.reference_start<read.next_reference_start-10:
			readdis=True
		elif read.mate_is_reverse and read.next_reference_start<read.reference_start-10:
			readdis=True
		if readdis and read.query_name not in Added:
			if read.query_name in Reads and Reads[read.query_name][0].is_read1!=read.is_read1:
				# Reads[read.query_name].append(read)
				thisReads=Reads[read.query_name]
				thisReads.append(read)
				# for r in Reads[read.query_name]:
				for r in thisReads:
					outfp.write(r)
					if r.is_read1:
						if r.is_reverse:
							outr1.write("@{}\n{}\n+\n{}\n".format(r.query_name, ReverseComp(r.query_sequence), "".join([chr(x+33) for x in r.query_qualities])))
						else:
							outr1.write("@{}\n{}\n+\n{}\n".format(r.query_name, r.query_sequence, "".join([chr(x+33) for x in r.query_qualities])))
					else:
						if r.is_reverse:
							outr2.write("@{}\n{}\n+\n{}\n".format(r.query_name, ReverseComp(r.query_sequence), "".join([chr(x+33) for x in r.query_qualities])))
						else:
							outr2.write("@{}\n{}\n+\n{}\n".format(r.query_name, r.query_sequence, "".join([chr(x+33) for x in r.query_qualities])))
				Added=Added|set(read.query_name)
				del Reads[read.query_name]
			else:
				Reads[read.query_name]=[read]
		
	infp.close()
	outfp.close()
	outr1.close()
	outr2.close()

def CommandDiscord(inbam, outprefix):
	global thresh
	count=0
	infp=pysam.AlignmentFile(inbam)
	outfp=pysam.AlignmentFile(outprefix+"_tmp.bam", "wb", template=infp)
	# Write raw discordant reads 
	for read in infp:
		count+=1
		if count%10000000==0:
			print(count)
		if read.is_secondary or read.is_supplementary or read.is_proper_pair:
			continue
		readdis=False
		if read.is_unmapped or read.mate_is_unmapped:
			readdis=True
		elif read.reference_id!=read.next_reference_id or read.is_reverse==read.mate_is_reverse or abs(read.reference_start-read.next_reference_start)>thresh:
			readdis=True
		elif read.is_reverse and read.reference_start<read.next_reference_start-10:
			readdis=True
		elif read.mate_is_reverse and read.next_reference_start<read.reference_start-10:
			readdis=True
		if readdis:
			outfp.write(read)
	infp.close()
	outfp.close()
	# Sort raw discordant reads by read name
	command="samtools sort -n "+outprefix+"_tmp.bam -o "+outprefix+"_tmp_sort.bam"
	p=subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
	(out, err)=p.communicate()
	if err:
		print(err)
		return
	# Select reads with both ends only once and Generate fastq files
	infp=pysam.AlignmentFile(outprefix+"_tmp_sort.bam")
	outfp=pysam.AlignmentFile(outprefix+"_original.bam", "wb", template=infp)
	outr1=open(outprefix+"_1.fastq", 'w')
	outr2=open(outprefix+"_2.fastq", "w")
	Reads=[]
	for read in infp:
		if len(Reads)==0 or Reads[0].query_name==read.query_name:
			Reads.append(read)
		else:
			hasread1=False
			hasread2=False
			for r in Reads:
				if r.is_read1:
					read1=r
					hasread1=True
				else:
					read2=r
					hasread2=True
			if hasread1 and hasread2:
				outfp.write(read1)
				outfp.write(read2)
				for r in [read1, read2]:
					if r.is_read1:
						if r.is_reverse:
							outr1.write("@{}\n{}\n+\n{}\n".format(r.query_name, ReverseComp(r.query_sequence), "".join([chr(x+33) for x in r.query_qualities])))
						else:
							outr1.write("@{}\n{}\n+\n{}\n".format(r.query_name, r.query_sequence, "".join([chr(x+33) for x in r.query_qualities])))
					else:
						if r.is_reverse:
							outr2.write("@{}\n{}\n+\n{}\n".format(r.query_name, ReverseComp(r.query_sequence), "".join([chr(x+33) for x in r.query_qualities])))
						else:
							outr2.write("@{}\n{}\n+\n{}\n".format(r.query_name, r.query_sequence, "".join([chr(x+33) for x in r.query_qualities])))
			Reads=[read]
	infp.close()
	outfp.close()
	outr1.close()
	outr2.close()
	# remove tmp files
	command="rm -fr "+outprefix+"_tmp.bam "+outprefix+"_tmp_sort.bam"
	p=subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
	(out, err)=p.communicate()
	if err:
		print(err)
		return

def Part2(inbam, outprefix):
	# Sort raw discordant reads by read name
	command="samtools sort -n "+outprefix+"_tmp.bam -o "+outprefix+"_tmp_sort.bam"
	p=subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
	(out, err)=p.communicate()
	if err:
		print(err)
		return

def Part3(inbam, outprefix):
	# Select reads with both ends only once and Generate fastq files
	infp=pysam.AlignmentFile(outprefix+"_tmp_sort.bam")
	outfp=pysam.AlignmentFile(outprefix+"_original.bam", "wb", template=infp)
	outr1=open(outprefix+"_1.fastq", 'w')
	outr2=open(outprefix+"_2.fastq", "w")
	Reads=[]
	for read in infp:
		if len(Reads)==0 or Reads[0].query_name==read.query_name:
			Reads.append(read)
		else:
			hasread1=False
			hasread2=False
			for r in Reads:
				if r.is_read1:
					read1=r
					hasread1=True
				else:
					read2=r
					hasread2=True
			if hasread1 and hasread2:
				outfp.write(read1)
				outfp.write(read2)
				for r in [read1, read2]:
					if r.is_read1:
						if r.is_reverse:
							outr1.write("@{}\n{}\n+\n{}\n".format(r.query_name, ReverseComp(r.query_sequence), "".join([chr(x+33) for x in r.query_qualities])))
						else:
							outr1.write("@{}\n{}\n+\n{}\n".format(r.query_name, r.query_sequence, "".join([chr(x+33) for x in r.query_qualities])))
					else:
						if r.is_reverse:
							outr2.write("@{}\n{}\n+\n{}\n".format(r.query_name, ReverseComp(r.query_sequence), "".join([chr(x+33) for x in r.query_qualities])))
						else:
							outr2.write("@{}\n{}\n+\n{}\n".format(r.query_name, r.query_sequence, "".join([chr(x+33) for x in r.query_qualities])))
			Reads=[read]
	infp.close()
	outfp.close()
	outr1.close()
	outr2.close()
	# remove tmp files
	command="rm -fr "+outprefix+"_tmp.bam "+outprefix+"_tmp_sort.bam"
	p=subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
	(out, err)=p.communicate()
	if err:
		print(err)
		return

if __name__=="__main__":
	if len(sys.argv)<3:
		print("python WGSLowMapqReads.py <InBam> <OutPrefix>")
	else:
		MinAlignLen=40
		thresh=5000
		InBam=sys.argv[1]
		OutPrefix=sys.argv[2]
		# WriteDiscordant(InBam, OutPrefix)
		CommandDiscord(InBam, OutPrefix)
		# Part3(InBam, OutPrefix)
