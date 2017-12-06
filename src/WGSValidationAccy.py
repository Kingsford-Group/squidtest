#!/bin/python

import pysam
import sys

def GetNonreadthruChr(squidfasta):
	global thresh
	fp=open(squidfasta, 'r')
	ChrNames=[]
	for line in fp:
		if line[0]=='>':
			strs=line.strip().split(" ")
			segs1=strs[2].split("|")[0].split(",")
			segs2=strs[2].split("|")[1].split(",")
			chr1=segs1[0]
			chr2=segs2[0]
			if segs1[3]=="False":
				bp1=int(segs1[2])
			else:
				bp1=int(segs1[1])
			if segs2[3]=="False":
				bp2=int(segs2[2])
			else:
				bp2=int(segs2[1])
			if chr1!=chr2:
				ChrNames.append(strs[0][1:])
			elif abs(bp1-bp2)>thresh:
				ChrNames.append(strs[0][1:])
			elif segs1[3]!="False" or segs2[3]!="True":
				ChrNames.append(strs[0][1:])
	fp.close()
	ValidateChrNames=[]
	for c in set(ChrNames):
		if ChrNames.count(c)>2:
			ValidateChrNames.append(c)
	return set(ChrNames)

def GetNonImmunoChr(squidfasta):
	fp=open(squidfasta, 'r')
	IgH=['14', 105300000, 106900000]
	IgL=['22',  22230000, 22930000]
	IgK=['2', 88900000, 90240000]
	ChrNames=[]
	for line in fp:
		if line[0]=='>':
			strs=line.strip().split(" ")
			segs1=strs[2].split("|")[0].split(",")
			segs2=strs[2].split("|")[1].split(",")
			chr1=segs1[0]
			chr2=segs2[0]
			if segs1[3]=="False":
				bp1=int(segs1[2])
			else:
				bp1=int(segs1[1])
			if segs2[3]=="False":
				bp2=int(segs2[2])
			else:
				bp2=int(segs2[1])
			if chr1==IgH[0] and chr2==IgH[0] and bp1>=IgH[1] and bp1<=IgH[2] and bp2>=IgH[1] and bp2<IgH[2] and segs1[3]=="False" and segs2[3]=="True":
				continue
			elif chr1==IgL[0] and chr2==IgL[0] and bp1>=IgL[1] and bp1<=IgL[2] and bp2>=IgL[1] and bp2<IgL[2] and segs1[3]=="False" and segs2[3]=="True":
				continue
			elif chr1==IgK[0] and chr2==IgK[0] and bp1>=IgK[1] and bp1<=IgK[2] and bp2>=Igk[1] and bp2<Igk[2] and segs1[3]=="False" and segs2[3]=="True":
				continue
			else:
				ChrNames.append(strs[0][1:])
	fp.close()
	return set(ChrNames)

def GetHitChr(bamfile):
	fp=pysam.AlignmentFile(bamfile)
	ChrNames=[]
	for read in fp:
		ChrNames.append(read.reference_name)
	fp.close()
	return set(ChrNames)

if __name__=="__main__":
	if len(sys.argv)<2:
		print("python3 WGSValidationAccy.py <squidfasta> <Better.bam>")
	else:
		thresh=5000

		SquidFasta=sys.argv[1]
		Bamfile=sys.argv[2]
		
		strs=Bamfile.split("/")

		# ChrNamesFasta=GetNonreadthruChr(SquidFasta)
		ChrNamesFasta=GetNonImmunoChr(SquidFasta)
		ChrNamesBam=GetHitChr(Bamfile)
		print("{}\t{}\t{}\t{}".format(len(ChrNamesBam&ChrNamesFasta), len(ChrNamesFasta), strs[-6], strs[-4]))
