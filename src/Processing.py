#!/bin/python3

import sys
import re

def isint(x):
	try:
		int(x)
		return True
	except ValueError:
		return False

def ProcessGalante(filename):
	fp=open(filename, 'r')
	SVcode=""
	Bedpe=[]
	for line in fp:
		if "interchromosomal" in line.lower():
			SVcode="TRA"
			continue
		elif "deletion" in line.lower():
			SVcode=""	
			continue
		elif "duplication" in line.lower():
			SVcode="DUP"
			continue
		elif "inversion" in line.lower():
			SVcode="INV"
			continue
		if SVcode=="" or len(line)<5:
			continue
		else:
			strs=re.split(r'[: \t]', line.strip())
			BP1=["", -1, -1, "."]
			BP2=["", -1, -1, "."]
			for e in strs:
				if "chr" in e.lower():
					if BP1[0]=="":
						BP1[0]=e.lower()
					else:
						BP2[0]=e.lower()
				elif "-" in e and len(e)>5:
					numbers=[int(e.split("-")[0]), int(e.split("-")[1])]
					if BP1[1]==-1:
						BP1[1:3]=numbers
					else:
						BP2[1:3]=numbers
				elif "(+)"==e or "(-)"==e:
					if BP1[3]==".":
						BP1[3]=e[1]
					else:
						BP2[3]=e[1]
			if BP1[0]=="chr24" or BP2[0]=="chr24":
				continue
			if SVcode=="TRA" and BP1[1]!=-1 and BP1[2]!=-1 and BP1[3]!="." and BP2[1]!=-1 and BP2[2]!=-1 and BP2[3]!=".":
				Bedpe.append([BP1[0], BP1[1], BP1[2], BP2[0], BP2[1], BP2[2], BP1[3], BP2[3], "TRA"])
			elif SVcode=="INV":
				Bedpe.append([BP1[0], BP1[1], BP1[1]+1, BP1[0], BP1[2], BP1[2]+1, "-", "-", "INV"])
			elif SVcode=="DUP":
				Bedpe.append([BP1[0], BP1[1], BP1[1]+1, BP1[0], BP1[2]-1, BP1[2], "-", "+", "DUP"])
	fp.close()
	return Bedpe

def ProcessBignell(filename): # XXX SV type need to be added
	fp=open(filename, 'r')
	Bedpe=[]
	for line in fp:
		strs=line.strip().split()
		if len(strs)>=8 and strs[1]=="HCC1954":
			index=[i for i in range(len(strs)) if strs[i]=="+" or strs[i]=="-"]
			BP1=[]
			if strs[index[0]]=="+" and isint(strs[index[0]+2]):
				BP1=["chr"+strs[index[0]-1], int(strs[index[0]+3]), int(strs[index[0]+3])+1, "-"]
			elif strs[index[0]]=="-" and isint(strs[index[0]+2]):
				BP1=["chr"+strs[index[0]-1], int(strs[index[0]+3])-1, int(strs[index[0]+3]), "+"]
			elif strs[index[0]]=="+" and not isint(strs[index[0]+2]):
				BP1=["chr"+strs[index[0]-1], int(strs[index[0]+4]), int(strs[index[0]+4])+1, "-"]
			elif strs[index[0]]=="-" and not isint(strs[index[0]+2]):
				BP1=["chr"+strs[index[0]-1], int(strs[index[0]+4])-1, int(strs[index[0]+4]), "+"]
			BP2=[]
			if strs[index[1]]=="+" and isint(strs[index[1]-2]):
				BP2=["chr"+strs[index[1]+1], int(strs[index[1]-3])-1, int(strs[index[1]-3]), "+"]
			elif strs[index[1]]=="-" and isint(strs[index[1]-2]):
				BP2=["chr"+strs[index[1]+1], int(strs[index[1]-3]), int(strs[index[1]-3])+1, "-"]
			elif strs[index[1]]=="+" and not isint(strs[index[1]-2]):
				BP2=["chr"+strs[index[1]+1], int(strs[index[1]-4])-1, int(strs[index[1]-4]), "+"]
			elif strs[index[1]]=="-" and not isint(strs[index[1]-2]):
				BP2=["chr"+strs[index[1]+1], int(strs[index[1]-4]), int(strs[index[1]-4])+1, "-"]
			if BP1[0]=="chr24" or BP2[0]=="chr24":
				continue
			if int(BP1[0][3:])>int(BP2[0][3:]) or (int(BP1[0][3:])==int(BP2[0][3:]) and int(BP1[1])>int(BP2[1])):
				Bedpe.append([BP2[0], BP2[1], BP2[2], BP1[0], BP1[1], BP1[2], BP2[3], BP1[3]])
			else:
				Bedpe.append([BP1[0], BP1[1], BP1[2], BP2[0], BP2[1], BP2[2], BP1[3], BP2[3]])
			if BP1[0]==BP2[0]:
				if BP1[3]==BP2[3]:
					Bedpe[-1].append("INV")
				elif Bedpe[-1][6]=="-" and Bedpe[-1][7]=="+":
					Bedpe[-1].append("DUP")
				elif abs(Bedpe[-1][1]-Bedpe[-1][4])<100000:
					Bedpe.pop()
				else:
					Bedpe[-1].append("INS")
			else:
				Bedpe[-1].append("TRA")
	fp.close()
	return Bedpe

def ProcessStephens(filename):
	fp=open(filename, 'r')
	Bedpe=[]
	for line in fp:
		strs=line.strip().split("\t")
		if len(strs)>1 and strs[0]=="HCC1954" and strs[2].lower()!="deletion":
			if strs[5]=="+":
				BP1=["chr"+strs[3], int(strs[4])-1, int(strs[4]), strs[5]]
			else:
				BP1=["chr"+strs[3], int(strs[4]), int(strs[4])+1, strs[5]]
			if strs[8]=="+":
				BP2=["chr"+strs[6], int(strs[7])-1, int(strs[7]), strs[8]]
			else:
				BP2=["chr"+strs[6], int(strs[7]), int(strs[7])+1, strs[8]]
			if BP1[0]==BP2[0] and BP1[3]=="+" and BP2[3]=="-" and abs(BP1[1]-BP2[1])<100000:
				continue
			if BP1[0]=="chr24" or BP2[0]=="chr24":
				continue
			if BP1[0]==BP2[0]:
				if BP1[3]==BP2[3]:
					Bedpe.append([BP1[0], BP1[1], BP1[2], BP2[0], BP2[1], BP2[2], BP1[3], BP2[3], "INV"])
				elif BP1[3]=="-" and BP2[3]=="+":
					Bedpe.append([BP1[0], BP1[1], BP1[2], BP2[0], BP2[1], BP2[2], BP1[3], BP2[3], "DUP"])
				else:
					Bedpe.append([BP1[0], BP1[1], BP1[2], BP2[0], BP2[1], BP2[2], BP1[3], BP2[3], "INS"]) # XXX do we have INS?
			else:
				Bedpe.append([BP1[0], BP1[1], BP1[2], BP2[0], BP2[1], BP2[2], BP1[3], BP2[3], "TRA"])
	fp.close()
	return Bedpe

def WritePlain(filename, Bedpe):
	fp=open(filename, 'w')
	for b in Bedpe:
		fp.write(b[0]+":"+str(b[1])+"-"+str(b[2])+"\n")
		fp.write(b[3]+":"+str(b[4])+"-"+str(b[5])+"\n")
	fp.close()

def SubstituteBedpe(filename, Bedpe):
	fp=open(filename, 'r')
	flag=True
	count=0
	for line in fp:
		newchr=line[:line.index(":")]
		startpos=line[line.index(":")+1:line.index("-")]
		endpos=line[line.index("-")+1:-1]
		if flag:
			Bedpe[count][0]=newchr
			Bedpe[count][1]=int(startpos)
			Bedpe[count][2]=int(endpos)
		else:
			Bedpe[count][3]=newchr
			Bedpe[count][4]=int(startpos)
			Bedpe[count][5]=int(endpos)
			count+=1
		flag=not flag
	fp.close()
	return Bedpe

def WriteBedpe(filename, Bedpe):
	fp=open(filename, 'w')
	for b in Bedpe:
		fp.write(b[0]+"\t"+str(b[1])+"\t"+str(b[2])+"\t"+b[3]+"\t"+str(b[4])+"\t"+str(b[5])+"\t.\t.\t"+b[6]+"\t"+b[7]+"\n")
	fp.close()

Bedpe=ProcessGalante("Galante_text.txt")
#x=ProcessBignell("Bignell_text.txt")
#for e in x:
#	Bedpe.append(e)
#x=ProcessStephens("Stephens.csv")
#for e in x:
#	Bedpe.append(e)
#WritePlain("SVbp.txt", Bedpe)
