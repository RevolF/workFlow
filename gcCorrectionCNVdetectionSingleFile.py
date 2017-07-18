# -*- coding: utf-8 -*-
"""
Created on Tue May 02 16:23:09 2017

@author: thor
"""

# single file processing mode of GC correction and CNV detection

####################### log #######################
# 2017.05.02:
#	mappable positions are indicated in file: DMD100_target.mappable.sites
#	 all positions are mappable in our case, thus we do not put this in our params
####################### log #######################

from __future__ import division
import os,sys
import subprocess as sp
from optparse import OptionParser
import traceback
import numpy as np

parser=OptionParser()

parser.add_option(
	'-D',
	'--raw-dir',
	dest='rawdir',
	help='rawdir containing target bam file, should be sorted and indexed'
	)
	
parser.add_option(
	'-T',
	'--tmp-dir',
	dest='tmpdir',
	help='temporary dir for keeping all intermedian files, should be given in form of full path'
	)
	
parser.add_option(
	'-B',
	'--bed-file',
	dest='bedFile',
	help='full path to bedFile, default set to: /media/disk2/ljzhang/data/DMD_plugin/DMD100_Targets.bed',
	default='/media/disk2/ljzhang/data/DMD_plugin/DMD100_Targets.bed'
	)
	
parser.add_option(
	'-P',
	'--plugin-dir',
	dest='pluginDir',
	help='plugin dir containing bed file, default set to: /media/disk2/ljzhang/data/DMD_plugin',
	default='/media/disk2/ljzhang/data/DMD_plugin'
	)
	
parser.add_option(
	'-F',
	'--ref-fasta',
	dest='refFa',
	help='reference fasta, default set to: /media/disk2/ljzhang/data/hg19/ucsc.hg19.chrX.fasta',
	default='/media/disk2/ljzhang/data/hg19/ucsc.hg19.chrX.fasta'
	)
	
parser.add_option(
	'-I',
	'--input-bam',
	dest='inBam',
	help='input bam, should be sorted and indexed'
	)
	
parser.add_option(
	'-A',
	'--a-value',
	dest='aval',
	help='a value, by default set to 70',
	default=70,
	type='int'
	)
	
parser.add_option(
	'-L',
	'--l-value',
	dest='lval',
	help='l value, by default set to 30',
	default=30,
	type='int'
	)

(options,args)=parser.parse_args()

if not options.rawdir or not options.tmpdir or not options.inBam:
	parser.print_help()
	sys.exit(1)
	
if not os.path.exists(options.rawdir) or not os.path.exists(options.inBam):
	parser.print_help()
	sys.exit(1)
	
if not os.path.exists(options.tmpdir):
	os.mkdir(options.tmpdir)
	
def main():
	global options
	mainExc(options.rawdir,options.tmpdir,options.pluginDir,options.inBam,options.bedFile,options.refFa,options.aval,options.lval)
	return
	
def mainExc(rawdir,tmpdir,pluginDir,inBam,bedFile,refFasta,a,l):
	percentile=[0.3,0.7]
	getBedReads(rawdir,tmpdir,pluginDir,inBam,bedFile)
	depthFile=tmpdir+'/'+inBam.rstrip('bam')+'targetBed.depth'
	bedIntvDct=getGcRef(rawdir,tmpdir,bedFile,depthFile,refFasta,pluginDir,a,l,percentile[0],percentile[1])
	bedDepthInfoFile=tmpdir+'/'+inBam.rstrip('bam')+'targetBed.crcted.depth.info'
	with open(bedDepthInfoFile,'w') as infoOutfh:
		srtedKeys=sorted(bedIntvDct.keys(),key=lambda x:int(x))
		for key in srtedKeys:
			infoOutfh.write(str(key)+'\t'+'\t'.join([str(x) for x in bedIntvDct[key]])+'\n')
	return

def getBedReads(rawdir,tmpdir,pluginDir,inBam,bedFile):
	'''
	get out all reads info foreach bedFile line
	note that in the temporary fileholder, there's filtered sam file and depth file according to start positions
	'''
	os.chdir(rawdir)
	outSam=tmpdir+'/'+inBam.rstrip('bam')+'targetBed.sam'
	try:
		sp.call('rm '+outSam,shell=True)
	except Exception:
		traceback.print_exc()
		pass
	with open(bedFile,'r') as bedInfh:
		for line in iter(bedInfh):
			linear=line.strip().split('\t')
			samCmd=pluginDir+'/samtools view '+inBam+' -q 20 chrX:'+linear[1]+'-'+linear[2]+' >> '+outSam
			sp.call(samCmd,shell=True)
	outBedDepth=outSam.rstrip('sam')+'depth'
	with open(outSam,'r') as samInfh,open(outBedDepth,'w') as depthOutFh:
		posDpthDct={}
		for line in iter(samInfh):
			linear=line.strip().split('\t')
			if posDpthDct.has_key(linear[3]):
				posDpthDct[linear[3]] += 1
			else:
				posDpthDct[linear[3]] = 1
		positions=sorted(posDpthDct.keys(),key=lambda x:int(x))
		for pos in positions:
			depthOutFh.write(pos+'\t'+str(posDpthDct[pos])+'\n')
	return
	
def getGcRef(rawdir,tmpdir,bedFile,depthFile,refFasta,pluginDir,a,l,*percentile):
	'''
	get gc reference data for calculating gc-fragment ratio: lambda
	the calculating method is based on the coverage info, use coverage quantile of 30-70% for position sampling
	use bed file, cut each interval to half and make sure that the cut tag length falls into 100-150
	'''
	os.chdir(tmpdir)
	'''
	depthList:
		[int(position),int(reads)]
	'''
	depthList=[]
	with open(depthFile,'r') as infh:
		for line in iter(infh):
			linear=line.strip().split('\t')
			depthList.append([int(linear[0]),int(linear[1])])
	'''
	bedIntv:
		[int(stt),int(intvEnd),int(intvDpth)]
	'''
	bedIntv=[]
	with open(bedFile,'r') as bedInfh:
		for line in iter(bedInfh):
			linear=line.strip().split('\t')
			stt=int(linear[1])
			end=int(linear[2])
			if abs(stt-end)>=150:
				medpos=(stt+end)//2
				bedIntv.append([stt,medpos])
				bedIntv.append([medpos,end])
			else:
				bedIntv.append([stt,end])
	for subBedIntv in bedIntv:
		stt,end=subBedIntv
		depth=sum([x[1] for x in depthList if stt<=x[0]<=end])
		subBedIntv.append(depth)
	univDpth=sorted([x[2] for x in bedIntv])
	print percentile
	intvGap=np.percentile(univDpth,percentile)
	'''
	quanGap:
		[int(stt),int(intvEnd),int(intvDpth)]
	'''
	quanGapPosList=[x for x in bedIntv if intvGap[0]<=x[2]<=intvGap[1]]
	gcCorDct=gcLambda(quanGapPosList,refFasta,pluginDir,depthList,a,l)
	bedIntvDct=crctDpthGc(refFasta,pluginDir,bedFile,depthList,gcCorDct,a,l)
	with open('gcCorDct.txt','w') as outfh:
		srtedKeys=sorted(gcCorDct.keys(),key=lambda x:int(x))
		for key in srtedKeys:
			outfh.write(str(key)+'\t'+str(gcCorDct[key])+'\n')
	return bedIntvDct

def gcLambda(quanGapPosList,refFasta,pluginDir,depthList,a,l):
	'''
	returns gc lambda for each gc content for correct depth file
	gcCorDct:
		gcContent => lambda
	'''
	sttPosList=[x[0] for x in quanGapPosList]
	gcDct={}
	for i in range(l+1):
		gcDct[i]=[0,0]
	depthDct=dict(depthList)
	for sttPos in sttPosList:
		samCmd=pluginDir+'/samtools faidx '+refFasta+' chrX:'+str(sttPos+a)+'-'+str(sttPos+a+l)
		faStrList=os.popen(samCmd).read().split('\n')
		faStr=''.join(faStrList[1:len(faStrList)]).upper()
		gcContent=faStr.count('G')+faStr.count('C')
		gcDct[gcContent][0]+=1
		if depthDct.has_key(sttPos):
			gcDct[gcContent][1]+=depthDct[sttPos]
	gcCorDct={}
	for key in gcDct.keys():
		if not gcDct[key][0] == 0:
			gcCorDct[key]=round(gcDct[key][1]/gcDct[key][0])
		else:
			gcCorDct[key] = 0
	return gcCorDct

def crctDpthGc(refFasta,pluginDir,bedFile,depthList,gcCorDct,a,l):
	'''
	depthList:
		[int(position),int(reads)]
	gcCorDct:
		gcContent => lambda
	bedIntvDct:
		start position => [int(original depth),int(corrected depth), int(dist betw crted depth and the original)]
	'''
	bedPos=[]
	with open(bedFile,'r') as bedInfh:
		for line in iter(bedInfh):
			linear=line.strip().split('\t')
			bedPos.append([int(linear[1]),int(linear[2])])
	bedIntvDct={}
	for subBedPos in bedPos:
		depth=sum([x[1] for x in depthList if subBedPos[0]<=x[0]<=subBedPos[1]])/abs(subBedPos[1]-subBedPos[0])
		bedIntvDct[subBedPos[0]]=[depth]
	'''
	calculating GC frame according to each 
	'''
	for subBedPos in bedPos:
		stt=subBedPos[0]
		end=subBedPos[1]
		next=stt+a+l
		gcIntvCont={}
		while next <= end:
			samCmd=pluginDir+'/samtools faidx '+refFasta+' chrX:'+str(stt+a)+'-'+str(stt+a+l)
			faStrList=os.popen(samCmd).read().split('\n')
			faStr=''.join(faStrList[1:len(faStrList)]).upper()
			gcCt=faStr.count('G')+faStr.count('C')
			if gcIntvCont.has_key(gcCt):
				gcIntvCont[gcCt] += 1
			else:
				gcIntvCont[gcCt] = 1
		gcCrctDepth=0
		for key in gcIntvCont.keys():
			if gcCorDct.has_key(key):
				gcCrctDepth += int(round(gcIntvCont[key] * gcCorDct[key],0))
		oriDepth=bedIntvDct[subBedPos[0]][0]
		bedIntvDct[subBedPos[0]].append(gcCrctDepth)
		bedIntvDct[subBedPos[0]].append(gcCrctDepth-oriDepth)
	bedIntvDct=pvalIntv(bedIntvDct)
	return bedIntvDct

def pvalIntv(bedIntvDct):
	'''
	bedIntvDct now in form of:
		bedSttPos => [int(ori depth),int(crcted depth), int(ori crcted gap), float(p value of gap)]
	'''
	allGap=[x[2] for x in bedIntvDct.values()]
	meanGap=np.mean(allGap)
	sdGap=np.std(allGap)
	for key in bedIntvDct.keys():
		pval = round((bedIntvDct[key][2]-meanGap)/sdGap,5)
		bedIntvDct[key].append(pval)
	return bedIntvDct

if __name__=='__main__':
	main()


