# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 10:57:28 2017

@author: thor
"""

#################### log ####################
# clarify the syntex and control for the bed input
# should include a plot procedure, with each interval depth scaled by the length of the region
# also note that in 2M pannal, all DMD region is covered by probes, use dmdExons as bedFile
# note that we have extend dmd region with 300 bp range
#############################################

from __future__ import division
import re,glob,sys,os
from optparse import OptionParser
import subprocess as sp
#import shutil
from multiprocessing import Pool


parser=OptionParser()

parser.add_option(
	'-D',
	'--sortedBam-dir',
	dest='sortedDir',
	help='directory containing mapped bams, if _sorted.bam exists the better'
	)
	
parser.add_option(
	'-P',
	'--plugin-dir',
	dest='pluginDir',
	help='plugin dir containing bed file, default set to: /media/disk2/ljzhang/data/DMD_plugin',
	default='/media/disk2/ljzhang/data/DMD_plugin'
	)
	
parser.add_option(
	'-B',
	'--bed-file',
	dest='bedFile',
	help='complete path to bed file, default set to: /media/disk2/ljzhang/data/dmdExons/dmd.exons.txt, in fmt of EXON57  31514904        31515061',
	default='/media/disk2/ljzhang/data/dmdExons/dmd.exons.txt'
	)
	
parser.add_option(
	'-R',
	'--res-dir',
	dest='resDir',
	help='results dir using for containing all bam stats results, full path should be given'
	)

parser.add_option(
	'-U',
	'--cpu-nbr',
	dest='cpuNbr',
	help='cpu nbr used for processing files, default set to 1',
	default='1'
	)
	
parser.add_option(
	'-M',
	'--processing-mode',
	dest='pmode',
	help='processing mode, 0 stands for multiprocessing, 1 stands for single file processing, default set to 0',
	default='0'
	)
	
parser.add_option(
	'-I',
	'--input-bam',
	dest='tgtBam',
	help='single processing bam input, name should be like: IonXpress_018_rawlib.bam'
	)
	
(options,args)=parser.parse_args()

if not options.sortedDir or not options.resDir:
	parser.print_help()
	sys.exit(1)
	
if options.pmode == '1' and not options.tgtBam:
	parser.print_help()
	sys.exit(1)

if not os.path.exists(options.resDir):
	os.mkdir(options.resDir)
	
def main():
	print 'main procedure starts ...'
	global options
	if not os.path.exists(options.resDir+'/targetSams'):
		os.mkdir(options.resDir+'/targetSams')
	bedDct=getBed(options.bedFile)
	checkPreFiles(options.sortedDir,options.resDir,options.pluginDir,options.cpuNbr)
	processingSam(options.resDir,options.cpuNbr,bedDct)
	return

def getBed(dmdExon):
	'''
	getting dmd bed File from dmdExon
	format like:
		EXON57  31514904        31515061
	should return in form like:
		bedDct['EXON57']=[31514904,31515061]
		in type of int
	'''
	bedDct={}
	with open(dmdExon,'r') as infh:
		intronLst=[]
		ind = 0
		for line in iter(infh):
			linear=line.strip().split('\t')
			exonNbr=linear[0].lstrip('EXON')
			bedDct[linear[0]]=[int(linear[1]),int(linear[2])]
			if ind == 0:
				intronLst.append(int(linear[2]))
				ind = 1
				continue
			intronLst.append(int(linear[1]))
			bedDct['INTRON'+exonNbr]=intronLst
			intronLst=[int(linear[2])]
	print 'bedDct generated .. '
	return bedDct

def checkPreFiles(rawDir,resDir,pluginDir,cpuNbr):
	'''
	check bam files mapped
	check mapped bams sorted
	check sorted bams transformed into sam format
	'''
	os.chdir(rawDir)
	mappedBams=glob.glob('*rawlib.bam')
	if mappedBams == []:
		rawBams=glob.glob('*basecaller.bam')
		if rawBams == []:
			print 'no bam file found at '+rawDir
		mapPool=Pool(cpuNbr)
		for rawBam in rawBams:
			mapPool.apply_async(mapBam,args=(rawDir,rawBam,pluginDir))
		mapPool.close()
		mapPool.join()

	mappedBams=glob.glob('*rawlib.bam')
	sortPool=Pool(cpuNbr)
	for mappedBam in mappedBams:
		sortPool.apply_async(sortBam,rawDir,mappedBam,pluginDir)
	sortPool.close()
	sortPool.join()
	
	samPool=Pool(cpuNbr)
	for mappedBam in mappedBams:
		samPool.apply_async(getSam,args=(resDir,mappedBam,pluginDir))
	samPool.close()
	samPool.join()
	
	print 'all bam mapped and sorted and samFiles are ready...'
	return

def mapBam(rawDir,rawBam,pluginDir):
	os.chdir(rawDir)
	print 'tmap for '+rawBam+' starts ...'
	ionXpress=re.findall(r'IonXpress\_\d+',rawBam)[0]
	tmapCmd='%s mapall -f %s -r %s -i bam -v -Y -o 2 -a 0 stage1 map1 map2 map3 map4 > %s' % (pluginDir+'/tmap',pluginDir+'/hg19/hg19.fasta',rawBam,rawDir+'/'+ionXpress+'_rawlib.bam')
	sp.call(tmapCmd,shell=True)
	return

def sortBam(rawDir,mappedBam,pluginDir):
	os.chdir(rawDir)
	ionXpress=re.findall(r'IonXpress\_\d+',mappedBam)[0]
	if not os.path.exists(ionXpress+'_rawlib.sorted'):
		sortCmd=pluginDir+'/samtools sort %s %s' % (mappedBam,ionXpress+'_rawlib.sorted')
		sp.call(sortCmd,shell=True)
		if not os.path.exists(ionXpress+'_rawlib.sorted.bam.bai'):
			sp.call(pluginDir+'/samtools index %s' % (ionXpress+'_rawlib.sorted.bam'),shell=True)
	return

def getSam(resDir,mappedBam,pluginDir):
	os.chdir(resDir+'/targetSams')
	ionXpress=re.findall(r'IonXpress\_\d+',mappedBam)[0]
	sortedBam=ionXpress+'_rawlib_sorted.bam'
	samFile=ionXpress+'_rawlib_sorted.sam'
	if not os.path.exists(samFile):
		samCmd=pluginDir+'/samtools view %s > %s' % (sortedBam,resDir+'/targetSams/'+samFile)
		sp.call(samCmd,shell=True)
	return

def processingSam(resDir,cpuNbr,bedDct):
	os.chdir(resDir)
	samFiles=glob.glob('*.sam')
	
	proPool=Pool(cpuNbr)
	for samFile in samFiles:
		proPool.apply_async(samInfoProc,args=(resDir,samFile,bedDct))
	proPool.close()
	proPool.join()
		
	return

def samInfoProc(resDir,samFile,bedDct):
	'''
	bedDct:
	format like:
		EXON57  31514904        31515061
	should return in form like:
		bedDct['EXON57']=[31514904,31515061]
		in type of int
	'''
	os.chdir(resDir)
	
	infh = open(resDir+'/targetSams/'+samFile,'r')
	
	totalReads=0
	totalBases=0
	
	basegr20=0
	basegr30=0
	
	mappedReads=0
	mappedBases=0
	
	onTgtReads=0
	onTgtBases=0
	
	dupReads=0
	
	x1bases=0
	x10bases=0
	x20bases=0
	x30bases=0
	x50bases=0
	x100bases=0
	x200bases=0
	
	mismReads=0
	mismBases=0
	
	insReads=0
	insBases=0
	delReads=0
	delBases=0
	
	mappedReadLen=0
	mappedReadTtlen=0
	
	unmappedReadLen=0
	unmappedReadTtlen=0
	unmappedReads=0
	
	tgtReadLen=0
	
	univRangDct={}
	
	for line in infh.xreadlines():
		if line.startswith('@'):
			continue
		
		linear=line.strip().split('\t')
		
#		get total reads and total bases
		totalReads += 1
		totalBases += len(linear[9])
		
		quals=linear[10]
		for qual in quals:
			score=ord(qual) - 33
#			getting base quality greater than 20 and 30
			if score >= 30:
				basegr30 += 1
			elif score >= 20:
				basegr20 += 1
		
		if linear[1] == '4':
#			get out unmapped reads
			unmappedReads += 1
			unmappedReadTtlen += len(linear[9])
			continue
		
#		reads that are not unmapped is considerred as mapped reads
		mappedReads += 1
		mappedReadTtlen += len(linear[9])
		mappedBases += len(linear[9])
		
#==============================================================================
# 		if linear[1] == '1024':
##			get out duplicated reads
# 			dupReads += 1
##		this returns dup reads of 0, using hash to calculate same start and ending position reads
#==============================================================================
		
		mdTag=linear[20]
#		remove deleted positions from mdTag
		mdNodel=re.sub('\^[AGCT]+','',mdTag)
		
#		mismatch info counts has error here
#		mismatch reads also returns 0
		mismAr=[i for i in mdNodel if i in 'TCGA']
		
		mismBases += len(mismAr)
		
		if len(mismAr) > 0:
			mismReads += 1
		'''
		DMD gene region on chrX:
			31137344,33146544
			extend head and tail by 300
		'''
		dmdRange=[31137044,33146844]
		
		if not linear[2] == 'chrX' or int(linear[3]) < dmdRange[0] or int(linear[3]) > dmdRange[-1]:
			continue
		
		insPat=re.compile(r'(\d+)I')
		
		if 'I' in linear[5]:
			insReads += 1
			insLst=insPat.findall(linear[5])
			insLst=[int(i) for i in insLst]
			insBases += sum(insLst)
		
		delPat=re.compile(r'(\d+)D')
		if 'D' in linear[5]:
			delReads += 1
			delLst=delPat.findall(linear[5])
			delLst=[int(i) for i in delLst]
			delBases += sum(delLst)
			
		sttP=int(linear[3])
		
		softTpat=re.compile(r'(\d+)S$')
		hardTpat=re.compile(r'(\d+)H$')
		sftTnbr=sum([int(i) for i in softTpat.findall(linear[5])])
		hardTnbr=sum([int(i) for i in hardTpat.findall(linear[5])])
		
		endP=sttP+len(linear[9])-sftTnbr-hardTnbr+1
		poskey=str(sttP)+'_'+str(endP)
		
#		univRangDct is used for storing read stt-end ranging mapping according to each read
		if univRangDct.has_key(poskey):
			univRangDct[poskey] += 1
		else:
			univRangDct[poskey] = 1
	
	dupLst=[i for i in univRangDct.values() if i>1]
	dupReads += len(dupLst)
	'''
	sortedKeys stores starting and ending postitions, value is appearing time
	while sortedBedKeys has starts points, value corresponds to ending positions
	these would be applied for calculating nX coverage
	'''
	sortedKeys = sorted(univRangDct,key=lambda x:int(x.split('_')[0]))
	sortedBedValues=sorted(bedDct.values(),key=lambda x:x[0])
	
#	univPosDct stores point mapping reads info
	univPosDct={}
#	onTgtReadsLst is used for calculating uniq on target reads
	onTgtReadsLst=[]
#	sortedKeys are range str with stt_end info
	for sortedKey in sortedKeys:
		startP,tailP=[int(subs) for subs in sortedKey.split('_')]
		
		for valuePair in sortedBedValues:
			valueStt=valuePair[0]
			valueEnd=valuePair[1]
			if startP > valueEnd or tailP < valueStt:
				continue
			for pointpos in range(valueStt,valueEnd+1):
				tgtReadLen += (tailP-startP+1)
				onTgtReadsLst.append(sortedKey)
				if univPosDct.has_key(pointpos):
					univPosDct[pointpos] += 1
				else:
					univPosDct[pointpos] = 1
	
	onTgtReadsLst = list(set(onTgtReadsLst))
	onTgtReads += len(onTgtReadsLst)
	
	tgtTotalLen = 0
	for subTgtRead in onTgtReadsLst:
		subsst,subend=subTgtRead.split('_')
		tgtTotalLen += (int(subend)-int(subsst)+1)
	
	tgtRegionLen=round(tgtTotalLen/onTgtReads)
	
	sortedPosKeys=sorted(univPosDct.keys())
	totalReadsDpt=0
	
	for subPosKey in sortedPosKeys:
		onTgtBases += univPosDct[subPosKey]
		
		if univPosDct[subPosKey] >= 200:
			x1bases += 1
			x10bases += 1
			x20bases += 1
			x30bases += 1
			x50bases += 1
			x100bases += 1
			x200bases += 1
		elif univPosDct[subPosKey] >= 100:
			x1bases += 1
			x10bases += 1
			x20bases += 1
			x30bases += 1
			x50bases += 1
			x100bases += 1
		elif univPosDct[subPosKey] >= 50:
			x1bases += 1
			x10bases += 1
			x20bases += 1
			x30bases += 1
			x50bases += 1
		elif univPosDct[subPosKey] >= 30:
			x1bases += 1
			x10bases += 1
			x20bases += 1
			x30bases += 1
		elif univPosDct[subPosKey] >= 20:
			x1bases += 1
			x10bases += 1
			x20bases += 1
		elif univPosDct[subPosKey] >=10:
			x1bases += 1
			x10bases += 1
		elif univPosDct[subPosKey] >= 1:
			x1bases += 1
			
		totalReadsDpt += univPosDct[subPosKey]
	
	tgtSize=dmdRange[1]-dmdRange[0]
	
	x1rate=round(x1bases/tgtSize,4)*100
	x10rate=round(x10bases/tgtSize,4)*100
	x20rate=round(x20bases/tgtSize,4)*100
	x30rate=round(x30bases/tgtSize,4)*100
	x50rate=round(x50bases/tgtSize,4)*100
	x100rate=round(x100bases/tgtSize,4)*100
	x200rate=round(x200bases/tgtSize,4)*100
	
	mappedReadLen += round(mappedReadTtlen/mappedReads)
	unmappedReadLen += round(unmappedReadTtlen/unmappedReads)
	
	onTgtReadRate=round(onTgtReads/totalReads,4)
	onTgtBaseRate=round(onTgtBases/totalBases,4)
	dupReadsRate=round(dupReads/totalReads,4)
	mismReadRate=round(mismReads/totalReads,4)
	mismBaseRate=round(mismBases/totalBases,4)
	insReadRate=round(insReads/totalReads,4)
	insBaseRate=round(insBases/totalBases,4)
	delReadRate=round(delReads/totalReads,4)
	delBaseRate=round(delBases/totalBases,4)
	mappedReadRate=round(mappedReads/totalReads,4)
	mappedBaseRate=round(mappedBases/totalBases,4)
	
	tgtMeanDpt=round(onTgtBases/tgtSize)
	
	basegr20pct=round(basegr20/totalBases,4)
	basegr30pct=round(basegr30/totalBases,4)
	meanReadLen=round((mappedReadTtlen+unmappedReadTtlen)/totalReads)
	evenScMe=evenScore(univPosDct,tgtSize,tgtMeanDpt)
	
	basicStatFile=re.findall(r'IonXpress_\d+',samFile)[0]+'_Basic_statistics_for_mapping_data.txt'
	
	outfh=open(resDir+'/'+basicStatFile,'w')
	
	outfh.write('target_region_size(bp)\t'+str(tgtSize)+'\n')
	outfh.write('total_reads\t'+str(totalReads)+'\n')
	outfh.write('total_bases\t'+str(totalBases)+'\n')

	outfh.write('percent_of_bases_with_quality >= 20\t'+str(100*basegr20pct)+'%\n')
	outfh.write('percent_of_bases_with_quality >= 30\t'+str(100*basegr30pct)+'%\n')

	outfh.write('mapped_reads\t'+str(mappedReads)+'\n')
	outfh.write('mapped_reads_rate\t'+str(100*mappedReadRate)+'%\n')
	outfh.write('mapped_bases\t'+str(mappedBases)+'\n')
	outfh.write('mapped_base_rate\t'+str(100*mappedBaseRate)+'%\n')

	outfh.write('on-target_reads\t'+str(onTgtReads)+'\n')
	outfh.write('on-target_read_rate\t'+str(100*onTgtReadRate)+'%\n')
	outfh.write('on-target_bases\t'+str(onTgtBases)+'\n')
	outfh.write('on-target_base_rate\t'+str(100*onTgtBaseRate)+'%\n')
	
	outfh.write('duplicated_reads\t'+str(dupReads)+'\n')
	outfh.write('duplicated_read_rate\t'+str(100*dupReadsRate)+'%\n')
	
	outfh.write('1X_coverage_of_target_region\t'+str(x1rate)+'%\n')
	outfh.write('10X_coverage_of_target_region\t'+str(x10rate)+'%\n')
	outfh.write('20X_coverage_of_target_region\t'+str(x20rate)+'%\n')
	outfh.write('30X_coverage_of_target_region\t'+str(x30rate)+'%\n')
	outfh.write('50X_coverage_of_target_region\t'+str(x50rate)+'%\n')
	outfh.write('100X_coverage_of_target_region\t'+str(x100rate)+'%\n')
	outfh.write('200X_coverage_of_target_region\t'+str(x200rate)+'%\n')
	
	outfh.write('mean_depth_of_target_region\t'+str(tgtMeanDpt)+'\n')
	outfh.write('mismatch_reads\t'+str(mismReads)+'\n')
	outfh.write('mismatch_read_rate\t'+str(mismReadRate*100)+'%\n')
#	outfh.write('mismatch bases\t'+str(mismBases)+'\n')
	outfh.write('mismatch_base_rate\t'+str(mismBaseRate*100)+'%\n')
	outfh.write('insertion_read_rate\t'+str(insReadRate*100)+'%\n')
	outfh.write('insertion_base_rate\t'+str(insBaseRate*100)+'%\n')
	outfh.write('deletion_read_rate\t'+str(delReadRate*100)+'%\n')
	outfh.write('deletion_base_rate\t'+str(delBaseRate*100)+'%\n')
	
	outfh.write('eveness_score\t'+str(evenScMe)+'\n')
	
	outfh.write('mean_read_length\t'+str(meanReadLen)+'\n')
	outfh.write('mapped_read_length\t'+str(mappedReadLen)+'\n')
	outfh.write('unmapped_read_length\t'+str(unmappedReadLen)+'\n')
	outfh.write('target_region_read_length\t'+str(tgtRegionLen)+'\n')
	
	infh.close()
	outfh.close()
	
	return


def evenScore(univPosDctEs,tgtSizeEs,meanCovEs):
	evenSc=0
	covLst=univPosDctEs.values()
	covLst.sort(reverse=True)
	
	for i in range(1,int(meanCovEs)+1):
		ct=countNbr(covLst,i)
		evenSc += round(ct/(tgtSizeEs*meanCovEs),3)
	
	return evenSc

def countNbr(covLstCn,candi):
	ctCn=0
	for j in covLstCn:
		if candi <= j:
			ctCn += 1
		else:
			return ctCn
	return ctCn

if __name__=='__main__':
	main()
