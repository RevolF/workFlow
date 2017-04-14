# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 09:39:17 2017

@author: thor
"""

#################### log ####################
#	this is for getting out bam stat file after given a bed file 
#	given bed file should be specified as:
#		segmentID     startPos    endPos
#	if not for DMD, specific bed file should be given in params
#	for easy computing burden, modify bedDct to list type
#		and this list is in right order sorted
#
#	20170413:
#		try to add a function for generating foreach bam a dataframe containning pverall depth for each bed intv
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
	help='complete path to bed file, default set to: /media/disk2/ljzhang/data/DMD_plugin/DMD100_Targets.bed',
	default='/media/disk2/ljzhang/data/DMD_plugin/DMD100_Targets.bed'
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
	type = 'int',
	default=1
	)
	
#==============================================================================
# parser.add_option(
# 	'-M',
# 	'--processing-mode',
# 	dest='pmode',
# 	help='processing mode, 0 stands for multiprocessing, 1 stands for single file processing, default set to 0',
# 	default='0'
# 	)
# 	
# parser.add_option(
# 	'-I',
# 	'--input-bam',
# 	dest='tgtBam',
# 	help='single processing bam input, name should be like: IonXpress_018_rawlib.bam'
# 	)
#==============================================================================
	
parser.add_option(
	'-X',
	'--remove-files',
	dest='rmFile',
	help='signal of if or not remove all intermediate files, y stands for yes and n for no, default set to n',
	default='n'
	)
	
(options,args)=parser.parse_args()

if not options.sortedDir or not options.resDir:
	parser.print_help()
	sys.exit(1)
	
#==============================================================================
# if options.pmode == '1' and not options.tgtBam:
# 	parser.print_help()
# 	sys.exit(1)
#==============================================================================

if not os.path.exists(options.resDir):
	os.mkdir(options.resDir)
	
def main():
	print 'main procedure starts ...'
	global options
	if not os.path.exists(options.resDir+'/targetSams'):
		os.mkdir(options.resDir+'/targetSams')
	bedDct=getBed(options.bedFile)
	checkPreFiles(options.sortedDir,options.resDir,options.pluginDir,options.cpuNbr)
	processingSam(options.sortedDir,options.resDir,options.cpuNbr,bedDct,options.pluginDir)
	if options.rmFile == 'y':
		sp.call('rm '+options.sortedDir+'/*sorted.bam',shell=True)
		sp.call('rm '+options.sortedDir+'/*sorted.bam.bai',shell=True)
		sp.call('rm '+options.sortedDir+'/*rawlib.bam',shell=True)
		sp.call('rm -rf '+options.resDir+'/targetSams')
	return

def getBed(dmdExon):
	'''
	getting dmd bed File from dmdExon
	format like:
		EXON57  31514904        31515061
	should return in form like:
		bedDct['EXON57']=[31514904,31515061]
		in type of int
	for other bam bed file, should make corresponding modifications
	'''
	'''
	bed briefed info
		X       10911688        10911888
		X       10914229        10914429
		X       28932233        28932483
		X       29003094        29003244
		X       30211769        30211913
		X       30885651        30885777
		X       30887994        30888198
		X       31137296        31140096
	returned bedDct is a list:
		[[31514904,31515061], ... ]
	'''
	bedDct=[]
	
	with open(dmdExon,'r') as infh:
		for line in iter(infh):
			linear=line.strip().split('\t')
			bedDct.append([int(linear[1]),int(linear[2])])
	print 'bedDct generated ...'
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
		sortPool.apply_async(sortBam,args=(rawDir,mappedBam,pluginDir))
	sortPool.close()
	sortPool.join()
	
	samPool=Pool(cpuNbr)
	for mappedBam in mappedBams:
		samPool.apply_async(getSam,args=(rawDir,resDir,mappedBam,pluginDir))
	samPool.close()
	samPool.join()
	
	print 'all bam mapped and sorted and samFiles are ready ...'
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
	if not os.path.exists(ionXpress+'_rawlib.sorted.bam'):
		print 'sort bam '+ mappedBam +' starts ...'
		sortCmd=pluginDir+'/samtools sort %s %s' % (mappedBam,ionXpress+'_rawlib.sorted')
		sp.call(sortCmd,shell=True)
		if not os.path.exists(ionXpress+'_rawlib.sorted.bam.bai'):
			sp.call(pluginDir+'/samtools index %s' % (ionXpress+'_rawlib.sorted.bam'),shell=True)
	return

def getSam(rawDir,resDir,mappedBam,pluginDir):
	os.chdir(rawDir)
	ionXpress=re.findall(r'IonXpress\_\d+',mappedBam)[0]
	sortedBam=ionXpress+'_rawlib.sorted.bam'
	
	samFile=ionXpress+'_rawlib.sorted.sam'
	if not os.path.exists(samFile):
		print 'samtools for getting sam file starts ...'
		samCmd=pluginDir+'/samtools view %s > %s' % (rawDir+'/'+sortedBam,resDir+'/targetSams/'+samFile)
		sp.call(samCmd,shell=True)
	return

def processingSam(rawDir,resDir,cpuNbr,bedDct,pluginDir):
	print 'processing all samFiles starts ...'
	os.chdir(resDir+'/targetSams')
	samFiles=glob.glob('*.sam')
	
	proPool=Pool(cpuNbr)
	for samFile in samFiles:
		proPool.apply_async(infoProc,args=(rawDir,resDir,samFile,bedDct,pluginDir))
	proPool.close()
	proPool.join()
	print 'all samFile processing completed ..'
	return

def infoProc(rawDir,resDir,samFile,bedDct,pluginDir):
	totalReads,totalBases,basegr20,basegr30,mappedReads,mappedBases,onTgtReads,dupReads,mismReads,mismBases,insReads,insBases,delReads,delBases,mappedReadTtlen,unmappedReadTtlen,unmappedReads,tgtTotalLen=samInfoProc(rawDir,resDir,samFile,bedDct)
	
	tgtRegionLen=round(tgtTotalLen/onTgtReads)
	
#	depthInfoProc(rawDir,sortedBam,pluginDir,bedDct)
	sortedBamDp=samFile.rstrip('sam')+'bam'
	x1bases,x10bases,x20bases,x30bases,x50bases,x100bases,x200bases,onTgtBases,univPosLst=depthInfoProc(rawDir,sortedBamDp,pluginDir,bedDct)
	
	tgtSize=0
	for subList in bedDct:
		tgtSize = tgtSize+abs(subList[1]-subList[0])+1
	
	x1rate=round(x1bases/tgtSize,4)*100
	x10rate=round(x10bases/tgtSize,4)*100
	x20rate=round(x20bases/tgtSize,4)*100
	x30rate=round(x30bases/tgtSize,4)*100
	x50rate=round(x50bases/tgtSize,4)*100
	x100rate=round(x100bases/tgtSize,4)*100
	x200rate=round(x200bases/tgtSize,4)*100
	
	mappedReadLen = round(mappedReadTtlen/mappedReads)
	unmappedReadLen = round(unmappedReadTtlen/unmappedReads)
	
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
#	evenScore(covLst,tgtSizeEs,meanCovEs)
	
	evenScMe=evenScore(univPosLst,tgtSize,tgtMeanDpt)
	
	basicStatFile=re.findall(r'IonXpress_\d+',samFile)[0]+'_Basic_statistics_for_mapping_data.txt'
	
	outfh=open(resDir+'/'+basicStatFile,'w')
	print 'starting writing info to '+resDir+'/'+basicStatFile
	
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
	
	outfh.close()
	
	return

def samInfoProc(rawDir,resDir,samFile,bedDct):
	'''
	bedDct:
	format like:
		[[31514904, 31515061], ... ]
	'''
	print 'samInfoProc for '+samFile+' starts ...'
	os.chdir(resDir)
	infh = open(resDir+'/targetSams/'+samFile,'r')
	
	totalReads=0
	totalBases=0
	
	basegr20=0
	basegr30=0
	
	mappedReads=0
	mappedBases=0
	
	onTgtReads=0
	
	dupReads=0
	
	mismReads=0
	mismBases=0
	
	insReads=0
	insBases=0
	delReads=0
	delBases=0
	
	mappedReadTtlen=0
	unmappedReadTtlen=0
	unmappedReads=0
	tgtTotalLen = 0
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
		
		sttP=int(linear[3])
		softTpat=re.compile(r'(\d+)S$')
		hardTpat=re.compile(r'(\d+)H$')
		sftTnbr=sum([int(i) for i in softTpat.findall(linear[5])])
		hardTnbr=sum([int(i) for i in hardTpat.findall(linear[5])])
		endP=sttP+len(linear[9])-sftTnbr-hardTnbr+1
		poskey=str(sttP)+'_'+str(endP)
		
#		univRangDct is used for storing read stt-end ranging mapping according to each read
		if univRangDct.has_key(poskey):
			dupReads += 1
		else:
			univRangDct[poskey] = 1
			
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
		
		if linear[2] != 'chrX' or endP < bedDct[0][0] or sttP > bedDct[-1][-1]:
			continue
		
		for subList in bedDct:
			if subList[0] <= sttP <= subList[1] or subList[0] <= endP <= subList[1]:
				onTgtReads += 1
				tgtTotalLen = tgtTotalLen + abs(endP - sttP) + 1
				break
	infh.close()
	print 'samInfoProc ends for '+samFile+' ...'
	return totalReads,totalBases,basegr20,basegr30,mappedReads,mappedBases,onTgtReads,dupReads,mismReads,mismBases,insReads,insBases,delReads,delBases,mappedReadTtlen,unmappedReadTtlen,unmappedReads,tgtTotalLen

def depthInfoProc(rawDir,sortedBam,pluginDir,bedDct):
	'''
	add processing depth info for each bed interval here
	'''
	os.chdir(rawDir)
	if not os.path.exists(sortedBam+'.dpth'):
		samDptCmd=pluginDir+'/samtools depth '+sortedBam+' > '+sortedBam+'.dpth'
		sp.call(samDptCmd,shell=True)
	
	print 'depthInfoProc starts ...'
	
	x1bases=0
	x10bases=0
	x20bases=0
	x30bases=0
	x50bases=0
	x100bases=0
	x200bases=0
	
	ontgtBases=0
	
	univPosDptLst=[]
	infoDct={}
	
	for item in bedDct:
		key=str(item[0])+'_'+str(item[1])
		infoDct[key]=0
	print infoDct
	with open(sortedBam+'.dpth','r') as infh, open(sortedBam+'.bedDpt.df','w') as outfh:
		print 'entered depth loop ...'
		for line in iter(infh):
			print line
			linear=line.strip().split('\t')
			if linear[0] != 'chrX' or int(linear[1]) < bedDct[0][0] or int(linear[1]) > bedDct[-1][-1]:
				continue
			
			for subList in bedDct:
				if subList[0] <= int(linear[1]) <= subList[1]:
					ontgtBases += int(linear[2])
					univPosDptLst.append(int(linear[2]))
					
					key=str(subList[0])+'_'+str(subList[1])
					infoDct[key] += int(linear[2])
					
					if int(linear[2]) >= 200:
						x1bases += 1
						x10bases += 1
						x20bases += 1
						x30bases += 1
						x50bases += 1
						x100bases += 1
						x200bases += 1
					elif int(linear[2]) >= 100:
						x1bases += 1
						x10bases += 1
						x20bases += 1
						x30bases += 1
						x50bases += 1
						x100bases += 1
					elif int(linear[2]) >= 50:
						x1bases += 1
						x10bases += 1
						x20bases += 1
						x30bases += 1
						x50bases += 1
					elif int(linear[2]) >= 30:
						x1bases += 1
						x10bases += 1
						x20bases += 1
						x30bases += 1
					elif int(linear[2]) >= 20:
						x1bases += 1
						x10bases += 1
						x20bases += 1
					elif int(linear[2]) >= 10:
						x1bases += 1
						x10bases += 1
					elif int(linear[2]) >= 1:
						x1bases += 1
					break
		for key in infoDct.keys():
			outfh.write(key+'\t'+str(infoDct[key])+'\n')
	print 'depthInfoProc ends for '+ sortedBam + ' ...'
	return x1bases,x10bases,x20bases,x30bases,x50bases,x100bases,x200bases,ontgtBases,univPosDptLst

def evenScore(covLst,tgtSizeEs,meanCovEs):
	evenSc=0
	covLst.sort(reverse=True)
	divbase=tgtSizeEs*meanCovEs

	for i in range(1,int(meanCovEs)+1):
		ct=countNbr(covLst,i)
		evenSc += ct/divbase
	evenSc = round(evenSc,3)
	return evenSc

def countNbr(covLstCn,candi):
	ctCn = 0
	for j in covLstCn:
		if j >= candi:
			ctCn += 1
		else:
			return ctCn
	return ctCn
	
if __name__=='__main__':
	main()

