# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 11:36:13 2017

@author: thor
"""


#==============================================================================
# this is for implement a pcf function using dynamic programming
# this procedure would be based on PCF function
# note that pcf matrix structure is modified
#==============================================================================

from __future__ import division
import os
import sys
#==============================================================================
# from copy import deepcopy
# from operator import attrgetter
#==============================================================================
from optparse import OptionParser
import itertools
import math

parser=OptionParser()
parser.add_option(
	'-D',
	'--raw-dir',
	dest='rawDir',
	help='rawDir containning target depth file and also as work dir'
	)
	
parser.add_option(
	'-B',
	'--bed-file',
	dest='bedFile',
	help='full path to bedFile, default set to: /media/disk2/ljzhang/data/DMD_plugin/DMD100_Targets.bed',
	default='/media/disk2/ljzhang/data/DMD_plugin/DMD100_Targets.bed'
	)
	
parser.add_option(
	'-G',
	'--gamma-value',
	dest='gamma',
	help='gamma value, default set to: 50',
	default=50,
	type='int'
	)
	
parser.add_option(
	'-E',
	'--gap-thres',
	dest='gapThres',
	help='gap threshold for defining window size, default set to 100',
	default=100,
	type='int'
	)
	
parser.add_option(
	'-I',
	'--input-depth',
	dest='dpthFile',
	help='depth file, ends with '
	)
	
(options,args)=parser.parse_args()

if not options.rawDir or not options.dpthFile:
	parser.print_help()
	sys.exit(1)
	
rawDir=options.rawDir
dpthFile=options.dpthFile
bedFile=options.bedFile
gamma=options.gamma
gapThres=options.gapThres


def main():
	global rawDir,dpthFile,bedFile,gapThres,gamma,gapThres
	os.chdir(rawDir)
	
	bedLst=getBed(bedFile,gapThres)
	'''
	deptObjList:
		[DepthIntv,DepthIntv, ...]
	DepthIntv:
		stt, end, meanDpt(mean depth over position)
	'''
	print('generating depthObjList starts ...')
	depthObjList=splitDpth(dpthFile,bedLst)
	depthLen=len(depthObjList)
#	print(depthLen)
	
	sdVal=calSig(depthObjList)
	
	print('using dynamic programming for getting list ...')
	gammaSpl=gamma*sdVal
	idxList=dynamicPcf(depthObjList,depthLen,gammaSpl)
	print(idxList)
	
	idxIntv=cutIdxToIntv(idxList)
	
	
#==============================================================================
# 	'''
# 	use a plot for plotting possible CNV region
# 	'''
# 	fileRawName=dpthFile.rstrip('dpth')
# 	print('data loading to file starts ...')
# 	writeDpthIndx(depthObjList,maxValPath,fileRawName+'depth.df',fileRawName+'.indx')
#==============================================================================
	
	return
	
def idxValIntv(depthObjList,idxIntv):
	'''
	this is for plotting method
	calculate winsorized depth for given intervals
	also note that in the final result, using plots for all depth and use line for intervals
	'''
	
	return
	
def cutIdxToIntv(idxList):
	idxList.sort()
	idxIntvList=[(idxList[i],idxList[i+1]) for i in range(len(idxList)-1)]
	return idxIntvList
	
def calSig(depthObjList):
	meanDptLst=[i.meanDpt for i in depthObjList]
	meanVal=sum(meanDptLst)/len(meanDptLst)
	sdLst=[(j-meanVal)**2 for j in meanDptLst]
	sdVal=sum(sdLst)/len(meanDptLst)
	sdVal=math.sqrt(sdVal)
	return sdVal
	
def dynamicPcf(depthObjList,depthLen,gamma):
	'''
	the key idea of this pcf function is to find the segment which minimize the variance
	'''
	indexList=[]
	i=1
	ak=[]
	ek=[0]
	while i<depthLen:
		ak.append(0)
		ak=[j+depthObjList[i-1].meanDpt for j in ak]
		dk=[-ak[k]*ak[k]/(len(ak)-k) for k in range(len(ak))]
#		print('print dk ...')
#		print(dk)
		ekMinList=[dk[id]+ek[id]+gamma for id in range(i)]
#		print('print ekMinList ...')
#		print(ekMinList)
#		print('\n\n')
		ekMin=min(ekMinList)
		indexList.append(ekMinList.index(ekMin))
		ek.append(ekMin)
		i += 1
	return indexList
	
def writeDpthIndx(depthObjList,maxValPath,dfFile,dfIndx):
	with open(dfFile,'w') as dffh, open(dfIndx,'w') as idxfh:
		for depthObj in depthObjList:
			dffh.write(str(depthObj.stt)+'\t'+str(depthObj.end)+'\t'+str(depthObj.meanDpt))
		for valPath in maxValPath:
			for position in valPath:
				idxfh.write(str(depthObjList[position].stt)+'\t')
			idxfh.write('\n')
	return
	
def rtrMaxPathMat(pathValueLst):
	maxIndxLst=[0]
	mark=pathValueLst[0]
	i=1
	pathLen=len(pathValueLst)
	while i < pathLen:
		if pathValueLst[i] == mark:
			maxIndxLst.append(i)
		elif pathValueLst[i] > mark:
			maxIndxLst=[i]
			mark = pathValueLst[i]
	return maxIndxLst

def valueForAllPath(pcfMat,allPath,gamma):
	pathValueLst=[]
	for index in allPath:
		valueList=[pcfMat[i] for i in index]
		segNbr=len(index)-1
		value=sum(valueList)-gamma*segNbr
		pathValueLst.append(value)
	return pathValueLst
	
def rtnAllPathIter(depthLen):
	allPathList=[]
	gap=depthLen - 2
	while gap >= 0:
		allComb=itertools.combinations(range(1,depthLen-1),gap)
		for comb in allComb:
			comb=list(comb)
			comb.append(depthLen-1)
			comb.insert(0,0)
			allPathList.append(comb)
		gap -= 1
	
	return allPathList
	
def getBed(bedFile,gapThres):
	bedLst=[]
	with open(bedFile,'r') as infh:
		for line in iter(infh):
			linear=line.strip().split('\t')
			stt = int(linear[1])
			end = int(linear[2])
			while (stt + gapThres) <= end:
				bedLst.append([stt,stt + gapThres])
				stt += gapThres
			if stt != end:
				bedLst.append([stt,end])
	return bedLst
	
def splitDpth(dpthFile,bedLst):
	depthDct=transBed2Dct(bedLst)
	with open(dpthFile,'r') as infh:
		for line in iter(infh):
			if not line.startswith('chrX'):
				continue
			linear=line.strip().split('\t')
			chrXpos=int(linear[1])
			for subPosLst in bedLst:
				if subPosLst[0] <= chrXpos <= subPosLst[1]:
					depthDct[str(subPosLst[0])+'_'+str(subPosLst[1])] += int(linear[2])
					break
	deptObjList=transDpthDct(depthDct)
	return deptObjList
	
def transBed2Dct(bedLst):
	bedDct={}
	for subLst in bedLst:
		bedDct[str(subLst[0])+'_'+str(subLst[1])]=0
	return bedDct
	
def transDpthDct(depthDct):
	'''
	returns a list of DepthIntv objects
	'''
	depthInfoList=[]
	for key in depthDct.keys():
		depObj=DepthIntv(key,depthDct[key])
		depthInfoList.append(depObj)
	depthInfoList=sorted(depthInfoList,key=lambda x:x.stt)
	return depthInfoList
	
class DepthIntv:
	'''
	this would be viewed as continuous windows
	'''
	def __init__(self,subPosStr,depth):
		subPosList=[int(i) for i in subPosStr.split('_')]
		self.stt=subPosList[0]
		self.end=subPosList[1]
		self.meanDpt=round(0.01*depth/abs(self.end-self.stt),3)

if __name__=='__main__':
	main()

