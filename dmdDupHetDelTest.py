# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 16:48:00 2017

@author: thor
"""

#==============================================================================
# use deep copy to generate all tree
# generate from the smallest
# use itertools to avoid the huge mem consumption
#==============================================================================

from __future__ import division
import os
import sys
from copy import deepcopy
from operator import attrgetter
from optparse import OptionParser
import itertools

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

#==============================================================================
# os.chdir('D:\\thor\\documents\\WorkProjects\\DMD_Project\\dmdDupHeterDelTest')
# dpthFile='IonXpress_001695_rawlib.sorted.bam.dpth'
# bedFile='DMD100_Targets.bed'
# gamma=50
# gapThres=100
#==============================================================================

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
	print(depthLen)
	print('generating PCF matrix starts ...')
	pcfMat=PcfMat(depthObjList)
    
	'''
	print('starting building index tree ...')
	idxTree=getIdxTree(depthLen)
	print('get all possible path starts ...')
	allPath=getAllPath(idxTree)
	'''
	
	print('get all possible path using iter starts ...')
	allPath=rtnAllPathIter(depthLen)
	
	'''
	below allPath is transformed into:
	[[(1, 2), (2, 3)]]
	'''
	print('transform all possible paths')
	allPathTup=[[(allPath[i][j],allPath[i][j+1]) for j in range(len(allPath[i])-1)] for i in range(len(allPath))]
	print('generating value for all possible paths')
	pathValueLst=valueForAllPath(pcfMat,allPathTup,gamma)
	'''
	note that we've applied a positive sign to the target equation, thus max should be applied instead of min
	'''
	print('returns maximum value index starts ...')
	maxIndxLst=rtrMaxPathMat(pathValueLst)
	maxValPath=[allPath[i] for i in maxIndxLst]
	'''
	use a plot for plotting possible CNV region
	'''
	fileRawName=dpthFile.rstrip('dpth')
	print('data loading to file starts ...')
	writeDpthIndx(depthObjList,maxValPath,fileRawName+'depth.df',fileRawName+'.indx')
	return
	
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


#==============================================================================
# def getAllPath(idxTree):
# 	pathLst=[]
# 	queue=[deepcopy(idxTree)]
# 	tmpLst=[]
# 	while len(queue)>0:
# 		curNode=queue.pop()
# 		if curNode.son != []:
# 			queue.extend(curNode.son)
# 			tmpLst.append(curNode.rowIdx)
# 		else:
# 			candPath=deepcopy(tmpLst)
# 			candPath.append(curNode.rowIdx)
# 			pathLst.append(candPath)
# 			if len(tmpLst)>1 and len(queue)>0:
# 				while tmpLst[-1] > queue[-1].rowIdx:
# 					tmpLst.pop()
# 	return pathLst
# 	
# def getIdxTree(depthLen):
# 	idxLst=sorted(range(depthLen),reverse=True)
# 	nodeLst=[]
# 	for i in idxLst:
# 		if nodeLst == []:
# 			nodeLst.append(Node(i))
# 		else:
# 			sonLst=deepcopy(nodeLst)
# #			sonLst=sorted(sonLst,key=lambda x:x.rowIdx)
# 			sonLst=sorted(sonLst,key=attrgetter('rowIdx'))
# 			newNode=Node(i,sonLst)
# 			nodeLst.append(newNode)
# 	return nodeLst[-1]
# 	
# class Node:
# 	'''
# 	use belinked and linked to construct a directed graph
# 	'''
# 	def __init__(self,rowIdx,son=[]):
# 		self.rowIdx=rowIdx
# 		self.son=son
#==============================================================================
		
	
	
class PcfMat:
	'''
	note that when using index in PcfMat, this follows that initialize index is 0
	also note that in PcfMat, all values are positive, thus when applying function, signs should be reversed
	'''
	def __init__(self,deptObjList):
		self.deptObjLst=deptObjList
		self.pcfmat=[]
		self.__generatePcfMat()
	def __generatePcfMat(self):
		for i in range(len(self.deptObjLst)-1):
			tmpList=[]
			for j in range(i+1,len(self.deptObjLst)):
				
				tmpObjList=self.deptObjLst[i:(j+1)]
				ni=j+1-i
				valueList=[round(k.meanDpt**2,3) for k in tmpObjList]
				value=round(sum(valueList)/ni,3)
				tmpList.append(value)
			self.pcfmat.append(tmpList)
	def __getitem__(self,idx):
		if isinstance(idx,tuple) and idx[0] < idx[1]:
			return self.pcfmat[idx[0]][idx[1]-1]
	
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

