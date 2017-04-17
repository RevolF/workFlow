# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 15:17:16 2017

@author: thor
"""

# annotate hgvs to genomic positions

import os,re
os.chdir(r'D:\\thor\\documents\\WorkProjects\\cDNA')

cdnaBed='hg19_Homo_sapiens.GRCh37.63.DMD.enst7033.cds.with.cDNA.numbering.bed'
hgvsFile='HGMD_profession_processed_splicing.deep.intronic.txt'

def main():
	global cdnaBed,hgvsFile
	cdnaLst=getCdna(cdnaBed)
	annoHgvs(hgvsFile,cdnaLst)
	return

def getCdna(cdnaBed):
	cdnaLst=[]
	with open(cdnaBed,'r') as infh:
		for line in iter(infh):
			linear=line.strip().split('\t')
			cdnaLst.append([int(i) for i in linear[1:5]])
	return cdnaLst
	
def annoHgvs(hgvsFile,cdnaLst):
	
	with open(hgvsFile,'r') as infh, open(hgvsFile.rstrip('txt')+'genoPos.txt','w') as outfh:
		for line in iter(infh):
			hgvs=re.findall(r'(c\..*?\>\S)',line)[0]
			print hgvs
			
			cdnaPos, intronSlide = [int(i) for i in re.findall(r'(\d+)',hgvs)]
			orient = re.findall(r'\d+(.*?)\d+',hgvs)[0]
			
			for item in cdnaLst:
				if item[2] <= cdnaPos <= item[3]:
					gcdsPos=item[0]+(item[3]-cdnaPos)
					genoPos=0
					if orient == '+':
						genoPos = gcdsPos-intronSlide
					elif orient == '-':
						genoPos =gcdsPos+intronSlide
					chrxGpos='chrX:'+str(genoPos)
					outfh.write(chrxGpos+'\t'+line)
					break
	return
	
if __name__=='__main__':
	main()


