#!/usr/bin/env python
#import numpy as np
import sys
import pandas as pd
import xmltodict
#from os import listdir
from os.path import join
from os import walk
import csv
#from numpy import matlib
#import pickle

####
# XML header example:
# 0      <Item Name="Caption" Type="String">NC_012207</Item>
# 1      <Item Name="Title" Type="String">Mycobacterium bovis BCG str. Tokyo 172 DNA, complete genome</Item>
# 2      <Item Name="Extra" Type="String">gi|224988383|ref|NC_012207.1|[224988383]</Item>
# 3      <Item Name="Gi" Type="Integer">224988383</Item>
# 4      <Item Name="CreateDate" Type="String">2009/03/13</Item>
# 5      <Item Name="UpdateDate" Type="String">2015/07/30</Item>
# 6      <Item Name="Flags" Type="Integer">800</Item>
# 7      <Item Name="TaxId" Type="Integer">561275</Item>
# 8      <Item Name="Length" Type="Integer">4371711</Item>
# 9      <Item Name="Status" Type="String">live</Item>
# 10     <Item Name="ReplacedBy" Type="String"></Item>
# 11     <Item Name="Comment" Type="String"><![CDATA[  ]]></Item>
####

with open(join(sys.argv[2],'taxIDs.txt'), mode='r') as infile:
	reader = csv.reader(infile,delimiter=';')
	bacteriaName={}
	bacNameAcc={}
	tiAcc={}
	genLenAcc={}
	for rows in reader:
		with open(join(sys.argv[2],rows[0]+".xml")) as fd:
			obj = xmltodict.parse(fd.read())
			genLenAcc[rows[0]]=str(obj['eSummaryResult']['DocSum']['Item'][8]['#text'])

		if (not rows[1]):
			with open(join(sys.argv[2],rows[0]+".xml")) as fd:
				obj = xmltodict.parse(fd.read())	
				bacteriaName[str(obj['eSummaryResult']['DocSum']['Item'][7]['#text'])]=str(obj['eSummaryResult']['DocSum']['Item'][1]['#text'])
				bacNameAcc[rows[0]]=str(obj['eSummaryResult']['DocSum']['Item'][1]['#text'])
				tiAcc[rows[0]]=str(obj['eSummaryResult']['DocSum']['Item'][7]['#text'])
		else:
			bacteriaName[rows[1]]=rows[2]
			bacNameAcc[rows[0]]=rows[2]
			tiAcc[rows[0]]=rows[1]


print("Traversing directories and generating coverage table, please wait.")

coverageTable=pd.DataFrame()

for (dirpath, dirnames, filenames) in walk(sys.argv[1]):
	folderName=dirpath[-9:]
	for fn in filenames:
		if (fn.endswith('coverage.csv')):
			print(folderName)
			tmpTable=pd.read_csv(join(dirpath,fn),delimiter=' ',index_col=0,usecols=[0, 1])
			for ACC in tmpTable.index.values:
				try:
					if (float(tmpTable.loc[ACC])/float(genLenAcc[ACC])>0.001):
						try:
							coverageTable.loc[ACC,'TaxID']
						except:
							coverageTable.loc[ACC,'TaxID']=tiAcc[ACC]
						coverageTable.loc[ACC,folderName]=float(tmpTable.loc[ACC])/float(genLenAcc[ACC])
				except KeyError, e:
					print 'IndexError - "%s"' % str(e)

maxInd=[];
resACC=[];
fo = open("referenceACC.txt", "wb")
for tid in list(set(coverageTable['TaxID'].values)):
	tmpTable=coverageTable.loc[coverageTable['TaxID'] == tid]
	tmpMean=list(tmpTable.mean(axis=1,numeric_only=True))
	maxInd=tmpMean.index(max(tmpMean))
	tmpACClist=tmpTable.index.values;
	print>>fo, tmpACClist[maxInd]
	print(str(bacNameAcc[tmpACClist[maxInd]]))

fo.close()
