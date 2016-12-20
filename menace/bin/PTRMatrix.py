#!/usr/bin/env python
#import matplotlib
#matplotlib.use('Agg')

import numpy as np
#from scipy.optimize import minimize
# from scipy import signal
# from lmfit.models import Model
# from lmfit import conf_interval
import sys
import pandas as pd
# from pandas import DataFrame
#import matplotlib.pyplot as plt
import xmltodict

#from os import listdir
from os.path import join

#onlyfiles = [ f for f in listdir(sys.argv[1]) if isfile(join(mypath,f)) ]

from os import walk
from os import makedirs
import csv
from numpy.matlib import repmat
import re

#import pickle

def circularMedian(p,g):
	def tmFun():	
		pm=np.asmatrix(p)
		pLen=p.shape[0]
		m1=np.max(np.mod(repmat(pm.T,1,pLen)-repmat(pm,pLen,1),g),0)
		m2=np.min(np.mod(repmat(pm.T,1,pLen)-repmat(pm,pLen,1),g),0)
		return pm[0,np.argmin(m1-m2)]
	tm=tmFun()
	return np.mod(np.median(np.mod(p-tm,g))+tm,g)

def growthRate(p,g):
	pc=np.float64(np.log2(p))
	gc=np.float64(g)
	return 2*60*np.log(2)*pc/gc

def doublingTime(p,g):
	a=-4.602311
	b=1050.197
	return float(b)/(np.log2(p)/float(g)-a)

remove=["python" "-W" "ignore"]
for r in remove:
	try:
		sys.argv.remove("python")
	except:
		pass

if(len(sys.argv) <= 1):
	print("Usage: ./PTRMatrix.py DATA_PATH REF_PATH OUT_PATH DORIC_PATH LOC_FILE")
	sys.exit(1)

if(len(sys.argv) == 5):
	try:
		oriData=pd.DataFrame.from_csv(join(sys.argv[5],"bacteria_record.dat"),  sep='	',index_col=1)
	except:
		koremOriData=pd.DataFrame.from_csv(sys.argv[5],  sep=';',index_col=0)

if(len(sys.argv) == 6):
	oriData=pd.DataFrame.from_csv(join(sys.argv[5],"bacteria_record.dat"),  sep='   ',index_col=1)
	koremOriData=pd.DataFrame.from_csv(sys.argv[6],  sep=';',index_col=0)

#print(koremOriData.loc['NC_016845.1','OriC'])

#oTable=pd.DataFrame(columns=[ f for f in listdir(sys.argv[1]) if isdir(join(sys.argv[1],f)) ])
oTable=pd.DataFrame()
tTable=pd.DataFrame()
ptrTable=pd.DataFrame()
headerTable=pd.DataFrame()
covTable=pd.DataFrame()
#headerTable.colums=['ACC','Name','PTR','OriC','TerC','Length']

abundanceTable=pd.DataFrame()
#abundanceTable.index.name = 'Name'
#abundanceTable.set_index('Name',inplace=True)
#taxIDs = pickle.load( open( "taxIds.p", "rb" ) )
#bacteriaName = dict (zip(taxIDs.values(),taxIDs.keys()))

with open(join(sys.argv[2],'taxIDs.txt'), mode='r') as infile:
	reader = csv.reader(infile,delimiter='	')
	bacteriaName={}
	bacNameAcc={}
	for rows in reader:
#		print(str(rows))
		if (not rows[1]):
			with open(join(sys.argv[2],'Headers',rows[0]+".xml")) as fd:
				obj = xmltodict.parse(fd.read())
				#print(str(obj['eSummaryResult']['DocSum']['Item'][7]['#text']))
				bacteriaName[str(obj['DocSum']['Item'][7]['#text'])]=str(obj['DocSum']['Item'][1]['#text'])
				bacNameAcc[rows[0]]=str(obj['DocSum']['Item'][1]['#text'])
		else:
			bacteriaName[rows[1]]=rows[2]
			bacNameAcc[rows[0]]=rows[2]

#	tmpList= set([rows[2],rows[1] for rows in reader])
#	bacteriaName = {rows[0]:rows[1] for rows in tmpList}

#bacteriaName = dict (zip(taxIDs.values(),taxIDs.keys()))


print("Traversing directories and generating tables, please wait.")
for (dirpath, dirnames, filenames) in walk(sys.argv[1]):
	folderName=dirpath[-9:]
	for fn in filenames:
		#print(fn)
		#if (fn.endswith('coverage.csv')):
		#	tmpTable=pd.read_csv(join(sys.argv[1],dirpath,fn),delimiter=' ',index_col=0,usecols=[0,1])
		#	for ACC in tmpTable.index.values:
		#		covTable.loc[ACC,folderName]=tmpTable.loc[ACC]

		if (fn.endswith('report.tsv')):
			print(folderName)
			tmpTable=pd.read_csv(join(sys.argv[1],dirpath,fn),delimiter='	',skiprows=1,index_col=0,usecols=[0, 1])
			
			for ti in tmpTable.index.values:
				
				try:
					tmp_ti=re.match(r".*\|([0-9]+)\|.*",ti)
					tmp_ti=tmp_ti.group(1)
					#print(bacteriaName[ti[3:]])
					abundanceTable.loc[tmp_ti,'Name']=bacteriaName[tmp_ti]
					abundanceTable.loc[tmp_ti,folderName]=float(tmpTable.loc[ti])
				except KeyError as e:
					print('IndexError - "%s"' % str(e))
					
					
		if (fn.endswith('.depth.best.npy')):
			ACC=fn[:-15]
			try:
				#data=pd.read_csv(join(sys.argv[1],folderName,fn),delimiter=" ")
				vec=np.load(join(sys.argv[1],dirpath,fn))
				#dataVec=np.load(join(sys.argv[1],dirpath,ACC+'.depth.npy'))
				#print("load")
				try:
					tmp=headerTable.loc[ACC]
				except:
					with open(join(sys.argv[2],'Headers',ACC+".xml")) as fd:
						obj = xmltodict.parse(fd.read())

						genomeLen=int(obj['DocSum']['Item'][8]['#text'])
						bacteriaName2=obj['DocSum']['Item'][1]['#text']
					headerTable.loc[ACC,'Name']=bacteriaName2
					headerTable.loc[ACC,'Species Name']=bacNameAcc[ACC]
					headerTable.loc[ACC,'Length']=genomeLen
					#covTable.loc[ACC]=covTable.loc[ACC]/genomeLen
					try:
						oriLoc=oriData.loc[ACC,'oriC location']
						oriLoc=oriLoc.split('..',1)
						oriLoc=int(oriLoc[0])
#						print(oriData.loc[sys.argv[1][:-6],'Organism'])
						headerTable.loc[ACC,'DoriC']=oriLoc
						headerTable.loc[ACC,'DterC']=np.mod(oriLoc+genomeLen/float(2),genomeLen)
#						headerTable.loc[ACC,'DterC']=koremOriData.loc[ACC,'TerC']
					except:
						try:
							headerTable.loc[ACC,'DoriC']=koremOriData.loc[ACC,'OriC']
							headerTable.loc[ACC,'DterC']=koremOriData.loc[ACC,'TerC']
						except:
							pass
#							headerTable.loc[ACC,'DterC']=np.mod(oriData.loc[ACC,'oriC location']+headerTable.loc[ACC,'Length']/float(2),headerTable.loc[ACC,'Length'])
			except IOError as e:
				print("I/O error({0}): {1}".format(e.errno, e.strerror))
				break
			#print(ACC+","+folderName+": "+str(vec[0]))
            
			oTable.loc[ACC,folderName]=vec[0]
			tTable.loc[ACC,folderName]=vec[1]

print(oTable.values)


lenVec=headerTable['Length'].values
noTable=oTable.values
ntTable=tTable.values

#headerTable.colums=['ACC','Name','PTR','Length']
#headerTable=headerTable.set_index('ACC')

#print(headerTable.to_string)

#print(lenVec)

#oLocs=np.apply_along_axis( myfunction, axis=1, arr=noTable )

oLoc=np.zeros(noTable.shape[0])
tLoc=np.zeros(ntTable.shape[0])

for i in range(noTable.shape[0]):
	idx=np.invert(np.isnan(noTable[i]))
	if (np.sum(idx)>=int(sys.argv[3])):
		oLoc[i]=circularMedian(noTable[i][idx],lenVec[i])
		headerTable.ix[i,'OriC']=oLoc[i]

	idx=np.invert(np.isnan(ntTable[i]))
	if (np.sum(idx)>=int(sys.argv[3])):
		tLoc[i]=circularMedian(ntTable[i][idx],lenVec[i])
		headerTable.ix[i,'TerC']=tLoc[i]



#	ptrTable.loc[headerTable.ix[i,'ACC']]=

#headerTable=headerTable.isnan()
#index = np.invert(headerTable['OriC'].index[headerTable['OriC'].apply(np.isnan)])
#headerTable['OriC']=headerTable['OriC'].ix[index]



#headerTable = headerTable[headerTable['OriC'] != np.nan]
res=pd.notnull(headerTable['OriC'])
headerTable=headerTable[res]
#res2=headerTable['ACC'].loc[res]
#print(res2.to_string)
#print(headerTable.ix[res].to_string)
#sys.quit()
print(oTable.to_string(max_rows=1000))
#print(headerTable.to_string)
headerTable=headerTable.sort('Name')
print(headerTable.to_string(max_rows=1000))


#x0=np.linspace(0,10**6,100)
#y=np.zeros_like(x0)
#index=np.isnan(noTable[0])
#for i in range(0,100-1):
#	y[i]=circularMedian(noTable[0][index],5000000,x0[i])
#	print(y[i])

#print(oTable.to_string)
#ot2=oTable[0].values
#ot2i=np.isnan(ot2)

ptrTable=pd.DataFrame()
tauTable=pd.DataFrame()
abTable=pd.DataFrame()

for (dirpath, dirnames, filenames) in walk(sys.argv[1]):
	folderName=dirpath[-9:]
	for fn in filenames:
		if (fn.endswith('.depth.npy')):
			cont=True
			ACC=fn[:-10]
			try:
				dataVec=np.power(2,np.load(join(sys.argv[1],dirpath,fn)))
			except IOError as e:
				print("I/O error({0}): {1}".format(e.errno, e.strerror))
				cont=False	
			try:
				tmp=headerTable.loc[ACC,'TerC']
			except:
				cont=False
			if(cont):
				OriC=headerTable.loc[ACC,'OriC']
				TerC=headerTable.loc[ACC,'TerC']
				Glen=headerTable.loc[ACC,'Length']
				#print('ACC: '+ACC)
				#print(Glen)
				#print(OriC)
				#print(dataVec.shape[0])
				#print(folderName)
				#print(np.round(OriC/Glen*dataVec.shape[0]))
				#print(np.round(TerC/Glen*dataVec.shape[0]))
				PTR=dataVec[np.floor((OriC/Glen)*dataVec.shape[0])]/dataVec[np.floor((TerC/Glen)*dataVec.shape[0])]
				ptrTable.loc[ACC,'Name']=headerTable.loc[ACC,'Species Name']
				ptrTable.loc[ACC,folderName]=PTR
				tauTable.loc[ACC,'Name']=headerTable.loc[ACC,'Species Name']
				tauTable.loc[ACC,folderName]=growthRate(PTR,Glen)
				#covTable.loc[ACC,'Name']=headerTable.loc[ACC,'Species Name']
#				tauTable.loc[ACC,folderName]=PTR


#ot2=oTable[oTable[[0]].apply(lambda x: x[0].isnan(), axis=0)]
#print(np.array_str(ot2[ot2i]))
#df = df[[]]

ptrTable=ptrTable.sort('Name')
ptrTable=ptrTable.sort_index(axis=1)
#ptrTable=tauTable.sort('Name')
#ptrTable=tauTable.sort_index(axis=1)

#print(str(doublingTime(1.3,6000000)))

#print(ptrTable.columns.values)


#for i in ptrTable.index:
#abTable.columns=ptrTable.columns

#abTable=abundanceTable.loc[ptrTable['Name'].isin(ptrTable['Name'])]
#print(ptrTable.to_string)
####
#try:
#	abTable=abundanceTable[abundanceTable['Name'].isin(ptrTable['Name'])]
#	abTable=abTable.sort('Name')
#	abTable=abTable.sort_index(axis=1)
#	abTable.to_csv('AbundancePTR.csv',sep=";")
#except:
#	pass
###
#try:
#	covTable=covTable[covTable['Name'].isin(ptrTable['Name'])]
#	covTable=covTable.sort('Name')
#	covTable=covTable.sort_index(axis=1)
#	covTable.to_csv('AnalysisCov.csv',sep=";")
#except:
#	pass

#for i, row in ptrTable.iterrows():
#	abTable.ix[i]=abundanceTable.loc[ptrTable.ix[i,'Name']]


print("\n"+ptrTable.to_string(index=False))

tauTable=tauTable.sort('Name')
tauTable=tauTable.sort_index(axis=1)
#covTable=covTable.sort_index(axis=1)

out_path=join(sys.argv[4],'Collect')

try:
	makedirs(out_path)
except:
	pass

headerTable.to_csv(join(out_path,'Header.csv'),sep=";")
ptrTable.to_csv(join(out_path,'PTR.csv'),sep=";",index=False)
tauTable.to_csv(join(out_path,'DoublingTime.csv'),sep=";",index=False)
abundanceTable.to_csv(join(out_path,'Abundance.csv'),sep=";",index=False)

print(" ".join(["\nOutput stored in"]+[join(sys.argv[4],'Collect')]))
#plt.plot(y, 'r-')
#plt.savefig("asd.png")

