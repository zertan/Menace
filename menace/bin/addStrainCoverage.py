#!/usr/bin/env python
import numpy as np
import pandas
import xmltodict
import os
import sys
#import time

def chunks(l, n):
	"""Yield successive n-sized chunks from l."""
	for i in range(0, len(l), n):
		yield l[i:i+n]

def binData(x,binSize):
	l=np.ceil((len(x)/binSize))
	out=np.zeros(l+1)
	tmp=chunks(x,binSize)
	for i,val in enumerate(tmp):
		out[i]=np.sum(val)
	return out

def interp(data,length):
	calls=0
	while len(data)!=length:
		calls+=1
		length1=len(data)
		cov=np.sum(data)/float(length1)
		diff=length-length1
		binSize=int(np.ceil(length1/float(diff)))
		if binSize==0 or diff>length1/float(2.3):
			break
		

		chunks1=chunks(data,int(binSize))
		#data2=data


		#print("nr " + str(diff))

		vals=[]
		insert_ind=[]
		for i,val in enumerate(chunks1):
			binSize2=min(len(val),binSize)
			insert_ind.append(i*binSize+np.round(binSize2/float(2)))
			vals.append(np.mean(val))

		vals=np.array(vals)
		insert_ind=np.array(insert_ind)
		#print(str(data.shape))
		#print(str(insert_ind.shape))
		#print(str(vals.shape))
		#print(np.array_str(vals))
		data=np.insert(data,insert_ind,vals)
		#print(str(data.shape))
		cov2=np.sum(data)/length
		data=cov/cov2*data
		#print("k")
		if calls>=10:
			return []
	#print(str(calls))
	return data

originFrame=pandas.read_csv(os.path.join(sys.argv[1],'extra','origins.txt'),delimiter=" ",index_col=0)
#print(originFrame.to_string())
#originFrame=originFrame.transpose()
#print(originFrame.to_string())

#print(sys.argv[3])
#print(originFrame.index.values)
#print(originFrame.loc[sys.argv[2]]['Origin'])

#for i,arg in enumerate(sys.argv[2:]):
#	print(originFrame.loc[:arg])

data=[]
origin=[]
genomeLen=[]
bacteriaName=[]
ind=0

for i,arg in enumerate(sys.argv[3:]):
	#print(arg+".depth")
	try:
		tmpData=np.array(pandas.read_csv(arg+".depth",delimiter=" ")).astype('float')
		data.append(tmpData)
		#ind=ind+1

		with open(os.path.join(sys.argv[2],'Headers',arg+".xml")) as fd:
			obj = xmltodict.parse(fd.read())

		#genomeLen.append(int(obj['DocSum']['Item'][8]['#text']))
		bacteriaName.append(obj['DocSum']['Item'][1]['#text'])
		
		genomeLen.append(int(len(tmpData)))

		origin.append(int(originFrame.loc[arg]['Origin']))
		
		data[ind]=np.roll(data[ind],-origin[ind])
		ind=ind+1
	except IOError as e:
		print(arg+" depth not found.")
	except ValueError as e:
		print(arg+" origin not found.")
	
if not data or len(data)==1:
	sys.exit()

#print("Len data: "+str(len(data)))

#ind=np.where(np.max(genomeLen))
#print(str(ind))
#print(np.array_str(ind))
#print(str(genomeLen))
ind=genomeLen.index(max(genomeLen))
genomeLen=max(genomeLen)#[int(ind[0])]
#print(str(genomeLen))

print("ind "+str(ind))# + "len "+str(len(genomeLen)))
#print(str(len(data[ind])))

dataLen=len(data[ind])
data2=np.array(data[ind]).flatten()
del data[ind]

#print("init shape")
bap=data2.shape
#print(str(genomeLen))

for i,val in enumerate(data):
	tmp=interp(val,dataLen)
	
	print(sys.argv[3+i]+" tmp shape: "+ str(tmp.shape) + " init shape: "+str(bap))

	if tmp.shape==data2.shape:
		data2=np.add(data2,tmp)

print(sys.argv[3])
np.save(sys.argv[3], data2) 
