#!/usr/bin/env python3
'''
Created on Sep 8, 2015

@author: Daniel Hermansson
'''
import argparse
import re
from Bio import Entrez
from os import listdir,environ,mkdir
from os.path import isfile, join, exists
from sys import exit
from time import sleep
from math import ceil

def chunks(l, n):
	"""Yield successive n-sized chunks from l."""
	for i in range(0, len(l), n):
		yield l[i:i+n]

parser = argparse.ArgumentParser(description='Download nucleotide sequences from NCBI in fasta format.')
parser.add_argument('-d', metavar="pathname", dest="dataPath", default="./", help="Path to directory containing fasta files. This will be the download directory. Already donloaded files in this directory are automatically left on the server. Default is './' (current directory).")
parser.add_argument("-s", metavar="filepath", dest="searchFile", default="searchStrings", help="Path to file containing search strings to pass to Entrez. It should contain one search string per line. Default is './searchStrings'.")
parser.add_argument("-e", metavar="adress", dest="email", default=environ.get('EMAIL', ""), help="Specify an email adress to use with Entrez. If not specified the environment variable EMAIL is used.")
parser.add_argument("-m", metavar="N", dest="fetchNr",type=int, default=40, help="Max number of sequences to fetch into memory.")
parser.add_argument("-t", metavar="bool",dest="taxBool",type=bool,default=False,help="Download taxonomy information.")
args = parser.parse_args()

# create directories
if not exists(args.dataPath):
	mkdir(args.dataPath)
if not exists(join(args.dataPath,"Headers")):
	mkdir(join(args.dataPath,"Headers"))
if not exists(join(args.dataPath,"Fasta")):
	mkdir(join(args.dataPath,"Fasta"))
	
if (isfile(args.searchFile)):
	with open(args.searchFile) as f:
		inputString = f.read().splitlines()
	
	# basic checks
	inputString = [a for a in inputString if a != ""]
	inputString=list(set(inputString))
	files = [ f[0:-6] for f in listdir(join(args.dataPath,"Fasta")) if f.endswith(".fasta") ]
	inputString = [string for string in inputString if string not in files]

	if (len(inputString)==0): 
		print("No search strings present, exiting.")
		exit()
else:
	print("No search strings present, exiting.")
	exit()

Entrez.email = args.email
if(args.email==""):
	print("No email adress specified, exiting.")
	exit()

# search Entrez
idStr=[]
matchNr=0;
inputString2=chunks(inputString,3);
strLen=int(float(len(inputString))/3+ceil(float(len(inputString)%3)/3))
errIndex=[];

print("Querying Entrez in chunks of 3.")
for i, searchStr in enumerate(inputString2):
	for j,searchStr2 in enumerate(searchStr):
		searchHandle = Entrez.esearch(db="nucleotide",term=searchStr2)
		record = Entrez.read(searchHandle)
		try:
			idStr.append(record["IdList"][0])
			matchNr+=int(record["Count"])
		except IndexError:
			errIndex.append(i*3+j)

	print(repr(strLen-i)+" ",end="")
	print(', '.join(searchStr))
	if(strLen-i-1>0):
		sleep(1)

if (matchNr==0):
	print("Found 0 matches. Exiting.")
	exit()
else: 
	print("Found "+repr(matchNr)+" matche(s). Downloading to "+repr(args.dataPath)+".")

# delete non found items from inputString
inputString = [i for j, i in enumerate(inputString) if j not in errIndex]
inputStringChunks=chunks(inputString,args.fetchNr)

# post matched IDs to Entrez
idResults = Entrez.read(Entrez.epost("nuccore", id=",".join(idStr)))
webenv = idResults["WebEnv"]
queryKey = idResults["QueryKey"]

# get sequences from Entrez in chunks of 40 (default)
for ind, searchStrings in enumerate(inputStringChunks):
	print("Fetching batch " +repr(ind+1)+ " to "+str(args.dataPath))
	fetchHandle = Entrez.efetch(db="nuccore", query_key=queryKey,WebEnv=webenv,rettype="fasta",retmode="text",retstart=repr(ind*args.fetchNr),retmax=repr(args.fetchNr))
	data = fetchHandle.read()
	fetchHandle.close()

	data=data.split('>')
	data=data[1:]

	print("Writing batch " +repr(ind+1)+ " to "+str(args.dataPath))
	for i, searchStr in enumerate(searchStrings):
		outHandle = open(join(args.dataPath,"Fasta",searchStr+".fasta"), "w")
		outHandle.write(">"+data[i])
		outHandle.close()

# get taxonomy info at strain and organism level
if (args.taxBool==True):
	fetchHandle = Entrez.efetch(db="nuccore", query_key=queryKey,WebEnv=webenv,rettype="docsum",retmode="xml")
	data = fetchHandle.read()
	fetchHandle.close()
	
	tIdArr=[]
	
	data=data.split('<DocSum>')
	data=data[1:]

	print("Writing headers " +repr(ind+1)+ " to "+str(args.dataPath))
	for i, searchStr in enumerate(searchStrings):
		outHandle = open(join(args.dataPath,"Headers",searchStr+".xml"), "w")
		outHandle.write("<DocSum>"+data[i])
		outHandle.close()
		
		tmp=re.search('<Item Name="TaxId" Type="Integer">([0-9]+)',data[i])
		tIdArr.append(tmp.group(1))

	print("Retrieving species taxonomic ids.")
	fetchHandle = Entrez.efetch(db='taxonomy',id=",".join(tIdArr),retmode='xml')
	records = Entrez.read(fetchHandle)
	fetchHandle.close()

	orgIdArr=[]

	for i,record in enumerate(records):
		lineage = record['LineageEx'] # get the entire lineage of the first record
		#assert len(records) == 1 # die if we get more than one record, unlikely?
		for entry in lineage:
			if entry['Rank'] == 'species':
				#print(entry) # prints: {u'ScientificName': 'Viridiplantae', u'TaxId': '33090', u'Rank': 'kingdom '}
				orgIdArr.append(entry['TaxId'])
	
	outHandle = open(join(args.dataPath,"taxIDs.txt"), "w")

	#outHandle.write("ACC"+"\tOrg"+"\tStr"+"\n")
	for i, searchStr in enumerate(searchStrings):	
		outHandle.write(searchStr+"\t"+orgIdArr[i]+"\t"+tIdArr[i]+"\n")	
				
	outHandle.close()

