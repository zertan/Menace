#!/usr/bin/env python
'''
Created on Sep 8, 2015

@author: Daniel Hermansson
'''
import argparse
import re
from Bio import Entrez
from os import listdir,environ,makedirs
from os.path import isfile, join, exists
from sys import exit
from time import sleep
from math import ceil

def chunks(l, n):
	"""Yield successive n-sized chunks from l."""
	for i in range(0, len(l), n):
		yield l[i:i+n]

def search_entrez(inputString,args):
	idStr=[]
	matchNr=0;
	inputString2=chunks(inputString,3);
	strLen=int(float(len(inputString))/3+ceil(float(len(inputString)%3)/3))
	errIndex=[];

	print("Querying Entrez in chunks of 3.")
	for i, searchStr in enumerate(inputString2):
		for j,searchStr2 in enumerate(searchStr):
			if (re.match('^N[CTZ]_([A-Z]{2})*[0-9]{6}\.[0-9]+',searchStr2)):
				extra="[ACCN] AND srcdb_refseq_known[PROP]"
			else:
				extra='[PORG] AND "complete genome"[title] AND srcdb_refseq_known[PROP] NOT plasmid[title]'
				org_string=True

			searchHandle = Entrez.esearch(db="nucleotide",term=searchStr2+extra,retmax=1000)
			record = Entrez.read(searchHandle)
			try:
				if(extra=='[PORG] AND "complete genome"[title] AND srcdb_refseq_known[PROP] NOT plasmid[title]'):
					for rec in record["IdList"]:
						idStr.append(rec)
				else:
					idStr.append(record["IdList"][0])

				matchNr+=int(record["Count"])
			except IndexError:
				errIndex.append(i*3+j)

		#print(repr(strLen-i)+" ",end="")
		print((', '.join(searchStr)))
		if(strLen-i-1>0):
			sleep(1)

	if (matchNr==0):
		print("Found 0 matches. Exiting.")
		exit()
		#inputString=files

	print(("Found "+repr(matchNr)+" matche(s). Downloading to "+repr(args.dataPath)+"."))

	# delete non found items from inputString
	errorString = [i for j, i in enumerate(inputString) if j in errIndex]
	inputString = [i for j, i in enumerate(inputString) if j not in errIndex]

	if errorString:
		print("\nNon found strings:")
		print(('\n'.join(errorString)))

	idResults = Entrez.read(Entrez.epost("nuccore", id=",".join(idStr)))
	webenv = idResults["WebEnv"]
	queryKey = idResults["QueryKey"]
	return [webenv,queryKey,inputString]

def download_tax(webenv,queryKey,inputString,args):
	fetchHandle = Entrez.efetch(db="nuccore", query_key=queryKey,WebEnv=webenv,rettype="docsum",retmode="xml")
	data = fetchHandle.read()
	fetchHandle.close()

	####
	tIdArr=[]

	data=data.split('<DocSum>')
	data=data[1:]

	print(("Writing headers to "+str(args.dataPath)))
	for i, searchStr in enumerate(inputString):
		data[i]=re.sub('</eSummaryResult>', '', data[i])
		outHandle = open(join(args.dataPath,"Headers",searchStr+".xml"), "w")
		outHandle.write("<DocSum>"+data[i])
		outHandle.close()

		tmp=re.search('<Item Name="TaxId" Type="Integer">([0-9]+)',data[i])
		tIdArr.append(tmp.group(1))
	###

	print("Retrieving species taxonomic ids.")
	fetchHandle = Entrez.efetch(db='taxonomy',id=','.join(tIdArr),retmode='xml')
	records = Entrez.read(fetchHandle)
	fetchHandle.close()

	orgIdArr=[]
	orgNameArr=[]
	tmpIdArr=[]

	for i,record in enumerate(records):
		#print(record['ParentTaxId'])
		lineage = record['LineageEx'] # get the entire lineage of the first record
		#assert len(records) == 1 # die if we get more than one record, unlikely?

		tmpId=record['TaxId']#lineage[-1]['TaxId']
		tmpName=record['ScientificName']
		tmpIdArr.append(tmpId)
		for entry in lineage:
			if entry['Rank'] == 'species':
				#print(entry) # prints: {u'ScientificName': 'Viridiplantae', u'TaxId': '33090', u'Rank': 'kingdom '}
				tmpId=entry['TaxId']
				tmpName=entry['ScientificName']

		orgIdArr.append(tmpId)
		orgNameArr.append(tmpName)

	outHandle = open(join(args.dataPath,"taxIDs.txt"), "a")

	for i, searchStr in enumerate(inputString):
		outHandle.write(searchStr+"\t"+orgIdArr[i]+"\t"+orgNameArr[i]+"\n")
	outHandle.close()

def main():
	parser = argparse.ArgumentParser(description='Download nucleotide sequences from NCBI in fasta format.')
	parser.add_argument('-d', metavar="pathname", dest="dataPath", default="./", help="Path to directory containing fasta files. This will be the download directory. Already donloaded files in this directory are automatically left on the server. Default is './' (current directory).")
	parser.add_argument("-s", metavar="filepath", dest="searchFile", default="searchStrings", help="Path to file containing search strings to pass to Entrez. It should contain one search string per line. Default is './searchStrings'.")
	parser.add_argument("-e", metavar="adress", dest="email", default=environ.get('EMAIL', ""), help="Specify an email adress to use with Entrez. If not specified the environment variable EMAIL is used.")
	parser.add_argument("-m", metavar="N", dest="fetchNr",type=int, default=100, help="Max number of sequences to fetch into memory.")
	parser.add_argument("-t", metavar="bool",dest="taxBool",type=bool,default=False,help="Download taxonomy information.")
	args = parser.parse_args()

	# create directories
	if not exists(args.dataPath):
		makedirs(args.dataPath)
	if (not exists(join(args.dataPath,"Headers")) and args.taxBool):
		makedirs(join(args.dataPath,"Headers"))
	if not exists(join(args.dataPath,"Fasta")):
		makedirs(join(args.dataPath,"Fasta"))
	
	Entrez.email = args.email
	if(args.email==""):
		print("No email adress specified, exiting.")
		exit()
	
	# read input strings and reorganize based on presence in folder
	if (isfile(args.searchFile)):
		print(("Using search file: " + args.searchFile))
		with open(args.searchFile) as f:
			inputStringInit = f.read().splitlines()
		# basic checks
		inputString = [a for a in inputStringInit if a != ""]
		inputString=list(set(inputString))
		print(("Length of input: "+str(len(inputString))))
		files = [ f[0:-6] for f in listdir(join(args.dataPath,"Fasta")) if f.endswith(".fasta") ]
		#files_xml = [ f[0:-4] for f in listdir(join(args.dataPath,"Headers")) if f.endswith(".xml") ]
		#headers = [ f[0:-4] for f in listdir(join(args.dataPath,"Headers")) if f.endswith(".xml") ]
		inputString = [string for string in inputString if string not in files]
		#inputString2 = [string for string in inputString if string not in files]

		if (not len(inputString)==0):# and len(files)==0):
			# search Entrez
			idStr=[]
			matchNr=0;
			inputString2=chunks(inputString,3);
			strLen=int(float(len(inputString))/3+ceil(float(len(inputString)%3)/3))
			errIndex=[];

			print("Querying Entrez in chunks of 3.")
			for i, searchStr in enumerate(inputString2):
				for j,searchStr2 in enumerate(searchStr):
					if (re.match('^N[CTZ]_([A-Z]{2})*[0-9]{6}\.[0-9]+',searchStr2)):
						extra="[ACCN] AND srcdb_refseq_known[PROP]"
					else:
						extra='[PORG] AND "complete genome"[title] AND srcdb_refseq_known[PROP] NOT plasmid[title]'
						org_string=True

					searchHandle = Entrez.esearch(db="nucleotide",term=searchStr2+extra,retmax=1000)
					record = Entrez.read(searchHandle)
					try:
						if(extra=='[PORG] AND "complete genome"[title] AND srcdb_refseq_known[PROP] NOT plasmid[title]'):
							for rec in record["IdList"]:
								idStr.append(rec)
						else:
							idStr.append(record["IdList"][0])

						matchNr+=int(record["Count"])
					except IndexError:
						errIndex.append(i*3+j)

				#print(repr(strLen-i)+" ",end="")
				print((', '.join(searchStr)))
				if(strLen-i-1>0):
					sleep(1)

			if (not matchNr==0):
				#print("Found 0 matches. Exiting.")
				#exit()
				#inputString=files

				print(("Found "+repr(matchNr)+" matche(s). Downloading to "+repr(args.dataPath)+"."))

				# delete non found items from inputString
				errorString = [i for j, i in enumerate(inputString) if j in errIndex]
				inputString = [i for j, i in enumerate(inputString) if j not in errIndex]

				if errorString:
					print("\nNon found strings:")
					print(('\n'.join(errorString)))

				# update inputString to matched accessions if searched using organism name
				if extra=='[PORG] AND "complete genome"[title] AND srcdb_refseq_known[PROP] NOT plasmid[title]':
					fetchHandle = Entrez.efetch(db="nuccore", id=",".join(idStr),rettype="acc",retmode="text")
					inputString = fetchHandle.read()
					fetchHandle.close()
					inputString=inputString.splitlines()

				inputStringChunks=chunks(inputString,args.fetchNr)

				#print(len(inputStringChunks))

				# post matched IDs to Entrez
				idResults = Entrez.read(Entrez.epost("nuccore", id=",".join(idStr)))
				webenv = idResults["WebEnv"]
				queryKey = idResults["QueryKey"]

				# get sequences from Entrez in chunks of 100 (default)
				for ind, searchStrings in enumerate(inputStringChunks):
					print(("Fetching batch " +repr(ind+1)))
					fetchHandle = Entrez.efetch(db="nuccore", query_key=queryKey,WebEnv=webenv,rettype="fasta",retmode="text",retstart=repr(ind*args.fetchNr),retmax=repr(args.fetchNr))
					data = fetchHandle.read()
					fetchHandle.close()

					data=data.split('>')
					data=data[1:]

					print(("Writing batch " +repr(ind+1)))
					print((len(data)))
					for i, searchStr in enumerate(searchStrings):
						outHandle = open(join(args.dataPath,"Fasta",searchStr+".fasta"), "w")

						# limit fasta headers to first part of accession number for downstream tools
						nl_ind = data[i].find('\n')+1
						tmp_header = searchStr.split(" ")[0]+'\n'

						outHandle.write(">"+tmp_header+data[i][nl_ind:])
						outHandle.close()

				# get taxonomy info at strain and organism level
				if (args.taxBool):
					download_tax(webenv,queryKey,inputString,args)

	# check files present in Fasta folder and download corresponding Headers
	files_xml = [ f[0:-4] for f in listdir(join(args.dataPath,"Headers")) if f.endswith(".xml") ]
	files = [ f[0:-6] for f in listdir(join(args.dataPath,"Fasta")) if f.endswith(".fasta") ]
	inputString_xml = [string for string in files if string not in files_xml]

	[webenv,queryKey,inputString_xml]=search_entrez(inputString_xml,args)
	download_tax(webenv,queryKey,inputString_xml,args)


	files_xml = [ f[0:-4] for f in listdir(join(args.dataPath,"Headers")) if f.endswith(".xml") ]
	files = [ f[0:-6] for f in listdir(join(args.dataPath,"Fasta")) if f.endswith(".fasta") ]
	inputString_xml = [string for string in files if string not in files_xml]
	if inputString_xml:	
		print("*Note* There is a mismatch between Fasta and Header files. Any downstream analyses may be affected.")

if "__name__"=="__main__":
	sys.exit(main())
