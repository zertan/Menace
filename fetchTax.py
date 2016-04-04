from Bio import Entrez
Entrez.email = "daniel@murrays.nu" # we're nice and tell our address
#handle = Entrez.esearch(db='taxonomy', term='NC_006643.1', retmax=2, usehistory='y') # return 2 entries at the most - I expect only 1
#search_results = Entrez.read(handle)
#handle.close()
#count = int(search_results['Count'])
#print('Got %s records' %count) # should be 1

# get the cookies from the search
#webenv = search_results["WebEnv"]
#query_key = search_results["QueryKey"]

# download results from the database using our cookies from above

fetch_handle = Entrez.efetch(db='taxonomy',id="NC_006643.1",retmode='xml')
#, webenv=webenv, query_key=query_key)
records = Entrez.read(fetch_handle)
#assert len(records) == 1 # die if we get more than one record, unlikely?
lineage = records[0]['LineageEx'] # get the entire lineage of the first record
for entry in lineage:
	if entry['Rank'] == 'species':
		print(entry) # prints: {u'ScientificName': 'Viridiplantae', u'TaxId': '33090', u'Rank': 'kingdom '}
		print(entry['ScientificName']) # prints: 'Viridiplantae'
		print(entry['TaxId']) # prints: 'Viridiplantae'
fetch_handle.close()
