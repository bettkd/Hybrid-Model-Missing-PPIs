import csv

data1 = 'Data/MCTResults.txt'
data2 = 'Data/DKTResults.txt'
resultFile = 'Results/truePositivesMissingPPIs.txt'

resultFile1 = 'Results/uniqueProteins.txt'


def readPPIsFromFile(data):
	ppi_set = set()
	ppis = list(csv.reader(open(data, 'rU'), dialect=csv.excel_tab))
	for ppi in ppis:
		ppi_set.add('\t'.join(sorted([ppi[0], ppi[1]])))
	return ppi_set

def saveResults(resultsFile, data):
	with open(resultsFile, 'w') as f:
		for ppi in sorted(data):
			f.write( ppi )
			#print ppi
			f.write( '\n' )
	f.close()

def intersect(a, b):
    return list(set(a) & set(b))

def getUniqueProteins(ppis):
	proteins = set()
	for ppi in ppis:
		ppi = ppi.split('\t')
		proteins.add(ppi[0])
		proteins.add(ppi[1])

	#print list(sorted(proteins))
	return list(sorted(proteins))

results1 = readPPIsFromFile(data1)
results2 = readPPIsFromFile(data2)
ppis =  intersect(results1, results2)
saveResults(resultFile, ppis)
proteins = getUniqueProteins(ppis)
saveResults(resultFile1, proteins)
print 'Done'

