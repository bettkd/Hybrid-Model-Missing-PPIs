import csv
import numpy as np
import threading
import time

class Enricher:

	inputFile = ""
	outputFile = ""
	ppis = [[]]
	ppisdict = dict()
	proteins = []
	proteinHubs = dict()
	branches = dict()

	discoveredEdges = set()

	DEPTH = 3
	MAX_SCORE = 1000.0
	ALPHA = 1.5
	confidencemean = 0
	confidencestddev = 0

	confidenceindex = 0
	branchcount = 0
	start_time = 0
	

	"""docstring for ProcessData"""
	def __init__(self, inputFile="Data/sample.txt", outputFile="Results/results.txt"):

		self.inputFile = inputFile
		self.outputFile = outputFile

	def execute(self):
		self.start_time = time.time()
		self.readPPIsFromFile()
		self.getUniqueProteins()
		self.getProteinHubs()
		self.getNetworkStats()
		self.getAllBranches()
		self.findMissingEdges()
		print "Total time to execute: ", time.time() - self.start_time, "seconds"

	def readPPIsFromFile(self):

		self.ppis = list(csv.reader(open(self.inputFile, 'rb'), delimiter='\t'))

		for ppi in self.ppis:

			self.ppisdict['\t'.join(sorted([ppi[0], ppi[1]]))] = int(ppi[2])

	def getUniqueProteins(self):

		proteins = set()

		for ppi in self.ppis:

			proteins.add(ppi[0])
			proteins.add(ppi[1])

		self.proteins =  list(sorted(proteins))

	def getProteinHubs(self):

		for protein in self.proteins:

			children = []
			for ppi in self.ppis:

				if protein in ppi:

					if protein == ppi[0]:
						children.append(ppi[1])
					else:
						children.append(ppi[0])

			self.proteinHubs[protein] = children

	def getChildren(self, tail, ancestors):

		children = list(self.proteinHubs[tail])

		try:
			children = [x for x in children if x not in ancestors]
		except Exception, e:
			pass

		return children


	def getInteractionBranches(self, protein, depth):

		branches = [[protein]]

		print "Working on branch for: ", protein

		for i in xrange(depth):

			_mybranches = []
			for branch in branches:

				tail = branch[-1]

				ancestors = branch[:-1]

				children = self.getChildren(tail, ancestors)

				_mybranch = []

				for child in children:

					n = list(branch)
					n.append(child)
					_mybranch.append(n)

				for item in _mybranch:

					_mybranches.append(item)

				self.branchcount += 1

				branches = _mybranches

		
		print 'Number of branches to process on this node = ', len(branches)

		'''
		t = threading.Thread(target=self.getMissingEdges, args=(protein, branches))
		t.start()
		t.join()
		'''

		self.branches[protein] = branches

	def getAllBranches(self):

		print "[", time.time() - self.start_time, "sec] Processing branches...."

		for protein in self.proteins:

			# multi-threading to do tasks concurrently hence faster
			t = threading.Thread(target=self.getInteractionBranches, args=(protein, self.DEPTH))
			t.start()
			t.join()
			#self.getInteractionBranches(protein, depth = self.DEPTH)

		print "Number of Branches Generated ", self.branchcount

		#for edge in self.discoveredEdges:
		#	print "Found: ", edge


	def getNetworkStats(self):

		scores = np.array([int(i[2])/self.MAX_SCORE for i in self.ppis])
		self.confidencemean = np.mean(scores)
		self.confidencestddev = np.std(scores)
		self.confidenceindex = self.ALPHA * self.confidencestddev + self.confidencemean



	def getPPIScore(self, protein1, protein2):

		score = 0
		try:
			score = self.ppisdict['\t'.join(sorted([protein1, protein2]))]
		except Exception, e:
			pass
		return score


	def getMissingEdges(self, root, branches):

		print "[", time.time() - self.start_time, "sec] Finding missing edges for: ", root
		for branch in branches:

				child = branch[1]
				descendants = branch[2:]
				confidence = self.getPPIScore(root, child) / self.MAX_SCORE


				for desc in descendants:

					descconfidence = self.getPPIScore(child, desc) / self.MAX_SCORE
					
					confidence = descconfidence * confidence

					if confidence >= self.confidenceindex:

						self.discoveredEdges.add('\t'.join(sorted([root, desc])))
					else:

						continue

					child = desc


	def removeDuplicates(self, rt, tmp):

		for brs in self.branches[rt]:

			if tmp.count(list(reversed(brs))) + tmp.count(brs) > 1:
				self.branches[rt].remove(brs)



	def findMissingEdges(self):

		################ Pre-processing ###############
		print "[", time.time() - self.start_time, "sec] PROCESSING...."
		tmp = []
		for rt in self.branches:
			for brs in self.branches[rt]:
				tmp.append(brs)

		print "Initial number of branches", len(tmp)
		'''
		print "Removing duplicate branches..."
		

		for rt in self.branches:
			# multi-threading to do tasks concurrently hence faster
			t = threading.Thread(target=self.removeDuplicates, args=(rt, tmp,))
			t.start()
			t.join()
		'''

		################ END Pre-processing ###############

		################ Finding missing edges ###############

		print ("Finding missing edges....")

		for root in self.branches:

			# multi-threading to do tasks concurrently hence faster
			t = threading.Thread(target=self.getMissingEdges, args=(root, self.branches[root]))
			t.start()
			t.join()
			#self.getMissingEdges(root)
			
		print ("Done finding mising edges....")

		for edge in self.discoveredEdges:
			if sorted(edge.split('\t')) in sorted(self.ppis):
				print edge, "----> in ppi ----"
				self.discoveredEdges.remove(edge)


		#for edge in self.discoveredEdges:
		#	print "Found: ", edge

		print "Total number found: ", len(self.discoveredEdges)

		#print self.discoveredEdges

		self.saveResults(self.outputFile)


	def saveResults(self, outputFile):

		#print len(self.discoveredEdges)

		print "Saving the results to ", outputFile

		with open(outputFile, 'w') as f:

			f.write( "Results for Asthma & Allergy\n\n" )
			f.write( "Score mean: %.2f\n" % self.confidencemean )
			f.write( "Score standard deviation: %.2f\n" % self.confidencestddev )
			f.write( "ALPHA (std devs from mean): %.2f\n" % self.ALPHA )
			f.write( "Confidence Index: %.2f\n" % self.confidenceindex )
			f.write( "Maximum graph depth: %d\n" % self.DEPTH )
			f.write( "Number of Proteins: %d\n" % len(self.proteins) )
			f.write( "Number of PPI: %d\n" % len(self.ppis) )
			f.write( "Number of interactions discovered: %d\n\n" % len(self.discoveredEdges) )
			for edge in sorted(self.discoveredEdges):
				f.write( edge )
				f.write( '\n' )
			f.close()

		print "DONE!!"

		################ END Finding missing edges ###############





inputFile = "Data/Allergy_and_Asthma.txt"
#inputFile = "Data/sample.txt"

#outputFile = "Results/results.txt"

#process = Enricher(inputFile=inputFile)
#process.execute()


########### Toogle comment to execute #######
#process.saveResults(outputFile)
