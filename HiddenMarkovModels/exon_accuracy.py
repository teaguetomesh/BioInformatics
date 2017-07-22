''' Compare the output of predict_exons.py with the TRUE values

As part of course work for CS 576: Introduction to Bioinformatics
11/22/2016
'''


import sys

__author__ = "Teague Tomesh"

#Load in files and store the sequences in lists
trueExonsFile = open(sys.argv[1], 'r')
predictExonsFile = open(sys.argv[2], 'r')

trueExonSeqs = []
predictExonSeqs = []

while True:
	newSeq = trueExonsFile.readline()
	if not newSeq: break
	trueExonSeqs.append(newSeq[:-1])

while True:
	newSeq = predictExonsFile.readline()
	if not newSeq: break
	predictExonSeqs.append(newSeq[:-1])
	
numSeqs = len(trueExonSeqs)

#Count and keep track of all the required parameters to be used in the equations
numCPrPos = 0
numCPrExonPos = 0
numPos = 0
numTrueExonPos = 0
numPredExonPos = 0
#Crawl through each character of each sequence and check if it is a match, if its an exon, etc
for i in range(numSeqs):
	currTrueSeq = trueExonSeqs[i]
	currPredSeq = predictExonSeqs[i]
	for j in range(len(currTrueSeq)):
		numPos += 1
		if(currTrueSeq[j].isupper()):
			numTrueExonPos += 1
			if(currPredSeq[j] is currTrueSeq[j]):
				numCPrExonPos += 1
		if(currPredSeq[j].isupper()):
			numPredExonPos += 1
		if(currPredSeq[j] is currTrueSeq[j]):
			numCPrPos += 1
			
#Calculate and print the results of the equations
print(format(numCPrPos/numPos, '.3f'))
print(format(numCPrExonPos/numTrueExonPos, '.3f'))
print(format(numCPrExonPos/numPredExonPos, '.3f'))