''' Use HMMs to predict where the Exons and Introns are in a sequence of DNA 

As part of course work for CS 576: Introduction to Bioinformatics
11/22/2016
'''

import sys
import math

__author__ = "Teague Tomesh"

#Load in the given files
trainingFile = open(sys.argv[1], 'r')
predictionFile = open(sys.argv[2], 'r')

trainSeqs = []
predictSeqs = []

#Put the sequences in the files into different lists
while True:
	newSeq = trainingFile.readline()
	if not newSeq: break
	trainSeqs.append(newSeq[:-1])

while True:
	newSeq = predictionFile.readline()
	if not newSeq: break
	predictSeqs.append(newSeq[:-1])
	
#This class will represent a state in the HMM
#It will keep track of which state it is, its emission and transition probabilites
# and it will also keep track of incoming and outgoing edges
class State:
	def __init__(self, label):
		self.stateLabel = label
		if(label is 0 or label is 3):
			self.silentState = True
		else:
			self.silentState = False
			
		self.emitProb = [0,0,0,0]
		
		if(label is 3):
			self.incomingStates = [1]
			self.outgoingStates = []
		elif(label is 2):
			self.incomingStates = [1,2]
			self.outgoingStates = [1,2]
		elif(label is 1):
			self.incomingStates = [0,1,2]
			self.outgoingStates = [1,2,3]
		else: 
			self.incomingStates = []
			self.outgoingStates = [1]
	
	def getLabel(self):
		return self.stateLabel
		
	def emission(self, char):
		char = char.upper()
		if char == 'A': return self.emitProb[0]
		elif char == 'C': return self.emitProb[1]
		elif char == 'G': return self.emitProb[2]
		elif char == 'T': return self.emitProb[3]
		
	def setEmitProb(self,char,n,n_tot):
		if char is 'A':
			self.emitProb[0] = float(n+1)/float(n_tot+4)
		if char is 'C':
			self.emitProb[1] = float(n+1)/float(n_tot+4)
		if char is 'G':
			self.emitProb[2] = float(n+1)/float(n_tot+4)
		if char is 'T':
			self.emitProb[3] = float(n+1)/float(n_tot+4)

#Initializing the transition matrix and the HMM which is stored in a list			
width = 4
height = 4
transitionProb = [[0 for x in range(width)] for y in range(height)]

hmmModel = [State(0), State(1), State(2), State(3)]

#Maximum Likelihood Estimation
#Attain the estimations for the parameters using Laplace smoothing
intronA, intronC, intronG, intronT = 0,0,0,0
exonA, exonC, exonG, exonT = 0,0,0,0
tran11, tran12, tran22, tran21, tran13 = 0,0,0,0,0
exonEmit, intronEmit = 0,0
tranOutIntron, tranOutExon = 0,0
for j in trainSeqs:
	pos = 0
	for i in j:
		if(i.isupper()):
			exonEmit += 1
			if i is 'A': exonA += 1
			elif i is 'C': exonC += 1
			elif i is 'G': exonG += 1
			else: exonT += 1
			if(pos is not 0):
				if(j[pos-1].isupper()):
					tran11 += 1
					tranOutExon += 1
				else:
					tran21 += 1
					tranOutIntron += 1
			if(pos is len(j)-1):
				tran13 += 1
				tranOutExon += 1
		else:
			intronEmit += 1
			if i is 'a': intronA += 1
			elif i is 'c': intronC += 1
			elif i is 'g': intronG += 1
			else: intronT += 1
			if(pos is not 0):
				if(j[pos-1].isupper()):
					tran12 += 1
					tranOutExon += 1
				else:
					tran22 += 1
					tranOutIntron += 1
		pos += 1

#Set the emission probabilities for the 2 states
hmmModel[1].setEmitProb('A',exonA,exonEmit)
hmmModel[1].setEmitProb('C',exonC,exonEmit)
hmmModel[1].setEmitProb('G',exonG,exonEmit)
hmmModel[1].setEmitProb('T',exonT,exonEmit)

hmmModel[2].setEmitProb('A',intronA,intronEmit)
hmmModel[2].setEmitProb('C',intronC,intronEmit)
hmmModel[2].setEmitProb('G',intronG,intronEmit)
hmmModel[2].setEmitProb('T',intronT,intronEmit)
	
#Set the transition probabilities for the 2 states
transitionProb[0][1] = 1
transitionProb[1][1] = float(tran11+1)/float(tranOutExon+3)
transitionProb[1][2] = float(tran12+1)/float(tranOutExon+3)
transitionProb[1][3] = float(tran13+1)/float(tranOutExon+3)
transitionProb[2][1] = float(tran21+1)/float(tranOutIntron+2)
transitionProb[2][2] = float(tran22+1)/float(tranOutIntron+2)

seqsToFile = []

#Viterbi Algorithm
#This function will compute the natural log of the given number
def ln(val):
	if val == 0:
		return -1*float("inf")
	else:
		return math.log(val)

#This function will return the index that gave the max value of the given list
def argmax(vals):
	#tempMax = -1*math.inf
	tempMax = -1*float("inf")
	#tempMaxIndex = -1*math.inf
	tempMaxIndex = -1*float("inf")
	for index in range(len(vals)):
		if(vals[index]>tempMax):
			tempMax = vals[index]
			tempMaxIndex = index
	for index in range(len(vals)):
		if index is not tempMaxIndex:
			if tempMax == vals[index]:
				return 2
	return tempMaxIndex

#For each of the sequences, run the Viterbi algorithm
for j in predictSeqs:
	#Initialize matrices
	numCols = len(j)+1
	numRows = 4
	vMatrix = [[ln(0) for col in range(numCols)] for row in range(numRows)]
	vMatrix[0][0] = ln(1)
	ptrMatrix = [[-1 for col in range(numCols)] for row in range(numRows)]
	#Matrix crawler recursion
	for i in range(numCols):
		for l in range(numRows):
			if l is not 0:
				if(hmmModel[l].silentState):
					#Update viterbi matrix and pointer matrix
					possibleVals = []
					for k in range(numRows):
						tranProb = ln(transitionProb[k][l])
						possibleVals.append(vMatrix[k][i] + tranProb)
					vMatrix[l][i] = max(possibleVals)
					ptrMatrix[l][i] = argmax(possibleVals)
				else:
					#Update viterbi matrix and pointer matrix
					emitVal = ln(hmmModel[l].emission(j[i-1]))
					possibleVals = []
					for k in range(numRows):
						tranProb = ln(transitionProb[k][l])
						possibleVals.append(vMatrix[k][i-1] + tranProb)
					vMatrix[l][i] = emitVal + max(possibleVals)
					ptrMatrix[l][i] = argmax(possibleVals)
	#Termination
	tempVals = []
	for k in range(numRows):
		tranProb = ln(transitionProb[k][numRows-1])
		tempVals.append(vMatrix[k][numCols-1] + tranProb)
	probability = max(tempVals)
	finalState = argmax(tempVals)
	
	#Traceback: follow pointers back starting at finalState
	bestPath = []
	n = list(range(0,numCols))
	bestPath.insert(0,finalState)
	nextState = ptrMatrix[finalState][numCols-1]
	for col in reversed(n):
		bestPath.insert(0,nextState)
		nextState = ptrMatrix[nextState][col-1]
	del bestPath[0]
	
	#Create the predicted sequence and print to standard output
	predictedSeq = ""
	for i in range(len(j)):
		curChar = j[i]
		if bestPath[i+1] is 2:
			predictedSeq += curChar.lower()
		else:
			predictedSeq += curChar
	print(predictedSeq)
	seqsToFile.append(predictedSeq)
	
outputFile = open('gene.out', 'w')
for seq in seqsToFile:
	outputFile.write(seq)
	outputFile.write('\n')
outputFile.close()

print("Exon:")
print(hmmModel[1].emission('A'), hmmModel[1].emission('C'), hmmModel[1].emission('G'), hmmModel[1].emission('T'),)
print(hmmModel[2].emission('A'), hmmModel[2].emission('C'), hmmModel[2].emission('G'), hmmModel[2].emission('T'),)












