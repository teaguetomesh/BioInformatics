import sys
import csv
import math
import numpy as np

dataFile = open(sys.argv[1], 'r')
clusterType = sys.argv[2]
numReturnCluster = int(sys.argv[3])

trueExonSeqs = []
with dataFile as f:
	reader = csv.reader(f, delimiter="\t")
	data = list(reader)

class Gene:
	def __init__(self, indicator, descrip, pVector):
		self.indicator = indicator
		self.descrip = descrip
		self.pVector = pVector

def initial_dist(c_u, c_v):
	sumSquares = 0
	for j in range(len(c_u)):
		sumSquares = sumSquares + (c_u[j]-c_v[j])**2
	return math.sqrt(sumSquares)
	
def updateDist(u, v, uLen, vLen):
	#print('distMat dimensions:',len(distMat),len(distMat[0]))
	#print(np.matrix(distMat))
	#print('u:',u,'v:',v)
	newMatrix = []
	for i in range(len(distMat)):
		newRow = []
		for k in range(len(distMat)):
			#print((i!=u and i!=v) and (k!=u and k!=v))
			if((i!=u and i!=v) and (k!=u and k!=v)):
				#print("adding to new row")
				newRow.append(distMat[i][k])
				#print("len new row inside:",len(newRow))
		if(len(newRow) != 0): 
			#print("inside adding to newMatrix")
			newRow.append(-1)
			newMatrix.append(newRow)		
	lastRow = [-1 for i in range(len(newMatrix[0]))]
	#print(len(lastRow),len(newRow))
	lastRow[len(lastRow)-1] = 0
	newMatrix.append(lastRow)
	#print(lastRow)
	#print(newMatrix)
	k = len(newMatrix)-1
	if(clusterType == 'S'):
		#print(len(newMatrix),len(C))
		myAdder = 0
		i=0
		#alreadyAdded = False
		while((i+myAdder)<len(newMatrix)+1):
			if((i+myAdder)==u or (i+myAdder)==v):
				nextIndex = i+myAdder+1
				if((nextIndex==u or nextIndex==v) and nextIndex==len(distMat)-1):
					break
				if((nextIndex==u or nextIndex==v) and nextIndex!=len(distMat)-1):
					myAdder = myAdder+2
				#print("in adder")
				if(i+myAdder == len(distMat)-1):
					myAdder = myAdder
				else: myAdder = myAdder + 1	
			#print("curr index:",i+myAdder)
			dist1 = distMat[i+myAdder][u]
			dist2 = distMat[i+myAdder][v]
			#print('[',i+myAdder,u,"]dist1:",dist1,"\t[",i+myAdder,v,"]dist2:",dist2)
			if(dist1 < dist2): 
				newMatrix[i][k] = dist1
				newMatrix[k][i] = dist1
			else: 
				newMatrix[i][k] = dist2
				newMatrix[k][i] = dist2
			#print("distance added:",newMatrix[i][k])
			i = i + 1
	elif(clusterType == 'C'):
		#print(len(newMatrix),len(C))
		myAdder = 0
		i=0
		#alreadyAdded = False
		while((i+myAdder)<len(newMatrix)+1):
			if((i+myAdder)==u or (i+myAdder)==v):
				nextIndex = i+myAdder+1
				if((nextIndex==u or nextIndex==v) and nextIndex==len(distMat)-1):
					break
				if((nextIndex==u or nextIndex==v) and nextIndex!=len(distMat)-1):
					myAdder = myAdder+2
				#print("in adder")
				if(i+myAdder == len(distMat)-1):
					myAdder = myAdder
				else: myAdder = myAdder + 1	
			#print("curr index:",i+myAdder)
			dist1 = distMat[i+myAdder][u]
			dist2 = distMat[i+myAdder][v]
			#print('[',i+myAdder,u,"]dist1:",dist1,"\t[",i+myAdder,v,"]dist2:",dist2)
			if(dist1 < dist2): 
				newMatrix[i][k] = dist2
				newMatrix[k][i] = dist2
			else: 
				newMatrix[i][k] = dist1
				newMatrix[k][i] = dist1
			#print("distance added:",newMatrix[i][k])
			i = i + 1
	else:
		#print(len(newMatrix),len(C))
		myAdder = 0
		i=0
		#alreadyAdded = False
		while((i+myAdder)<len(newMatrix)+1):
			if((i+myAdder)==u or (i+myAdder)==v):
				nextIndex = i+myAdder+1
				if((nextIndex==u or nextIndex==v) and nextIndex==len(distMat)-1):
					break
				if((nextIndex==u or nextIndex==v) and nextIndex!=len(distMat)-1):
					myAdder = myAdder+2
				#print("in adder")
				if(i+myAdder == len(distMat)-1):
					myAdder = myAdder
				else: myAdder = myAdder + 1	
			#print("curr index:",i+myAdder)
			dist1 = distMat[i+myAdder][u]
			dist2 = distMat[i+myAdder][v]
			#print('[',i+myAdder,u,"]dist1:",dist1,"\t[",i+myAdder,v,"]dist2:",dist2)
			newMatrix[i][k] = float(uLen*dist1 + vLen*dist2)/float(uLen+vLen)
			newMatrix[k][i] = float(uLen*dist1 + vLen*dist2)/float(uLen+vLen)
			i = i + 1
	
	#print(np.matrix(newMatrix))
	#print("newMat dimensions:",len(newMatrix),len(newMatrix[0]))
	return newMatrix

def argmin():
	currMin = 100000000
	for u in range(len(distMat)):
		for v in range(len(distMat)):
			if(u < v):
				tempVal = distMat[u][v]
				if(tempVal < currMin):
					currMin = tempVal
					#print("currMin:", currMin)
					minU = u
					minV = v
	return minU,minV

allGenes = []
setX = []
for i in range(len(data)):
	vecData = [float(j) for j in data[i][2:]]
	allGenes.append(Gene(data[i][0], data[i][1], vecData))
	setX.append(vecData)	
	
C = []
for i in range(len(setX)):
	c_i = setX[i]
	C.append(c_i)
	
#initialize distance matrix
distMat = [[-1 for x in range(len(C))] for y in range(len(C))]
for i in range(len(C)):
	for k in range(len(C)):
		if(i == k): distMat[i][k] = 0
		else:
			distMat[i][k] = initial_dist(C[i], C[k])
#print(distMat)
	
j = len(setX)
while len(C)>numReturnCluster:
	#print('\n')
	j = j+1
	#Find the least distant pair in C
	[minU, minV] = argmin()
	#print("minU:",minU)
	#print("minV:",minV)
	#Create a new cluster for pair
	#print("len U:", len(C[minU]), "len V:", len(C[minV]))
	c_j = []
	clust_u = C[minU]
	clust_v = C[minV]
	if(isinstance(clust_u[0],list)):
		for g in range(len(clust_u)):
			c_j.append(clust_u[g])
	else:
		c_j.append(clust_u)
	if(isinstance(clust_v[0],list)):
		for g in range(len(clust_v)):
			c_j.append(clust_v[g])
	else:
		c_j.append(clust_v)
	#print(clust_u)
	#print('\n',clust_v)
	#print('\n',c_j)
	#print("Length before merge:",len(C))
	#Add new cluster to list of clusters to be joined in the tree
	if(isinstance(C[minU][0],list)): uLen = len(C[minU])
	else: uLen = 1
	if(isinstance(C[minV][0],list)): vLen = len(C[minV])
	else: vLen = 1
	if(minU < minV):
		del C[minU]
		del C[minV-1]
	else:
		del C[minV]
		del C[minU-1]
	C.append(c_j)
	distMat = updateDist(minU, minV, uLen, vLen)
	#print("length after:",len(C))

#print("num clusters:", len(C))
#print(len(C[0]),len(C[1]))

geneClusters = []
for cluster in C:
	currCluster = []
	if(isinstance(cluster[0],list)):
		for pVec in cluster:
			for gene in allGenes:
				if(gene.pVector == pVec):
					currCluster.append(gene)
	else:
		for gene in allGenes:
			if(gene.pVector == cluster):
				currCluster.append(gene)
	geneClusters.append(currCluster)

#throwAwayClusters = copy.deepcopy(geneClusters)
OrderedgeneClusters = []
allAvgs = []
clusDict = {}
for val in range(len(geneClusters)):
	cluster = geneClusters[val]
	avgs = []
	for i in range(len(cluster)):
		curAvg = sum(cluster[i].pVector)/len(cluster[i].pVector)
		avgs.append(curAvg)
	clusAvg = sum(avgs)/len(avgs)
	clusDict[clusAvg] = val
	allAvgs.append(clusAvg)
	
allAvgs.sort()
for avg in allAvgs:
	OrderedgeneClusters.append(clusDict[avg])

for val in OrderedgeneClusters:
	cluster = geneClusters[val]
	avgs = []
	while(len(cluster)>0):
		minAvg = 1000
		bestGene = -1
		for i in range(len(cluster)):
			curAvg = sum(cluster[i].pVector)/len(cluster[i].pVector)
			if(curAvg < minAvg):
				minAvg = curAvg
				bestGene = i
		print(cluster[bestGene].indicator,cluster[bestGene].descrip,"{0:.3f}".format(minAvg))
		avgs.append(minAvg)
		del cluster[bestGene]
	print("{0:.3f}".format(sum(avgs)/len(avgs)))
	print()


























