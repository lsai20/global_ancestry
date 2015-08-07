### Measure variance explained and plot PC1 vs PC2 ###

import PCA_nocluster
#from sklearn.decomposition import PCA
import kmeans
import parseHapmap
import random
import numpy as np
import copy

indivs_, genoArr_ = parseHapmap.runParse()

def makeInputCopies(N=250, M=9000):
	indices = random.sample(range(len(genoArr_)), N)
	indivs_copy = np.array([copy.deepcopy(indivs_[i]) for i in indices])
	for i in range(N):
		indivs_copy[i].geno = np.array(indivs_copy[i].geno[:M])
		indivs_copy[i].j = i 		# also update position in new indiv list

	genoArr_copy = np.array([genoArr_[i][:M] for i in indices])
	return indivs_copy, genoArr_copy	

def findVarExplained(N=250, M=9000, n_components = 100):
	'''print list of var explained for each component'''
	indivs_copy, genoArr_copy = makeInputCopies(N,M)
	pcaObj = PCA_nocluster.PCA(n_components = n_components)
	pcaObj.fit(genoArr_copy)
	#print(pcaObj.explained_variance_ratio_)
	for var in pcaObj.explained_variance_ratio_:
		print(var)

	return pcaObj


def findPC1_PC2(N=250, M=9000):
	'''print variance explained, and print tab-sep of
	(PC1, PC2, truePop, assignedPop) for each indiv in a run'''
	n_components = 2
	K = 3
	indivs_copy, genoArr_copy = makeInputCopies(N,M)

	# get PC1 and PC2
	pcaObj, genoArr_copy = PCA_nocluster.pca_transform(indivs_copy, genoArr_copy, n_components)

	# get assignedPops
	centers = kmeans.kmeans(indivs_copy, genoArr_copy, K, maxIter = 1000, verbose = False)

	# indivs_copy: 
	# 	.geno now transformed to have two components
	# 	has both assigned and true pop
	indivs_copy = sorted(indivs_copy, key=lambda indiv: indiv.truePop)
	for indiv in indivs_copy:
		s = '\t'.join([str(val) for val in [indiv.geno[0], indiv.geno[1], indiv.truePop, indiv.assignedPop]])
		print(s)

	return


findVarExplained(N=250,M=400)
print("\n\n\n")
findPC1_PC2(N=250,M=400)


