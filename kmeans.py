import numpy as np
import scipy.stats as stat 
import scipy.spatial.distance as sc.dist
from copy import copy, deepcopy

# lsgai
# k-means clustering (lloyd's algorithm)
# objective: minimize (total within cluster sum of squares)
# Lloyds k-means
# 0. pick k centers at random
# repeat 1-3 until convergence (new centroids same as points)
# 	1. each data point finds closest center (center owns that point)
# 	2. each center finds centroid of its cluster
# 	3. update centroids to be new centers. 

# single linkage hierarchical clustering (either here or another file)

# TODO implement k-means in numpy/scipy (sklearn is low priority)
# TODO implement hierachical clustering in numpy/scipy
# TODO implement EM (at least the E step, think about the M step)

### distance options

# EUCLIDEAN
'''square of Euclidean distance = sum (g_i - g_j)^2.
uses squared difference at each snp, but don't bother sqrting sum. 
penalizes difference between 0 and 2 more strongly than manDist.'''


# MANHATTAN/TAXICAB/CITYBLOCK
'''Manhattan distance = sum (abs(g_i - g_j)). uses unsquared difference. '''

# HAMMING
'''Hamming dist = sum (g_i == g_j). exact matches only, so 0 is equally different 
 from 1 and from 2'''



class Center:
	'''center of cluster, used in k-means'''

	def __init__(self, label, geno):
		self.label = label 	# label is int 0, ..., k-1
		self.geno = geno	# ex/ np.arr([0.7, 1.5, ...])
		self.members = [] # list of member indices [j_1, j_2, ...]


def kmeans(indivs, genos, k, maxIter = 10000):
	'''run lloyds on a list of Individuals'''
	# note: can't use only array of genos as input - need to track labels for indivs for visualizing
	# return: indivs list updated with assigned pop, or dict of assigned pop, or something

	# initial step - pick k genos to be centers at random
	centers = []
	label = 0 	# clusterid ranges from 0 to k-1
	randomJs = random.sample(range(N), k):
	for j in randomJs: 	# j = which indiv's geno is selected
		center = Center(label, indivs[j].geno)
		centers.append(center)
		label += 1


	for t in range(maxIter):

		# each Individual finds its closest center, which then tracks that Individual
		center.members = [] # clear slate from any prev run
		for indiv in indivs:
			closestLabel = np.argmin([sc.dist(indiv.geno, center.geno) for center in centers])
			centers[closestLabel].members.append(indiv.j)

		# each center finds the centroid of its own cluster
		haveAnyUpdated = False 	# if it converges (no update to centers), bail

		for center in centers:
			a = genos[center.members,] # rows that are member of this center
			centroidGeno = np.mean(a, axis = 0) # row of col avgs
			if center.geno != centroidGeno:
				haveAnyUpdated = True
				center.geno = centroidGeno

		if !haveAnyUpdated: # if there were no changes
			break

		centers = newCenters

	# update indivs with final cluster assignments
	for center in centers:
		for j in center.members:
			indivs[j].assignedPop = center.label

	return


truePopsL = [CHB, MKK]

def countClustering(indivs, truePopsL = , k):
	'''create count dict of clustered indivs, given number of clusters'''


	for indiv in indivs:
		indiv.truePop




