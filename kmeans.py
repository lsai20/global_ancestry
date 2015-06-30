import random
import numpy as np
import scipy.stats as stat 
import scipy.spatial.distance as spdist
from copy import copy, deepcopy
from parseHapmap import *

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


def pickInitial(indivs, k):
	'''randomly pick centers. return the centers if successful, 
	return None if any centers are duplicates (identical genotype/same position)
	'''
	# NOTE: must have at least k unique values for genotype
	# if you have multiple centers in exact same position, 
	#	other points won't have a unique closest center. may end up with some centers empty depending on tie break.
	N = len(indivs)
	M = len(indivs[0].geno)

	alreadyUsedGenos = [] # can't use set bc genos are np array/list (unhashable)
	numTries = 0
	label = 0 	# clusterid ranges from 0 to k-1
	centers = []

	randomJs = random.sample(range(N), k)

	for j in randomJs: 	# j = which indiv's geno is selected
		# (have to convert to list for checking membership)
		if list(indivs[j].geno) in alreadyUsedGenos:
			return None

		center = Center(label, indivs[j].geno)
		centers.append(center)
		alreadyUsedGenos.append(list(indivs[j].geno))
		label += 1

	return centers


def kmeans(indivs, genos, k, maxIter = 10000, verbose = False):
	'''run lloyds on a list of Individuals'''

	# note: can't use only array of genos as input - need to track labels for indivs for visualizing
	# return: indivs list updated with assigned pop, or dict of assigned pop, or something

	N = len(indivs)
	M = len(indivs[0].geno)

	# initial step - pick k genos to be centers at random

	# initialize randomly up to 3+ times (can't allow duplicate values in centers)
	# (allow more chances if k small, i.e. k = 2 for 50/50 split of identical values is ok, but may need more chances)
	# after that, warn user to pick smaller k
	# this is kind of dumb, but k << N, M so w/e
	numTries = 3 + max(0,5-k)**2
	for i in range(numTries):
		centers = pickInitial(indivs, k)
		if centers: # if we successfully picked centers, use them
			break

	# if fail to pick within allowed tries
	if centers == None:
		print("WARNING: unable to pick k = %d unique initial center values after 3+ attempts" % k)
		print("Your data may not contain enough unique values for k = %d clusters, or the cluster sizes may be very uneven." % k)
		print("Consider using a smaller k, or using an alternative to k-means clustering.")
		print("Now quitting.")
		#print("ERROR: failed to pick unique initial center values after 3 attempts. Quitting.")
		return

	# the loop
	for t in range(maxIter):
		# clear centers for upcoming round
		for center in centers:
			center.members = []

		# each Individual finds its closest center, which then tracks that Individual
		for indiv in indivs:
			closestLabel = np.argmin([spdist.euclidean(indiv.geno, center.geno) for center in centers])
			# NOTE: this indiv.j step requires indiv.j to correspond to the j which is indiv's actual position in list
			centers[closestLabel].members.append(indiv.j)

		# each center finds the centroid of its own cluster
		noneUpdated = True 	# if it converges (no update to centers), bail

		for center in centers:
			a = genos[center.members,] # rows that are member of this center

			centroidGeno = np.mean(a, axis = 0) # row of col avgs
			if not np.array_equal(centroidGeno, center.geno): # if new centroid != current center
				noneUpdated = False
				center.geno = centroidGeno

		###
		if verbose:
			# update indivs with intermediate cluster assignments
			'''
			for center in centers:
				for j in center.members:
					indivs[j].assignedPop = center.label

			print('\nK-MEANS')
			countClustering(indivs, k) # TODO verbose, see how clusters change over time
			print("\n")

			print("\ndistances in each cluster")
			print([int(
					sum([spdist.euclidean(center.geno, indivs[j].geno) for j in center.members])) 
				for center in centers])
			print("members in each cluster")
			print([[indivs[j].geno[0] for j in center.members]
				for center in centers ]
			)
			'''
			pass
			#print(kmeansObj(indivs, centers))
			#print(int( kmeansObj(indivs, centers) ))  # round to int
		###

		if noneUpdated: # if there were no changes
			if verbose:
				print("Stopped at iteration %d" % t)
			break


	# update indivs with final cluster assignments
	for center in centers:
		for j in center.members:
			indivs[j].assignedPop = center.label

	return centers


def countClustering(indivs, k):
	'''create count dict of clustered indivs, given number of clusters'''

	clusterL = [[] for i in range(k)]
	for indiv in indivs:
		clusterL[indiv.assignedPop].append(indiv.truePop)
		#clusterL[indiv.assignedPop].append((indiv.truePop, indiv.famID))
		#clusterL[indiv.assignedPop].append(indiv.famID)

	#clusterL[0].sort()
	#clusterL[1].sort()
	#famLL = [[indiv.fam for indiv in cluster] for cluster in clusterL]

	for i in range(k):
		print("\ncluster %d" % i)
		print(sorted(clusterL[i]))

	return clusterL

def kmeansObj(indivs, centers):
	'''measure cluster homogenity with sum of intercluster distance'''
	# requires indivs to be in order of indiv.j
	return sum([spdist.euclidean(center.geno, indivs[j].geno) 
		for center in centers 
			for j in center.members])
	
'''
indivs, genoArr = runParse()
k = 2
centers = kmeans(indivs, genoArr, k, maxIter = 2, verbose = True)
clusterL = countClustering(indivs, k)
'''





