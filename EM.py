import random
import numpy as np
import scipy.stats as stat 
import scipy.spatial.distance as spdist
from copy import copy, deepcopy
from parseHapmap import *

import mixture # python mixed model pkg

# Expectation maximization
indivs, genoArr = runParse()
K = 2

N = len(indivs)
M = len(indivs[0].geno)

thetas = np.zeros(k, N) # each row is pop, each col is snp

Ps = np.empty(N,k) # each row is indiv, each col is pop
Ps.fill(1.0/k)     # assume P_js are uniform to start



def EMobj(indivs, thetas, Ps):
	'''log likelihood version of EM objective function. Want to maximize this likelihood.'''

	SNPs, clusters, indivs

	obj = 0

	# dumb for loops, can replace with numpy element wise stuff
	# have to draw out arrays to figure out how to do it element-wise
	# may have to rearrange array/summation order (ex/ which is col, which is outer loop)
	for j in range(N):
		for s in range(K):
			for i in range(M):
				g = indivs[j].geno[i]
				theta_si = thetas[s][i]
				P_js = Ps[j][s]
				obj += g*log(thetas_si) + (2-g)*log(1-thetas_si) + log(P_js)
	
	# 	numpy.log(arr)

	return obj


def myEM(indivs, thetas, Ps):
	# could move the globals in here

	# 1) Expectation
	# find best theta_si given P_js's
	# since (theta^g)(1-theta)^(2-g) follows binomial dist, 
	#	best theta_si = (# in S with geno i)/|S|
	popSizes = np.empty(K) # how many have been assigned to pop S
	for s in range(K):


	for s in range(K):
		for i in range(M):
			thetas[s][i] = 1.0*sum([indivs[j].geno[i] for j in range(N)])/popSizes[s]


	# 2) Maximization
	# find best P_js given thetas
	# P_js - maybe restrict them to {0,1/normalized, 1} or {0,0.33, 0.66, 1} to make problem easier?
	#		would be an exponential number of combos, and still need to normalize too
	# 
	# note: theta_si's aren't independent since some snps correlated but may assume otherwise
	for j in range(N):
		for s in range(K):
			oldP = Ps[j][s]
			thetas[s]
			newP = oldP * theta_si

			print)









'''
 The EM is called with a maximum number of steps of 40 and convergence criterion (difference in log likelihood of subsequent iterations) of 0.1
     >>> DNA = mixture.Alphabet(['A','C','G','T']) 
    >>> comp = mixture.MultinomialDistribution(6,4,[0.25,0.25,0.25,0.25],DNA)


genoAlphabet = mixture.Alphabet(['0', '1', '2'])
model = mixture.MultinomialDistribution(M, 3, [theta_si] , DNA)
m = some_model
m.EM(data,40,0.1)
print m
clust = m.classify(data)
# clust is now an array of cluster labels


'''

















