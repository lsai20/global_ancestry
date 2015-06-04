import random
import numpy as np
import scipy.stats as stat 
import scipy.spatial.distance as spdist
from copy import copy, deepcopy
from parseHapmap import *

# Expectation maximization!!
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

def EM(indivs, thetas, Ps):
	# could move the globals in here

	# 1) Expectation
	# find best theta_si given P_js's
	# since (theta^g)(1-theta)^(2-g) follows binomial dist, 
	#	best theta_si = (# in S with geno i)/|S|
	popSizes = np.empty(K) # how many have been assigned to pop S

	for s in range(K):
		for i in range(M):
			thetas[s][i] = 1.0*sum([indivs[j].geno[i] for j in range(N)])/popSizes[s]

	# 2) Maximization
	# find best P_js given thetas











