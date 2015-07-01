'''
TODO come up with some functions to measure cluster quality 
	(k means objective function, inter-cluster distance)
TODO decide final pops in parseHap and hardcode appropriate labels
'''


import kmeans
#import EM # TODO finish EM

import parseHapmap as parseHap
import copy
import timeit


def majorityPop(indivs, k):
	'''return majority element and size of majority in each cluster.
	ex/ if cluster0 members are 0.8 CHB and cluster1 members are 0.8 MKK, 
	return [CHB, MKK],[0.8, 0.9]'''

	clusters = [[] for i in range(k)] # true pops in each cluster 
		# e.g. [['mkk', 'mkk'], ['mkk', 'chb', 'chb']]
	clusterSizes = [0 for i in range(k)]

	for indiv in indivs:
		clusters[indiv.assignedPop].append(indiv.truePop)
		clusterSizes[indiv.assignedPop] += 1

	# find proportion made by majority element
	# (indivs in mostCommonPop)/(total cluster size)
	majPops = []
	majFracs = []
	for cluster in clusters:
		mostCommonPop = max(cluster, key=cluster.count)
		mostCommonCount = cluster.count(mostCommonPop)
		frac = 1.0*mostCommonCount/len(cluster) 
		majPops.append(mostCommonPop)
		majFracs.append(frac)

	return majPops, majFracs, clusterSizes





indivs, genoArr = parseHap.runParse()
print("Using DEU, CHB and MKK") # TODO update if i change pops
# though uneven pop size may be good to show weakness of kmeans


# TODO refactor so for N in range(200) you could run benchmark and append its return values to appropriate lists
# actual N = 94 + 171 + 165 for CHB + MKK + CEU
# actual M = 20,000+
# true K = 3 (CHB, MKK, CEU)

def runBenchmark(N=200, M=10000, K=3):
	# TODO change N. 
	# varying N
	kmeansTimes = []
	kmeansMajPops = [] 		# kmeansMaxPops[s] is majority pop of cluster s
	kmeansMajFracs = [] 	# kmeansMaxFracs[s] = % cluster s that is max 
	# note that order of clusters will vary from run to run, so track which true pop with each cluster
	#	along with true frac count


	# sample N indivs without replacement
	indices = random.choice(range(len(genoArr), N))
	indivs_copy = indivs[indices]
	genoArr_copy = genoArr[indices]

	# old version: just take first N since already shuffled
	#indivs_copy = copy.deepcopy(indivs[:N])
	#genoArr_copy = copy.deepcopy(genoArr[:N])

	k = K

	def kmeans_i():
		'''zero input fxn for timeit'''
		return kmeans.kmeans(indivs_copy, genoArr_copy, k)

	# timing kmeans
	time = timeit.timeit(kmeans_i, number = 1) # only run once per data pt
	kmeansTimes.append(time)

	# measure with kmeans objective fxn
	centers = kmeans.kmeans(indivs_copy, genoArr_copy, k, maxIter = 1000, verbose = True)
	kmeansObj = kmeans.kmeansObj(indivs_copy, centers)
	majPops, majFracs, clusterSizes = majorityPop(indivs_copy, k)

	kmeansMajPops.append(majPops)
	kmeansMajFracs.append(majFracs)


	# TODO should weight mean by size of cluster, not useful atm
	# TODO don't weight homogeneity equal for cluster with 1 indiv and cluster with 1000
	kmeansMajFrac_avg = 1.0*sum(kmeansMajFracs*kmeans)/k


	return kmeansObj, kmeansMajFrac_avg, #kmeansTimes, kmeansMajPops, kmeansMajFracs


# TODO print, or work with the variables interactively to make graphs
# K-means,varying total number of individuals from 2 populations
for N in range(20, 200, 10): # could change skip size to get more data
	print("Now testing N = %d" % N)
	runBenchmark(N, M, K)
		



# TODO also vary M

'''
# repeat tests for EM
for N in range(10, 150, 10):

	# EM
	indivs_copy = copy.deepcopy(indivs[:N])
	genoArr_copy = copy.deepcopy(genoArr_copy[:N])
	# timing EM

	# measure EM with kmeans obj
	# measure EM with EM obj

'''
