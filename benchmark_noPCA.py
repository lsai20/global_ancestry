# now redundant w usePCA=False on benchmark.py
import kmeans
#import EM # TODO finish EM

import parseHapmap as parseHap
import copy
import timeit
import random
import numpy as np


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



def runBenchmark(N=200, M=10000, K=3):

	# note that order of clusters will vary from run to run, so track which true pop with each cluster
	#	along with true frac count

	# sample N indivs without replacement, and also get first M snps for each geno
	indices = random.sample(range(len(genoArr)), N)
	indivs_copy = np.array([copy.deepcopy(indivs[i]) for i in indices])
	for i in range(N):
		indivs_copy[i].geno = np.array(indivs_copy[i].geno[:M])
		indivs_copy[i].j = i 		# also update position in new indiv list

	genoArr_copy = np.array([genoArr[i][:M] for i in indices])

	# old version: just take first N since already shuffled
	#indivs_copy = copy.deepcopy(indivs[:N])
	#genoArr_copy = copy.deepcopy(genoArr[:N])


	def kmeans_i():
		'''zero input fxn for timeit'''
		return kmeans.kmeans(indivs_copy, genoArr_copy, K)

	# timing kmeans, # run 10x per data pt, avg
	kmeansTime = timeit.timeit(kmeans_i, number = 10)/10.0 

	# measure with kmeans objective fxn
	centers = kmeans.kmeans(indivs_copy, genoArr_copy, K, maxIter = 1000, verbose = True)
	kmeansObj = kmeans.kmeansObj(indivs_copy, centers)
	majPops, majFracs, clusterSizes = majorityPop(indivs_copy, K)

	# majPops[s] = majority pop of cluster s
	# majFracs[s] = % cluster s that is max 

	# TODO should weight mean by size of cluster?
	# TODO is it ok to weight homogeneity equal for cluster with 1 indiv and cluster with 1000
	kmeansMajFrac_avg = 1.0*sum(majFracs)/K


	return kmeansObj, kmeansMajFrac_avg, kmeansTime, majPops, majFracs



		

# actual N = 94 + 171 + 165 for CHB + MKK + CEU
# actual M = 20109
# true K = 3 (CHB, MKK, CEU)

indivs, genoArr = parseHap.runParse()
print("Using DEU, CHB and MKK") # TODO update if i change pops
# though uneven pop size may be good to show weakness of kmeans

M = 500
K = 3
print("Benchmarking K-means, varying N with M = %d, K = %d" % (M, K))
# lists of results
kmeans_N = [] # the value of N used for the test
kmeansObj_N = []
kmeansMajFrac_avg_N = []
kmeansTime_N = []
kmeansMajPops_N = []
kmeansMajFracs_N = []
for N in range(20, 200, 10): # could change skip size to get more data
	print("Now testing N = %d" % N)
	kmeansObj, kmeansMajFrac_avg, kmeansTime, kmeansMajPops, kmeansMajFracs = runBenchmark(N, M, K)
	kmeans_N.append(N)
	kmeansObj_N.append(kmeansObj)
	kmeansMajFrac_avg_N.append(kmeansMajFrac_avg)
	kmeansTime_N.append(kmeansTime)
	kmeansMajPops_N.append(kmeansMajPops)
	kmeansMajFracs_N.append(kmeansMajFracs)
# output file:
print("writing kmeans_N results to file")
with open("results_kmeans_N_M500.txt", 'w') as outf:
	outf.write('\t'.join(['N', 'kmeans_obj', 'kmeans_majFrac_avg', 'kmeans_time']) + '\n')
	for i in range(len(kmeans_N)):
		outf.write("\t".join(
			[str(var) for var in [ kmeans_N[i], kmeansObj_N[i], kmeansMajFrac_avg_N[i], kmeansTime_N[i] ]]
		) + "\n")
with open("results_kmeans_N_M500_majorities.txt", 'w') as outf:
	for i in range(len(kmeans_N)):
		outf.write("\n".join(str(var) for var in kmeansMajPops_N[i] + kmeansMajFracs_N[i]) + "\n")
	outf.write("\n")



'''
# repeat for varying M
N = 150
K = 3
print("Benchmarking K-means, varying M with N = %d, K = %d" % (N, K))

# lists of results
kmeans_M = range(200, 20000, 200) # the value of N used for the test
kmeansObj_M = []
kmeansMajFrac_avg_M = []
kmeansTime_M = []
kmeansMajPops_M = []
kmeansMajFracs_M = []

for M in kmeans_M: 
	print("Now testing M = %d" % M)
	kmeansObj, kmeansMajFrac_avg, kmeansTime, kmeansMajPops, kmeansMajFracs = runBenchmark(N, M, K)
	kmeansObj_M.append(kmeansObj)
	kmeansMajFrac_avg_M.append(kmeansMajFrac_avg)
	kmeansTime_M.append(kmeansTime)
	kmeansMajPops_M.append(kmeansMajPops)
	kmeansMajFracs_M.append(kmeansMajFracs)

# output file:
print("writing kmeans_M results to file")
with open("results_kmeans_M.txt", 'w') as outf:
	outf.write('\t'.join(['M', 'kmeans_obj', 'kmeans_majFrac_avg', 'kmeans_time']) + '\n')
	for i in range(len(kmeans_M)):
		outf.write("\t".join(
			[str(var) for var in [ kmeans_M[i], kmeansObj_M[i], kmeansMajFrac_avg_M[i], kmeansTime_M[i] ]]
		) + "\n")

with open("results_kmeans_M_majorities.txt", 'w') as outf:
	for i in range(len(kmeans_M)):
		outf.write("\n".join(str(var) for var in kmeansMajPops_M[i] + kmeansMajFracs_M[i]) + "\n")
	outf.write("\n")
'''
