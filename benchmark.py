'''
TODO come up with some functions to measure cluster quality 
	(k means objective function, inter-cluster distance)
TODO decide final pops in parseHapmap and hardcode appropriate labels
'''

# TODO why is sklearn's pca so much slower here than on iris dataset (N=150, 4 features)???

import kmeans
#import EM # TODO finish EM

import parseHapmap
#import PCA_nocluster
import copy
import timeit
import random
import numpy as np

from sklearn.decomposition import PCA, TruncatedSVD, SparsePCA
#import sklearn.decomposition.IncrementalPCA



def pca_transform(indivs, genoArr, n_components):
	'''transforms genotypes into components and returns sklearn PCA obj. 
	note: alters input indivs and genoArr, so pass in copy if needed'''

	# run PCA on genotypes, transform genotypes into components
	pcaObj = SparsePCA(n_components = n_components)
	pcaObj.fit_transform(genoArr)

	# also transform genotypes of indivs into components
	for i in range(len(indivs)):
		indivs[i].geno = genoArr[i]

	# can then run kmeans on resulting components, updating assigned cluster for each
	#centers = kmeans(indivs, genos, k, maxIter = 10000, verbose = False)

	return pcaObj, indivs, genoArr
	# not sure why indivs and genoArr aren't getting updated after being passed in



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



def runBenchmark(N=200, M=10000, K=3, usePCA=False, n_components=2):

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

	pcaTime = 0
	def pca_i():
		'''zero input fxn for timeit'''
		pca_transform(indivs_copy, genoArr_copy, n_components)

	if usePCA:
		print("timing pca...")
		pcaTime = timeit.timeit(pca_i) # only run once since it changes size of genos
		#print(genoArr_copy)			# check that genoArr_copy was transformed properly

		# have to run PCA outside of timeit to get actual components
		# note: pca_transform(...) returns the pcaObj, can be used for analysis elsewhere
		# 		or visualizing a single run
		#pcaObj, indivs_copy, genoArr_copy = PCA.pca_transform(indivs_copy, genoArr_copy, n_components)
		print(indivs_copy[0].geno)

	def kmeans_i():
		'''zero input fxn for timeit'''
		return kmeans.kmeans(indivs_copy, genoArr_copy, K)

	# timing kmeans, # run 10x per data pt, avg
	print("timing k-means 10x...")
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

	if usePCA: # return types w/e
		return kmeansObj, kmeansMajFrac_avg, kmeansTime, majPops, majFracs, pcaTime

	return kmeansObj, kmeansMajFrac_avg, kmeansTime, majPops, majFracs



		
# repeat tests with and without PCA

# actual N = 94 + 171 + 165 for CHB + MKK + CEU
# actual M = 20109
# true K = 3 (CHB, MKK, CEU)

indivs, genoArr = parseHapmap.runParse()
print("Using DEU, CHB and MKK") # TODO update if i change pops
# though uneven pop size may be good to show weakness of kmeans

def runBenchmark_multiNM(n_components = 10):
	# TODO currently uses PCA all the time

	M = 500 #M = 10000
	K = 3
	#print("Benchmarking K-means, varying N with M = %d, K = %d" % (M, K))
	print("Benchmarking PCA + K-means, varying N with M = %d, K = %d, n_components = %d" % (M, K, n_components))

	# lists of results
	kmeans_N = [] # the value of N used for the test
	kmeansObj_N = []
	kmeansMajFrac_avg_N = []
	kmeansTime_N = []
	kmeansMajPops_N = []
	kmeansMajFracs_N = []
	pcaTime_N = []

	for N in range(20, 200, 10): # could change skip size to get more data
		print("Now testing N = %d" % N)
		#kmeansObj, kmeansMajFrac_avg, kmeansTime, kmeansMajPops, kmeansMajFracs = runBenchmark(N, M, K)
		kmeansObj, kmeansMajFrac_avg, kmeansTime, kmeansMajPops, kmeansMajFracs, pcaTime = runBenchmark(N, M, K, usePCA=True, n_components=n_components)

		kmeans_N.append(N)
		kmeansObj_N.append(kmeansObj)
		kmeansMajFrac_avg_N.append(kmeansMajFrac_avg)
		kmeansTime_N.append(kmeansTime)
		kmeansMajPops_N.append(kmeansMajPops)
		kmeansMajFracs_N.append(kmeansMajFracs)
		pcaTime_N.append(pcaTime)


	# output file:
	print("writing kmeans_N results to file")
	#print("writing kmeans_N results to file")
	#outfName = "results_kmeans_N.txt"
	outfBase = 'results_pca%d_kmeans_N_M%d' % (n_components,M)
	with open(outfBase + '.txt', 'w') as outf:
		outf.write('\t'.join(['N', 'kmeans_obj', 'kmeans_majFrac_avg', 'kmeans_time', 'pca_time']) + '\n')
		for i in range(len(kmeans_N)):
			outf.write("\t".join(
				[str(var) for var in [ kmeans_N[i], kmeansObj_N[i], kmeansMajFrac_avg_N[i], kmeansTime_N[i], pcaTime_N[i] ]]
			) + "\n")

	with open(outfBase + '_majorities.txt', 'w') as outf:
		for i in range(len(kmeans_N)):
			outf.write("\n".join(str(var) for var in kmeansMajPops_N[i] + kmeansMajFracs_N[i]) + "\n")
		outf.write("\n")


	# repeat for varying M
	N = 150
	K = 3
	#print("Benchmarking K-means, varying M with N = %d, K = %d" % (N, K))
	print("Benchmarking PCA + K-means, varying M with N = %d, K = %d, n_components = %d" % (N, K, n_components))

	# lists of results
	#kmeans_M = range(200, 20000, 200) # the value of M used for the test
	kmeans_M = range(200, 800, 200) + range(1, 200, 20) # the value of M used for the test
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
	outfBase = 'results_pca%d_kmeans_M_N%d' % (n_components,N)

	with open(outfBase + ".txt", 'w') as outf:
		outf.write('\t'.join(['M', 'kmeans_obj', 'kmeans_majFrac_avg', 'kmeans_time', 'pca_time']) + '\n')
		for i in range(len(kmeans_M)):
			outf.write("\t".join(
				[str(var) for var in [ kmeans_M[i], kmeansObj_M[i], kmeansMajFrac_avg_M[i], kmeansTime_M[i], pcaTime_N[i] ]]
			) + "\n")

	with open(outfBase + "_majorities.txt", 'w') as outf:
		for i in range(len(kmeans_M)):
			outf.write("\n".join(str(var) for var in kmeansMajPops_M[i] + kmeansMajFracs_M[i]) + "\n")
		outf.write("\n")

	return



runBenchmark(N=20, M=10, K=3, usePCA=True, n_components=2)


# run benchmarks with PCA for varying number of components
'''
for n_components in [1,2,10,50]:
	print("running n_components = %d" % n_components)
	runBenchmark_multiNM(n_components)
'''


