'''
# TODO weirdass bug array containing nan spdist

'''
import kmeans
#import EM # TODO finish EM

import parseHapmap
import PCA_nocluster
import copy
import timeit
import random
import numpy as np

import time

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
	indices = random.sample(range(len(genoArr_)), N)
	indivs_copy = np.array([copy.deepcopy(indivs_[i]) for i in indices])
	for i in range(N):
		indivs_copy[i].geno = np.array(indivs_copy[i].geno[:M])
		indivs_copy[i].j = i 		# also update position in new indiv list

	genoArr_copy = np.array([genoArr_[i][:M] for i in indices])
	

	### TIMING ###
	pcaTime = 0
	def pca_i(): # zero input fxn for timeit
		return PCA_nocluster.pca_transform(indivs_copy, genoArr_copy, n_components)

	if usePCA:
		print("timing pca...") # avg 2 runs. genoArr_copy isn't changed from run-to-run
		#(indivs_copy does get changed, but shouldn't affect run since it restarts each time)
		# unlike kmeans, pca is deterministic so runtime shouldn't vary. also pca step is slower, bottleneck.
		genoArr_copy_geno = genoArr_copy  # make a copy of genotype data first (as opposed to components)

		pcaTime = timeit.timeit(pca_i, number = 2)/2.0 
		# also update genoArr_copy for kmeans 
		pcaObj, genoArr_copy = PCA_nocluster.pca_transform(indivs_copy, genoArr_copy, n_components)


	def kmeans_i(): # zero input fxn for timeit
		return kmeans.kmeans(indivs_copy, genoArr_copy, K)

	# timing kmeans, # run 10x per data pt, avg
	print("timing k-means 10x...")
	kmeansTime = timeit.timeit(kmeans_i, number = 10)/10.0 

	### QUALITY ###
	if usePCA: # reset indivs_copy before running pca
		for i in range(len(indivs_copy)):
			indivs_copy[i].geno = genoArr_copy_geno[i]


	# TODO pca deterministic so only need to run once for all 10 quality trials
	majFracAvgsByRun = np.zeros(10)

	print("----")
	print("genoArr_copy",genoArr_copy[0])
	print("indivs_copy[0].geno", indivs_copy[0].geno)

	if usePCA:
		# make fresh copy of genoArr to run PCA on
		genoArr_copy = genoArr_copy_geno
		print("genoArr_copy_geno", genoArr_copy_geno[0])
		pcaObj, genoArr_copy = PCA_nocluster.pca_transform(indivs_copy, genoArr_copy, n_components)}


	for run in range(10):


		print("genoArr_copy",genoArr_copy[0], len(genoArr_copy[0]))
		print("indivs_copy[0].geno", indivs_copy[0].geno, len(indivs_copy[0].geno))
		centers = kmeans.kmeans(indivs_copy, genoArr_copy, K, maxIter = 1000, verbose = False)
		kmeansObj = kmeans.kmeansObj(indivs_copy, centers)
		majPops, majFracs, clusterSizes = majorityPop(indivs_copy, K)

		# majPops[s] = majority pop of cluster s
		# majFracs[s] = % cluster s that is max 

		# TODO should weight mean by size of cluster?
		# TODO is it ok to weight homogeneity equal for cluster with 1 indiv and cluster with 1000
		majFracAvgsByRun[run] = 1.0*sum(majFracs)/K

		if usePCA: # reset from components to genotypes for next run
			genoArr_copy = genoArr_copy_geno
			for i in range(len(indivs_copy)):
				indivs_copy[i].geno = genoArr_copy_geno[i]


	majFrac_avg = np.mean(majFracAvgsByRun)
	majFrac_std = np.std(majFracAvgsByRun)

	# pcaTime is 0 if pca isn't used
	return kmeansObj, majFrac_avg, majFrac_std, kmeansTime, majPops, majFracs, pcaTime



		
# repeat tests with and without PCA

# actual N = 94 + 171 + 165 = 430 for CHB + MKK + CEU
# actual M = 20109
# true K = 3 (CHB, MKK, CEU)

indivs_, genoArr_ = parseHapmap.runParse()
print("Using DEU, CHB and MKK") # TODO update if i change pops
# though uneven pop size may be good to show weakness of kmeans

def runBenchmark_multiNM(Nd, Md, K, usePCA = False, n_components = 10):
	'''run tests for given number of components and write results to file. Use input as default N,M,k'''

	#print("Benchmarking K-means, varying N with M = %d, K = %d" % (M, K))
	print("Benchmarking PCA + K-means, varying N with M = %d, K = %d, n_components = %d" % (Md, K, n_components))

	# lists of results
	kmeans_N = [] # the value of N used for the test
	kmeansObj_N = []
	majFrac_avg_N = []
	majFrac_std_N = []
	kmeansTime_N = []
	majPops_N = []
	majFracs_N = []
	pcaTime_N = []

	for N in range(20, 410, 10): # sample N indivs randomly
		print("Now testing N = %d" % N)
		#kmeansObj, majFrac_avg, kmeansTime, majPops, majFracs = runBenchmark(N, M, K)
		kmeansObj, majFrac_avg, majFrac_std, kmeansTime, majPops, majFracs, pcaTime = runBenchmark(N, Md, K, usePCA=usePCA, n_components=n_components)

		kmeans_N.append(N)
		kmeansObj_N.append(kmeansObj)
		majFrac_avg_N.append(majFrac_avg)
		majFrac_std_N.append(majFrac_std)
		kmeansTime_N.append(kmeansTime)
		majPops_N.append(majPops)
		majFracs_N.append(majFracs)
		pcaTime_N.append(pcaTime)



	outfBase = 'results_pca%d_kmeans_N_M%d' % (n_components, Md)
	if not usePCA:
		outfBase = 'results_kmeans_N_M%d' % (Md)

	print("writing %s ..." % outfBase)

	with open(outfBase + '.txt', 'w') as outf:
		outf.write('\t'.join(['N', 'kmeans_obj', 'majFrac_avg', 'majFrac_std', 'kmeans_time', 'pca_time']) + '\n')
		for i in range(len(kmeans_N)):
			outf.write("\t".join(
				[str(var) for var in [ kmeans_N[i], kmeansObj_N[i], majFrac_avg_N[i], majFrac_avg_N[i], kmeansTime_N[i], pcaTime_N[i] ]]
			) + "\n")


	# repeat for varying M, use first M snps
	#print("Benchmarking K-means, varying M with N = %d, K = %d" % (N, K))
	print("Benchmarking PCA + K-means, varying M with N = %d, K = %d, n_components = %d" % (Nd, K, n_components))

	# lists of results
	#kmeans_M = range(200, 20000, 200) # the value of M used for the test
	kmeans_M = range(200, 20200, 200) # the value of M used for the test
	kmeansObj_M = []
	majFrac_avg_M = []
	majFrac_std_M = []
	kmeansTime_M = []
	majPops_M = []
	majFracs_M = []
	pcaTime_M = []


	for M in kmeans_M: 
		print("Now testing M = %d" % M)
		kmeansObj, majFrac_avg, majFrac_std, kmeansTime, majPops, majFracs, pcaTime = runBenchmark(Nd, M, K)
		kmeansObj_M.append(kmeansObj)
		majFrac_avg_M.append(majFrac_avg)
		majFrac_std_M.append(majFrac_std)
		kmeansTime_M.append(kmeansTime)
		majPops_M.append(majPops)
		majFracs_M.append(majFracs)

	# output file:
	outfBase = 'results_pca%d_kmeans_M_N%d' % (n_components, Nd)
	if not usePCA:
		outfBase = 'results_kmeans_M_N%d' % (Nd)
	print("writing %s ..." % outfBase)
	with open(outfBase + ".txt", 'w') as outf:
		outf.write('\t'.join(['M', 'kmeans_obj', 'majFrac_avg', 'majFrac_std','kmeans_time', 'pca_time']) + '\n')
		for i in range(len(kmeans_M)):
			outf.write("\t".join(
				[str(var) for var in [ kmeans_M[i], kmeansObj_M[i], majFrac_avg_M[i], majFrac_std_M[i], kmeansTime_M[i], pcaTime_N[i] ]]
			) + "\n")


	return


''' # testing if copy is successful
N = 100
M = 20000
n_components = 5
indices = random.sample(range(len(genoArr_)), N)
indivs_copy = np.array([copy.deepcopy(indivs_[i]) for i in indices])
for i in range(N):
	indivs_copy[i].geno = np.array(indivs_copy[i].geno[:M])
	indivs_copy[i].j = i 		# also update position in new indiv list

genoArr_copy = np.array([genoArr_[i][:M] for i in indices])
print("starting pca_transform")
pcaObj, genoArr_copy = PCA_nocluster.pca_transform(indivs_copy, genoArr_copy, n_components)
print("finished pca_transform")
'''

# testing with one benchmark
runBenchmark(N=200, M=100, K=3, usePCA=True, n_components=2)

# run benchmarks with PCA for varying number of components

# default values of N,M,K (when not varied for testing)
N = 300
M = 10000
K = 3

'''
for n_components in [2, 10]: #[2,10,50,100]:
	print("running n_components = %d" % n_components)
	runBenchmark_multiNM(N,M,K,usePCA=True,n_components = n_components)

runBenchmark_multiNM(N,M,K,usePCA=False, n_components = 0)
'''



