# alternate clustering:
# perform PCA and then apply kmeans to components

import parseHapmap
import kmeans

import copy
import numpy as np
from sklearn.decomposition import PCA


def examples():
	# example
	X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
	pca = PCA(n_components=2)
	pca.fit(X)
	#print(pca.explained_variance_ratio_) 
	#[ 0.99244...  0.00755...]


	indivs, genoArr = parseHapmap.runParse()
	genoArr_copy = copy.deepcopy(genoArr)

	# with 2 components
	pca2 = PCA(n_components=2)
	pca.fit(genoArr_copy)
	print(pca2)
	#print(pca2.explained_variance_ratio_)
	print(genoArr_copy)
	print('\n\n\n')

	# with 2 components and transform data in place
	pca2_trans = PCA(n_components=2)
	genoArr_trans = pca2_trans.fit_transform(genoArr_copy)
	print(pca2_trans)
	#print(pca2_trans.explained_variance_ratio_) 
	print(genoArr_trans)

	# with 10 components, first 2 components are same as before
	pca10_trans = PCA(n_components = 10)
	genoArr_trans10 = pca10_trans.fit_transform(genoArr_copy)
	print(genoArr_trans10)

	
# TODO convert both indivs_copy[i].geno and genoArr_copy to components, then run kmeans
# TODO cluster on 2, 10, 100, and M components, see if diff
# TODO plot on PC1 and PC2, color-coded by true population
# TODO test whether clusters on same data are identical w/wo pca


def pca_and_kmeans(indivs, genoArr, n_components, k):
	'''transforms genotypes into components and runs kmeans. 
	note: alters input indivs and genoArr, so pass in copy if needed'''

	# run PCA on genotypes, transform genotypes into components
	pcaObj = PCA(n_components = n_components)
	pcaObj.fit_transform(genoArr)

	# also transform genotypes of indivs into components
	for i in range(len(indivs)):
		indivs[i].geno = genoArr[i]

	# run kmeans on resulting components, updating assigned cluster for each
	centers = kmeans(indivs, genos, k, maxIter = 10000, verbose = False)

	return pcaObj

