import timeit
import numpy as np
import scipy.spatial.distance as dist

# NOTE: this file is tests to help make design choices. this is not overall benchmarking.




'''
# EUCLIDEAN DISTANCE
# what's the faster way to compute euclidean distance?
# 	np.linalg.norm(a-b) vs scipy.distance.euclidean
# 	int arrays 			vs default arrays

# results of test
# random int arrays take longer to generate than default arrays
# after subtracting out array creation time, 
#	euclidean distance is faster on ints for both numpy and scipy (like a third faster)
# for ints, numpy is tiny bit faster 
#	0.0014 sec difference between scipy and numpy
#	the overall numpy function takes 0.0255 sec (excluding array creation), so 5% difference
# since scipy.distance comes with a lot of other distance functions 
# 	and marginal speed difference, gonna use that with int arrays.

'''
def f1():
	a = np.random.randint(0, 10, size=100000)
	b = np.random.randint(0, 10, size=100000)
	dist1 = np.linalg.norm(a-b)

def f2():
	a = np.random.randint(0, 10, size=100000)
	b = np.random.randint(0, 10, size=100000)
	dist2 = dist.euclidean(a,b)

def f3():
	c = np.random.rand(100000)
	d = np.random.rand(100000)
	dist3 = np.linalg.norm(c-d)

def f4():
	c = np.random.rand(100000)
	d = np.random.rand(100000)
	dist4 = dist.euclidean(c,d)


def f_init12():
	a = np.random.randint(0, 10, size=100000)
	b = np.random.randint(0, 10, size=100000)

def f_init34():
	c = np.random.rand(100000)
	d = np.random.rand(100000)


# results
'''
>>> timeit.timeit(f_init12, number = 100)
0.4179811477661133
>>> timeit.timeit(f1, number = 100)
0.4435720443725586
>>> timeit.timeit(f2, number = 100)
0.44500303268432617
>>> 
>>> 
>>> timeit.timeit(f_init34, number = 100)
0.25023698806762695
>>> timeit.timeit(f3, number = 100)
0.28459787368774414
>>> timeit.timeit(f4, number = 100)
0.30179691314697266
>>> 
'''