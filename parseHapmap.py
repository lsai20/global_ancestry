​​​​​import numpy as np
import scipy.stats as stat 
from copy import copy, deepcopy

# Lisa Gai

### Globals ###
GENO_LEN = 10 	# num SNPs in geno
NOT_REP = 9    	# '-' in orig hapmap data
MISSING = -9  	# 'N'in orig hapmap data

# population codes (strings for now, could be ints)
ASW = 'ASW' # african americans in SW usa
CEU = 'CEU' # northern european-amer in utah
CHB = 'CHB' # han chinese in beijing
CHD = 'CHD' # chinese in denver, colorado
GIH = 'GIH' # gujarati in houston, texas
JPT = 'JPT' # jp in tokyo
LWK = 'LWK' # luhya in webuye, kenya
MEX = 'MEX' # mexican-amer in LA
MKK = 'MKK' # maasai in kinyawa, kenya
TSI = 'TSI' # toscan in italy
YRI = 'YRI' # yoruba in italy

### Classes ###
class FullGenotype:
	'''reference and alternate allele(s)'''

	def __init__(self):
		# ref allele (a where aa = 0) and alt allele (A where AA = 2)
		# for this chr, sites with more than two variants are filtered out
		alleles = np.chararray([2, GENO_LEN])

		# list allowing multiple alt alleles (i.e. 1 to 3 alt alleles)
		# measure how many sites actually have multiple alts 
		# (see if i'm losing sig amount of info)
		numAlleles = np.empty([1, GENO_LEN]) 	# how many variants at each site
		allelesMulti = np.chararray([4, GENO_LEN])
		# TODO will probably discard sites with multiple variants 
		# (assuming very few and lots of two-variant sites to use), 
		#	unless can think of a way to calculate the likelihood 
		# 	p^g + (1-p)^(2-g) won't cut it

class Individual:
	'''Genotype 0/1/2 for each site (with special values) and label'''

	# shared variables here: (not sure if I need any)

	# initiator
	def __init__(self, 
				j, 
				pop, 
				indivID, 
				assignedPop = None, 
				geno = np.empty([1, GENO_LEN]), 
				momID = None, 
				dadID = None):
		self.j = j  	# individuals have j = 0, 1, 2, ..., N-1 (N total across all populations)
		self.truePop = pop 	# keep track of pop for later graph coloring/labeling
		self.assignedPop = assignedPop
		self.indivID = indivID
		self.geno = geno
		self.momID = momID
		self.dadID = dadID





