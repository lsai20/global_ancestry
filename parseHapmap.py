​​​​​import numpy as np
import scipy.stats as stat 
from copy import copy, deepcopy

# Lisa Gai

### Globals ###
GENO_LEN = 10 	# num SNPs in geno
NOT_REP = 9    	# '-' in orig hapmap data
MISSING = -9  	# 'N'in orig hapmap data
NO_ID = "NO_ID"
NO_POP = "NNN" 	# population not yet assigned

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

	# ref allele (a where aa = 0) and alt allele (A where AA = 2)
	alleles = np.chararray([2, GENO_LEN])

	# if allowing multiple alt alleles (at least one alt)
	numAlleles = np.empty([1, GENO_LEN]) 	# how many variants at each site
	allelesMulti = np.chararray([4, GENO_LEN])
	# TODO will probably discard sites with multiple variants 
	# (assuming very few and lots of two-variant sites to use), 
	#	unless can think of a way to calculate the likelihood 
	# 	p^g + (1-p)^(2-g) won't cut it



class Individual:
	'''Genotype 0/1/2 for each site (with special values) and label'''

	# keep track of pop for later graph coloring/labeling
	truePop = NO_POP
	assignedPop = NO_POP
	indivID = NO_ID
	momID = NO_ID # if their mom not in sample, it's no id
	dadID = NO_ID # else it's 
	geno = np.empty([1, GENO_LEN])	# [0, 1, 2, 1, ...], np array


	def __init__(self, pop, indivID):
		self.truePop = pop
		self.indivID = indivID


