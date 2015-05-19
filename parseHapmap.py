​​​​​import numpy as np
import scipy.stats as stat 
from copy import copy, deepcopy

# Lisa Gai

### Globals ###
GENO_LEN = 10 	# num SNPs in geno # may not know ahead of time bc filter sites
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
		refAlleles = []
		altAlleles = []

		# list allowing multiple alt alleles (i.e. 1 to 3 alt alleles)
		# measure how many sites actually have multiple alts 
		# (see if i'm losing sig amount of info)
		variantCounts = []
		allelesMulti = [] # tuples of (ref, alt1, <alt2>, <alt3>)
		# TODO will probably discard sites with multiple variants 
		# (assuming very few such sites and lots of two-variant sites left to use), 
		#	unless can think of a way to calculate the likelihood 
		# 	p^g + (1-p)^(2-g) won't cut it
		# TODO undecided on sites where no variation in the populations used

class Individual:
	'''Genotype 0/1/2 for each site (with special values) and label'''

	# shared variables here: (not sure if I need any)

	# initiator
	def __init__(self, 
				pop, 
				indivID):
		self.truePop = pop 	# keep track of pop for later graph coloring/labeling
		self.indivID = indivID

		self.assignedPop = None
		self.geno = []
		self.momID = None
		self.dadID = None




def parseFile(fileName, pop):
	'''given hapmap formatted data, get 0/1/2 genotypes for individuals'''

	# 0		1		2	3	4		5		6		7			8		9		  10	11
	#rs# alleles chrom pos strand assembly# center protLSID assayLSID panelLSID QCcode NA21741 NA21740
	
	variantCounts = []
	allelesMulti = [] # tuples of (ref, <alt1>, <alt2>, <alt3>) (alts are optional)

	with open(fileName, 'r') as inf:

		# read header
		indivIDs = inf.readline().rstrip().split()[11:]
		indivs = [] # list of Individuals

		for indivID in indivIDs:
			indivs.append(Individual(pop, indivID))

		N = len(indivs) # number of indivs
		
		# read info for each snp
		# NOTE: snps must be in same order in all files, same number of snps
		#		i.e. the ith snp in file has to be same site as in another file
		for line in inf:
			lineL = line.strip().split()

			ref, alt = allelesMulti.append(lineL[1].split("/")) # e.g. A/C

			for indi in indivs:
				g = ()
				indiv.geno.append()


