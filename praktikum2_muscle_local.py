rom Bio import SeqIO
from Bio import AlignIO
from Bio.Cluster import cluster
from Bio.Align.Applications import MuscleCommandline
from TreeConstruction import DistanceCalculator
from StringIO import StringIO
import sys
import urllib 
import urllib2 
import csv 

# the function read a .fasta-file and perform a multiple alignment
def readFastaAndAlign (filename):
	#create a muscle statement with the sequences within the .fasta-file
	muscle_client = MuscleCommandline(input=filename)
	#write the muscle output (default .fasta format) into a standart output string
	stdout, stderr = muscle_client()
	#safe the alignment result into a handle
	align = AlignIO.read(StringIO(stdout), "fasta")
	return align

#upgma() perform the upgma algorithmn based on a given multiple alignment
def upgma (multalign):
	calculator = DistanceCalculator('identity')
	dm = calculator.get_distance(multalign)
	print dm

#def bootstrapping ():
	


if __name__ == '__main__':
	multalign = readFastaAndAlign(sys.argv[1])
	upgma(multalign)
