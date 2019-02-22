import re 
import argparse
import seq_utils as sequ
import Bio
from Bio import SeqIO
from Bio import Restriction

#======================================================================#
#								ARGUMENTS
#======================================================================#

parser = argparse.ArgumentParser()

parser.add_argument("-a", "--assembly", help="assembled contigs")
parser.add_argument("-e", "--enzyme", help="restriction enzyme")

#parser.add_argument("-m", "--mapping", help="mapping of read to contigs in bam format")
#parser.add_argument("-d", "--dir", help="output directory for results",default='out')

args = parser.parse_args()

#======================================================================#
#								MODULES
#======================================================================#

# PARSE_ENZYME_INPUT
# input :
#	enzyme_input : str, one or several enzymes separated by commas
# output :
#	enzyme_list : list, list of enzymes
def parse_enzyme_input(enzyme_input: str) :
	"""
	Parse the list of enzymes given as input.
	"""
	enzyme_list = args.enzyme.strip(" ").split(',')
	return Restriction.RestrictionBatch(enzyme_list)
	
def main():
	seq_data = SeqIO.parse(args.assembly, "fasta")
	
	# Parse one or several enzymes given as input
	enzymes_input = parse_enzyme_input(args.enzyme)

	for record in seq_data :
		enzymes_input.search(record.seq)
	final_enzymes = []

	for seq in seq_data :
		
		id_seq, seq = key, f[key]
		left_count = 0
		right_count = 0
		for enzyme in final_enzymes :
			pos  = [m.start(0) for m in re.finditer(enzyme,seq)]
		 
			length = len(seq)	
			for each in pos:
				if each < length/2:
					left_count += 1
				else:
					right_count += 1

		print("{0}\t{1}\t{2}\n".format(id_seq, left_count, right_count))

	

if __name__ == '__main__':
	main()
