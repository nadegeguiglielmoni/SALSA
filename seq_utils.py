#======================================================================#
#	Functions to handle sequences for people who hold a grudge against 
#	Biopython: 
#		- load fasta file
#		- generate reverse complement
#======================================================================#

def parse_fasta(fasta_file) -> dict :
	"""
	Reads a fasta file.
	"""
	fasta = open(fasta_file, "r")
	seq_record = {}
	current_short_name = ""
	# Part 1: compile list of lines per sequence
	for line in fasta :
		if ">" in line :
			# new name line; remember current sequence's short name
			short_name = long_name.strip().split()[0]
			seq_record[short_name] = ""
		else:
			# append nucleotides to current sequence
			seq_record[short_name] = seq_record[short_name] \
									+ line.strip()
	return fa

def rev_comp( seq: str ) -> str :
	"""
	Generates the reverse complement of a sequence.
	"""
	comp = {	'A':'T', 'C':'G', 'G':'C', 'T':'A', 
				'B':'N', 'N':'N', 'R':'N', 'M':'N', 
				'Y':'N', 'S':'N', 'W':'N', 'K':'N',
				'a':'t', 'c':'g', 'g':'c', 't':'a', 
				'n':'n',' ':'' 
				}
	
	rev_seq = "".join(comp.get(base, base) for base in reversed(seq))
						
	return rev_seq

def seq_length(seq_record: dict) -> dict :
	"""
	Compute sequences length.
	"""
	seq_length = {}
	
	for seq in seq_record.keys() :
		seq_length[seq] = len(seq_record[seq])
		
	return seq_length 
