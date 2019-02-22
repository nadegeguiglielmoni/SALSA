# ======================================================================#
# Digest genome
# Last edition : 2019/02/29
# ======================================================================#

import argparse
import numpy as np
import Bio
from Bio import SeqIO
from Bio import Restriction

# ======================================================================#
#                                ARGUMENTS
# ======================================================================#

# assembly : assembled contigs
# enzyme : restriction enzyme(s)
# outname : path to output file

# ======================================================================#
#                                MODULES
# ======================================================================#

# PARSE_ENZYME_INPUT
# input :
#    enzyme_input : str, one or several enzymes separated by commas
# output :
# 	enzyme_list : list, list of enzymes
def parse_enzyme_input(enzyme_input: str):
    """
	Parse the list of enzymes given as input.
	"""
    enzyme_list = enzyme_input.strip(" ").split(",")
    return enzyme_list


# FIND_RESTRICTION_SITES
# input :
# 	fasta_file : fasta file of genome to digest
# 	enz_list : list of enzymes to use for digestion
# output :
# 	restrict_dict : dictionnary of cutting sites for each sequence
def find_restriction_sites(fasta_file, enzyme_list: list) -> dict:
    # Create Restriction enzyme object
    enzymes = Restriction.RestrictionBatch(enzyme_list)

    # Load fasta file
    seq_data = SeqIO.parse(fasta_file, "fasta")

    restrict_dict = {}

    for record in seq_data:
        rest_sites_dict = enzymes.search(record.seq)
        pos = [0, len(record.seq)]
        for enz in enzymes:
            pos.extend(rest_sites_dict[enz])
        restrict_dict[record.id] = list(np.unique(sorted(pos)))

    return restrict_dict


# WRITE_FRAGMENTS_FILE
# input :
# 	restrict_dict : dictionnary, list of restriction sites found for
# 					each sequence
# 	outfile : path to output file
# output :
# 	output file written with all the digested fragments in the sequence
# 	'seqid	frag_start	frag_end'
def write_fragments_file(restrict_dict: dict, outfile):
    """
	Write digested fragments coordinates.
	"""
    for seq in restrict_dict:
        short_name = seq.split(" ")
        for pos in range(0, (len(restrict_dict[seq]) - 1)):
            outfile.write(
                "{0}\t{1}\t{2}\n".format(
                    short_name, restrict_dict[seq][pos], restrict_dict[seq][pos + 1]
                )
            )


# ======================================================================#
#                                MAIN
# ======================================================================#

# GENERATE_DIGESTED_FRAGMENTS
# input :
# 	assembly : fasta file with genome
# 	enzyme : string of comma separated enzymes
# 	outname : path to output file
# output :
# 	output file written with all the digested fragments in the sequence
# 	'seqid	frag_start	frag_end'
def generate_digested_fragments(assembly, enzyme: str, outname):
    """
	Takes a genome and one or several enzymes to digest the genome,
	and returns the digested fragments.
	"""
    # Parse one or several enzymes given as input
    enzymes_input = parse_enzyme_input(enzyme)

    # Find restriction sites for all sequences
    restrict_sites = find_restriction_sites(assembly, enzymes_input)

    # Output file
    write_fragments_file(restrict_sites, outname)
