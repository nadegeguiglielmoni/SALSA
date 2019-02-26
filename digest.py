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
# 	restrict_dict : dictionnary with length and cutting sites for each sequence of your genome after digestion by the restriction enzyme
def find_restriction_sites(fasta_file, enzyme_list: list) -> dict:
    # Create Restriction enzyme object
    enzymes = Restriction.RestrictionBatch(enzyme_list)

    # Load fasta file
    seq_data = SeqIO.parse(fasta_file, "fasta")

    restrict_dict = {}

    for record in seq_data:
        rest_sites_dict = enzymes.search(record.seq)
        pos = []
        for enz in enzymes:
            pos.extend(rest_sites_dict[enz])
        restrict_dict[record.id] = [len(record.seq), list(np.unique(sorted(pos)))]

    return restrict_dict


# COUNT_CUT_SITES
# input :
#   restrict_dict : dictionnary of cutting sites for each sequence fragments of your genome after digestion by the restriction enzyme
# output :
#   count_dict : dictionnary of counts of cutting sites for each half of a sequence
def count_cut_sites(restrict_dict: dict) -> dict:
    count_cut_sites = {}

    for contig in restrict_dict.keys():

        mid = restrict_dict[contig][0] / 2
        count_cut_sites[contig] = [0, 0]

        for cut in restrict_dict[contig][1]:

            if cut < mid:
                count_cut_sites[contig][0] += 1
            else:
                count_cut_sites[contig][1] += 1

    return count_cut_sites


# WRITE_FRAGMENTS_FILE
# input :
#   count_dict : dictionnary of counts of cutting sites for each half of a sequence
# 	outfile : path to output file
# output :
# 	output file written with number of cut sites on each half of each contig
def write_count_cuts_file(count_dict: dict, outfile):
    """
	Write number of restriction enzyme sites on each half of each sequence.
	"""
    for seq in count_dict:
        short_name = seq.split(" ")
        outfile.write(
            "{0}\t{1}\t{2}\n".format(short_name, count_dict[seq][0], count_dict[seq][1])
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

    # Count restriction sites on each half of each contig
    count_restrict_sites = count_cut_sites(restrict_sites)

    # Output file
    write_count_cuts_file(count_restrict_sites, outname)
