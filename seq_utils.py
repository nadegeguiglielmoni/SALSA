# ======================================================================#
#     Functions to handle sequences for people who hold a grudge against
#     Biopython:
#         - load fasta file
#         - generate reverse complement
# ======================================================================#

import Bio
from Bio import SeqIO


def parse_fasta(fasta_file) -> dict:
    """
    Reads a fasta file.
    """
    fasta = open(fasta_file, "r")
    seq_record = {}
    # Part 1: compile list of lines per sequence
    for line in fasta:
        if ">" in line:
            # new name line; remember current sequence's short name
            short_name = line.strip().split()[0]
            seq_record[short_name] = ""
        else:
            # append nucleotides to current sequence
            seq_record[short_name] = seq_record[short_name] + line.strip()
    return seq_record


def rev_comp(seq: str) -> str:
    """
    Generates the reverse complement of a sequence.
    """
    comp = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "B": "N",
        "N": "N",
        "R": "N",
        "M": "N",
        "Y": "N",
        "S": "N",
        "W": "N",
        "K": "N",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a",
        "n": "n",
        " ": "",
    }

    rev_seq = "".join(comp.get(base, base) for base in reversed(seq))

    return rev_seq


def make_seq_length_file(assembly, output_file, namelist_file="abc"):
    """
    Compute sequences length.
    """
    seq_data = SeqIO.parse(assembly, "fasta")

    outlength = open(output_file, "w")
    if namelist_file != "abc":
        namelist = open(namelist_file, "w")

    for record in seq_data:
        outlength.write("{0}\t{1}\n".format(record.id, len(record.seq)))
        if namelist_file != "abc":
            namelist.write("{0}\n".format(record.id))

    outlength.close()
    namelist.close()
