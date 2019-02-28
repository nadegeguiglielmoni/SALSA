# ======================================================================#
# Process Hi-C data to make links
# Last edition : 2019/02/28
# ======================================================================#

import sys
import numpy as np

# ======================================================================#
#                                ARGUMENTS
# ======================================================================#

# mapping_data : path to Hi-C data
# restrict_data : path to restriction sites counts per half of contigs
# contig_len_data : path to file with contigs length
# iteration : integer, iteration number
# contig_links_file : path to file where links should be written
# dup_data : path to duplicated links file
# avoid_data : path to links to avoid

# ======================================================================#
#                                MODULES
# ======================================================================#

# MAKE_LINK_KEY
# input :
#   line : string with a link's contig names separated by a delimiter
#   delimiter : character to separate the name of two contigs, by default tab
# output :
#   contig_str : string of contigs that make up a link
def make_link_key(line: str, delimiter="\t") -> str:
    """
    Make a formatted string containing both contigs that make up a link.
    """
    attrs = line.strip().split(delimiter)
    attrs = sorted(attrs)
    contig_str = attrs[0] + "$" + attrs[1]
    return contig_str


# LOAD_LINKS
# input :
#   infile : path to file with a bunch of links
# output :
#   link_list : list of links loaded from the input file
def load_links(infile) -> list:
    """
    Load links from a file where a link is described by the names of the two
    contigs linked, separated by a tab by default.
    """
    link_list = []
    with open(infile, "r") as f:
        for line in f:
            link_list.append(make_link_key(line))
    link_list = np.unique(link_list)
    return link_list


# COMPUTE_LINKS
# input :
#   mapping_file : path to file with processed Hi-C reads
#   contig_len : dictionnary with length of each contig
#   restrict_dict : dictionnary with restriction sites count for each half of
#                   each contig
#   dup_link : list of duplicated links
#   avoid_link : list of links to avoid
# output :
#   links_dict : dictionnary where keys are link names and values are the
#               number of chimer mapped on the link
#   norm_dict : dictionnary where keys are link names and values are the
#               normalization score
def compute_links(
    mapping_file, contig_len: dict, restrict_dict: dict, dup_link=[], avoid_link=[]
):
    """
    Process links, compute the number of mapped chimeric reads and the 
    normalization score.
    """

    links_dict = {}
    norm_dict = {}

    prev_line = ""

    with open(mapping_file, "r") as f:
        for line in f:

            # At first line
            if prev_line == "":
                prev_line = line
                continue

            attrs = line.strip().split()
            prev_attrs = prev_line.strip().split()

            prev_read = prev_attrs[3].split("/")[0]
            curr_read = attrs[3].split("/")[0]

            if prev_read == curr_read and prev_attrs[0] != attrs[0]:

                # Check that the link is not among links duplicated/to avoid
                link = sorted([prev_attrs[0], attrs[0]])
                ks = link[0] + "$" + link[1]
                if ks in dup_link or ks in avoid_link:
                    continue

                # Set the middle of the mapped read as its position
                pos1 = int((int(prev_attrs[1]) + int(prev_attrs[2])) / 2)
                pos2 = int((int(attrs[1]) + int(attrs[2])) / 2)

                key = ""

                reads_contigs = sorted([attrs[0], prev_attrs[0]])
                if attrs[0] == reads_contigs[0]:
                    temp = pos1
                    pos1 = pos2
                    pos2 = temp

                len1 = contig_len[reads_contigs[0]]
                len2 = contig_len[reads_contigs[1]]

                r1 = 1 / len1 * 2
                r2 = 1 / len2 * 2

                # The ends of each contig are annotated by two tags, B and E
                # B : beginning of the contig
                # E : end of the contig
                # BE : junction beginning contig1 - end contig2
                # EB : junction end contig1 - beginning contig2
                # BB : junction beginning contig1 - beginning contig2
                # EE : junction end contig1 - end contig2

                mid1 = int(len1 / 2)
                mid2 = int(len2 / 2)

                # BB : junction beginning contig1 - beginning contig2
                if pos1 <= mid1 and pos2 <= mid2:
                    key = reads_contigs[0] + ":B$" + reads_contigs[1] + ":B"
                    norm_dict[key] = (
                        restrict_dict[reads_contigs[0]][0] * r1
                        + restrict_dict[reads_contigs[1]][0] * r2
                    )
                # BE : junction beginning contig1 - end contig2
                elif pos1 <= mid1 and pos2 > mid2:
                    key = reads_contigs[0] + ":B$" + reads_contigs[1] + ":E"
                    norm_dict[key] = (
                        restrict_dict[reads_contigs[0]][0] * r1
                        + restrict_dict[reads_contigs[1]][1] * r2
                    )
                # EB : junction end contig1 - beginning contig2
                elif pos1 > mid1 and pos2 <= mid2:
                    key = reads_contigs[0] + ":E$" + reads_contigs[1] + ":B"
                    norm_dict[key] = (
                        restrict_dict[reads_contigs[0]][1] * r1
                        + restrict_dict[reads_contigs[1]][0] * r2
                    )
                # EE : junction end contig1 - end contig2
                elif pos1 > mid1 and pos2 > mid2:
                    key = reads_contigs[0] + ":E$" + reads_contigs[1] + ":E"
                    norm_dict[key] = (
                        restrict_dict[reads_contigs[0]][1] * r1
                        + restrict_dict[reads_contigs[1]][1] * r2
                    )
                else:
                    print("ERROR : Unexpected length in contig links attribution.")
                    sys.exit()

                if key != "":
                    if key not in links_dict.keys():
                        links_dict[key] = 0
                    else:
                        links_dict[key] += 1

            prev_line = line

    return links_dict, norm_dict


# WRITE_LINKS_COUNT
# input :
#   links_dict : dictionnary where keys are link names and values are the
#               number of chimer mapped on the link
#   norm_dict : dictionnary where keys are link names and values are the
#               normalization score
#   outfile : path to output
# output :
#   output file written with the names of the linked contigs, the link score
#   and the number of mapped chimer
def write_links_count(links_dict: dict, norm_dict: dict, outfile):
    """
    Write contig links file.
    """
    output = open(outfile, "w")

    for key in links_dict:
        edge = key.split("$")
        if norm_dict[key] == 0:
            score = 0
        else:
            score = links_dict[key] / norm_dict[key]
            output.write(
                edge[0]
                + "\t"
                + edge[1]
                + "\t"
                + str(score)
                + "\t"
                + str(links_dict[key])
                + "\n"
            )

    output.close()


# ======================================================================#
#                                MAIN
# ======================================================================#

# MAKE_CONTIG_LINKS
# input :
#   mapping_data : path to Hi-C data
#   restrict_data : path to restriction sites counts per half of contigs
#   contig_len_data : path to file with contigs length
#   iteration : integer, iteration number
#   contig_links_file : path to file where links should be written
#   dup_data : path to duplicated links file
#   avoid_data : path to links to avoid
# output :
#   output written to contig_links_file with links
def make_contig_links(
    mapping_data,
    restrict_data,
    contig_len_data,
    iteration: int,
    contig_links_file,
    dup_data="abc",
    avoid_data="abc",
):
    """
    Process links from Hi-C data file.
    """

    # Load duplicate links
    dup_links = []
    if dup_data != "abc" and iteration == 1:
        dup_links = load_links(dup_data)

    # Load links to avoid
    avoid_links = []
    if int(iteration) > 1:
        avoid_links = load_links(avoid_data)

    # Read contig length file
    contig_lengths = {}
    with open(contig_len_data, "r") as cfile:
        for line in cfile:
            attrs = line.strip().split()
            contig_lengths[attrs[0]] = float(attrs[1])

    # Load restriction sites counts file
    restrict_sites_counts = {}
    with open(restrict_data, "r") as f:
        for line in f:
            attrs = line.strip().split()
            restrict_sites_counts[attrs[0]] = [int(attrs[1]), int(attrs[2])]

    contig_links = {}
    norm_score = {}

    print("Loading bedfile...")

    contig_links, norm_score = compute_links(
        mapping_file=mapping_data,
        contig_len=contig_lengths,
        restrict_dict=restrict_sites_counts,
        dup_link=dup_links,
        avoid_link=avoid_links,
    )

    print("Bedfile loaded.")

    write_links_count(contig_links, norm_score, contig_links_file)
