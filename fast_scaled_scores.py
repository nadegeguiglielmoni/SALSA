import networkx as nx
import sys
import argparse

# ======================================================================#
#                                ARGUMENTS
# ======================================================================#

# ======================================================================#
#                                MODULES
# ======================================================================#

# LOAD_CONTIG_LINKS
# input :
#
# output :
#
def load_contig_links(filename):
    G = nx.Graph()

    links = {}

    # map for node to max edge weight for it
    node2weight = {}

    f = open(filename, "r")
    for line in f:
        line = line.strip().split()

        key = line[0] + "$" + line[1]
        links[key] = float(line[2])

        # add first node to node2weight
        if line[0] not in node2weight:
            node2weight[line[0]] = float(line[2])
        else:
            if float(line[2]) >= node2weight[line[0]]:
                node2weight[line[0]] = float(line[2])

        # add second node to node2weight
        if line[1] not in node2weight:
            node2weight[line[1]] = float(line[2])
        else:
            if float(line[2]) >= node2weight[line[1]]:
                node2weight[line[1]] = float(line[2])

        # add nodes to graph
        G.add_node(line[0])
        G.add_node(line[1])

        # add edge to graph
        G.add_edge(
            line[0],
            line[1],
            weight=float(line[2]),
            isNeighbor="False",
            isGood="False",
            links=int(line[3]),
        )
    f.close()

    return G, node2weight, links


# GET_MAX_INCIDENT
# input :
#   start
#   end
#   G
# output :
#   max
def get_max_incident(start, end, G):
    max = 0
    try:
        for e in G.successors(start):
            if e != end and G[start][e]["weight"] > max:
                max = G[start][e]["weight"]
        for e in G.predecessors(start):
            if e != end and G[e][start]["weight"] > max:
                max = G[start][e]["weight"]
    except:
        for e in G.neighbors(start):
            if e != end and G[start][e]["weight"] > max:
                max = G[start][e]["weight"]

    return max


# GET_MAX_WEIGHT
# input :
#   start
#   end
#   G
# output :
#   max
def get_max_weight(start, end, G):
    return max(get_max_incident(start, end, G), get_max_incident(end, start, G))


# PROCESS_LINKS
# input :
#   G
#   node2weight
#   links
#   outfile
# output :
#   scaled scores written to outfile
def process_links(G, node2weight, links, outfile):

    ofile = open(outfile, "w")

    for key in links:
        attrs = key.split("$")
        u = attrs[0]
        v = attrs[1]
        w = links[key]

        max_u = node2weight[u]
        max_v = node2weight[v]

        bestAlt = max(max_u, max_v)
        if bestAlt == w:
            bestAlt = get_max_weight(u, v, G)
        if bestAlt == 0:
            bestAlt = 1

        ofile.write(
            str(u)
            + "\t"
            + str(v)
            + "\t"
            + str(w)
            + "\t"
            + str(bestAlt)
            + "\t"
            + str(w / bestAlt)
            + "\t"
            + str(G[u][v]["links"])
            + "\t"
            + str(G[u][v]["isNeighbor"])
            + "\t"
            + str(G[u][v]["isGood"])
            + "\n"
        )

    ofile.close()


# ======================================================================#
#                                MAIN
# ======================================================================#


def fast_scaled_scores(infile, outfile):
    G, node_to_weight, links_data = load_contig_links(infile)
    process_links(G, node_to_weight, links_data, outfile)

