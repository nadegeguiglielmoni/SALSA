import pickle
import argparse
import seq_utils as sequ

# ======================= ARGUMENTS ====================================#

parser = argparse.ArgumentParser()

parser.add_argument("-a", "--cleaned", help="cleaned assembly")
parser.add_argument("-f", "--scaffold", help="final scaffold file")
parser.add_argument("-g", "--agp", help="agp file")
parser.add_argument("-p", "--map", help="pickle map of scaffolds")

args = parser.parse_args()

# ======================================================================#
scaff_map = pickle.load(open(args.map, "r"))

contig_length = {}
id2seq = {}

id2seq = sequ.parse_fasta(args.cleaned)

for key in id2seq:
    id2seq[key] = id2seq[key]
    contig_length[key] = len(id2seq[key])

# first sort scaffolds in decreasing order of length

scaff2length = {}
for scaffold in scaff_map:
    path = scaff_map[scaffold]
    length = 0
    for i in xrange(0, len(path) - 1, 2):
        length += contig_length[path[i].split(":")[0]]
    scaff2length[scaffold] = length

sorted_scaffolds = sorted(scaff2length.items(), key=lambda x: x[1], reverse=True)

c_id = 1
line = ""
agp_output = open(args.agp, "w")
ofile = open(args.scaffold, "w")

for key in sorted_scaffolds:
    # print 'scaffold_'+str(c_id) + '\t' + key
    key = key[0]
    start = 1
    local_comp = 1
    # if len(scaff_map[key]) >= 4:
    path = scaff_map[key]
    scaff_len = 0
    curr_contig = ""
    # print c_id
    line = ""
    for i in range(0, len(path) - 1, 2):
        line += "scaffold_" + str(c_id) + "\t" + str(start) + str("\t")
        curr = path[i]
        nextitem = path[i + 1]
        curr = curr.split(":")
        nextitem = nextitem.split(":")
        curr_len = contig_length[curr[0]]
        scaff_len += curr_len
        end = curr_len + start - 1
        line += str(end) + "\t"
        start = end + 1
        line += str(local_comp)
        local_comp += 1
        line += "\tW\t" + curr[0] + "\t" + "1\t"
        line += str(curr_len) + "\t"
        # print curr
        if curr[1] == "B" and nextitem[1] == "E":
            curr_contig += id2seq[curr[0]]
            line += "+\t"
        if curr[1] == "E" and nextitem[1] == "B":
            line += "-\t"
            curr_contig += sequ.rev_comp(id2seq[curr[0]])

        agp_output.write(line + "\n")
        if i != len(path) - 2:
            for j in range(0, 500):
                curr_contig += "N"
            line = "scaffold_" + str(c_id) + "\t" + str(start) + "\t"
            end = 500 + start - 1
            line += str(end) + "\t"
            start = end + 1
            line += str(local_comp) + "\t"
            local_comp += 1
            line += "N\t500\tscaffold\tyes\tna"
            agp_output.write(line + "\n")
            line = ""

    chunks = [curr_contig[i : i + 80] for i in xrange(0, len(curr_contig), 80)]
    ofile.write(">scaffold_" + str(c_id) + "\n")
    for chunk in chunks:
        ofile.write(chunk + "\n")
    c_id += 1
