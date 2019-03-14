import os
import argparse
import sys
import subprocess
from subprocess import Popen, PIPE
import seq_utils as sequ
import digest
import make_links
import fast_scaled_scores as fss


def check(path):
    with open(path, "r") as f:
        for line in f:
            attrs = line.split()
            if float(attrs[4]) >= 1:
                return False
            else:
                return True


ng50 = []


def NG50(lengths, sz=0):
    genome_size = sz
    if genome_size == 0:
        genome_size = sum(lengths.values())
    contig_lengths = sorted(lengths.values(), reverse=True)
    lensum = 0
    for i in range(len(contig_lengths)):
        lensum += contig_lengths[i]
        if lensum >= genome_size / 2:
            return contig_lengths[i]


def main():

    workdir = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(description="SALSA Iterative Pipeline")

    parser.add_argument(
        "-a", "--assembly", help="Path to initial assembly", required=True
    )
    parser.add_argument(
        "-b", "--bed", help="Bed file of alignments sorted by read names", required=True
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output directory to put results",
        required=False,
        default="SALSA_output",
    )
    parser.add_argument(
        "-c",
        "--cutoff",
        help="Minimum contig length to scaffold, default=1000",
        required=False,
        default=1000,
    )
    parser.add_argument(
        "-g", "--gfa", help="GFA file for assembly", required=False, default="abc"
    )
    parser.add_argument(
        "-e", "--enzyme", help="Restriction Enzyme used for experiment", required=True
    )
    parser.add_argument(
        "-i",
        "--iter",
        help="Number of iterations to run, default = 3",
        required=False,
        default=3,
    )
    parser.add_argument(
        "-x",
        "--dup",
        help="File containing duplicated contig information",
        required=False,
        default="abc",
    )
    parser.add_argument(
        "-s",
        "--exp",
        help="Expected Genome size of the assembled genome",
        required=False,
        default=0,
    )
    parser.add_argument(
        "-m",
        "--clean",
        help='Set this option to "yes" if you want to find '
        "misassemblies in input assembly",
        required=False,
        default="yes",
    )
    parser.add_argument(
        "-f",
        "--filter",
        help="Filter bed file for contigs present in the assembly",
        required=False,
        default="no",
    )
    parser.add_argument(
        "-p",
        "--prnt",
        help="Set this option to 'yes' if you want to output the "
        "scaffolds sequence and gap file for each iteration",
        required=False,
        default="no",
    )

    # CURRENTLY NOT USED
    # parser.add_argument('-u', '--unitigs',
    #     help='The tiling of unitigs to contigs in bed format',
    #     required=False,
    #     default='abc')
    # parser.add_argument('-t', '--tenx',
    #     help='10x links tab separated file, sorted by last columnls',
    #     required=False,
    #     default='abc')
    # parser.add_argument("-d", "--dist",
    #     help="Maximum distance between pairs to consider for misassembly detection",
    #     required=False,
    #     default=2000000)

    args = parser.parse_args()

    iter_max = int(args.iter)

    # iteration counter
    iter_num = 1

    genome_size = int(args.exp)

    # ----------------------- Handle files and directories -----------------------#

    if not os.path.exists(args.assembly):
        print("ERROR : Could not find the assembly file {0}".format(args.assembly))
        sys.exit(1)
    if not os.path.exists(args.bed):
        print("ERROR : Could not find the Hi-C data bed file {0}".format(args.bed))
        sys.exit(1)

    # create output directory if it does not exist
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    log = open(args.output + "/commands.log", "w", 1)

    # ----------------------- Generate length file -------------------------------#

    sequ.make_seq_length_file(
        assembly=args.assembly,
        output_file=args.output + "/scaffold_length_iteration_1",
        namelist_file=args.output + "/contig_names.txt",
    )

    # ----------------------- Filtering bed file ---------------------------------#

    if args.filter == "yes":
        # filter Hi-C reads for contigs present in the assembly
        os.system(
            "grep -f "
            + args.output
            + "/contig_names.txt -w  "
            + args.bed
            + " > "
            + args.output
            + "/alignment_iteration_1.bed"
        )
    else:
        os.symlink(
            os.path.abspath(args.bed), args.output + "/alignment_iteration_1.bed"
        )

    os.system(
        "ln -s "
        + os.path.abspath(args.assembly)
        + " "
        + args.output
        + "/assembly.cleaned.fasta"
    )

    # ----------------------- FIND MISASSEMBLIES ---------------------------------#

    if args.clean == "yes":

        cmd = (
            workdir
            + "/break_contigs_start -a "
            + args.output
            + "/alignment_iteration_1.bed -l "
            + args.output
            + "/scaffold_length_iteration_1 > "
            + args.output
            + "/input_breaks -s 100"
        )
        log.write(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as err:
            print(
                "ERROR : Could not run break_contigs_start to detect misassemblies."
                " Will continue without detecting misassemblies."
            )

        cmd = (
            "python "
            + workdir
            + "/correct.py  "
            + args.assembly
            + " "
            + args.output
            + "/input_breaks "
            + args.output
            + "/alignment_iteration_1.bed "
            + args.output
        )
        log.write(cmd)
        try:
            p = subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as err:
            print("Could not run module 'correct.'")
            sys.exit(1)

        os.system(
            "mv "
            + args.output
            + "/alignment_iteration_1.tmp.bed "
            + args.output
            + "/alignment_iteration_1.bed"
        )
        os.system(
            "mv "
            + args.output
            + "/asm.cleaned.fasta "
            + args.output
            + "/assembly.cleaned.fasta"
        )

    print("Digesting genome.")
    # ------------- GET RESTRICTION ENZYME CUTTING SITES ---------------#

    digest.generate_digested_fragments(
        assembly=args.output + "/assembly.cleaned.fasta",
        enzyme=args.enzyme,
        outname=args.output + "/re_counts_iteration_" + str(iter_num),
    )

    # ------------- MAKE LINKS -----------------------------------------#

    print("Starting iteration {0}".format(iter_num))

    make_links.make_contig_links(
        mapping_data=args.output + "/alignment_iteration_" + str(iter_num) + ".bed",
        restrict_data=args.output + "/re_counts_iteration_" + str(iter_num),
        contig_len_data=args.output + "/scaffold_length_iteration_" + str(iter_num),
        iteration=iter_num,
        contig_links_file=args.output + "/contig_links_iteration_" + str(iter_num),
        dup_data=args.dup,
    )

    # ------------- FAST SCALED SCORES ---------------------------------#

    # now use Serge's code to calculate

    fss.fast_scaled_scores(
        infile=args.output + "/contig_links_iteration_" + str(iter_num),
        outfile=args.output + "/contig_links_scaled_iteration_" + str(iter_num),
    )

    # ------------- SORT LINKS -----------------------------------------#
    # Sort the links by column 5
    if not os.path.isfile(
        args.output + "/contig_links_scaled_sorted_iteration_" + str(iter_num)
    ):
        try:
            cmd = (
                "sort -k 5 -gr "
                + args.output
                + "/contig_links_scaled_iteration_"
                + str(iter_num)
                + " > "
                + args.output
                + "/contig_links_scaled_sorted_iteration_"
                + str(iter_num)
            )
            log.write(cmd + "\n")

            p = subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as err:
            if os.path.isfile(
                args.output + "/contig_links_scaled_sorted_iteration_" + str(iter_num)
            ):
                os.system(
                    "rm "
                    + args.output
                    + "/contig_links_scaled_sorted_iteration_"
                    + str(iter_num)
                )
            print("ERROR : Could not run sort.")
            sys.exit(1)

            # ------------- LOAD GFA -------------------------------------------#
    if args.gfa != "abc" and not os.path.isfile(args.output + "/tmp.links"):
        try:
            cmd = (
                workdir
                + "/correct_links -g "
                + args.gfa
                + " -l "
                + args.output
                + "/contig_links_scaled_sorted_iteration_1 > "
                + args.output
                + "/tmp.links"
            )
            log.write(cmd + "\n")
            p = subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as err:
            print("ERROR : Could not run correct_links.")
            sys.exit(1)

        os.system(
            "mv "
            + args.output
            + "/tmp.links "
            + args.output
            + "/contig_links_scaled_sorted_iteration_1"
        )

        # ------------- LAYOUT UNITIGS -------------------------------------#
    if not os.path.isfile(args.output + "/scaffolds_iteration_1.p"):
        try:
            cmd = (
                "python "
                + workdir
                + "/layout_unitigs.py -x "
                + args.gfa
                + " -l "
                + args.output
                + "/contig_links_scaled_sorted_iteration_"
                + str(iter_num)
                + " -c "
                + str(args.cutoff)
                + " -i 1 -d "
                + args.output
            )
            log.write(cmd + "\n")

            p = subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as err:
            if os.path.isfile(args.output + "scaffolds_iteration_1.p"):
                os.system("rm " + args.output + "scaffolds_iteration_1.p")
            print("ERROR : Could not run module layout_unitigs.")
            sys.exit(1)

            # ------------- BREAK CONTIGS --------------------------------------#
    if not os.path.isfile(
        args.output + "/misasm_iteration_" + str(iter_num + 1) + ".report"
    ):
        try:
            cmd = (
                workdir
                + "/break_contigs -a "
                + args.output
                + "/alignment_iteration_"
                + str(iter_num + 1)
                + ".bed -b "
                + args.output
                + "/breakpoints_iteration_"
                + str(iter_num + 1)
                + ".txt -l "
                + args.output
                + "/scaffold_length_iteration_"
                + str(iter_num + 1)
                + " -i "
                + str(iter_num + 1)
                + " -s 100   > "
                + args.output
                + "/misasm_iteration_"
                + str(iter_num + 1)
                + ".report"
            )
            p = subprocess.check_output(cmd, shell=True)

            log.write(cmd + "\n")
        except subprocess.CalledProcessError as err:
            print("ERROR : Could not run break_contigs.")
            sys.exit(1)

            # ------------- REFACTOR BREAKS ------------------------------------#
    if not os.path.isfile(args.output + "/misasm_" + str(iter_num + 1) + ".DONE"):
        cmd = (
            "python "
            + workdir
            + "/refactor_breaks.py -d "
            + args.output
            + " -i "
            + str(iter_num + 1)
        )
        try:
            p = subprocess.check_output(cmd, shell=True)

            log.write(cmd + "\n")
        except subprocess.CalledProcessError as err:
            print("ERROR : Could not run module refactor_breaks.")
            sys.exit(1)

            # if it is required to output scaffolds and gaps at each iteration
    if args.prnt == "yes":
        cmd = (
            "python "
            + workdir
            + "/get_seq.py -a "
            + args.output
            + "/assembly.cleaned.fasta -f "
            + args.output
            + "/scaffolds_ITERATION_"
            + str(iter_num)
            + ".fasta -g "
            + args.output
            + "/scaffolds_ITERATION_"
            + str(iter_num)
            + ".agp -p "
            + args.output
            + "/scaffolds_iteration_"
            + str(iter_num)
            + ".p"
        )
        log.write(cmd + "\n")
        try:
            p = subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError as err:
            print(
                "ERROR : Could not run module get_seq to output " "scaffolds and gaps."
            )

    iter_num += 1

    scaffold_length = {}

    with open(args.output + "/scaffold_length_iteration_" + str(iter_num), "r") as f:
        for line in f:
            attrs = line.split()
            scaffold_length[attrs[0]] = int(attrs[1])
        ng50.append(NG50(scaffold_length, genome_size))

    if iter_num - 1 == int(args.iter):
        cmd = (
            "python "
            + workdir
            + "/get_seq.py -a "
            + args.output
            + "/assembly.cleaned.fasta -f "
            + args.output
            + "/scaffolds_FINAL.fasta -g "
            + args.output
            + "/scaffolds_FINAL.agp -p "
            + args.output
            + "/scaffolds_iteration_"
            + str(args.iter)
            + ".p"
        )
        log.write(cmd + "\n")
        os.system(cmd)
        sys.exit(0)

        # now do iterative
    while iter_num <= iter_max:

        print("... Starting iteration {0}".format(iter_num))

        make_links.make_contig_links(
            mapping_data=args.output + "/alignment_iteration_" + str(iter_num) + ".bed",
            restrict_data=args.output + "/re_counts_iteration_" + str(iter_num),
            contig_len_data=args.output + "/scaffold_length_iteration_" + str(iter_num),
            iteration=iter_num,
            contig_links_file=args.output + "/contig_links_iteration_" + str(iter_num),
            dup_data=args.dup,
        )

        fss.fast_scaled_scores(
            infile=args.output + "/contig_links_iteration_" + str(iter_num),
            outfile=args.output + "/contig_links_scaled_iteration_" + str(iter_num),
        )

        if not os.path.isfile(
            args.output + "/contig_links_scaled_sorted_iteration_" + str(iter_num)
        ):
            try:
                cmd = (
                    "sort -k 5 -gr "
                    + args.output
                    + "/contig_links_scaled_iteration_"
                    + str(iter_num)
                    + " > "
                    + args.output
                    + "/contig_links_scaled_sorted_iteration_"
                    + str(iter_num)
                )
                log.write(cmd + "\n")
                p = subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError as err:
                if os.path.isfile(
                    args.output + "/new_links_scaled_sorted_iteration_" + str(iter_num)
                ):
                    os.system(
                        "rm "
                        + args.output
                        + "/new_links_scaled_sorted_iteration_"
                        + str(iter_num)
                    )
                print("ERROR : Could not run sort.")
                sys.exit(1)

                # NOW check if any useful link here
        if check(
            args.output + "/contig_links_scaled_sorted_iteration_" + str(iter_num)
        ):
            # FOR FUCK SAKE, WHO THOUGHT THIS WAS A GOOD IDEA ???
            print("No more useful links to use :/ ")
            break

        if not os.path.isfile(
            args.output + "/scaffolds_iteration_" + str(iter_num) + ".p"
        ):
            try:
                cmd = (
                    "python "
                    + workdir
                    + "/layout_unitigs.py -x abc -l "
                    + args.output
                    + "/contig_links_scaled_sorted_iteration_"
                    + str(iter_num)
                    + " -c "
                    + str(args.cutoff)
                    + " -i "
                    + str(iter_num)
                    + " -d "
                    + args.output
                )
                log.write(cmd + "\n")
                p = subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError as err:
                if os.path.isfile(
                    args.output + "scaffolds_iteration_" + str(iter_num) + ".p"
                ):
                    os.system(
                        "rm "
                        + args.output
                        + "scaffolds_iteration_"
                        + str(iter_num)
                        + ".p"
                    )
                print("ERROR : Could not run module layout_unitigs.")
                sys.exit(1)

        if not os.path.isfile(
            args.output + "/misasm_iteration_" + str(iter_num + 1) + ".report"
        ):
            try:
                cmd = (
                    workdir
                    + "/break_contigs -a "
                    + args.output
                    + "/alignment_iteration_"
                    + str(iter_num + 1)
                    + ".bed -b "
                    + args.output
                    + "/breakpoints_iteration_"
                    + str(iter_num + 1)
                    + ".txt -l "
                    + args.output
                    + "/scaffold_length_iteration_"
                    + str(iter_num + 1)
                    + " -i "
                    + str(iter_num + 1)
                    + " -s 100  > "
                    + args.output
                    + "/misasm_iteration_"
                    + str(iter_num + 1)
                    + ".report"
                )
                log.write(cmd + "\n")
                p = subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError as err:
                print("ERROR : Could not run break_contigs.")
                sys.exit(1)

        if not os.path.isfile(args.output + "/misasm_" + str(iter_num + 1) + ".DONE"):
            cmd = (
                "python "
                + workdir
                + "/refactor_breaks.py -d "
                + args.output
                + " -i "
                + str(iter_num + 1)
                + " > "
                + args.output
                + "/misasm_"
                + str(iter_num + 1)
                + ".log"
            )
            try:
                log.write(cmd + "\n")
                p = subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError as err:
                print("ERROR : Could not run module refactor_breaks.")
                sys.exit(1)

        if args.prnt == "yes":
            cmd = (
                "python "
                + workdir
                + "/get_seq.py -a "
                + args.output
                + "/assembly.cleaned.fasta -f "
                + args.output
                + "/scaffolds_ITERATION_"
                + str(iter_num)
                + ".fasta -g "
                + args.output
                + "/scaffolds_ITERATION_"
                + str(iter_num)
                + ".agp -p "
                + args.output
                + "/scaffolds_iteration_"
                + str(iter_num)
                + ".p"
            )
            log.write(cmd + "\n")
            try:
                p = subprocess.check_output(cmd, shell=True)
            except subprocess.CalledProcessError as err:
                print(
                    "ERROR : Could not run module get_seq to output "
                    "scaffolds and gaps."
                )

        scaffold_length = {}

        with open(
            args.output + "/scaffold_length_iteration_" + str(iter_num + 1), "r"
        ) as f:
            for line in f:
                attrs = line.split()
                scaffold_length[attrs[0]] = int(attrs[1])
            ng50.append(NG50(scaffold_length, genome_size))
            curr_sz = len(ng50)

            if ng50[curr_sz - 1] == ng50[curr_sz - 2]:
                cmd = (
                    "python "
                    + workdir
                    + "/get_seq.py -a "
                    + args.output
                    + "/assembly.cleaned.fasta -f "
                    + args.output
                    + "/scaffolds_FINAL.fasta -g "
                    + args.output
                    + "/scaffolds_FINAL.agp -p "
                    + args.output
                    + "/scaffolds_iteration_"
                    + str(iter_num - 1)
                    + ".p"
                )
                log.write(cmd + "\n")
                os.system(cmd)
                sys.exit(0)

        if iter_num - 1 == int(args.iter):
            cmd = (
                "python "
                + workdir
                + "/get_seq.py -a "
                + args.output
                + "/assembly.cleaned.fasta -f "
                + args.output
                + "/scaffolds_FINAL.fasta -g "
                + args.output
                + "/scaffolds_FINAL.agp -p "
                + args.output
                + "/scaffolds_iteration_"
                + str(args.iter)
                + ".p"
            )
            log.write(cmd + "\n")
            os.system(cmd)
            sys.exit(0)

        iter_num += 1

    log.close()


if __name__ == "__main__":
    main()
