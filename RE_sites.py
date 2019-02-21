
import re 
import argparse
import get_seq

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--assembly", help="assembled contigs")
    parser.add_argument("-e", "--enzyme", help="restriction enzyme")
    #parser.add_argument("-m", "--mapping", help="mapping of read to contigs in bam format")
    #parser.add_argument("-d", "--dir", help="output directory for results",default='out')
    args = parser.parse_args()
    f = get_seq.parse_fasta(args.assembly)
    enzymes_input = args.enzyme.replace(' ','').split(',')
    final_enzymes = []
    for each in enzymes_input:
        if 'N' in each:
            final_enzymes.append(each.replace('N','G'))
            final_enzymes.append(each.replace('N','A'))
            final_enzymes.append(each.replace('N','T'))
            final_enzymes.append(each.replace('N','C'))
        else:
            final_enzymes.append(each)

    for key in f:
        
        id,seq = key, f[key]
        left_count = 0
        rigt_count = 0
        for enzyme in final_enzymes:
            pos  = [m.start(0) for m in re.finditer(enzyme,seq)]
         
            length = len(seq)    
            for each in pos:
            	if each < length/2:
            		left_count += 1
            	else:
            		rigt_count += 1

        print("{0}\t{1}\t{2}\n".format(id, left_count, rigt_count))

    

if __name__ == '__main__':
    main()
