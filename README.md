# SALSA: A tool to scaffold long read assemblies with Hi-C 


This code is used to scaffold your assemblies using Hi-C data. This version implements some improvements in the original SALSA algorithm. If you want to use the old version, it can be found in the `old_salsa` branch. 

To use the latest version, first run the following commands:
```
  cd SALSA
  make
```

To run the code, you will need Python 2.7, [BOOST](http://www.boost.org/) libraries and [Networkx](https://networkx.github.io/)(version lower than 1.2).


If you consider using this tool, please cite our publication which describes the methods used for scaffolding.

Ghurye, J., Pop, M., Koren, S., Bickhart, D., & Chin, C. S. (2017). Scaffolding of long read assemblies using long range contact information. BMC genomics, 18(1), 527. [Link](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3879-z)

Ghurye, J., Rhie, A., Walenz, B.P., Schmitt, A., Selvaraj, S., Pop, M., Phillippy, A.M. and Koren, S., 2018. Integrating Hi-C links with assembly graphs for chromosome-scale assembly. bioRxiv, p.261149 [Link](https://www.biorxiv.org/content/early/2018/02/07/261149)

For any queries, please either ask on github issue page or send an email to Jay Ghurye (jayg@cs.umd.edu).
## How to run the code?

The new version of SALSA has been designed to consider several use cases depending on the input. Some assemblers output assembly graph as well along with the contig sequences. We provide options to use different information provided by the assembly to use for the scaffolding. Here is the what input options look like

```
python run_pipeline.py -h
usage: run_pipeline.py [-h] -a ASSEMBLY -l LENGTH -b BED [-o OUTPUT]
                       [-c CUTOFF] [-g GFA] -e ENZYME [-i ITER] [-x DUP]
                       [-s EXP] [-m CLEAN] [-f FILTER] [-p PRNT]

SALSA Iterative Pipeline

optional arguments:
  -h, --help            show this help message and exit
  -a ASSEMBLY, --assembly ASSEMBLY
                        Path to initial assembly
  -l LENGTH, --length LENGTH
                        Length of contigs at start
  -b BED, --bed BED     Bed file of alignments sorted by read names
  -o OUTPUT, --output OUTPUT
                        Output directory to put results
  -c CUTOFF, --cutoff CUTOFF
                        Minimum contig length to scaffold, default=1000
  -g GFA, --gfa GFA     GFA file for assembly
  -e ENZYME, --enzyme ENZYME
                        Restriction Enzyme used for experiment
  -i ITER, --iter ITER  Number of iterations to run, default = 3
  -x DUP, --dup DUP     File containing duplicated contig information
  -s EXP, --exp EXP     Expected Genome size of the assembled genome
  -m CLEAN, --clean CLEAN
                        Set this option to "yes" if you want to find
                        misassemblies in input assembly
  -f FILTER, --filter FILTER
                        Filter bed file for contigs present in the assembly
  -p PRNT, --prnt PRNT  Set this option to "yes" if you want to output the
                        scaffolds sequence and agp file for each iteration
```

### Mapping Reads

To start the scaffolding, first step is to map reads to the assembly. We recommend using BWA or BOWTIE2 aligner to map reads. The read mapping generates a `bam` file. SALSA requires `bed` file as the input. This can be done using the `bamToBed` command from the [Bedtools](http://bedtools.readthedocs.io/en/latest/) package. Also, SALSA requires bed file to be sored by the read name, rather than the alignment coordinates. Once you have bam file, you can run following commands to get the bam
file needed as an input to SALSA.

Since Hi-C reads and alignments contain experimental artifacts, the alignments needs some postprocessing. To align and postprocess the alignments, you can use the pipeline released by Arima Genomics which can be found here [ https://github.com/ArimaGenomics](https://github.com/ArimaGenomics).

```
bamToBed -i alignment.bam > alignment.bed
sort -k 4 alignment.bed > tmp && mv tmp alignment.bed
```

### Generating contig lengths file

SALSA requires contig lengths as an input. You can generate this file using this command on your contig sequence file.
```
samtools faidx contigs.fasta
```

This will generate `contigs.fasta.fai` as an input for `-l` option. You can use this file for contig lengths while running SALSA.

### Restriction Enzyme input

Hi-C experiments can use different restriction enzymes. We use the enzyme frequency in contigs to normalize the Hi-C interaction frequency. You will need to specify which enzyme was used for Hi-C experiment while running SALSA in `-e` option. If multiple enzymes were used, they can specified by separating with comma without space, like `-e GATC,AAGCTT`. Note that you need to specify the actual sequence of the cutting site for a restriction enzyme and not the enzyme name. For example, if you use `MboI` in the Hi-C protocol ,then you would specify it as `-e GATC`.


### 1) I have contig sequences and the alignment bam file
This is the minimum input you will require Suppose you only have contig sequences generated. Once you prepare the bed file as described above, the code can be run as follows:
```
python run_pipeline.py -a contigs.fasta -l contigs.fasta.fai -b alignment.bed -e {Your Enzyme} -o scaffolds 
```

### 2) I have contig sequences and the alignment bam file but also want to use Hi-C data to correct input assembly errors

We also implemented a method in SALSA that can correct some of the errors in the assembly with Hi-C data. To use this method, you need to run following
```
python run_pipeline.py -a contigs.fasta -l contigs.fasta.fai -b alignment.bed -e {Your Enzyme} -o scaffolds -m yes
```

If you want to know what were the locations in the contigs where SALSA found errors, you can look at the `input_breaks` file in the output directory.

### 3) I want to use the assembly graph to assist scaffolding 

Some assembles output gfa file for the assembly graph. You can use that as an input for SALSA as follows

```
python run_pipeline.py -a contigs.fasta -l contigs.fasta.fai -b alignment.bed -e {Your Enzyme} -o scaffolds -m yes -g contigs_graph.gfa
```

We utilize graph to guide the scaffolding, which in turn reduces the errors.

## How to interpret the output?

SALSA generates a bunch of files in the output folder. SALSA is an iterative algorithm, so it generates files for each iteration. The files you will be interested in are `scaffolds_FINAL.fasta`, which contains the sequences for the scaffolds generated by the algorithm. Another file which is of interest is `scaffolds_FINAL.agp`, which is the [agp](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/) style output for the scaffolds describing the assignment, orientation and ordering of contigs along the scaffolds.

These are the postprocessing options which SALSA provides.

### I want scaffolds sequences for all iterations rather than just final scaffolds

To dig more into the output, we have added an option that can output scaffolds along with the agp file for all intermediate iterations. This option is usually helpful in debugging and exploring the errors. Here is how you can run it:

```
 python run_pipeline.py -a unitigs.fasta -l unitigs.fasta.fai -b alignment.bed -e {Your Enzyme} -o scaffolds -m yes -g unitigs_graph.gfa -u    unitigs_tiling.bed -p yes
```

