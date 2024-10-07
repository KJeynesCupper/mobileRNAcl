mobileRNA command-line <a ><img src="docs/logo.png" align="right" height="138" /></a>
======================================================================
***UNDER DEVELOPMENT*** <br>
This is a command-line package to undertake the merging of genome reference, alignment and clustering dircecting from the command-line instead of
using R as an interface. Please note that this is still under development. 

This draft works on a Mac M1 chip computer, but has not been tested on other systems. 

<br> 

This is an extention of the R package, mobileRNA, where here we offer
a command-line package that undertakes the pre-processing steps via 
the command-line rather than via the R. This includes:

1. Merging genome references and/or genome annotations
2. Undertaking the alignement step & clustering for mRNA or sRNA sequencing reads

<br>

Author
======================================================================
Katie Jeynes-Cupper, University of Illinois Urbana-Champaign, kejc@illinois.edu
<br>
Marco Catoni, University of Birmingham

<br>


Table of Contents
======================================================================

-   [Installation](#installation)
-   [Installation of OS dependencies](#Installation-of-OS-dependencies)
-   [Basic Usage](#Basic-Usage)
-   [Merging genome reference assemblies and/or genome annotations](#Merging-genome-reference-assemblies-and/or-genome-annotations)
-   [Small RNA sequencing read alignment and clustering](#Small-RNA-sequencing-read-alignment-and-clustering)
-   [Messenger RNA sequencing read alignment and clustering  ](#Messenger-RNA-sequencing-read-alignment-and-clustering)
-   [Issues](#Issues)

  
<br>


Installation 
======================================================================
1. mobileRNAcl can be cloned from github:

``` bash
https://github.com/KJeynesCupper/mobileRNAcl

```

2. Change directory to where mobileRNAcl was cloned, and make the script executable by running:

``` bash
chmod +x mobileRNAcl.sh
mv mobileRNAcl.sh mobileRNAcl

```

3. Place the script in a directory that is included in your PATH environment variable for easy access, or use the full path to call the script.

``` bash
export PATH="/path/to/mobileRNAcl/folder:$PATH"
```
<br>

4. Re-open command line

``` bash
source ~/.bashrc   # or ~/.zshrc if you're using Zsh
```

5. Run package, to help page

``` bash
mobileRNAcl --help
```

<br>

Installation of OS dependencies 
======================================================================
Please install the following OS dependencies within a Conda environment: 

For sRNA data, `ShortStack` (>= 4.0.0) (Axtell 2013) is requires. Please consider 
that `ShortStack` is not available for Windows, hence, Windows users will either 
need to opt to use a virtual machine or 
[Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10)
`ShortStack` will need to be installed and used on 
the Linux side. Please head to [ShortStack](https://github.com/MikeAxtell/ShortStack#install-using-conda-recommended) 
to see the recommended installation instructions with Conda.
This will ensure all dependencies are available within the same environment. 

For mRNA data, [HISAT2](https://anaconda.org/bioconda/hisat2) (Kim 2015), 
[HTSeq](https://htseq.readthedocs.io/en/master/install.html) 
(Anders, Pyl, and Huber 2014), 
[SAMtools](https://anaconda.org/bioconda/samtools) (Danecek P 2021) are required 
within the same Conda environment [@Anaconda].

<br>

Basic Usage
======================================================================
The general syntax for using mobileRNAcl is:

```bash
mobileRNAcl [function] [options...]

    function:
	RNAmergeGenomes			merge genome references
	RNAmergeAnnotations	    	merge genome annoations
	map_sRNA			Map and cluster small RNA sequencing reads
	map_mRNA			Map and cluster messenger RNA sequencing reads
    	--help                  	Display this help message
```


<br>

Merging genome reference assemblies and/or genome annotations
======================================================================
The two function to generated either a merged genome reference or annotation function in a similar way. 
They can handle any number of input files and added the prefixes in the order supplied, in an alphabetical order.
For instance, this means that for the first input file the identifying prefix is "A", the second input file is "B" 
and so on. 


<br>

## 1. RNAmergeGenomes

Merge the any number of FASTA genome assembly files into a final merged genome reference. 
Each individual genome reference will have an ID attached to the chromosome names to 
make it identifable. 

```bash

mobileRNAcl RNAmergeGenomes -i ./path/to/genome1.fasta ./path/to/genome2.fasta -o ./path/to/output.fasta

```
note: users can compress the output (gz) but simply adding the extention to the output file name. 
(for example, ./path/to/output.fasta.gz )
<br>

## 2. RNAmergeAnnotations

Merge the any number of GFF genome annotation files into a final merged annotation. 
Each individual genome annotation will have an ID attached to the chromosome names to 
make it identifable. 

```bash

mobileRNAcl RNAmergeAnnotations -i ./path/to/genome1.gff ./path/to/genome2.gff -o ./path/to/output.gff

```

<br>


Small RNA sequencing read alignment and clustering  
======================================================================
## method 
The pipeline undertakes de novo detection of sRNA-producing loci and alignment, where the output of each are stored in their respective folders in 
the users desired location. The de novo detection of sRNA-producing loci analyses each sample to identify de novo sRNA-producing loci (ie. sRNA
clusters), and joins these results into a single file called "locifile.txt". The alignment step aligns and clusters each sample to the genome reference
long with the file containing the de novo sRNA clusters. The final reports are imported into R using [RNAimport()].
Input RNA sequencing files can be gzipped. 


## usage
Here is a basic usage with default settings:

```bash
mobileRNAcl map_sRNA -FASTA ./path/to/fasta -index ./path/to/index -i ./path/to/fastq -o ./path/to/output 

```

<br>

## other options
- `-f`      merged genome reference (FASTA), path to a FASTA genome reference file. 
- `-x`      bowtie genome reference index, path to save index
- `-i`          input files (FASTQ), path to sequencing reads  
- `-o`          output location of results, directory to store output. 
- `-t`    threads (default = 4), set the number of threads to use where more threads means a faster completion time. 
- `-p`        pad (default = 200), initial peaks are merged if they are this distance or less from each other. Must >= 1
- `-m`     mincov (default = 0.5),minimum alignment depth, in units of reads per million, required to nucleate a small RNA cluster during de novo cluster search. Must be a number > 0.
- `-d`   dicermin (default = 20), the minimum size in nucleotides of a valid small RNA. This option sets the bounds to discriminate dicer-derived small RNA loci from other loci.
- `-c`   dicermax (default = 24), the minimum size in nucleotides of a valid small RNA. This option sets the bounds to discriminate dicer-derived small RNA loci from other loci.
- `-n`   dn_mirna (defalt = FALSE), activates a de novo comprehensive genome-wide search for miRNA loci


<br>


## output files 

Wihtin the desired output directory, the function generates:

- 1_de_novo_detection: Stores output from the detection of de novo sRNA-producing loci
- locifile.txt
- 2_sRNA_results: Stores results 


In the frist folder, `1_de_novo_detection`, there is a folder called `1_alignment` which stored the 
alignment file (BAM) for each sample. The `1_de_novo_detection` folder then stores
one folder per sample which holds the results for the respective 
samples *de novo* sRNA analysis.

The `locifile.txt` contains all detected sRNA-producing genes from the experimental design. 

Lastly, the `2_sRNA_results` folder stores final clustering results for each sample,
and as before the results of each sample are stored within it's respective folder. 
The 'Results.txt' file imported into R using using [RNAimport()] in the `mobileRNA` R package. 



Messenger RNA sequencing read alignment and clustering  
======================================================================
## method 
The function invokes a number of OS commands, and is dependent on the installation of `HISAT2`,`HTSeq` and `SAMtools` with `Conda.` 
The pipeline can undertake single- or pair-end analysis, and to do so requires a data frame stating the sample information where each row 
represents a sample. The reads are mapped using `HISAT` and then the raw counts are estimated by `htseq-count`. The output alignment file (BAM) and 
raw counts file for each sample are stored within the samples own folder within the desired directory. 



## usage
Here is a basic usage with default settings:

```bash
mobileRNAcl map_mRNA -FASTA ./path/to/fasta -index ./path/to/index -GFF ./path/to/gff -i ./path/to/fastq -o ./path/to/output 

```

<br>

## other options
- `-f`      merged genome reference (FASTA), path to a FASTA genome reference file. 
- `-x`      bowtie genome reference index, path to save index
- `-g`      merged genome annotation (GFF), path to a GFF genome reference file. 
- `-i`      input files (FASTQ), path to sequencing reads  
- `-o`     	output location of results, directory to store output. 
- `-t`    	threads (default = 6), set the number of threads to use where more threads means a faster completion time. 
- `-p`     	paired, is the data pair-end (default = FALSE)
- `-f`     	format, format of alignment files (default = bam)
- `-a`      a, minaqual, skips all reads with a MAPQ alignment quality lower than the given value (default: 0).
- `-o`      order, the alignment file is sorted by read name or by alignment position (default = pos).
- `-s`  	stranded, define whether the RNAseq data is strand-specific. Choose from yes, no, reverse (default = no)
- `-m`      mode, how to handle reads that overlap with more than one feature. Choose from union, intersection-strict or intersection-nonempty (default = union)
- `-n`  	nonunique, how to handle reads which aligned to or are assigned to more than one feature (default = none)
- `-j`      type, the feature type defined by the 3rd column in GFF file (default = mRNA)
- `-d`     	idattr, the attribute to be used as feature ID from 9th column in GFF (default = Name)
- `-y`		python, the python version installed, for example python or python3. (default = python)

**IMPORTANT:** This function is relies on using the following naming system for your files. 
If the data is single-ended, ensure that bother the forward and reverse files have the same name and end with "_1" or "_2". 
For example, sampleA_1.fastq & sampleA_2.fastq 
Whereas, if you have pair-ended dated, please ensure you are labeling as such:
- sampleA_L1_1.fastq
- sampleA_L1_2.fastq
- sampleA_L2_1.fastq
- sampleA_L2_1.fastq
Please notice the use of *_L1_1*, *_L1_2* and *_L2_1*, *_L2_2* - these are the key labels to utilise before the file extention. 
Input RNA sequencing files can be gzipped. 

<br>

## output files 
For mRNA analysis, the fucntion generates one folder for each sample in the desired location. 
For each sample, their folder stores several files:

- [samplename]_AllReads.bam - The alignemnt file for all reads
- [samplename]_unique_sorted.bam - The sorted, uniquely aligned reads
- [samplename]_unique_index.bam - The sorted and indexed, uniquely aligned reads
- Results.txt - Raw count results to be imported into R using using [RNAimport()] in the `mobileRNA` R package. 

<br>


Testing data   
======================================================================
You can test the method using the test data (examples below), however know that these 
are truncated files and are not a true representation of real data. In fact expect that no biological results 
from the sRNA or mRNA data - whereby in the Results.txt file show results but the values are all zero. 

Make sure to direct the path to the mobileRNAcl location, for example


```bash

location=./Users/KJeynesCupper/mobileRNAcl

```


## Test help 
```bash
mobileRNAcl RNAmergeGenomes --help
```
## Test merging genome
```bash
mobileRNAcl RNAmergeGenomes -i $location/tests/Genomes/reduced_chr2_Tomato.fa.gz $location/tests/Genomes/reduced_chr12_Eggplant.fa.gz -o $location/tests/Genomes/mergedgenometest.fa
mobileRNAcl RNAmergeAnnotations -i $location/tests/Genomes/reduced_chr2_Tomato.gff.gz $location/tests/Genomes/reduced_chr12_Eggplant.gff.gz -o $location/tests/Genomes/mergedannotationtest.gff

``` 

## Test mapping sRNA -- data taken from mobileRNA package
```bash
mobileRNAcl map_sRNA -f $location/tests/Genomes/mergedgenometest3.fa \
-x $location/tests/Genomes/mergedgenometest3 \
-i $location/tests/sRNAtestdata/ \
-o $location/tests/sRNAtestoutput/ 
``` 

## Test mapping mRNA 
```bash
mobileRNAcl map_mRNA -f $location/tests/Genomes/mergedgenometest3.fa \
-x $location/tests/Genomes/mergedgenometest3 \
-g $location/tests/Genomes/mergedgenometest3.gff \
-i $location/tests/mRNAtestdata/ \
-o $location/tests/mRNAtestoutput/ 
```

Issues
======================================================================
Please post issues, comments, bug reports, questions, etc. to the project github page at <https://github.com/KJeynesCupper/mobileRNAcl>.
