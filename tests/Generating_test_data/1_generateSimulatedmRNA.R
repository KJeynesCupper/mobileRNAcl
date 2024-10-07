## make a simated reads from my merged fasta file..
## 07th October 2024
library(polyester)
library(Biostrings)




# FASTA annotation
fasta_file = "~/PACKAGES/mobileRNAcl/tests/Genomes/mergedgenometest3.fa"
small_fasta = readDNAStringSet(fasta_file)


# set up transcript-by-timepoint matrix:
num_timepoints = 4
countmat = matrix(readspertx, nrow=length(small_fasta), ncol=num_timepoints)

# add spikes in expression at certain timepoints to certain transcripts:
up_early = c(1,2)
up_late = c(3,4)
countmat[up_early, 2] = 3*countmat[up_early, 2]
countmat[up_early, 3] = round(1.5*countmat[up_early, 3])

# simulate reads:
simulate_experiment_countmat(fasta_file, readmat=countmat,
                             outdir='~/PACKAGES/mobileRNAcl/tests/mRNAtestdata')
