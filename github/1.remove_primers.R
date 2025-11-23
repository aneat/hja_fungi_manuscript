# -----------------------------------------------------------------------------#
# Microbiome analysis workshop
# Processing raw amplicon reads
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     dada2 v 1.14.1
#                     decontam v 1.6.0
#                     phyloseq v 1.30.0
#                     purrr v 0.3.4
#                     Biostrings v 2.54.0
#                     patchwork v 1.0.1
# -----------------------------------------------------------------------------#

#############################################################
#### This script calls cutadapt to remove any primers    ####
#### You must have cutadapt installed on your system     ####
#### and present in your PATH. See cutadapt installation ####
#### documents for instructions.                         ####
#############################################################

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")

# PARSE FILE PATHS ####

# File parsing. In this step you will name and sort your forward reads/reverse reads. 
# helpful definitions:
# demultiplexed: Your reads are considered 'demultiplexed' when each read has been associated with a sample. Usually this is done by the folks at the sequencing facility.
# fasta: this is a file containing all of your nucleotides. This type of file does not have quality scores.
# fastq: this is the same file format as fasta files, but also contains quality scores associated with each read.
path <- "/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324" # CHANGE to the subdirectory containing your demultiplexed fastq files
fnFs <- sort(list.files(path, pattern = "_R1", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2", full.names = TRUE))

# CHECK FOR AND REMOVE PRIMER SITES WITH CUTADAPT ####

############# You will need to change these two values to match your data #############

# Here, you should supply the primer sequences used during PCR
FWD <- "CTTGGTCATTTAGAGGAAGTAA" # Sequence of FWD primer
REV <- "GCTGCGTTCTTCATCGATGC"  # Sequence of REV primer
# you can find the absolute file path to your cutadapt program in your terminal by entering "which cutadapt"
# and then copying the output here
######################################################################################################

# this function searches through both your forward and reverse reads and identified the primer input sequence across all of the reads.
allOrients <- function(primer) {
  
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients; REV.orients

# Prefilter to remove reads with ambiguous (N) bases ####
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# Discover primer matches, regardless of orientation ####
primerHits <- function(primer, fn) {

# Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# Run cutadapt ####
# If the following command returns an error, you do not have cutadapt installed correctly
system2("/home/busbystudents/miniconda3/envs/cutadaptenv/bin/cutadapt", args = "--version")
system2("/home/busbystudents/miniconda3/pkgs/itsx-1.1.3-hdfd78af_1/bin/ITSx", args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

nthreads <- 12
# Run Cutadapt. This will remove the primers from your sequences.
for(i in seq_along(fnFs)) {
  system2("/home/busbystudents/miniconda3/envs/cutadaptenv/bin/cutadapt", args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut[i], "-p",fnRs.cut[i],"-l", 210,  # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
}
# sanity check
# This should show no occurences in any of the orientations now
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
