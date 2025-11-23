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

#################################################################################
#                               Main workflow                                   #
# Filter and trim, denoise, sample inferrence, chimera and contaminant removal, # 
# taxonomic assignment, combine sequence table and metadata                     #
#################################################################################

# PACKAGES, SCRIPTS, AND SETUP ####

# why each package (put in onboarding document)
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(patchwork); packageVersion("patchwork")
library(ShortRead)


# PARSE FILE PATHS ####

# File parsing - 
path <- "/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/cutadapt" # CHANGE to the directory containing your demultiplexed fastq files when using your own data
filtpath <- file.path(path, "filtered.dada.1124") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present

nthreads <- 12

# Here, the pattern is set to match the forward reads tagged as "...pass_1.fastq.gz" in their filenames; 
# "...pass_2.fastq.gz" for reverse reads

# Your data may differ, using "F" and "R" in the filenames, or something similar..
# Be sure to change that pattern to fit your own files when using your own data
fns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "R1.fastq.gz")) # make pattern match your FWD reads
rns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "R2.fastq.gz")) # make pattern match your REV reads

# This may have to be modified for your own data if your files follow a differnet naming convention
# The following line splits the filename on "_" and keeps the first element
sample.names <- unlist(map(strsplit(basename(fns), "_"), 1)) # this pulls out just the basename

# visualize a couple of fwd read quality profiles to help select reasonable filtration parameters
# you can select any number of files here...
# as-is, this just shows the Fwd and Rev quality profiles for the 1st and 2nd files
p1 <- plotQualityProfile(fns[2]) + ggtitle("Example forward reads")
p2 <- plotQualityProfile(rns[6]) + ggtitle("Example reverse reads")
p2
# display and save the plots
p1 / p2
ggsave("./output/figs/unfiltered_quality_plots.png",dpi=500,height = 6,width = 6)

# FILTER AND TRIM ####

# here, we decide what filenames we want our filtered reads to be called
# in this case, we use the base name of the sample and save the filtered files into the "filtered" subdirectory
filts_f <- file.path(path, "filtered.dada.1124", paste0(sample.names, "_R1.fastq.gz"))
filts_r <- file.path(path, "filtered.dada.1124", paste0(sample.names, "_R2.fastq.gz"))

# this is the actual quality control step
# These values are informed by our quality plots
out <- filterAndTrim(fns, filts_f, rns, filts_r, # input and output file names as denoted above
                     maxN=0, # uncalled bases are currently not supported in dada2
                     maxEE=c(2,4), # refers to the maximum expected errors allowed
                     rm.phix=TRUE, # automatically remove PhiX spike-in reads from sequencing center
                     compress=TRUE, # compress output files with gzip
                     multithread=4) # On Windows set multithread=FALSE

# save the filtration info in case we need it later
saveRDS(out, "/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/dada_out.1124.RDS")
out <- readRDS("/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/dada_out.1124.RDS")

# Did any samples have NO reads pass the filtration step?

length(fns);length(filts_f) # will be the same length if all samples had some passing reads
length(rns);length(filts_r) # will be the same length if all samples had some passing reads

# In case some samples may have had zero reads pass QC, reassign filts
filts_f <- sort(list.files(filtpath, full.names = TRUE,pattern = "R1"))
filts_r <- sort(list.files(filtpath, full.names = TRUE,pattern = "R2"))

# sanity check comparison of before and after filtration
p3 <- plotQualityProfile(rns[3]) + ggtitle("Unfiltered")
p4 <- plotQualityProfile(filts_r[3])+ ggtitle("Filtered")
p3 / p4
##ggsave("/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/filtered_quality_comparison.png",dpi=500,height = 6,width = 6)

# LEARN ERROR RATES ####

# learn errors
set.seed(123) # "random" seed for reproducibility
errF <- learnErrors(filts_f, multithread=TRUE, MAX_CONSIST = 20,verbose = 2,randomize = TRUE) # set multithread = FALSE on Windows
errF

set.seed(123)
errR <- learnErrors(filts_r, multithread=TRUE, MAX_CONSIST = 20,verbose = 2,randomize = TRUE) # set multithread = FALSE on Windows
errR

##save error files
saveRDS(errF,"/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/errF.1124.RDS")
errF <- readRDS("/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/errF.1124.RDS")
saveRDS(errR,"/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/errR.1124.RDS")
errR <- readRDS("/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/errR.1124.RDS")

# sanity check for error model
plotErrors(errF, nominalQ=FALSE)
##ggsave("/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/error_model.png",dpi=500,height = 6,width = 6)
plotErrors(errR, nominalQ=FALSE)

# dereplication
derepF <- derepFastq(filts_f, verbose=TRUE)
derepR <- derepFastq(filts_r, verbose=TRUE)

##save dereplication files
saveRDS(derepF,"/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/derepF.1124.RDS")
derepF <- readRDS("/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/derepF.1124.RDS")
saveRDS(derepR,"/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/derepR.1124.RDS")
derepR <- readRDS("/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/derepR.1124.RDS")

# Name the derep-class objects by the sample names
# If some samples were removed (no reads passed QC), reassign sample.names
if(length(derepF) != length(sample.names)){
  sample.names <- unlist(map(strsplit(basename(filts_f), "_filt"), 1))
}


if(identical(unlist(map(strsplit(basename(filts_f), "FWD_filt"), 1)),unlist(map(strsplit(basename(filts_r), "REV_filt"), 1)))){
  names(derepF) <- sample.names
  names(derepR) <- sample.names
} else {
  stop("Make sure fwd and rev files are in same order!")
}  

# SAMPLE INFERRENCE ####
set.seed(123)
dadaFs <- dada(derepF, err=errF, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo") # set multithread = FALSE on Windows
set.seed(123)
dadaRs <- dada(derepR, err=errR, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo") # set multithread = FALSE on Windows

##save sample inference files
saveRDS(dadaFs,"/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/dadaFs.1124.RDS")
dadaFs <- readRDS("/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/dadaFs.1124.RDS")
saveRDS(dadaRs,"/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/dadaRs.1124.RDS")
dadaRs <- readRDS("/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/dadaRs.1124.RDS")

# MERGE FWD and REV READS ####
mergers <- mergePairs(dadaFs, filts_f, dadaRs, filts_r, verbose=TRUE)

##save mergers files
saveRDS(mergers,"/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/mergers.1124.RDS")
mergers <- readRDS("/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/mergers.1124.RDS")

# MAKE SEQUENCE TABLE ####
seqtab <- makeSequenceTable(mergers)

# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim,"/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/seqtab.nochim.1124.RDS")
seqtab.nochim <- readRDS("/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/seqtab.nochim.1124.RDS")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

itsx <- colnames(seqtab.nochim)
names(itsx) <- paste("otu", c(1:length(itsx)), sep = ".")
dput(itsx)
itsx_dnastring <- DNAStringSet(dput(itsx))
writeXStringSet(itsx_dnastring, '/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/filtered_and_trimmed.1124.fasta')

ITSx.args <- paste("-i /home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/filtered_and_trimmed.1124.fasta",
                   "-t 'fungi'",
                   "--preserve T",
                   "--complement T", 
                   "--summary T",
                   "-o /home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/output/itsx/itsx_output.1124.fasta",
                   "--only_full T",
                   "--cpu 2",
                   "-E 1e-2")
system2("/home/busbystudents/miniconda3/pkgs/itsx-1.1.3-hdfd78af_1/bin/ITSx", args = ITSx.args)

# reassign "out" to remove any missing reads
out = out[as.data.frame(out)$reads.out > 0,]

# TRACK READS THROUGH PIPELINE ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$total.loss.proportion = (track[,1]-track$nonchim)/track[,1]
head(track)

write.csv(track, file = "/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/read_counts_at_each_step.1124.csv", row.names = TRUE)


##add itsx and hmmer3 to path
old_path <- Sys.getenv("PATH")
old_path #does it show your new stuff?
#you can add your new locations like this:

Sys.setenv(PATH = paste(old_path, "/home/busbystudents/miniconda3/pkgs/hmmer-3.4-hdbdd923_2/bin", sep = ":"))
Sys.setenv(PATH = paste(old_path, "/home/busbystudents/miniconda3/pkgs/itsx-1.1.3-hdfd78af_1/bin/ITSx", sep = ":"))
