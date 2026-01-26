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
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(MiscMetabar)
library(DECIPHER)
library(dada2)

##load data
meta <- read.csv("/home/busbystudents/Downloads/meta_data_raw.csv")
seqs <- readDNAStringSet("/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/output/itsx/itsx_output.1124.fasta.ITS1.fasta")
seqtab.nochim <- readRDS("/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/seqtab.nochim.1124.RDS")
row.names(seqtab.nochim) <- gsub("_R1.fastq.gz", "", rownames(seqtab.nochim))
colnames(seqtab.nochim) <- paste("otu", c(1:ncol(seqtab.nochim)), sep = ".")

##Create a phyloseq object with taxa table, otu table, and sequences.
rownames(meta) <- meta[,1]
seqtab.nochim <- subset(seqtab.nochim, row.names(seqtab.nochim) %in% row.names(meta))
otu.names <- names(seqs)
seqtab.nochim <- seqtab.nochim[, names(seqs)]
otu <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)
samples <- sample_data(meta)
ps<- merge_phyloseq(otu, seqs, samples)
head(ps@otu_table)

##Cluster with miscmetabarcoding.
otu_cluster <- asv2otu(ps, method = "clusterize", id = 0.95)
saveRDS(otu_cluster, file = "/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/asv2otu.1124.RDS")
otu_cluster <- readRDS("/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/asv2otu.1124.RDS")

# remove singleton taxa
otu_cluster@otu_table <- otu_cluster@otu_table[,colSums(otu_cluster@otu_table) > 1]

##subset seqs
clustered.otu.names <- colnames(otu_cluster@otu_table)
otu_cluster@refseq <- otu_cluster@refseq[clustered.otu.names]

# ASSIGN TAXONOMY ####

#First look at database with ALL eukaryotes then look at just fungal database
taxa_euk <- assignTaxonomy(otu_cluster@refseq, "/home/busbystudents/Downloads/sh_general_release_dynamic_all_18.07.2023.fasta", multithread=4)
saveRDS(taxa_euk, file = "/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/taxa_euk.1124.RDS")
taxa_euk <- readRDS("/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/taxa_euk.1124.RDS")

rownames(taxa_euk) <- names(otu_cluster@refseq)
taxa_euk <- as.data.frame(taxa_euk)
taxa_fungi <- subset(taxa_euk, Kingdom == "k__Fungi")
fungi_list <- row.names(taxa_fungi)
otu_subset <- otu_cluster@otu_table[,fungi_list]
seqs_fungi <- otu_cluster@refseq[fungi_list]

##remove taxa with less than 50 total reads
otu_fungi <- otu_subset
otu_colsums <- as.data.frame(colSums(otu_fungi))
otu_colsums$count <- otu_colsums$`colSums(otu_fungi)`
otu_colsums <- subset(otu_colsums, select = -c(1))
ggplot(otu_colsums, aes(factor(count))) +
  geom_bar(stat="count", position = "dodge") + coord_cartesian(xlim=c(0, 100)) 
colsums_50 <- subset(otu_colsums, count > 50)
colsums_50 <- row.names(colsums_50)
otu_subset <- otu_subset[,colsums_50]
seqs_fungi_50 <- seqs_fungi[colsums_50]
saveRDS(seqs_fungi_50, "/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/seqs_fungi_50_1124.RDS")

##with default bootstrap support
taxa_fungi_assign <- assignTaxonomy(seqs_fungi_50, "/home/busbystudents/Downloads/if-unite-insdc(2).fa", multithread=FALSE)


##with all assignments
taxa_fungi_assign <- assignTaxonomy(seqs_fungi, minBoot = 0, outputBootstraps = TRUE, "/home/busbystudents/Downloads/if-unite-insdc(2).fa")

# Save intermediate taxonomy file
saveRDS(taxa_fungi_assign, file = "/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/taxa_from_dada2_1124.RDS")
taxa_fungi_assign <- readRDS("/home/busbystudents/Desktop/busbyp/Projects/abbey/jones.hja/fwd+rev_11324/rds files/taxa_from_dada2_1124.RDS")


# inspect taxonomy
taxa_table <- taxa_fungi_assign

colnames(otu_subset)<-c(1:ncol(otu_subset))
colnames(otu_subset) <- paste("otu", colnames(otu_subset), sep = ".")

rownames(taxa_table)<-c(1:nrow(taxa_table))
rownames(taxa_table) <- paste("otu", rownames(taxa_table), sep = ".")

rownames(boot)<-c(1:nrow(boot))
rownames(boot) <- paste("otu", rownames(boot), sep = ".")

identical(colnames(otu_subset), rownames(taxa_table))
identical(colnames(otu_subset), rownames(boot))

write.csv(taxa_table, "/home/busbystudents/taxa_112424.csv")
write.csv(otu_subset, "/home/busbystudents/otu_112424.csv")
write.csv(boot, "/home/busbystudents/boot_112424.csv")