#!/usr/bin/env Rscript


# program setup -----------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(dada2)
  library(phyloseq)
})

max_treads <- 8


# load pre-trimmed reads --------------------------------------------------
trimmed_reads_path <- "./reads_filtered"
fwd_reads <- sort(list.files(trimmed_reads_path, pattern = "18S_R1.fq.gz", full.names = T))
rev_reads <- sort(list.files(trimmed_reads_path, pattern = "18S_R2.fq.gz", full.names = T))

sample_names <- substring(basename(fwd_reads), first = 1, last = nchar(basename(fwd_reads)) - 13)

filtered_reads_path <- "./reads_dada2"
fwd_filtered <- file.path(filtered_reads_path, basename(fwd_reads))
rev_filtered <- file.path(filtered_reads_path, basename(rev_reads))

names(fwd_filtered) <- sample_names
names(rev_filtered) <- sample_names

# Baseclear: PhiX removed, len<50 dropped, adapters removed, Phred check
filterAndTrim(fwd_reads, fwd_filtered, rev_reads, rev_filtered, compress = T, maxN = 0, rm.phix = T, minLen = 50, trimLeft = 0, trimRight = 0, multithread = max_treads, verbose = T)


# learn error rates -------------------------------------------------------
fwd_errors <- learnErrors(fwd_filtered, multithread = max_treads)
rev_errors <- learnErrors(rev_filtered, multithread = max_treads)

plotErrors(fwd_errors, nominalQ = T)
plotErrors(rev_errors, nominalQ = T)


# sample inference --------------------------------------------------------
fwd_dada <- dada(fwd_filtered, err = fwd_errors, pool = F, multithread = max_treads)
rev_dada <- dada(rev_filtered, err = rev_errors, pool = F, multithread = max_treads)

merged_reads <- mergePairs(fwd_dada, fwd_filtered, rev_dada, rev_filtered)


# create count table ------------------------------------------------------
sequence_table <- makeSequenceTable(merged_reads, orderBy = "abundance")
sequence_table_nochimeras <- removeBimeraDenovo(sequence_table, method = "consensus", multithread = max_treads)
cat(sprintf("Chimeric reads: %.1f%%", 100 - sum(sequence_table_nochimeras) / sum(sequence_table) * 100))


# add taxonomy ------------------------------------------------------------
taxonomies <- assignTaxonomy(sequence_table_nochimeras, "~/processed_16s18s/silva_nr_v132_train_set.fa.gz", multithread = max_treads)
taxonomies <- addSpecies(taxonomies, "~/processed_16s18s/silva_species_assignment_v132.fa.gz")


# export ASV table --------------------------------------------------------
taxonomy_table <- as.data.frame(taxonomies) %>%
  rownames_to_column("asv_seq")
colnames(taxonomy_table) <- tolower(colnames(taxonomy_table))

asv_output <- as.data.frame(t(sequence_table_nochimeras)) %>%
  rownames_to_column("asv_seq") %>%
  arrange(desc(rowSums(.[-1]))) %>%
  mutate(asv = row_number()) %>%
  left_join(taxonomy_table, by = "asv_seq") %>%
  select(-asv_seq)
# set asv and taxonomy as the first columns
asv_output <- asv_output[, c((ncol(asv_output)-7):ncol(asv_output), 1:(ncol(asv_output)-8))]

write.table(asv_output, file = "asv_18s.csv", quote = F, sep = "\t", row.names = F, col.names = T, fileEncoding = "utf-8")
