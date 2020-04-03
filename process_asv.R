#!/usr/bin/env Rscript

# program setup -----------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(RColorBrewer)
  library(cowplot)
  library(metagenomeSeq)
})

theme_set(theme_pubr(border = T, legend = "right"))
set.seed(56753)


# data loading ------------------------------------------------------------
asv_output <- read.delim("asv_16s.csv")
asv_18s_output <- read.delim("asv_18s.csv")

asv_taxonomy <- asv_output[, 1:8] %>%
  mutate(genus = sub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Rhizobium", genus, fixed = T))
asv_counts <- as.matrix(asv_output[, -c(1:8)])
accession_metadata <- data.frame(accession = colnames(asv_counts)) %>%
  mutate(
    accession = as.character(accession),
    group = as.factor(case_when(
      grepl("_P$", accession, ignore.case = T) ~ "Pink",
      grepl("_Y$", accession, ignore.case = T) ~ "Yellow",
      grepl("_BULK", accession, ignore.case = T) ~ "Bulk",
      grepl("_MM[\\d]$", accession, ignore.case = T, perl = T) ~ "Modern",
      grepl("_P[\\d]$", accession, ignore.case = T, perl = T) ~ "Wild",
      TRUE ~ "Other"
    )),
    sample_id = as.numeric(sub("^[X]?([\\d]+).*", "\\1", accession, perl = T)),
    ril_id = as.numeric(ifelse(group %in% c("Pink", "Yellow"), sub(".*_([\\d]+).*", "\\1", accession, perl = T), NA))
  )
row.names(accession_metadata) <- accession_metadata$accession


# ASV pre-filtering -------------------------------------------------------
# remove other genotypes; Y:253,226,259,275; P:217,237,249 (too long to germinate)
genotypes_to_keep <- !(
  accession_metadata$group == "Other" |
  (accession_metadata$group == "Yellow" & accession_metadata$ril_id %in% c(253, 226, 259, 275)) |
  (accession_metadata$group == "Pink" & accession_metadata$ril_id %in% c(217, 237, 249))
)


asv_counts <- asv_counts[, genotypes_to_keep]
accession_metadata <- accession_metadata[genotypes_to_keep, ] %>%
  mutate(group = droplevels(group))


# remove ASVs that are chloroplasts, mitochondria, or not bacteria
asvs_to_keep <- !(
  is.na(asv_taxonomy$kingdom) |
  (asv_taxonomy$kingdom != "Bacteria") |
  (!is.na(asv_taxonomy$order) & asv_taxonomy$order == "Chloroplast") |
  (!is.na(asv_taxonomy$family) & asv_taxonomy$family == "Mitochondria")
)
# remove zero-count ASVs
asvs_to_keep <- asvs_to_keep & (rowSums(asv_counts) > 0)

# filter and renumber ASV table
asv_taxonomy <- asv_taxonomy[asvs_to_keep, ]
asv_taxonomy$asv <- 1:nrow(asv_taxonomy)
saveRDS(asv_taxonomy, file = "asv_taxonomy.rds", compress = "xz")


# ASV CSS normalisation ---------------------------------------------------
# ADFs needed by metagenomeSeq MR experiments
phenotype_data <- AnnotatedDataFrame(data.frame(group = accession_metadata$group, row.names = accession_metadata$accession))
feature_data <- AnnotatedDataFrame(asv_taxonomy[asvs_to_keep, ])
feature_data@data$asv <- 1:nrow(feature_data)
row.names(feature_data@data) <- feature_data@data$asv
# create the MR exp with remaining ASVs
asv_mr_exp <- newMRexperiment(
  asv_counts[asvs_to_keep, ],
  phenoData = phenotype_data,
  featureData = feature_data
)

# metagenomeSeq cumulative sum scaling (CSS)
css_quantile <- cumNormStat(asv_mr_exp)
asv_css <- cumNorm(asv_mr_exp, css_quantile)

# zero-inflated gaussian mixture model
accession_group <- pData(asv_css)$group
normalisation_factor <- normFactors(asv_css)
normalisation_factor <- log2(normalisation_factor / median(normalisation_factor) + 1)
model_mtx <- model.matrix(~ accession_group + normalisation_factor)
model_control <- zigControl(maxit = 20, verbose = T)
model_fit <- fitZig(asv_css, model_mtx, control = model_control, useCSSoffset = T)
effective_samples <- calculateEffectiveSamples(model_fit)
effective_samples_minimum <- ave(effective_samples)[[1]]
# filter rare ASVs based on the effective sample size
rare_asvs <- unname(which(rowSums(MRcounts(asv_css)) < effective_samples_minimum))
normalised_counts <- as.data.frame(MRcounts(asv_css, norm = T))[-rare_asvs,]

# ZIG on rare-filtered ASVs to find rhizosphere ASVs
asv_css_common <- asv_css[-rare_asvs, ]
fz_model <- model.matrix(~ accession_group)
colnames(fz_model) <- levels(accession_group)
fz <- fitZig(asv_css_common, fz_model, control = model_control)

contrast_matrix <- makeContrasts(Bulk - (Yellow + Pink + Modern + Wild), levels = fz@fit$design)
fz_contrast <- contrasts.fit(fz@fit, contrast_matrix)
fz_bayes <- eBayes(fz_contrast)
fz_padj <- p.adjust(fz_bayes$p.value, method = "hochberg")

# final normalised and filtered ASV counts
asv_css_filtered <- asv_css_common[fz_padj <= 0.01, ]
asv_css_filtered_normalised <- MRcounts(asv_css_filtered, norm = TRUE, log = TRUE) %>%
  na.omit()

cat(sprintf(
  "%.1f%% ASVs filtered",
  100 - nrow(asv_css_filtered_normalised) / nrow(asv_output) * 100
))


plotOrd(asv_css_filtered_normalised, usePCA = F, useDist = T, bg = accession_group, pch = 21)
legend("bottomright", levels(accession_group), text.col = 1:5)


# make cut-offs for: stochastic, core, and variable -----------------------
core_cutoff <- 0.8
rare_cutoff <- 0.1

css_statistics_rhiz <- as.data.frame(t(apply(asv_css_filtered_normalised[, accession_metadata$group != "Bulk"], 1, function(x){c(
  "zeroes" = sum(x == 0),
  "nz" = sum(x > 0),
  "mean" = mean(x, na.rm = T),
  "nz_mean" = mean(x[x > 0]),
  "sum" = sum(x, na.rm = T),
  "n" = length(x),
  "min" = min(x, na.rm = T),
  "max" = max(x, na.rm = T),
  "norm_sd" = sd(x, na.rm = T) / max(x, na.rm = T),
  "norm_sd_nz" = sd(x[x > 0]) / max(x, na.rm = T)
)})), row.names = row.names(asv_css_filtered_normalised)) %>%
  rownames_to_column("asv") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(
    asv = as.numeric(asv),
    ecology = as.factor(case_when(
      nz / n > core_cutoff ~ "Core",
      nz / n < rare_cutoff ~ "Stochastic",
      TRUE ~ "Flexible"
    ))
  ) #%>% left_join(asv_taxonomy, by = "asv")


# visualise groups
ggplot(css_statistics_rhiz, aes(x = norm_sd_nz, y = nz / n, size = (sum/max(sum))^2, fill = ecology)) +
  geom_point(alpha = 0.25, shape = 21, color = "black") +
  geom_hline(yintercept = core_cutoff, linetype = "dotted") +
  annotate("text", x = 0, y = core_cutoff, label = sprintf("%.0f%%", core_cutoff*100), hjust = -0.2, vjust = -0.4) +
  geom_hline(yintercept = rare_cutoff, linetype = "dotted") +
  annotate("text", x = 0, y = rare_cutoff, label = sprintf("%.0f%%", rare_cutoff*100), hjust = -0.2, vjust = -0.4) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0, 0, 0.05), labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "ASV separation based on presence/absence",
    x = "Normalised standard deviation",
    y = "% non-zero counts",
    size = "norm.sum^2"
  )


# Exporting of tables -----------------------------------------------------
as.data.frame(asv_mr_exp@assayData$counts) %>%
  rowid_to_column("aasv") %>%
  write.table(., file = "TableS1_asv_raw.csv", sep = "\t", row.names = F, col.names = substring(colnames(.), 2), quote = F)
as.data.frame(asv_css_filtered_normalised) %>%
  rownames_to_column("aasv") %>%
  write.table(., file = "TableS2_asv_css.csv", sep = "\t", row.names = F, col.names = substring(colnames(.), 2), quote = F)
write.table(asv_mr_exp@featureData@data, file = "TableS3_asv_taxonomy.csv", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(css_statistics_rhiz, file = "TableS4_asv_css_statistics_rhizosphere.csv", sep = "\t", row.names = F, col.names = T, quote = F)


# for exporting both the normalised ASV and randomised accession ASV table
export_css <- function(asv_mtx, selected_accessions, selected_asvs) {
  # create final table
  asv_df <- as.data.frame(t(asv_mtx[css_statistics_rhiz$ecology == selected_asvs, accession_metadata$group == selected_accessions]))
  colnames(asv_df) <- paste0("ASV", rownames(asv_mtx)[css_statistics_rhiz$ecology == selected_asvs])
  asv_df <- cbind(ril_id = subset(accession_metadata, group == selected_accessions)$ril_id, asv_df) %>%
    group_by(ril_id) %>%
    summarise_all(~mean(.)) %>%
    arrange(ril_id)
  # export final table
  asv_df_fname <- sprintf("%s_%s.csv", str_to_lower(selected_accessions), str_to_lower(selected_asvs))
  write.table(asv_df, file = asv_df_fname, sep = ",", quote = F, row.names = F, col.names = T)
  cat(sprintf("Exported %s\n", asv_df_fname))
  # create genotype-randomised table and export
  # asv_df_randomised <- as.data.frame(t(apply(asv_df, 1, function(x) { c(x[1], sample(x[-1])) })))
  # colnames(asv_df_randomised) <- colnames(asv_df)
  asv_df_randomised <- cbind(asv_df[, 1], as.data.frame(apply(asv_df[, -1], 2, sample)))
  asv_df_randomised_fname <- sprintf("%s_%s_randomised.csv", str_to_lower(selected_accessions), str_to_lower(selected_asvs))
  write.table(asv_df_randomised, file = asv_df_randomised_fname, sep = ",", quote = F, row.names = F, col.names = T)
  cat(sprintf("Exported %s\n", asv_df_randomised_fname))
  # end
  return(invisible(NULL))
}

# export the tables
export_css(asv_css_filtered_normalised, "Pink", "Flexible")
export_css(asv_css_filtered_normalised, "Pink", "Core")
export_css(asv_css_filtered_normalised, "Yellow", "Flexible")
export_css(asv_css_filtered_normalised, "Yellow", "Core")

## Column-randomization à la Ben
# mtx <- as.matrix(asv_df[, -1])
# rownames(mtx) <- asv_df$ril_id
# mtx_rnd <- mtx
# for (i in 1:dim(mtx_rnd)[2]) {
#   random_vector <- sample(mtx_rnd[1:dim(mtx_rnd)[1],i])
#   mtx_rnd[1:dim(mtx_rnd)[1], i] <- random_vector
# }
