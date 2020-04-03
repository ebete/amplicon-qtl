#!/usr/bin/env Rscript

# program setup -----------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(RColorBrewer)
  library(qtl2)
  library(cowplot)
  library(yaml)
})

argv <- commandArgs(trailingOnly = T)
theme_set(theme_pubr(border = T, legend = "right"))
max_treads <- 4


# load R/qtl2 data --------------------------------------------------------
rqtl_phenofile <- argv[1]
rqtl_control <- ifelse(length(argv) > 1, argv[2], "./sl23_control.yaml")
asv_taxonomy <- readRDS("asv_taxonomy.rds")

output_basename <- tools::file_path_sans_ext(rqtl_phenofile)
output_lodscores <- sprintf("%s_lod.csv", output_basename)
output_peaks <- sprintf("%s_peaks.csv", output_basename)
output_heatmap <-sprintf("%s_lod.pdf", output_basename)

# dirty hack to change input phenotype file
rqtl_yaml <- read_yaml(rqtl_control)
rqtl_new_control <- tempfile(pattern = "rqtl2_", fileext = ".yaml", tmpdir = dirname(rqtl_control))
rqtl_yaml$pheno <- rqtl_phenofile
write_yaml(rqtl_yaml, file = rqtl_new_control)
in_cross <- read_cross2(rqtl_new_control, quiet = F)
unlink(rqtl_new_control)


# make metadata -----------------------------------------------------------
cat("Reading metadata from", rqtl_yaml$gmap, "...\n")
chr_metadata <- read.csv(rqtl_yaml$gmap) %>%
  arrange(chr, pos) %>%
  rowid_to_column("index") %>%
  group_by(chr) %>%
  summarise(num_markers = n(), start = min(index), end = max(index), width = end - start) %>%
  remove_missing() %>%
  arrange(chr) %>%
  mutate(chr = factor(chr, ordered = T))


# run R/qtl2 --------------------------------------------------------------
cat("Performing QTL analysis on", rqtl_phenofile, "...\n")
map <- insert_pseudomarkers(in_cross$gmap, step = 1, cores = 1)
pr <- calc_genoprob(in_cross, in_cross$gmap, error_prob = 1e-4, cores = max_treads)
out <- scan1(pr, in_cross$pheno)

# make LOD score dataframe
lod_scores <- as.data.frame(out) %>%
  rownames_to_column("marker") %>%
  rowid_to_column("index") %>%
  separate(marker, into = c("pos", "chr"), remove = F, sep = "-") %>%
  pivot_longer(-c(1:4), names_to = "phenotype", values_to = "lod_score", names_ptypes = list(phenotype = factor())) %>%
  mutate(asv_id = as.numeric(substring(phenotype, 4)))
lod_scores[!(lod_scores$chr %in% 1:12), c("pos", "chr")] <- NA

cat("Writing LOD scores to", output_lodscores, "...\n")
write.table(lod_scores, file = output_lodscores, sep = "\t", quote = F, row.names = F, col.names = T)

# infer LOD peak threshold using permutation testing
# out_perm <- scan1perm(pr, in_cross$pheno, n_perm = 1000, cores = max_treads)
# perm_thresh <- as.data.frame(summary(out_perm, alpha = 0.05)) %>%
#   gather("phenotype", "threshold")

lod_inclusion_threshold <- 2.5

found_peaks <- find_peaks(out, map, threshold = lod_inclusion_threshold) %>%
  mutate(asv = as.numeric(substring(lodcolumn, 4))) %>%
  left_join(asv_taxonomy, by = "asv")

cat("Writing QTL peaks to", output_peaks, "...\n")
write.table(found_peaks, file = output_peaks, sep = "\t", quote = F, row.names = F, col.names = T)


# plotting ----------------------------------------------------------------
y_limit <- max(lod_scores$lod_score, lod_inclusion_threshold) * 1.05

# # plot all columns separately
# ggplot(lod_scores, aes(x = index)) +
#   geom_rect(data = chr_metadata, inherit.aes = F, aes(xmin = start, xmax = end+1, ymin = 0, ymax = y_limit, fill = chr), alpha = 0.25, color="white") +
#   geom_area(aes(y = lod_score), color = "black") +
#   geom_hline(data = perm_thresh, aes(yintercept = threshold)) +
#   scale_x_continuous(expand = c(0, 0), breaks = chr_metadata$start + chr_metadata$width/2, labels = sprintf("%d\n[%d]", chr_metadata$chr, chr_metadata$num_markers)) +
#   scale_y_continuous(expand = c(0, 0), limits = c(0, y_limit)) +
#   theme_pubr(border = T, legend = "none") +
#   scale_fill_manual(values = rep(c("#f1a340", "#998ec3"), 6)) +
#   facet_wrap(vars(phenotype), strip.position = "bottom") +
#   labs(
#     title = "Quantitative Trait Locus (SL2.40)",
#     x = "SNP markers",
#     y = "LOD score"
#   )
#
# # overlay all columns
# ggplot(lod_scores, aes(x = index)) +
#   geom_rect(data = chr_metadata, inherit.aes = F, aes(xmin = start, xmax = end+1, ymin = 0, ymax = y_limit, fill = chr), alpha = 0.25, color="white") +
#   geom_area(aes(y = lod_score, group = phenotype), color = NA, position = "identity", alpha = 0.1) +
#   scale_x_continuous(expand = c(0, 0), breaks = chr_metadata$start + chr_metadata$width/2, labels = sprintf("%d\n[%d]", chr_metadata$chr, chr_metadata$num_markers)) +
#   scale_y_continuous(expand = c(0, 0), limits = c(0, y_limit)) +
#   theme_pubr(border = T, legend = "none") +
#   scale_fill_manual(values = rep(c("#f1a340", "#998ec3"), 6)) +
#   labs(
#     title = "Quantitative Trait Locus (SL2.40)",
#     x = "SNP markers",
#     y = "LOD score"
#   )

# heatmap
qtl_peak_heatmap <- lod_scores %>%
  group_by(phenotype) %>%
  filter(max(lod_score) >= lod_inclusion_threshold) %>%
  ungroup() %>%
  left_join(asv_taxonomy, by = c("asv_id" = "asv")) %>%
  ggplot(aes(x = index, y = phenotype, fill = lod_score)) +
  geom_raster() +
  geom_vline(data = chr_metadata, aes(xintercept = start), color = "darkgreen", size = 0.2) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0 ,0)) +
  scale_fill_gradient2(low = "#ffffff", mid = "#aaaaaa", high = "#aa0000", midpoint = lod_inclusion_threshold) +
  # scale_fill_gradient2(low = "#000000", mid = "#555555", high = "#ff0000", midpoint = lod_inclusion_threshold) +
  facet_wrap(vars(genus), ncol = 2, scales = "free_y", strip.position = "left") +
  theme(
    axis.text.y = element_blank(),
    panel.spacing = unit(0, "lines"),
    panel.border = element_rect(color = "white", size = 0.5)
  ) +
  labs(
    title = sprintf("QTL peaks %s", gsub("_", " ", basename(output_basename))),
    subtitle = rqtl_phenofile,
    x = "Marker",
    y = "",
    fill = "LOD"
  )
cat("Writing LOD heatmap to", output_heatmap, "...\n")
ggsave2(output_heatmap, plot = qtl_peak_heatmap, units = "mm", width = 297, height = 210)

# # peak hit frequency
# ggplot(found_peaks, aes(x = pos)) +
#   geom_bar() +
#   facet_wrap(vars(chr), scales = "free_x", nrow = 2, strip.position = "bottom") +
#   theme(panel.spacing = unit(0, "lines"))
