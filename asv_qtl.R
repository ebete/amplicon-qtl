#!/usr/bin/env Rscript

# program setup -----------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(RColorBrewer)
  library(cowplot)
})

theme_set(theme_cowplot())


# ASVs per group ----------------------------------------------------------
asv_output <- read.delim("asv_16s.csv")
asv_taxonomy <- asv_output[, 1:8]

asv_p <- asv_output[, grepl("_P$", colnames(asv_output))]
asv_y <- asv_output[, grepl("_Y$", colnames(asv_output))]
asv_money <- asv_output[, grepl("_MM[\\d]$", perl = T, colnames(asv_output))]
asv_pimpi <- asv_output[, grepl("_P[\\d]$", perl = T, colnames(asv_output))]
asv_bulk <- asv_output[, grepl("_BULK", colnames(asv_output))]

group_asv_means <- data.frame(
  p = rowMeans(asv_p),
  y = rowMeans(asv_y),
  money = rowMeans(asv_money),
  pimpi = rowMeans(asv_pimpi),
  bulk = rowMeans(asv_bulk)
)

group_asv_means[1:100, ] %>%
  rowid_to_column("asv") %>%
  pivot_longer(-asv, names_to = "group", values_to = "mean_count") %>%
  ggplot(aes(x = asv, y = mean_count, fill = group, color = group)) +
  geom_area(alpha = 0.8) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1")


asv_p.t <- as.data.frame(t(asv_p))[, 1:10] %>%
  rownames_to_column("id") %>%
  mutate(id = as.numeric(sub("^.*_([\\d]+).*$", "\\1", id, perl = T))) %>%
  arrange(id)
colnames(asv_p.t) <- c("id", paste0("ASV", 1:(ncol(asv_p.t)-1)))
write.csv(asv_p.t, file = "p_asv.csv", row.names = F)
