# Amplicon QTL
ASV calling from raw 16S reads in FASTQ format.


## 1. [asv_16s.R](asv_16s.R)
1. Put raw files in the directory `./reads_trimmed/`, named like `*_16S_R1.fq.gz` and `*_16S_R2.fq.gz` for the forward and reverse reads.
2. Execute script like `Rscript asv_16s.R`.
3. Output files are:
   * `asv_16s.csv`, containing the raw ASVs, taxonomy, and counts per sample.
   * `asv_16s_statistics.csv`, containing read filtering statistics.

## 2. [process_asv.R](process_asv.R)
1. Execute after running [asv_16s.R](asv_16s.R), as it uses `asv_16s.csv` as an input file.
2. Execute script like `Rscript asv_16s.R`.
3. Output files are:
   * Figure showing separation between the core, flexible, and stochastic microbiome.
   * `asv_taxonomy.rds`, an R data file containing renumbered ASV-IDs and taxonomy.
   * `TableS1_asv_raw.csv`, raw counts ASV table with Chloroplast and non-bacteria taxa removed. ASV-IDs were renumbered.
   * `TableS2_asv_css.csv`, CSS-normalised and fitZIG filtered ASVs.
   * `TableS3_asv_taxonomy.csv`, ASV-IDs and related taxonomical lineages.
   * `TableS4_asv_css_statistics_rhizosphere.csv`, statistics calculated on TableS2 counts.
   * `(yellow|pink)_(core|flexible).csv` and `(yellow|pink)_(core|flexible)_randomised.csv`, containing ASV counts from a subset of samples and one where the order of counts in each column are randomised.

## 3. [asv_qtl.R](asv_qtl.R)
1. Execute after running [process_asv.R](process_asv.R), as it uses `asv_taxonomy.rds` and an ASV table as an input file.
2. Also uses an `sl23_control.yaml` R/qtl2 control file. Is currently set up to use the supplied SL2.40 SNP map.
3. Execute script like `Rscript asv_qtl.R INPUT_ASV_TABLE [RQTL_CONTROL_FILE]`.
4. Output files are:
   * `*_lod.csv`, containing the LOD scores of each supplied ASV.
   * `*_peaks.csv`, peaks identified using a cut-off score of `LOD >= 2.5`
   * `*_lod.pdf`, visualisation of LOD scores across the genome.
