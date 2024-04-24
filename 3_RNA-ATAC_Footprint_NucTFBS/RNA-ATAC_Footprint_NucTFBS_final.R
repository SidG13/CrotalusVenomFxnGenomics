#### This script is to combine various TF explanatory figures together ####

library(cowplot)
library(tidyverse)
library(matrixStats)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(GenomicRanges)
library(ggsignif)

#### 1: RNA/ATAC PCAs ####
setwd('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Final_Figure_Assembly/3_RNA-ATAC_Footprint_NucTFBS')

# RNA PCA #
VST_count_mat_orig <- read.table('../../Figure_and_Scripts/z_data_files/RNAseq_VSTNormalizeCounts_noOutliers_05.23.23.txt', header=T)
VST_count_mat_orig <- VST_count_mat_orig %>% 
  dplyr::select(contains('LVG'), contains('RVG')) %>% 
  filter(grepl("Venom_", rownames(.)))

column_numbers <- 1:13
target_names <- paste0('VG_', column_numbers)

for (i in column_numbers) {
  target_columns <- grep(paste0('_', i, '$'), colnames(VST_count_mat_orig))
  
  if (length(target_columns) > 0) {
    if (length(target_columns) == 1) {
      # If only one column is found, copy it to become VG_i
      VST_count_mat_orig[target_names[i]] <- VST_count_mat_orig[, target_columns]
    } else {
      # Calculate row mean for multiple columns
      VST_count_mat_orig[target_names[i]] <- rowMeans(VST_count_mat_orig[, target_columns], na.rm = TRUE)
    }
  }
}

VST_count_mat_orig <- VST_count_mat_orig %>% 
  dplyr::select(-contains('LVG'),-contains('RVG'))

# Perform PCA
pca_vg <- prcomp(VST_count_mat_orig, scale. = TRUE)

# Create a data frame with sample names and corresponding PCA coordinates
pca_data_vg <- data.frame(Sample = colnames(VST_count_mat_orig),
                          PC1 = pca_vg$rotation[, 1],
                          PC2 = pca_vg$rotation[, 2])

sample_info <- read.table('../../Figure_and_Scripts/z_data_files/rnaSampleInfo_06.29.21.txt', header = T)
sample_info <- sample_info %>% 
  filter(str_detect(Sample_ID, 'RVG')) %>% 
  mutate(Sample_ID = gsub('R', '', Sample_ID)) %>% 
  mutate(Population = fct_relevel(Population, 
                                  "South", "Mid", "North", "Other")) %>% 
  arrange(Lineage, Population, Sample_ID)

pca_data_vg <- pca_data_vg %>% 
  left_join(sample_info, by = c('Sample'='Sample_ID')) %>% 
  mutate(LineageSpecific = ifelse(Population != 'Other', paste0(Lineage, Population), Lineage))

p_pca1 <- ggplot(pca_data_vg, aes(x = PC1, y = PC2, color = LineageSpecific)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(x = "PC1", y = "PC2", title = "Venom gene expression") +
  scale_color_manual(values = c("lutosus"="#8B481F", 
                                "concolor"= "#D9A528",
                                "cerberus"="black",
                                "viridisMid"="#D81B60",
                                "viridisSouth"="#1E88E5",
                                "viridisNorth"="#07856f"),
                     labels = (c(expression(italic("C. lutosus")), 
                                 expression(italic("C. concolor")), 
                                 expression(italic("C. cerberus")), 
                                 expression(paste(italic("C. viridis"), "-Mid")), 
                                 expression(paste(italic("C. viridis"), "-South")), 
                                 expression(paste(italic("C. viridis"), "-North")))),
                     breaks = c("lutosus", "concolor", "cerberus", "viridisMid", "viridisSouth", "viridisNorth")) +
  #geom_text(aes(label = Sample), vjust = -0.8) +
  guides(color = guide_legend(title="Lineage",
                              override.aes = list(shape = 16),
                              nrow = 2)) +
  theme(legend.direction = 'horizontal',
        legend.position = 'bottom',
        plot.margin = margin(20,30,0,30),
        legend.margin = margin(-10,50,0,10))

# ATAC PCA #
ATAC_mat <- read.table('scaleFactorNormCounts.June2023.txt', header = T)
ATAC_mat <- ATAC_mat %>%
  rename_with(~sub("\\..*", "", .), everything())

# Filter for venom regions
enh_pr <- read.table('../../Figure_and_Scripts/z_data_files/venom_enhancers_promoters-1000+100_09.12.23.txt', header = T)

ATAC_mat_gr <- with(ATAC_mat, GRanges(chr, IRanges(start, end)))
enh_pr_gr <- with(enh_pr, GRanges(seqnames, IRanges(start, end)))

# Find overlaps between atac_peaks scores and enhancers/promoters
overlaps <- findOverlaps(ATAC_mat_gr, enh_pr_gr) # motifs in merged peak set
ATAC_mat_venom_regs <- cbind(as.data.frame(ATAC_mat[queryHits(overlaps),]), enh_pr[subjectHits(overlaps),])
ATAC_mat_venom_regs <- ATAC_mat_venom_regs %>% 
  dplyr::select(-c(chr, start, end, seqnames, name, type)) # leave as matrix

other_meta <- readxl::read_xlsx('../../Figure_and_Scripts/z_data_files/fixed_12snake_metadata_05.21.23.xlsx', sheet = 2, col_names = T) %>% 
  janitor::clean_names() %>% 
  dplyr::select(corrected_name, pop_id, lineage, sex, rna_12_id) %>% 
  dplyr::slice(1:13) %>% 
  dplyr::slice(-10) %>% 
  mutate(full_name = paste(corrected_name, lineage, pop_id, sex, sep = '_'))

# Perform PCA
pca_atac <- prcomp(ATAC_mat_venom_regs, scale. = TRUE)

# Create a data frame with sample names and corresponding PCA coordinates
pca_data_atac <- data.frame(Sample = colnames(ATAC_mat_venom_regs),
                            PC1 = pca_atac$rotation[, 1],
                            PC2 = -pca_atac$rotation[, 2])
pca_data_atac <- pca_data_atac %>% 
  left_join(other_meta, by = c('Sample'='full_name')) %>% 
  mutate(LineageSpecific = ifelse(pop_id != 'Other', paste0(lineage, pop_id), lineage))

p_pca2 <- ggplot(pca_data_atac, aes(x = PC1, y = PC2, color = LineageSpecific)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(x = "PC1", y = "PC2", title = "ATAC-seq peak scores\nat venom regions") +
  scale_color_manual(values = c("lutosus"="#8B481F", 
                                "concolor"= "#D9A528",
                                "cerberus"="black",
                                "viridisMid"="#D81B60",
                                "viridisSouth"="#1E88E5",
                                "viridisNorth"="#07856f"),
                     labels = (c(expression(italic("C. lutosus")), 
                                 expression(italic("C. concolor")), 
                                 expression(italic("C. cerberus")), 
                                 expression(paste(italic("C. viridis"), "-Mid")), 
                                 expression(paste(italic("C. viridis"), "-South")), 
                                 expression(paste(italic("C. viridis"), "-North")))),
                     breaks = c("lutosus", "concolor", "cerberus", "viridisMid", "viridisSouth", "viridisNorth")) +
  #geom_text(aes(label = Sample), vjust = -0.8) +
  guides(color = 'none') +
  theme(plot.margin = margin(20,20,0,30))

p1 <- plot_grid(p_pca1, p_pca2, align = 'hv', axis = 'tb',  nrow = 1)

#### 2: Enhancer/Promoter Footprint/Chromatin accessibility ####
enh_pr <- read.table('../../Figure_and_Scripts/z_data_files/venom_enhancers_promoters-1000+100_08.08.23.txt', header = T)

# Read in footprint data
fp_thresh <- read.table('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/3_ATACseq_comparisons/4_Plot_FootPrint_bw/allBoundThresholds.txt', col.names = c('SampleID', 'threshold'))
fp_mat1 <- read.table('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/3_ATACseq_comparisons/4_Plot_FootPrint_bw/AvgCounts.FootprintMatrix.enhancerTFBS_Sept2023.txt', header = T)
fp_mat2 <- read.table('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/3_ATACseq_comparisons/4_Plot_FootPrint_bw/AvgCounts.FootprintMatrix.promoterTFBS_Nov2023.txt', header = T)
fp_mat1$Feature <- 'enhancer'
fp_mat2$Feature <- 'promoter'
fp_mat <- rbind(fp_mat1, fp_mat2)
rm(fp_mat1, fp_mat2)
fp_mat[fp_mat == 'NaN'] <- NA
# fp_mat[fp_mat == 'NaN'] <- 0 # change NaNs to 0, footprints of 0 should be treated as nothing binds
# fp_mat <- fp_mat %>% na.omit() # remove NAs entirely, leaving fp matrix as continuous

fp_mat_gr <- with(fp_mat, GRanges(chr, IRanges(start, end)))
enh_pr_gr <- with(enh_pr, GRanges(seqnames, IRanges(start, end)))

# Find overlaps between atac_peaks scores and enhancers/promoters
overlaps <- findOverlaps(fp_mat_gr, enh_pr_gr) # motifs in merged peak set
fp_mat_enh_pr <- cbind(as.data.frame(fp_mat[queryHits(overlaps),]), enh_pr[subjectHits(overlaps),])
rm(fp_mat_gr, enh_pr_gr, overlaps)

fp_mat_enh_pr <- fp_mat_enh_pr %>% 
  pivot_longer(cols = starts_with('CV'), values_to = 'fp_peak_score', names_to = 'SampleID') %>% 
  mutate(SampleID = gsub('_footprints', '', SampleID)) %>% 
  left_join(fp_thresh, by = 'SampleID') %>% 
  # filter(fp_peak_score > threshold) %>% # filter by sample-specific thresholds, remove if using as binary matrix
  mutate(fp_peak_score = if_else(fp_peak_score > threshold | fp_peak_score == 0, fp_peak_score, 0)) # filter by sample-specific thresholds, set NAs to 0 and anything below threshold to 0

fp_variances <- fp_mat_enh_pr %>%
  group_by(type, name) %>%
  summarise(variance = var(fp_peak_score, na.rm = T))

### Optional filter ###
# Calculate the interquartile range
Q1 <- quantile(fp_variances$variance, 0.25, na.rm = T)
Q3 <- quantile(fp_variances$variance, 0.75, na.rm = T)
IQR <- Q3 - Q1
# Define the threshold for outliers
outlier_threshold <- 1.5 * IQR

# Filter the dataframe to remove outliers
fp_variances <- fp_variances %>%
  filter(variance >= Q1 - outlier_threshold & variance <= Q3 + outlier_threshold)
### End Optional filter ###

kruskal.test(variance ~ type, data = fp_variances) # significant

p_fp <- ggplot(fp_variances, aes(x = type, y = variance, fill = type)) +
  geom_boxplot() +
  theme_bw() + 
  xlab(NULL) +
  ylab("Footprint score\nvariation") +
  scale_fill_manual(values = c('enhancer' = '#1B9E77', 'promoter' = '#D95F02')) +
  guides(fill = 'none') +
  theme(plot.margin = margin(15,5,15,5),
        axis.text.x = element_text(angle = 50, hjust = 1))


# Read ATAC-seq peak scores
ATAC_mat <- read.table('~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_12SnakeVenomRegulation/data/atacseq/scaleFactorNormCounts.June2023.txt', header = T)
ATAC_mat <- ATAC_mat %>%
  rename_with(~sub("\\..*", "", .), everything())

enh_pr <- read.table('../../Figure_and_Scripts/z_data_files/venom_enhancers_promoters-1000+100_09.12.23.txt', header = T)

ATAC_mat_gr <- with(ATAC_mat, GRanges(chr, IRanges(start, end)))
enh_pr_gr <- with(enh_pr, GRanges(seqnames, IRanges(start, end)))

# Find overlaps between atac_peaks scores and enhancers/promoters
overlaps <- findOverlaps(ATAC_mat_gr, enh_pr_gr) # motifs in merged peak set
ATAC_peak_scores_enh_pr <- cbind(as.data.frame(ATAC_mat[queryHits(overlaps),]), enh_pr[subjectHits(overlaps),])
ATAC_peak_scores_enh_pr <- ATAC_peak_scores_enh_pr %>% 
  dplyr::select(-c(chr, start, end, seqnames))
 ####
ATAC_peak_scores_enh_pr <- ATAC_peak_scores_enh_pr %>% 
  pivot_longer(cols = starts_with('CV'), values_to = 'ATAC_peak_score', names_to = 'SampleID') # %>% 
#filter(grepl('viridis', SampleID))

ATAC_variances <- ATAC_peak_scores_enh_pr %>%
  group_by(type, name) %>%
  summarise(variance = var(ATAC_peak_score))

### Optional filter ###
# Calculate the interquartile range
Q1 <- quantile(ATAC_variances$variance, 0.25)
Q3 <- quantile(ATAC_variances$variance, 0.75)
IQR <- Q3 - Q1
# Define the threshold for outliers
outlier_threshold <- 1.5 * IQR

# Filter the dataframe to remove outliers
ATAC_variances <- ATAC_variances %>%
  filter(variance >= Q1 - outlier_threshold & variance <= Q3 + outlier_threshold)
### End Optional filter ###

kruskal.test(variance ~ type, data = ATAC_variances) # significant
wilcox.test(variance ~ type, data = ATAC_variances) # significant

p_atac <- ggplot(ATAC_variances, aes(x = type, y = variance, fill = type)) +
  geom_boxplot() +
  theme_bw() + 
  xlab(NULL) +
  ylab("Accessibility\nvariation") +
  scale_fill_manual(values = c('enhancer' = '#1B9E77', 'promoter' = '#D95F02')) +
  guides(fill = 'none') +
  theme(plot.margin = margin(15,5,15,5),
        axis.text.x = element_text(angle = 50, hjust = 1))

p2 <- plot_grid(p_atac, p_fp, align = 'tb', nrow = 1)

#### 3. Gene family ATAC boxplots #### 
ATAC_variances <- ATAC_variances %>% 
  mutate(family = case_when(
    grepl('CRISP', name) ~ 'CRISP',
    grepl('CTL', name) ~ 'CTL', 
    grepl('LAAO', name) ~ 'LAAO',
    grepl('SVSP', name) ~ 'SVSP',
    grepl('SVMP', name) ~ 'SVMP',
    grepl('PLA2', name) ~ 'PLA2',
    grepl('vQC', name) ~ 'vQC',
    grepl('VEGF', name) ~ 'VEGF',
    grepl('Vespryn', name) ~ 'Vespryn',
    grepl('EXO', name) ~ 'EXO', 
    grepl('myo', name) ~ 'myotoxin',
    TRUE ~ 'others'
  ))

# Reorder levels of `family` based on mean variance
ordered_levels_atac <- ATAC_variances %>%
  group_by(family) %>%
  summarise(mean_variance = mean(variance)) %>%
  arrange(desc(mean_variance)) %>%
  pull(family)

# Convert `family` to a factor with the reordered levels
ATAC_variances$family <- factor(ATAC_variances$family, levels = ordered_levels_atac)

p_box_atac <- ggplot(ATAC_variances, aes(x = family, y = variance, fill = family)) +
  geom_boxplot(color="#0a0dc4", fill="#03dffc", alpha=0.2) + # or geom_boxplot
  theme_bw() + 
  xlab(NULL) +
  ylab("Accessibility\nvariation") +
  guides(fill = 'none') +
  theme(axis.text.x = element_text(angle = 70, hjust = 1),
        plot.margin = margin(-10,5,10,5))

fp_variances <- fp_variances %>% 
  mutate(family = case_when(
    grepl('CRISP', name) ~ 'CRISP',
    grepl('CTL', name) ~ 'CTL', 
    grepl('LAAO', name) ~ 'LAAO',
    grepl('SVSP', name) ~ 'SVSP',
    grepl('SVMP', name) ~ 'SVMP',
    grepl('PLA2', name) ~ 'PLA2',
    grepl('vQC', name) ~ 'vQC',
    grepl('VEGF', name) ~ 'VEGF',
    grepl('Vespryn', name) ~ 'Vespryn',
    grepl('EXO', name) ~ 'EXO', 
    grepl('myo', name) ~ 'myotoxin',
    TRUE ~ 'others'
  ))

# Reorder levels of `family` based on mean variance
ordered_levels_fp <- fp_variances %>%
  group_by(family) %>%
  summarise(mean_variance = mean(variance)) %>%
  arrange(desc(mean_variance)) %>%
  pull(family)

# Convert `family` to a factor with the reordered levels
fp_variances$family <- factor(fp_variances$family, levels = ordered_levels_fp)

p_box_fp <- ggplot(fp_variances, aes(x = family, y = variance, fill = family)) +
  geom_boxplot(color="#0a0dc4", fill="#03dffc", alpha=0.2) + # or geom_boxplot
  theme_bw() + 
  xlab(NULL) +
  ylab("Footprint score\nvariation") +
  guides(fill = 'none') +
  theme(axis.text.x = element_text(angle = 70, hjust = 1),
        plot.margin = margin(-10,5,5,5))

p3 <- plot_grid(p_box_atac, p_box_fp, align='hv', axis = 'lr', ncol=1)

# 3.52 x 4.35
p4 <- plot_grid(p2, p3, align='hv', axis = 'lr', ncol=1, rel_heights = c(0.8,1))

#### 4: ATACxFP GEX correlation plot ####
ATAC_variances <- ATAC_variances %>% 
  mutate(gene = str_split(name, "\\.")) %>% 
  unnest(gene) %>%
  relocate(gene, .after = name) %>% 
  mutate(gene = gsub("^PER\\d+_|^Promoter_", "", gene)) %>% 
  mutate(gene = paste0('Venom_', gene)) %>% 
  dplyr::rename(ATAC_variance = variance)

fp_variances <- fp_variances %>% 
  mutate(gene = str_split(name, "\\.")) %>% 
  unnest(gene) %>%
  relocate(gene, .after = name) %>% 
  mutate(gene = gsub("^PER\\d+_|^Promoter_", "", gene)) %>% 
  mutate(gene = paste0('Venom_', gene)) %>% 
  dplyr::rename(fp_variance = variance)

# Add VG expression information
VST_count_mat <- read.table('../../Figure_and_Scripts/z_data_files/RNAseq_VSTNormalizeCounts_noOutliers_05.23.23.txt', header=T)

# Filter for Venom genes, also combine LVG and RVG (take the average). Note that LVG_13 failed and is missing.
VST_count_mat <- VST_count_mat %>% 
  dplyr::select(matches('LVG|RVG')) %>% 
  filter(grepl('Venom', row.names(.)))

column_numbers <- 1:13
target_names <- paste0('VG_', column_numbers)

for (i in column_numbers) {
  target_columns <- grep(paste0('_', i, '$'), colnames(VST_count_mat))
  
  if (length(target_columns) > 0) {
    if (length(target_columns) == 1) {
      # If only one column is found, copy it to become VG_i
      VST_count_mat[target_names[i]] <- VST_count_mat[, target_columns]
    } else {
      # Calculate row mean for multiple columns
      VST_count_mat[target_names[i]] <- rowMeans(VST_count_mat[, target_columns], na.rm = TRUE)
    }
  }
}

VST_count_mat_filtered <- VST_count_mat %>% 
  dplyr::select(-contains('LVG'),-contains('RVG'))

VST_count_mat_filtered <- VST_count_mat_filtered %>% 
  rownames_to_column('Gene') %>% 
  pivot_longer(cols = contains('VG'), values_to = 'Expression', names_to = 'rna_12_id')

# Add metadata
other_meta <- readxl::read_xlsx('../../Figure_and_Scripts/z_data_files/fixed_12snake_metadata_05.21.23.xlsx', sheet = 2, col_names = T) %>% 
  janitor::clean_names() %>% 
  dplyr::select(corrected_name, pop_id, lineage, sex, rna_12_id) %>% 
  dplyr::slice(1:13) %>% 
  dplyr::slice(-10) %>% 
  mutate(full_name = paste(corrected_name, lineage, pop_id, sex, sep = '_'))

other_meta <- other_meta %>% 
  # filter(lineage == 'viridis') %>% # keep only viridis, make sure you do this for the ATAC too.
  mutate(rna_12_id = paste0('VG_', rna_12_id))

GEX_variances <- VST_count_mat_filtered %>% 
  filter(rna_12_id %in% other_meta$rna_12_id) %>% # keep only information for viridis
  group_by(Gene) %>%
  summarise(GEX_Variance = var(Expression)) # calculate GEX variation

getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
p5 <- fp_variances %>% 
  left_join(ATAC_variances %>% ungroup() %>% dplyr::select(gene, ATAC_variance), by = 'gene') %>% 
  mutate(fpxATAC = fp_variance * ATAC_variance) %>% 
  left_join(GEX_variances, by = c('gene' = 'Gene')) %>% 
  na.omit() %>% 
  ggplot() +
  theme_bw() +
  geom_point(aes(x = ATAC_variance, y = GEX_Variance, colour = family), shape = 17, size = 2) +
  geom_smooth(aes(x = ATAC_variance, y = GEX_Variance), method = "lm", se = FALSE, color = "red", linewidth = 0.5, linetype = 'dashed') +
  ggpubr::stat_cor(aes(x = ATAC_variance, y = GEX_Variance, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), method = "pearson", label.x = 1, label.y = 9, size = 3.2) +
  labs(x = 'Composite accessibility variance', y = 'Expression variance') +
  scale_color_manual(values = getPalette(12)) +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal', 
        plot.margin = margin(5,15,5,5),
        legend.margin = margin(5,20,5,5),
        legend.spacing.x = unit(0, 'cm')) +
  guides(color = guide_legend(nrow = 3, title = NULL))

p6 <- plot_grid(p1, plot_grid(p4, p5, nrow = 1, align = 'tb', rel_widths = c(1, 0.8)), nrow = 2, align = 'lr', rel_heights = c(0.7,1))

#### 5: Pi and TFBS plots ####
pi_df <- read.table('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/3_ATACseq_comparisons/5_Enh_Pr_ATACseq_comparisons/Enh_Pr_pi_calculations/enh_pr_alignment_pi_Nov2023.txt', header = T)

### Optional filter ###
# Calculate the interquartile range
Q1 <- quantile(pi_df$Pi_corr, 0.25)
Q3 <- quantile(pi_df$Pi_corr, 0.75)
IQR <- Q3 - Q1
# Define the threshold for outliers
outlier_threshold <- 1.5 * IQR

# Filter the dataframe to remove outliers
pi_df <- pi_df %>%
  filter(Pi_corr >= Q1 - outlier_threshold & Pi_corr <= Q3 + outlier_threshold)
### End Optional filter ###

kruskal.test(Pi_corr ~ feature, data = pi_df) # significant
wilcox.test(Pi_corr ~ feature, data = pi_df) # significant

p_pi <- ggplot(pi_df, aes(x = feature, y = Pi_corr, fill = feature)) +
  geom_boxplot() +
  theme_bw() + 
  xlab(NULL) +
  ylab(expression(pi)) +
  scale_fill_manual(values = c('enhancer' = '#1B9E77', 'promoter' = '#D95F02')) +
  guides(fill = 'none') +
  theme(plot.margin = margin(15,5,15,5),
        axis.text.x = element_text(angle = 50, hjust = 1))

## TFBS chart ##
TFBS_mat_prom <- read.table('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/2_TFBS_Presence_Absence/4_Formatted_Outs_For_Rich/unformatted_TFBS_variant_mat_promoters_June2023.txt', header = T)
TFBS_mat_enh <- read.table('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/2_TFBS_Presence_Absence/4_Formatted_Outs_For_Rich/unformatted_TFBS_variant_mat_enhancers_June2023.txt', header = T)

TFBS_PA_sum_p <- as.data.frame(TFBS_mat_prom) %>% colSums() %>% table()
TFBS_PA_sum_e <- as.data.frame(TFBS_mat_enh) %>% colSums() %>% table()

# 3.56 x 2.82
p_tfbs <- as.data.frame(cbind(TFBS_PA_sum_p, TFBS_PA_sum_e)) %>% 
  mutate(NumIndiv = row_number()) %>% 
  dplyr::rename(promoter = TFBS_PA_sum_p) %>% 
  dplyr::rename(enhancer = TFBS_PA_sum_e) %>% 
  pivot_longer(cols = -NumIndiv, names_to = "feature", values_to = "frequency") %>% 
  group_by(feature) %>%
  mutate(relative_frequency = frequency / sum(frequency)) %>%
  ungroup() %>% 
  ggplot() +
  theme_bw() +
  geom_col(aes(x = NumIndiv, y = relative_frequency, fill = feature), position = 'dodge') +
  labs(x = "Number of individuals", y = "Relative\nfrequency") +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        plot.margin = margin(5,10,5,5),
        legend.margin = margin(-10,0,0,0)) +
  scale_x_continuous(breaks = 1:10, labels = 1:10) +
  scale_fill_manual(values = c('promoter' = '#D95F02', 'enhancer' = '#7570B3')) +
  ggtitle('Relative frequency of variable TFBS\nin venom promoters and enhancers')

p7 <- plot_grid(p_pi, p_tfbs, nrow = 1, rel_widths = c(0.4,1))

#### Final ####
plot_grid(p6, p7, nrow = 2, rel_heights = c(1,0.3))

plot_grid(p1, plot_grid(p4, p5, nrow = 1, align = 'tb', rel_widths = c(1, 0.8)), nrow = 2, align = 'lr', rel_heights = c(0.7,1))

plot_grid(p1, p4, p7, nrow = 3, align = 'tb', )

# Supplement
pi_df %>% 
  mutate(gene = str_split(names, "\\.")) %>% 
  unnest(gene) %>%
  relocate(gene, .after = names) %>% 
  mutate(gene = gsub("^PER\\d+_|^Promoter_", "", gene)) %>% 
  mutate(gene = paste0('Venom_', gene)) %>% 
  left_join(GEX_variances, by = c('gene' = 'Gene')) %>% 
  mutate(family = case_when(
    grepl('CRISP', gene) ~ 'CRISP',
    grepl('CTL', gene) ~ 'CTL', 
    grepl('LAAO', gene) ~ 'LAAO',
    grepl('SVSP', gene) ~ 'SVSP',
    grepl('SVMP', gene) ~ 'SVMP',
    grepl('PLA2', gene) ~ 'PLA2',
    grepl('vQC', gene) ~ 'vQC',
    grepl('VEGF', gene) ~ 'VEGF',
    grepl('Vespryn', gene) ~ 'Vespryn',
    grepl('EXO', gene) ~ 'EXO', 
    grepl('myo', gene) ~ 'myotoxin',
    TRUE ~ 'others'
  )) %>% 
  na.omit() %>% 
  ggplot() +
  theme_bw() +
  geom_point(aes(x = Pi_corr, y = GEX_Variance, colour = family), shape = 17, size = 2) +
  geom_smooth(aes(x = Pi_corr, y = GEX_Variance), method = "lm", se = FALSE, color = "red", linewidth = 0.5, linetype = 'dashed') +
  ggpubr::stat_cor(aes(x = Pi_corr, y = GEX_Variance, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), method = "pearson", label.x = 0.005, label.y = 3, size = 3.2) +
  labs(x = expression(pi), y = 'Expression variance') +
  scale_color_manual(values = getPalette(2)) +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal', 
        plot.margin = margin(5,15,5,5),
        legend.margin = margin(5,20,5,5),
        legend.spacing.x = unit(0, 'cm')) +
  guides(color = guide_legend(nrow = 1, title = 'Feature'))