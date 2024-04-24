library(BSgenome.Cviridis.custom.CroVir.noSeqNamesMods)
library(tidyverse)
library(rtracklayer)
library(ggcoverage)
library(scales)
library(gggenes)
library(cowplot)
library(chromVAR)
library(ggpmisc)

setwd('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/3_ATACseq_comparisons/6_GeneVignettes')
#### Linear modelling ####
SVSP9_peak_df <- read.table('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/5_Tables_for_Linear_Modelling/multiFeature_Gene_specific_tables/multifeature_tables_Oct2023_for_Rich_modelling/SVSP9_multifeature_table.txt', header = T, check.names = F)

# Initialize a list to store the model results
model_results_Top3 <- list()

# Loop through each peak column and fit a linear regression model
for (col in grep('Top3', colnames(SVSP9_peak_df), value = T)) {
  model <- lm(GEX_Venom_SVSP9 ~ get(col), data = SVSP9_peak_df)
  model_results_Top3[[col]] <- summary(model)
}

model_results_CTCF <- list()
for (col in grep('^CTCF', colnames(SVSP9_peak_df), value = T)) {
  model <- lm(GEX_Venom_SVSP9 ~ get(col), data = SVSP9_peak_df)
  model_results_CTCF[[col]] <- summary(model)
}


# Extract p-values and R-squared values
p_values_Top3 <- sapply(model_results_Top3, function(x) x$coefficients[2, "Pr(>|t|)"]) # Top3VarPeak_peak_254::scaffold-mi2:8537054-8537582 Top3VarPeak_peak_359::scaffold-mi2:8852477-8852853
p_values_CTCF <- sapply(model_results_CTCF, function(x) x$coefficients[2, "Pr(>|t|)"]) # CTCF_11973::scaffold-mi2:8923875-8924426

# Plot the results
p4 <- ggplot(SVSP9_peak_df, aes(x = get(names(p_values_Top3)[1]), y = GEX_Venom_SVSP9)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = 'dashed', formula = y~x) +
  stat_poly_eq(aes(x =  get(names(p_values_Top3)[1]), y = GEX_Venom_SVSP9),
               formula = y ~ x) +
  labs(
    x = "Accessibility at\nde-novo peak 1",
    y = "SVSP9 mRNA expression"
  ) +
  theme_bw()

p5 <- ggplot(SVSP9_peak_df, aes(x = get(names(p_values_Top3)[2]), y = GEX_Venom_SVSP9)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = 'dashed', formula = y~x) +
  stat_poly_eq(aes(x =  get(names(p_values_Top3)[2]), y = GEX_Venom_SVSP9),
               formula = y ~ x) +
  labs(
    x = "Accessibility at\nde-novo peak 2",
    y = "SVSP9 mRNA expression"
  ) +
  theme_bw()

p6 <- ggplot(SVSP9_peak_df, aes(x = get(names(p_values_CTCF)[6]), y = GEX_Venom_SVSP9)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = 'dashed', formula = y~x) +
  stat_poly_eq(aes(x =  get(names(p_values_CTCF)[6]), y = GEX_Venom_SVSP9),
               formula = y ~ x) +
  labs(
    x = "Accessibility at\nCTCF peak",
    y = "SVSP9 mRNA expression"
  ) +
  theme_bw()

#### Look at SVSP9 specifically ####
peakfile <- paste('../1_chromVAR/merged_peaks/allSample_mergedPeaks.bed', sep = '/')
peaks <- getPeaks(peakfile, sort_peaks = FALSE)
peaks <- resize(peaks, width = 500, fix = "center")

# Load data and plot
track.folder = '/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/3_ATACseq_comparisons/1_chromVAR/bams/'

# Add BAM info
sample.meta = data.frame(SampleName = gsub('.bam', '' , grep('bam$', list.files(track.folder), value = T))) # change from bw/bam
sample.meta <- sample.meta %>% 
  mutate(Type = str_extract(SampleName, "^[^_]+")) %>% 
  mutate(Group = str_extract(SampleName, 'viridis|concolor|lutosus|cerberus')) %>% 
  filter(Group != 'lutosus') %>% # remove lutosus for now because its ATAC-seq quality is bad
  #filter(Type != 'CV1096') %>% # remove CV1096 for now because we don't have a genome
  arrange(match(Group, c('viridis', 'concolor', 'cerberus')))

# load regions of interest from bam files
chrom = "scaffold-mi2"
start_pos = floor(8923875 / 100 ) * 100 # denovopeak 1 8537054, denovopeak2 8852477,  CTCF peak 8923875
end_pos = ceiling(8924426 / 500 ) * 500 # denovopeak 1 8537582, denovopeak2 8852853,  CTCF peak 8924426
track.df = LoadTrackFile(track.folder = track.folder, format = "bam", # change to bam/bw
                         bamcoverage.path = '/Users/sidgopalan/miniconda3/bin/bamCoverage',
                         meta.info = sample.meta,
                         single.nuc = T, # change to T for BAM, F for BigWig
                         single.nuc.region = paste0(chrom, ':', start_pos, '-', end_pos) # short hand for getting a region, change "region" to single.nuc.region for BAM, "region" for
)


peak_shift <- 0
rect <- as.data.frame(peaks) %>% 
  filter(seqnames=='scaffold-mi2') %>% 
  select(start, end) %>% 
  mutate(start = start + peak_shift, end = end + peak_shift) %>% 
  mutate(Type = 'CV1090') %>% 
  slice(612, 713, 747) # 612 = de novo peak 1, 713 = de novo peak 2, 747 = CTCF peak


p3 <-
ggplot() +
  theme_minimal() +
  geom_col(data = track.df, aes(x = start, y = score, fill = Group)) +
  annotate(geom = "rect",
           xmin = rect$start,
           xmax = rect$end,
           ymin = -Inf,
           ymax = +Inf,
           alpha = 0.2) +
  geom_gene_arrow(data = example_genes, inherit.aes = F, aes(xmin = start, xmax = end, y = 40, fill="blue"), arrowhead_height = grid::unit(2, "mm"), arrow_body_height = grid::unit(1, "mm")) +
  facet_wrap(~factor(Type, levels = sample.meta$Type), ncol = 1, strip.position = "right") + # scales = "free_y"
  scale_x_continuous(limits = c(start_pos,end_pos), expand = c(0,0), labels = label_number(scale = 1e-6), breaks = pretty_breaks(n=3)) +
  scale_y_continuous(breaks = pretty_breaks(n=2)) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        plot.margin = margin(10,5,5,5),
        #strip.text.y = element_blank()
        ) +
  guides(fill = guide_legend(nrow = 1, title = NULL)) +
  scale_fill_manual(values = c("viridis" = "#2E8B58",
                               "concolor" = "#D9A528",
                               "cerberus" = "black"),
                    labels = c(expression(italic("C. viridis")),
                               expression(italic("C. concolor")),
                               expression(italic("C. cerberus")),
                               breaks = c("viridis", "concolor", "cerberus"))) +
  labs(x = paste0("pos. on ", chrom, ' (Mb)'), y = "ATAC-seq read density") +
  ggtitle("CTCF peak")


### Add Expression Barplot ###
VST_count_mat <- read.table('../../z_data_files/RNAseq_NormalizedCounts_noOutliers_06.29.23.txt', header=T) # use either the VST normalized mat or the normalized count mat

# Log Normalize counts if using normalized counts
# VST_count_mat =  log(VST_count_mat + 1)

# Filter for Venom genes, also combine LVG and RVG (take the average). Note that LVG_13 failed and is missing.
VST_count_mat <- VST_count_mat %>% 
  select(matches('LVG|RVG')) %>% 
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
  select(-contains('LVG'),-contains('RVG'))

samples <- read.table('../../z_data_files/rnaSampleInfo_06.29.21.txt', header = T)
samples <- samples %>% 
  mutate(rna_12_id = str_extract(Sample_ID, str_extract(Sample_ID, "\\d+")))

other_meta <- readxl::read_xlsx('../../z_data_files/fixed_12snake_metadata_05.21.23.xlsx', sheet = 2, col_names = T) %>% 
  select(Corrected_name, 7) %>% 
  dplyr::slice(1:13) %>% 
  janitor::clean_names() %>% 
  dplyr::slice(-10)

samples <- samples %>% 
  left_join(other_meta, by = 'rna_12_id') %>% 
  filter(grepl('RVG', Sample_ID)) %>% 
  select(Lineage, corrected_name, rna_12_id) %>% 
  dplyr::rename(SampleID = corrected_name) %>% 
  relocate(SampleID, .before=Lineage) %>% 
  filter(!grepl('CV1084', SampleID)) %>% # remove bad ATACseq viridis sample
  #filter(!grepl('CV1096', SampleID)) %>% # remove missing genome viridis sample
  filter(Lineage != 'lutosus') %>% # remove lutosus
  mutate(rna_12_id = paste0('VG_', rna_12_id)) %>% 
  arrange(match(SampleID, sample.meta$Type))

VST_count_mat_filtered <- VST_count_mat_filtered %>% 
  select(any_of(samples$rna_12_id))
colnames(VST_count_mat_filtered) <- samples$SampleID[match(colnames(VST_count_mat_filtered), samples$rna_12_id)]
VST_count_mat_filtered %>% 
  select(one_of(samples$SampleID), everything())

VST_count_mat_filtered <- VST_count_mat_filtered %>% 
  filter(rownames(.) == 'Venom_SVSP9') %>% 
  t() %>% 
  as.data.frame()
VST_count_mat_filtered$SampleID <- rownames(VST_count_mat_filtered)
p7 <- VST_count_mat_filtered %>% 
  left_join(samples %>% select(SampleID, Lineage), by = 'SampleID') %>%
  ggplot(., aes(x = SampleID, y = Venom_SVSP9, fill = Lineage)) +
  theme_bw() +
  geom_col() +
  labs(y = "SVSP9\nmRNA expression", x = NULL) +
  coord_flip() +
  scale_x_discrete(limits = rev(samples$SampleID),
                   labels = rev(samples$SampleID)) +
  scale_y_continuous(expand = c(0,0), labels = label_number(scale = 1e-3)) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_fill_manual(values =c("viridis" = "#2E8B58",
                              "concolor" = "#D9A528",
                              "cerberus" = "black")) +
  guides(fill = 'none')

# 6.69 x 
p_1x <- plot_grid(p4, p4, p5, p6, nrow = 1, align = 'hv', axis = 'tb') # first p4 is just to make space
p_2x <- plot_grid(p1, p2, p3, p7, nrow = 1, align = 'hv', axis = 'tb', rel_widths = c(2,2,2,1))
plot_grid(p_1x, p_2x, nrow = 2, rel_heights = c(1,2))

