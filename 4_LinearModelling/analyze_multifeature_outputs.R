# Analyze/visualize outputs of feature importance and P-values from Rich #
library(tidyverse)
library(cowplot)
library(colorspace)
library(viridis)

setwd('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/5_Tables_for_Linear_Modelling/ModelingResults_19Oct23')

file_paths <- list.files('.', pattern = '.txt', full.names = TRUE)
data_list <- list()

for (file_name in file_paths) {
  # Read the file into a data frame
  data <- read.table(file_name, header = TRUE, sep = "\t")  # Modify the separator if needed
  
  # Add a "File" column with the file name
  data$File <- file_name
  
  # Append the data frame to the list
  data_list <- c(data_list, list(data))
}

# Step 3: Combine data frames into one
combined_data <- do.call(rbind, data_list)

# Optional: If you want to reset row names, use the following line
rownames(combined_data) <- NULL

combined_data$File <- gsub('\\./RESULTS_', '', gsub('_multifeature_table.txt', '', combined_data$File))
colnames(combined_data)[4] <- 'Gene'

# All gene data
#p1 <- 
combined_data %>% 
  # filter(FeaturePvalue < 0.05) %>% 
  mutate(AbsFeatureCoefficients = abs(FeatureCoefficients)) %>% 
  ggplot(aes(x=Gene, y=FeatureName)) +
  geom_point(aes(size = AbsFeatureCoefficients, color = -log10(FeaturePvalue))) +
  #scale_color_viridis_c(option = 'magma') +
  scale_color_gradient2(low = 'grey50', high = 'tomato2', midpoint = -log10(0.05), breaks = c(0.5,1,1.5,2)) +
  labs(
    title = "Gene-feature dotplot of p-value and effect size",
    x = NULL,
    y = NULL,
    color = expression(-log[10]("p-value")),
    size = "abs(coefficient)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  scale_y_discrete(labels = c('Accessibility at\nCTCF binding sites', 'CRE TF occupancy\n(binarized footprints)', 'CRE genotype', 
                              'Enhancer accessibility', 'Promoter accessibility', 'TF expression',
                              'Accessibility at 3\nmost-variable peaks')) +
  scale_x_discrete(limits = venom_metrics$Toxin, labels = venom_metrics$Toxin)


# Plot
plot_grid(p2,p1,ncol=1)

# Plot individual gene tracks
p1 + scale_x_discrete(limits = 'SVSP2', breaks = 'SVSP2')
p1 + scale_x_discrete(limits = 'SVMP6', breaks = 'SVMP6')
p1 + scale_x_discrete(limits = 'SVMP4', breaks = 'SVMP4')
p1 + scale_x_discrete(limits = 'SVSP9', breaks = 'SVSP9')


#### Plot gene variance and mean GEX single rows
VST_count_mat <- read.table('../../z_data_files/RNAseq_VSTNormalizeCounts_noOutliers_05.23.23.txt', header=T)

# Filter for TFs and Venom genes, also combine LVG and RVG (take the average). Note that LVG_13 failed and is missing.
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

sample_info <- read.table('../../z_data_files/rnaSampleInfo_06.29.21.txt', header = T)
sample_info <- sample_info %>% 
  filter(str_detect(Sample_ID, 'RVG')) %>% 
  mutate(Sample_ID = gsub('R', '', Sample_ID)) %>% 
  mutate(Population = fct_relevel(Population, 
                                  "South", "Mid", "North", "Other")) %>% 
  arrange(Lineage, Population, Sample_ID)

Venom_VST_mat <- VST_count_mat_filtered %>% 
  filter(str_detect(row.names(.), 'Venom')) %>% 
  rownames_to_column(var = "Toxin") %>% 
  mutate(Toxin = gsub("Venom_", '', Toxin)) %>%
  mutate(Toxin = gsub("CRISP_", 'CRISP', Toxin)) %>% 
  mutate(Toxin = gsub("CTL_", 'CTL', Toxin)) %>% 
  #filter(str_detect(Toxin, 'ADAM|CRISP3|CRISP4|CTL1|CTL4|CTL5|CTL6|EXO|LAAO1|LAAO2|VEGF2|vQC|SVMP11', negate = T)) %>% # remove lowly expressed genes in all samples
  pivot_longer(cols = VG_1:VG_13, 
               names_to = "VG_Sample",
               values_to = "Average_Exp")
#toxin_order <- c('CRISP1','CRISP2','CTL2','LAAO3','PLA2A1','PLA2B1','PLA2C1','PLA2K', paste0('SVMP', seq(1,10)), paste0('SVSP', seq(1,11)), 'VEGF1')

venom_metrics <- Venom_VST_mat %>% 
  filter(Toxin %in% toxin_order) %>% 
  group_by(Toxin) %>% 
  summarize(Variance = var(Average_Exp),
            Mean = mean(Average_Exp)) %>% 
  arrange(desc(Variance))

# Variance
venom_metrics %>% 
  ggplot(., aes(x=Toxin,y=0,fill=Variance)) +
  geom_tile() +
  scale_fill_viridis() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.margin = margin(5,5,5,5)) +
  scale_x_discrete(labels=venom_metrics$Toxin,
                   limits=venom_metrics$Toxin)

# Average Expr
venom_metrics %>% 
  ggplot(., aes(x=Toxin,y=0,fill=Mean)) +
  geom_tile() +
  scale_fill_viridis() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.margin = margin(5,5,5,5)) +
  scale_x_discrete(labels=venom_metrics$Toxin,
                   limits=venom_metrics$Toxin)
