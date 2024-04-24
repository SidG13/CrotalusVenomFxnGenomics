#### This script is to understand differences in RNA and protein expression across samples ####

library(cowplot)
library(tidyverse)
library(matrixStats)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(readxl)
library(ggpmisc)
library(compositions)

setwd('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/1_RNAseq_and_TF_expression_variation')

#### RNA heatmap ####
# RNAseq VST normalized count matrix
VST_count_mat <- read.table('../z_data_files/RNAseq_VSTNormalizeCounts_noOutliers_05.23.23.txt', header=T)

# Filter for TFs and Venom genes, also combine LVG and RVG (take the average). Note that LVG_13 failed and is missing.
VST_count_mat <- VST_count_mat %>% 
  select(matches('LVG|RVG')) # %>% 
  # filter(grepl('Venom', row.names(.)))


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

# Make boxplots of NVP expression vs Venom expression
NVPs_mat <- rbind(VST_count_mat_filtered[grep("nbisL1-mrna-", rownames(VST_count_mat_filtered)),], # bdefensins
              VST_count_mat_filtered[grep("^ADAM", rownames(VST_count_mat_filtered)),], # ADAMs
              VST_count_mat_filtered[grep("^PLA2(?!R)", rownames(VST_count_mat_filtered), perl = TRUE), ], # PLA2s that aren't PLA2 receptors
              VST_count_mat_filtered[grep("^PRSS", rownames(VST_count_mat_filtered)),]) # serine proteases

d1 <- as.data.frame(rowVars(as.matrix(NVPs_mat))) # across samples
d2 <- as.data.frame(rowVars(as.matrix(NVPs_mat[, !(colnames(NVPs_mat) %in% c('VG_6','VG_7', 'VG_8'))]))) # within C. viridis
colnames(d1) <- "Variances"
d1$Group <- "Across"
colnames(d2) <- "Variances"
d2$Group <- "Within"

Venoms_mat <- VST_count_mat_filtered %>% 
  filter(str_detect(row.names(.), 'Venom'))
d3 <- as.data.frame(rowVars(as.matrix(Venoms_mat))) # across samples
d4 <- as.data.frame(rowVars(as.matrix(Venoms_mat[, !(colnames(Venoms_mat) %in% c('VG_6','VG_7', 'VG_8'))]))) # within C. viridis
colnames(d3) <- "Variances"
d3$Group <- "Across"
colnames(d4) <- "Variances"
d4$Group <- "Within"

variances_df <- rbind(d1, d2, d3, d4)
variances_df <- variances_df %>% 
  rownames_to_column("gene") %>% 
  mutate(Group = if_else(grepl("Venom", gene), paste0("Venom_",Group), paste0("NVP_",Group))) 

variances_df %>% 
  ggplot(aes(x = Group, y = Variances, fill = Group)) +
  geom_boxplot() +
  # stat_compare_means(method = "t.test", label = "p.format") +  # Display the p-value on the plot
  theme_bw() +
  scale_fill_manual(values = c("grey60", "grey60", "forestgreen", "forestgreen")) +
  scale_x_discrete(labels = c('NVPs (across)', 'NVPs (within)', 'Venoms (across)', 'Venoms (within)')) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# T tests
var.test(Variances ~ Group, data = variances_df %>% filter(Group %in% c("Venom_Across", "NVP_Across")))
var.test(Variances ~ Group, data = variances_df %>% filter(Group %in% c("Venom_Within", "NVP_Within")))
var.test(Variances ~ Group, data = variances_df %>% filter(Group %in% c("Venom_Within", "Venom_Across")))
t.test(Variances ~ Group, data = variances_df %>% filter(Group %in% c("Venom_Across", "NVP_Across")), var.equal = F) # < 0.01
t.test(Variances ~ Group, data = variances_df %>% filter(Group %in% c("Venom_Within", "NVP_Within")), var.equal = F) # < 0.01
t.test(Variances ~ Group, data = variances_df %>% filter(Group %in% c("Venom_Within", "Venom_Across")), var.equal = F) # < n.s.
#

# Heatmap of NVPs
NVPs_mat %>% 
  rownames_to_column(var = "NVP") %>% 
  pivot_longer(cols = VG_1:VG_13, 
               names_to = "VG_Sample",
               values_to = "Average_Exp") %>% 
  ggplot(., aes(y=NVP, x=VG_Sample, fill=Average_Exp)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  scale_x_discrete(labels=species_labs,
                   limits=c(sample_info$Sample_ID)) +
  theme(
        #legend.direction = 'horizontal',
        #legend.position = 'bottom',
        legend.text = element_text(size=10),
        legend.title = element_text(size=10, angle = 90),
        legend.title.align = 0.5,
        legend.direction = "vertical",
        axis.text.x = ggtext::element_markdown(size=10, color = species_cols, hjust = 1, vjust = 1, angle = 35),
        legend.margin = margin(-30,5,0,5)) +
  guides(fill=guide_colorbar(barwidth = 1, title.position = "right")) +
  labs(fill = "VST norm. exp.")

sample_info <- read.table('../z_data_files/rnaSampleInfo_06.29.21.txt', header = T)
sample_info <- sample_info %>% 
  filter(str_detect(Sample_ID, 'RVG')) %>% 
  mutate(Sample_ID = gsub('R', '', Sample_ID)) %>% 
  mutate(Population = fct_relevel(Population, 
                                  "South", "Mid", "North", "Other")) %>% 
  arrange(Lineage, Population, Sample_ID)

# Format colours and italics, and order
species_cols <- ifelse(sample_info$Lineage == 'viridis', "#2E8B58",
                       ifelse(sample_info$Lineage == 'cerberus', "black", 
                              ifelse(sample_info$Lineage == 'lutosus', "#8B481F", "#D9A528")))
species_labs <- paste0('C. ', gsub('-Other', '', paste(sample_info$Lineage, sample_info$Population, sep = '-')))
species_labs <- gsub('C','*C', gsub('viridis','viridis*', gsub('concolor','concolor*', gsub('lutosus','lutosus*', gsub('cerberus','cerberus*', species_labs)))))

Venom_VST_mat <- VST_count_mat_filtered %>% 
  filter(str_detect(row.names(.), 'Venom')) %>% 
  rownames_to_column(var = "Toxin") %>% 
  mutate(Toxin = gsub("Venom_", '', Toxin)) %>%
  mutate(Toxin = gsub("CRISP_", 'CRISP', Toxin)) %>% 
  mutate(Toxin = gsub("CTL_", 'CTL', Toxin)) %>% 
  filter(str_detect(Toxin, 'ADAM|CRISP3|CRISP4|CTL1|CTL4|CTL5|CTL6|EXO|LAAO1|LAAO2|VEGF2|vQC|SVMP11', negate = T)) %>% # remove lowly expressed genes in all samples
  pivot_longer(cols = VG_1:VG_13, 
               names_to = "VG_Sample",
               values_to = "Average_Exp")
toxin_order <- c('BPP', 'CRISP1', 'CRISP2', 'CTL2', 'LAAO3', 'myotoxin', 'ohanin', 'VEGF1', 'PLA2A1', 'PLA2B1', 'PLA2C1', 'PLA2K', paste0('SVMP',c(1:10)), paste0('SVSP',1:11))

Venom_VST_matvars <- Venom_VST_mat %>% 
  group_by(Toxin) %>% 
  summarize(Variance = var(Average_Exp)) %>% 
  arrange(desc(Variance))

p1 <- Venom_VST_mat %>% 
  ggplot(., aes(x=Toxin, y=VG_Sample, fill=Average_Exp)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  scale_y_discrete(labels=species_labs,
                   limits=c(sample_info$Sample_ID)) +
  scale_x_discrete(labels=Venom_VST_matvars$Toxin,
                   limits=Venom_VST_matvars$Toxin) +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #legend.direction = 'horizontal',
        #legend.position = 'bottom',
        legend.text = element_text(size=10),
        legend.title = element_text(size=10, angle = 90),
        legend.title.align = 0.5,
        legend.direction = "vertical",
        axis.text.y = ggtext::element_markdown(size=10, color = species_cols),
        legend.margin = margin(-30,5,0,5)) +
  guides(fill=guide_colorbar(barwidth = 1, title.position = "right")) +
  labs(fill = "VST norm. exp.")


Venom_VST_matvars_withinViridis <- Venom_VST_mat %>% 
  filter(!VG_Sample %in% c('VG_6', 'VG_7', 'VG_8')) %>% 
  group_by(Toxin) %>% 
  summarize(Variance = var(Average_Exp))


p_vars <- Venom_VST_matvars_withinViridis %>% 
  ggplot(., aes(x=Toxin,y=0,fill=Variance)) +
  geom_tile() +
  scale_fill_viridis(trans = 'sqrt') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.margin = margin(5,5,-30,5)) +
  scale_x_discrete(labels=Venom_VST_matvars$Toxin,
                   limits=Venom_VST_matvars$Toxin)

p1_x2 <- plot_grid(p_vars, p1, nrow = 2, align = 'hv', axis = 'lr', rel_heights = c(0.4,1))

#### Protein figures ####
prot <- read_excel('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/z_data_files/2023_Aug9_QEHF_Castoe.xlsx', sheet = 1) # QEHF protein

# Convert AS names to standard names (version 1)
name_conversion <- data.frame(
  stringsAsFactors = FALSE,
  Anthony_name = c("C_o_cerberus_AZ_Kingman",
                   "C_o_cerberus_AZ_NW_Kingman","C_o_concolor_UT_Duchesne",
                   "C_o_lutosus_528_UT_Beaver_Co","C_o_lutosus_UT_Beaver_Co",
                   "C_o_tigris_AZ_Pima_Co","C_v_viridis_CO_Hwy93",
                   "C_v_viridis_CO_Otero_Co","C_v_viridis_CO_Weld_Co",
                   "C_v_viridis_MT","C_v_viridis_NM_6","C_v_viridis_NM_Luna_Co",
                   "C_v_viridis_NM_Socorro","C_v_viridis_CO_San_Miguel_Co"),
  New_name = c("NONE",
               "CV1090_cerberus_Other_M","CV0985_concolor_Other_F","CV0987_lutosus_Other_F",
               "NONE","NONE","CV1096_viridis_North_F",
               "CV1084_viridis_Mid_M","CV1087_viridis_North_F",
               "CV0857_viridis_North_M","NONE","CV1086_viridis_South_M",
               "CV1089_viridis_South_M","CV1081_viridis_Mid_M")
)

# # Convert anthony names to standard names (version 2)
# name_conversion <- data.frame(
#   stringsAsFactors = FALSE,
#                          Anthony_name = c("20230809_SCP_CFS_cerberus_11_6ul_Slot2_12_1_4851",
#                                           "20230809_SCP_CFS_cerberus_9_6ul_Slot2_10_1_4849",
#                                           "20230809_SCP_CFS_concolor_13_6ul_Slot2_14_1_4853",
#                                           "20230809_SCP_CFS_cvv_1_6ul_Slot2_1_1_4839",
#                                           "20230809_SCP_CFS_cvv_2L_6ul_Slot2_2_1_4840",
#                                           "20230809_SCP_CFS_cvv_2R_6ul_Slot2_3_1_4842",
#                                           "20230809_SCP_CFS_cvv_3_6ul_Slot2_4_1_4843",
#                                           "20230809_SCP_CFS_cvv_4_6ul_Slot2_5_1_4844",
#                                           "20230809_SCP_CFS_cvv_5_6ul_Slot2_6_1_4845",
#                                           "20230809_SCP_CFS_cvv_6_6ul_Slot2_7_1_4846",
#                                           "20230809_SCP_CFS_cvv_7_6ul_Slot2_8_1_4847",
#                                           "20230809_SCP_CFS_cvv_8_6ul_Slot2_9_1_4848",
#                                           "20230809_SCP_CFS_lutosus_10_6ul_Slot2_11_1_4850",
#                                           "20230809_SCP_CFS_lutosus_12_6ul_Slot2_13_1_4852",
#                                           "20230809_SCP_CFS_tigris_14_6ul_Slot2_15_1_4854"),
#                              New_name = c("NONE",
#                                           "CV1090_cerberus_Other_M","CV0985_concolor_Other_F",
#                                           "CV1096_viridis_North_F",
#                                           "CV1084_viridis_Mid_M_L","CV1084_viridis_Mid_M_R",
#                                           "CV1081_viridis_Mid_M",
#                                           "CV1087_viridis_North_F","CV1089_viridis_South_M",
#                                           "NONE","CV0857_viridis_North_M",
#                                           "CV1086_viridis_South_M",
#                                           "CV0987_lutosus_Other_F","NONE","NONE")
#                    )

prot <- prot %>% 
  select(-51) %>% # remove indistinguishable proteins column # 51
  pivot_longer(cols = 9:50, names_sep = '\\s', names_to = c("SampleID", "Feature"), values_to = "Value") %>% #9:50
  left_join(name_conversion, by = c('SampleID'='Anthony_name')) %>% 
  filter(!New_name == 'NONE') %>% # remove samples for which we have no study samples
  select(-SampleID) %>% 
  rename(SampleID = New_name) %>% 
  filter(!is.na(Protein))

protein_to_orthos <- read.table('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/z_data_files/Cvv_GTF_to_converted_names_05.22.23.txt', header=T)  
protein_to_orthos <- protein_to_orthos %>% 
  filter(!grepl('fgenesh', gtf_gene)) %>% # remove fgenesh annotations
  select(crovir_transcript, converted_id_no_dups) %>% 
  mutate(crovir_transcript = gsub('transcript', 'protein', crovir_transcript)) %>% 
  mutate(crovir_transcript = ifelse(converted_id_no_dups == 'Venom_myotoxin', 'crovir-protein-myotoxin', crovir_transcript))

prot <- prot %>% 
  left_join(protein_to_orthos, by = c('Protein' = 'crovir_transcript')) %>% 
  select(-Protein) %>% 
  rename(Protein = converted_id_no_dups) %>% 
  filter(grepl('Venom', Protein)) # keep only Venoms 

#### RNA ####
# RNAseq VST normalized count matrix
VST_count_mat <- read.table('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/z_data_files/RNAseq_VSTNormalizeCounts_noOutliers_05.23.23.txt', header=T)

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

new_names <- data.frame(
  stringsAsFactors = FALSE,
  VG_name = c("VG_1","VG_2","VG_3","VG_4","VG_5","VG_6",
              "VG_7","VG_8","VG_9","VG_11","VG_12",
              "VG_13"),
  New_name = c("CV1084_viridis_Mid_M","CV1081_viridis_Mid_M",
               "CV1089_viridis_South_M",
               "CV0857_viridis_North_M","CV1087_viridis_North_F",
               "CV0987_lutosus_Other_F","CV0985_concolor_Other_F",
               "CV1090_cerberus_Other_M",
               "CV1086_viridis_South_M","CV1095_viridis_North_M",
               "CV1082_viridis_South_M","CV1096_viridis_North_F")
)


VST_count_mat_filtered <- VST_count_mat_filtered %>% 
  filter(str_detect(row.names(.), 'Venom')) %>% 
  rownames_to_column(var = "Toxin") %>% 
  mutate(Toxin = gsub("CRISP_", 'CRISP', Toxin)) %>% 
  mutate(Toxin = gsub("CTL_", 'CTL', Toxin)) %>% 
  filter(str_detect(Toxin, 'ADAM|CRISP3|CRISP4|CTL1|CTL4|CTL5|CTL6|EXO|LAAO1|LAAO2|VEGF2|SVMP11', negate = T)) %>% # remove lowly expressed genes in all samples
  pivot_longer(cols = VG_1:VG_13, 
               names_to = "VG_Sample",
               values_to = "Average_Exp") %>% 
  left_join(new_names, by = c('VG_Sample'='VG_name')) %>% 
  select(-VG_Sample) %>% 
  rename(VG_Sample = New_name)

# Add GEX to prot
prot_GEX <- prot %>% 
  left_join(VST_count_mat_filtered, by = c('SampleID'='VG_Sample', 'Protein'='Toxin')) %>% 
  separate(SampleID, into = c('CV_ID','Lineage', 'Population','Sex'), sep = '_', remove = F) %>% 
  mutate(Venom_family = case_when(grepl('SVMP', Protein) ~ 'SVMP',
                                  grepl('SVSP', Protein) ~ 'SVSP',
                                  grepl('PLA2', Protein) ~ 'PLA2',
                                  grepl('ADAM', Protein) ~ 'ADAM',
                                  grepl('CRISP', Protein) ~ 'CRISP',
                                  #grepl('CTL', Protein) ~ 'CTL',
                                  grepl('EXO', Protein) ~ 'EXO',
                                  grepl('LAAO', Protein) ~ 'LAAO',
                                  grepl('myotoxin', Protein) ~ 'myotoxin',
                                  grepl('BPP', Protein) ~ 'BPP',
                                  TRUE ~ 'others'))

# Intensity for quants

prot_GEX_transformed <- prot_GEX %>% 
  filter(Feature == 'Intensity') %>% 
  filter(Average_Exp > 10) %>% 
  select(CV_ID, Lineage, Protein, Feature, Value, Average_Exp, Venom_family) %>% 
  group_by(CV_ID) %>% 
  mutate(scaled_prot_exp = Value / sum(Value)) %>%# sum protein exps to 1
  mutate(scaled_gene_exp = Average_Exp / sum(Average_Exp)) %>% # sum gene exps to 1
  ungroup() %>% 
  mutate(clr_scaled_prot_exp = compositions::clr(scaled_prot_exp)) %>% 
  mutate(clr_scaled_gene_exp = compositions::clr(scaled_gene_exp))

prot_GEX_transformed %>% 
  ggplot() +
  geom_point(aes(x=clr_scaled_gene_exp, y=clr_scaled_prot_exp, color = Venom_family), size = 2.5, alpha = 0.4) +
  geom_smooth(
    aes(x = clr_scaled_gene_exp, y = clr_scaled_prot_exp),
    method = "lm",
    se = FALSE,  # Do not display confidence intervals
    color = "black",
    linetype = "dashed",
    # formula = y ~ poly(x, 3, raw = TRUE) # 3rd degree polynomial
    formula = y ~ x # linear
  ) +
  stat_poly_eq(aes(x = clr_scaled_gene_exp, y = clr_scaled_prot_exp),
                   formula = y ~ x) +
  scale_color_brewer(palette = 'Dark2') +
  xlab('clr(gene expression)') +
  ylab('clr(protein expression)') +
  theme_bw() +
  theme(legend.position = 'bottom', 
        legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 3))

g1 <- glm(formula = Value~Average_Exp+Lineage+Sex+Venom_family,
          data = subset(prot_GEX, Feature == 'Intensity'),
          family = 'poisson')
summary(g1)


#### Protein PCA ####
prot_PCA <- prot_GEX %>% 
  select(SampleID, Protein, Feature, Value) %>% 
  filter(Feature == 'Intensity') %>% 
  select(-Feature) %>% 
  pivot_wider(names_from = SampleID, values_from = Value) %>% 
  column_to_rownames(var = 'Protein')

pca <- prcomp(prot_PCA, scale. = TRUE)

# Create a data frame with sample names and corresponding PCA coordinates
pca_data <- data.frame(Sample = colnames(prot_PCA),
                       PC1 = pca$rotation[, 1],
                       PC2 = pca$rotation[, 2])
pca_data <- pca_data %>% 
  separate(Sample, into = c('CV_ID','Lineage', 'Population','Sex'), sep = '_', remove = F) %>% 
  mutate(LineageSpecific = ifelse(Population != 'Other', paste0(Lineage, Population), Lineage))

p3 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = LineageSpecific)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(x = paste0("PC1 ", '(', summary(pca)$importance[2,1]*100, '%)'), y = paste0("PC2 ", '(', summary(pca)$importance[2,2]*100, '%)')) +
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
  guides(color = guide_legend(
                              override.aes = list(shape = 16),
                              nrow = 3)) +
  theme(legend.direction = 'horizontal',
        legend.position = 'bottom',
        plot.margin = margin(20,30,0,30),
        legend.margin = margin(-10,0,0,0),
        legend.title = element_blank())

#### Protein Pie charts ####
pie1 <- prot_GEX %>% 
  filter(Feature == 'Intensity') %>% 
  filter(Lineage == 'viridis') %>% 
  filter(Population == 'North') %>% 
  group_by(Venom_family) %>% 
  summarize(Mean_Abundance = mean(Value)) %>% 
  ggplot(., aes(x = "", y = Mean_Abundance, fill = Venom_family)) +
  geom_col(color = "black") +
  #geom_text(aes(label = Venom_family),
  #          position = position_stack(vjust = 0.5), size = 6) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette="Dark2") +
  theme_void() +
  theme(legend.position = 'none')

pie2 <- prot_GEX %>% 
  filter(Feature == 'Intensity') %>% 
  filter(Lineage == 'viridis') %>% 
  filter(Population == 'Mid') %>% 
  group_by(Venom_family) %>% 
  summarize(Mean_Abundance = mean(Value)) %>% 
  ggplot(., aes(x = "", y = Mean_Abundance, fill = Venom_family)) +
  geom_col(color = "black") +
  #geom_text(aes(label = Venom_family),
  #          position = position_stack(vjust = 0.5), size = 6) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette="Dark2") +
  theme_void() +
  theme(legend.position = 'none')

pie3 <- prot_GEX %>% 
  filter(Feature == 'Intensity') %>% 
  filter(Lineage == 'viridis') %>% 
  filter(Population == 'South') %>% 
  group_by(Venom_family) %>% 
  summarize(Mean_Abundance = mean(Value)) %>% 
  ggplot(., aes(x = "", y = Mean_Abundance, fill = Venom_family)) +
  geom_col(color = "black") +
  #geom_text(aes(label = Venom_family),
  #          position = position_stack(vjust = 0.5), size = 6) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette="Dark2") +
  theme_void() +
  theme(legend.position = 'none')

pie4 <- prot_GEX %>% 
  filter(Feature == 'Intensity') %>% 
  filter(Lineage == 'lutosus') %>% 
  filter(Population == 'Other') %>% 
  group_by(Venom_family) %>% 
  summarize(Mean_Abundance = mean(Value)) %>% 
  ggplot(., aes(x = "", y = Mean_Abundance, fill = Venom_family)) +
  geom_col(color = "black") +
  #geom_text(aes(label = Venom_family),
  #          position = position_stack(vjust = 0.5), size = 6) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette="Dark2") +
  theme_void() +
  theme(legend.position = 'none')

pie5 <- prot_GEX %>% 
  filter(Feature == 'Intensity') %>% 
  filter(Lineage == 'concolor') %>% 
  filter(Population == 'Other') %>% 
  group_by(Venom_family) %>% 
  summarize(Mean_Abundance = mean(Value)) %>% 
  ggplot(., aes(x = "", y = Mean_Abundance, fill = Venom_family)) +
  geom_col(color = "black") +
  #geom_text(aes(label = Venom_family),
  #          position = position_stack(vjust = 0.5), size = 6) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette="Dark2") +
  theme_void() +
  theme(legend.position = 'none')

pie6 <- prot_GEX %>% 
  filter(Feature == 'Intensity') %>% 
  filter(Lineage == 'cerberus') %>% 
  filter(Population == 'Other') %>% 
  group_by(Venom_family) %>% 
  summarize(Mean_Abundance = mean(Value)) %>% 
  ggplot(., aes(x = "", y = Mean_Abundance, fill = Venom_family)) +
  geom_col(color = "black") +
  #geom_text(aes(label = Venom_family),
  #          position = position_stack(vjust = 0.5), size = 6) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette="Dark2") +
  theme_void() +
  theme(legend.position = 'none')


# Individual pies
pie7 <- prot_GEX %>% 
  filter(Feature == 'Intensity') %>% 
  filter(CV_ID == 'CV1089') %>% # South
  group_by(Venom_family) %>% 
  summarize(Mean_Abundance = mean(Value)) %>% 
  ggplot(., aes(x = "", y = Mean_Abundance, fill = Venom_family)) +
  geom_col(color = "black") +
  #geom_text(aes(label = Venom_family),
  #          position = position_stack(vjust = 0.5), size = 6) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette="Dark2") +
  theme_void() +
  theme(legend.position = 'none')

pie8 <- prot_GEX %>% 
  filter(Feature == 'Intensity') %>% 
  filter(CV_ID == 'CV1081') %>% # Mid
  group_by(Venom_family) %>% 
  summarize(Mean_Abundance = mean(Value)) %>% 
  ggplot(., aes(x = "", y = Mean_Abundance, fill = Venom_family)) +
  geom_col(color = "black") +
  #geom_text(aes(label = Venom_family),
  #          position = position_stack(vjust = 0.5), size = 6) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette="Dark2") +
  theme_void() +
  theme(legend.position = 'none')

pie9 <- prot_GEX %>% 
  filter(Feature == 'Intensity') %>% 
  filter(CV_ID == 'CV1087') %>% # North
  group_by(Venom_family) %>% 
  summarize(Mean_Abundance = mean(Value)) %>% 
  ggplot(., aes(x = "", y = Mean_Abundance, fill = Venom_family)) +
  geom_col(color = "black") +
  #geom_text(aes(label = Venom_family),
  #          position = position_stack(vjust = 0.5), size = 6) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette="Dark2") +
  theme_void() +
  theme(legend.position = 'none')

plot_grid(pie7,pie8,pie9,nrow = 1)

#### 6.69 x 6.71 ####
p_1x <- plot_grid(p1_x2, p_pies, nrow=2, axis = 'lr', rel_heights = c(1.8,1))
plot_grid(p_1x, plot_grid(p3,p2, align = 'hv'), axis = 'lr', nrow = 2, rel_heights = c(1.3,1))

