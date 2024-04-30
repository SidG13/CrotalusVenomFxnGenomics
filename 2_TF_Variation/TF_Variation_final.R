#### This script is to combine various TF explanatory figures together ####
## Genome Biology Figures ##
# width of 85 mm (3.34") for half page width figure
#width of 170 mm (6.69") for full page width figure
# maximum height of 225 mm (8.85") for figure and legend

library(cowplot)
library(tidyverse)
library(matrixStats)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(ggrepel)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggpubr)
library(ggvenn)
library(Hmisc)
library(corrplot)

#### 1: TF Expression Heatmap ####
setwd('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Final_Figure_Assembly/2_TF_Variation')

# List of 161 candidate TFs from Perry 2022 GR
TF_list <- c("SIRT1", "RERE", "NOC2L", "IKZF5", "SUFU", "PPRC1", "HES6", "PIAS2", "BMPR1A", "TP73", "RORA", "KLF13", "TSC22D4", "ZNF592", "CAMTA1", "PRDM2", "AEBP2", "TLE1", "ZNF710", "FOXC2", "CTCF", "ZBED4", "ZNF800", "SREBF2", "ATF4", "SOX10", "TDG", "CREB3L2", "CREB3", "OVOL1", "FOSB", "GATAD2B", "PPP1R13L", "RORC", "ARNT", "AEBP1", "ARHGAP35", "PURA", "PURB", "KAT6A", "BARX2", "MXD1", "FIGLA", "BOLA3", "ATF6B", "WNT5A", "CREBRF", "FOXI1", "ZNF750", "SOX9", "BUD31", "NR4A1", "DDIT3", "ZNF740", "CSRNP2", "POU6F1", "MBD1", "KDM5C", "HCFC1", "SMARCA2", "SMARCA4", "CARM1", "NFIX", "ZBTB6", "SP1", "ZBTB26", "ZNF462", "MLLT3", "NFIB", "GTF2H2", "JMY", "CDK7", "SREBF1", "MBNL3", "FOXO4", "FOXO3", "TSC22D3", "TERF2", "HIF3A", "LITAF", "XBP1", "TFAP4", "GLIS2", "CSRNP1", "CREM", "DNAJC1", "ZBTB33", "NKRF", "SKI", "RARA", "NFE2L1", "SP6", "ZNF652", "DLX4", "BHLHA15", "CREG1", "TAF13", "BLZF1", "MKL2", "PITX2", "SMARCA5", "CLOCK", "NFXL1", "TAF2", "MED30", "TRPS1", "GRHL2", "MTDH", "NCOA2", "PLAG1", "KCTD1", "NFATC1", "HIVEP1", "TFAP2A", "RREB1", "IRX2", "IRX1", "HDAC1", "HDAC3", "ATF6", "CDC73", "DR1", "BCL10", "ZZZ3", "NFIA", "SPDEF", "PPARD", "ELK4", "KDM5B", "ZNF341", "NCOA3", "ZNF217", "UHRF1", "KLF11", "KLF16", "GATAD2A", "JUN", "ARID3A", "ZNF511", "ERF", "FOS", "CIC", "ZNF410", "ZFP36L1", "CREB3L1", "NR1H3", "CCND1", "EHF", "ELF5", "SSRP1", "KDM2A", "TFCP2L1", "NR4A2", "ATF2", "IKZF2", "MEIS1", "GRHL1", "ZNF451", "TBX3", "GRHL3", "SUPT4H1")
regulons <- c("REL", "ETV6", "PSMD12", "KLF16", "HMBOX1", "BDP1", "TCF3", "STAT1", "STAT2", "IRF2", "NPDC1", "TBX3", "MGA", "KLF13", "MAX", "HLF", "PBX3", "RELA", "ZNF143", "GTF3C2", "RBPJ", "ENO1", "ZBTB7B", "MEIS1", "TRIM69", "NR2C2", "PPARD", "SPDEF", "NR2C1", "XBP1", "MEF2D", "SIN3A", "ZXDC", "SP1", "NFKB1", "YY1", "ATF4", "REST", "HDAC6", "ARNTL", "SP2", "KDM5B", "CCNT2", "FOXP1", "XRCC4", "NR1D2", "TCF4", "KLF3", "ZBED1", "HDAC2", "RAD21", "RUNX1", "BRCA1", "TAF1", "PHF8", "KDM5A", "TAF7", "ZBTB2", "BACH1", "TFEB", "NFATC3", "ZFX", "SP3", "NR3C1", "RCOR1", "NFAT5", "GABPB1", "FOXP4", "SRF", "RFX1", "CEBPZ", "SREBF2", "ELF5", "ERF", "ETS2", "EHF", "ELF3", "EP300", "ELF2", "ATF6B", "CLOCK", "ATF2", "CREB3L2", "ATF6", "CREB3", "ELK1", "CTCF", "BCLAF1", "BRF1", "ELF1", "HCFC1", "SMARCA4", "SREBF1", "GTF2F1", "RBBP5", "UBTF")

VST_count_mat <- read.table('../../Figure_and_Scripts/z_data_files/RNAseq_VSTNormalizeCounts_noOutliers_05.23.23.txt', header=T)

# Filter for candidate TFs and Venom genes, also combine LVG and RVG (take the average). Note that LVG_13 failed and is missing.
VST_count_mat <- VST_count_mat %>% 
  dplyr::select(matches('LVG|RVG')) # %>% 
  #filter(row.names(.) %in% TF_list | grepl('Venom', row.names(.)))

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

sample_info <- read.table('../../Figure_and_Scripts/z_data_files/rnaSampleInfo_06.29.21.txt', header = T)
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

# Make boxplots of Venom TF and nonvenom TF expression
Venom_TFs_mat <- VST_count_mat_filtered[row.names(VST_count_mat_filtered) %in% TF_list,]
d1 <- as.data.frame(rowVars(as.matrix(Venom_TFs_mat))) # across samples
d2 <- as.data.frame(rowVars(as.matrix(Venom_TFs_mat[, !(colnames(Venom_TFs_mat) %in% c('VG_6','VG_7', 'VG_8'))]))) # within C. viridis
colnames(d1) <- "Variances"
d1$Group <- "VenomTF_Across"
colnames(d2) <- "Variances"
d2$Group <- "VenomTF_Within"

human_tfs <- readLines('/Volumes/SeagatePortableDrive/sc_multiome/Crotalus-Dendroaspis_joint_multiome_analyses/STRINGdb_Networks_Cytoscape/Human_TFs.txt')
Non_Venom_TFs_mat <- na.omit(VST_count_mat_filtered[human_tfs[!human_tfs %in% TF_list],])
d3 <- as.data.frame(rowVars(as.matrix(Non_Venom_TFs_mat))) # across samples
d4 <- as.data.frame(rowVars(as.matrix(Non_Venom_TFs_mat[, !(colnames(Non_Venom_TFs_mat) %in% c('VG_6','VG_7', 'VG_8'))]))) # within C. viridis
colnames(d3) <- "Variances"
d3$Group <- "NVTF_Across"
colnames(d4) <- "Variances"
d4$Group <- "NVTF_Within"

variances_df <- rbind(d1, d2, d3, d4)
variances_df %>% 
  rownames_to_column("gene") %>% 
  ggplot(aes(x = Group, y = Variances, fill = Group)) +
  geom_boxplot() +
  # stat_compare_means(method = "t.test", label = "p.format") +  # Display the p-value on the plot
  theme_bw() +
  scale_fill_manual(values = c("grey60", "grey60", "forestgreen", "forestgreen")) +
  scale_x_discrete(labels = c('NV TFs (across)', 'NV TFs (within)', 'Venom TFs (across)', 'Venom TFs (within)')) +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# T tests
var.test(Variances ~ Group, data = variances_df %>% filter(Group %in% c("VenomTF_Across", "NVTF_Across")))
var.test(Variances ~ Group, data = variances_df %>% filter(Group %in% c("VenomTF_Within", "NVTF_Within")))
t.test(Variances ~ Group, data = variances_df %>% filter(Group %in% c("VenomTF_Across", "NVTF_Across")), var.equal = F) # < 0.01
t.test(Variances ~ Group, data = variances_df %>% filter(Group %in% c("VenomTF_Within", "NVTF_Within")), var.equal = F) # < 0.01
#

# Make TF Heatmap dataframe
TF_VST_mat <- VST_count_mat_filtered %>% 
  filter(!str_detect(row.names(.), 'Venom')) %>%
  mutate(Expr_Variance = rowVars(as.matrix(.))) %>%
  arrange(desc(Expr_Variance)) %>% 
  slice_head(n = 25) %>% # keep the top 25 TFs sorted by variance
  rownames_to_column(var = "TF") %>% 
  # select(-Expr_Variance) %>% 
  pivot_longer(cols = VG_1:VG_13, 
               names_to = "VG_Sample",
               values_to = "Average_Exp")

p_heatmap <- TF_VST_mat %>% 
  ggplot(., aes(x=TF, y=VG_Sample, fill=Average_Exp)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma", breaks = c(6,9,12)) +
  scale_y_discrete(labels=species_labs,
                   limits=c(sample_info$Sample_ID)) +
  scale_x_discrete(labels=unique(TF_VST_mat$TF),
                   limits=unique(TF_VST_mat$TF)) +
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
        legend.margin = margin(-10,5,0,5),
        plot.margin = margin(0,0,10,0)) +
  guides(fill=guide_colorbar(barwidth = 1, title.position = "right")) +
  labs(fill = "VST norm. exp.") +
  ggtitle('TF expression heatmap, sorted by variance')

gene_lookup <- read.csv('../../Figure_and_Scripts/z_data_files/gene_annotated_lookup.csv', header=T)
gene_lookup <- gene_lookup %>% 
  filter(Gene_ID %in% unique(TF_VST_mat$TF)) %>% 
  arrange(match(Gene_ID, unique(TF_VST_mat$TF))) %>% 
  dplyr::select(-c(TF,Blair_candidate, venom_family,ASP_related, SP1_complex)) %>% 
  pivot_longer(cols = -Gene_ID, 
               names_to = "Characteristic",
               values_to = "Yes_No") %>% 
  mutate(Yes_No = ifelse(Yes_No == 0, NA, 1))

p_balloon <- gene_lookup %>% 
  ggballoonplot(x = 'Gene_ID', y = 'Characteristic', size = 'Yes_No', fill = 'Characteristic', shape = 23, size.range = 4) +
  guides(size = 'none', fill = 'none') +
  scale_x_discrete(labels=unique(TF_VST_mat$TF),
                   limits=unique(TF_VST_mat$TF)) +
  scale_y_discrete(limits = c('ERK_related', 'UPR_related', 'pioneer', 'Westfall_regulon'), 
                   labels=c('ERK related', 'UPR related', 'Pioneer', 'Venom regulon')) +
  theme(plot.margin = margin(-10,15,0,5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values=brewer.pal(5, 'Dark2'))

p1 <- plot_grid(p_heatmap, p_balloon, ncol = 1, align = 'v', rel_heights = c(4.5,1), axis = 'lr')

#### 2: WGCNA Module-Gene Significance ####
TF_df <- read.csv('../../Figure_and_Scripts/z_data_files/gene_annotated_lookup.csv', header = T)
TF_df <- TF_df %>% 
  filter(TF == 1)

WGCNA_info_viridis <- read.csv('../../Figure_and_Scripts/1_RNAseq_and_TF_expression_variation/3_WGCNA/geneInfo_Viridis.csv', header = T)
colnames(WGCNA_info_viridis)[1] <- 'Gene'
WGCNA_info_lutosus <- read.csv('../../Figure_and_Scripts/1_RNAseq_and_TF_expression_variation/3_WGCNA/geneInfo_Lutosus.csv', header = T)
colnames(WGCNA_info_lutosus)[1] <- 'Gene'
WGCNA_info_concolor <- read.csv('../../Figure_and_Scripts/1_RNAseq_and_TF_expression_variation/3_WGCNA/geneInfo_Concolor.csv', header = T)
colnames(WGCNA_info_concolor)[1] <- 'Gene'
WGCNA_info_cerberus <- read.csv('../../Figure_and_Scripts/1_RNAseq_and_TF_expression_variation/3_WGCNA/geneInfo_Cerberus.csv', header = T)
colnames(WGCNA_info_cerberus)[1] <- 'Gene'

WGCNA_info_joined <- left_join(WGCNA_info_viridis, WGCNA_info_lutosus, by = 'Gene') %>% 
  left_join(WGCNA_info_concolor, by = 'Gene') %>% 
  left_join(WGCNA_info_cerberus, by = 'Gene') %>% 
  dplyr::select(-matches('MM|p.')) %>% 
  dplyr::select(-matches('\\.y|\\.x\\.x')) %>% 
  dplyr::rename(moduleColor = moduleColor.x) %>% 
  filter(Gene %in% TF_df$Gene_ID) %>% 
  # filter(Gene %in% TF_df$Gene_ID[TF_df$Blair_candidate == 1]) %>% 
  filter(moduleColor == 'pink' | moduleColor == 'lightyellow' | moduleColor == 'magenta' | moduleColor == 'royalblue') %>% 
  mutate(moduleColor = case_when(moduleColor == 'pink' ~ 'viridis',
                                 moduleColor == 'lightyellow' ~ 'lutosus',
                                 moduleColor == 'magenta' ~ 'concolor', 
                                 moduleColor == 'royalblue' ~ 'cerberus')) %>% 
  mutate(GS_cumulative = case_when(moduleColor == 'lutosus' ~ GS.Lutosus,
                                   moduleColor == 'viridis' ~ GS.Viridis,
                                   moduleColor == 'concolor' ~ GS.Concolor,
                                   moduleColor == 'cerberus' ~ GS.Cerberus)) %>% 
  dplyr::select(Gene, moduleColor, GS_cumulative)

WGCNA_info_joined <- WGCNA_info_joined %>%
  group_by(moduleColor) %>%
  arrange(moduleColor, desc(GS_cumulative), .by_group = TRUE) %>% 
  slice_head(n=10) %>% 
  filter(GS_cumulative > 0)

p_modulegenesignificance <- WGCNA_info_joined %>% 
  ggplot(aes(x = Gene, y = GS_cumulative, fill = moduleColor)) +
  geom_col() +
  theme_bw() +
  labs(x = "Transcription factors", y = "Module gene\nsignificance") +
  scale_x_discrete(labels = WGCNA_info_joined$Gene,
                   limits = WGCNA_info_joined$Gene) +
  scale_fill_manual(values = c('cerberus' = "black", 
                               'viridis' = "#2E8B58",
                               'concolor' = "#D9A528",
                               'lutosus' = "#8B481F")) +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=9),
        legend.title=element_blank(),
        legend.position = 'none',
        plot.margin = margin(5,5,0,5))

# Balloon plot
gene_lookup <- read.csv('../../Figure_and_Scripts/z_data_files/gene_annotated_lookup.csv', header=T)
gene_lookup <- gene_lookup %>% 
  filter(Gene_ID %in% WGCNA_info_joined$Gene) %>% 
  arrange(match(Gene_ID, WGCNA_info_joined$Gene)) %>% 
  dplyr::select(Gene_ID, Blair_candidate) %>% # none of these genes are UPR related, ERK related, SP1_complex or pioneer, or regulon
  pivot_longer(cols = -Gene_ID, 
               names_to = "Characteristic",
               values_to = "Yes_No") %>% 
  mutate(Yes_No = ifelse(Yes_No == 0, NA, 1))

p_balloon2 <- gene_lookup %>% 
  ggballoonplot(x = 'Gene_ID', y = 'Characteristic', size = 'Yes_No', fill = 'Characteristic', shape = 23, size.range = 3.5) +
  guides(size = 'none', fill = 'none') +
  scale_x_discrete(labels=WGCNA_info_joined$Gene,
                   limits=WGCNA_info_joined$Gene) +
  scale_y_discrete(labels='Perry (2022)\ncandidate') +
  theme(plot.margin = margin(5,15,2,5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values='#2E8B58')

p2 <- plot_grid(p_modulegenesignificance, p_balloon2, ncol = 1, rel_heights = c(7,1), axis = 'lr', align = 'v')

#### 3: Toxin - TF Correlation Heatmap ####
# Sample information
sample_info <- read.table('../../Figure_and_Scripts/z_data_files/rnaSampleInfo_06.29.21.txt', header = T)
sample_info <- sample_info %>% 
  filter(str_detect(Sample_ID, 'RVG')) %>% 
  mutate(Sample_ID = gsub('R', '', Sample_ID)) %>% 
  mutate(Population = fct_relevel(Population, 
                                  "South", "Mid", "North", "Other")) %>% 
  arrange(Lineage, Population, Sample_ID)

# RNAseq VST normalized count matrix
VST_count_mat <- read.table('../../Figure_and_Scripts/z_data_files/RNAseq_VSTNormalizeCounts_noOutliers_05.23.23.txt', header=T)
# Filter for TFs and Venom genes, also combine LVG and RVG (take the average). Note that LVG_13 failed and is missing.
VST_count_mat <- VST_count_mat %>% 
  dplyr::select(matches('LVG|RVG')) %>% 
  filter(row.names(.) %in% TF_list | grepl('Venom', row.names(.)))

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

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

sample_info$Sample_ID[(sample_info$Lineage=='viridis')] # viridis
sample_info$Sample_ID[!(sample_info$Lineage=='viridis')] # other snakes


VST_count_mat_filtered_t <- t(VST_count_mat_filtered)[sample_info$Sample_ID,] # All snakes samples
# res <- cor(VST_count_mat_filtered_t, method = "pearson")
res2 <- rcorr(as.matrix(VST_count_mat_filtered_t), type = 'pearson')

corrMatres <- flattenCorrMatrix(res2$r, res2$P)
corrMatres$fdr_p <- p.adjust(corrMatres$p, method = "fdr")

corrMatres_filtered <- corrMatres %>%
  filter(p < 0.05) %>%
  filter(fdr_p < 0.1) %>%
  mutate(
    row2 = ifelse(grepl("Venom", column), column, row),
    column2 = ifelse(grepl("Venom", column), row, column)) %>% 
  dplyr::select(-c(row, column)) %>% 
  filter((grepl("Venom", row2) & !grepl("Venom", column2)) | (!grepl("Venom", row2) & grepl("Venom", column2))) %>% 
  mutate(row2 = gsub("Venom_", '', row2))

toxin_order <- c('BPP', 'myotoxin', 'ohanin', paste0('ADAM28_',c(1:2)), paste0('CRISP',c(1:4)), paste0('CTL',c(1:6)), paste0('EXO',c(1:3)), paste0('LAAO',c(1:3)),  paste0('VEGF',c(1:2)), paste0('vQC',c(1:2)),
                 'PLA2A1', 'PLA2B1', 'PLA2C1', 'PLA2K', paste0('SVMP',c(1:10)), paste0('SVSP',1:11))

toxin_order[sort(match(unique(corrMatres_filtered$row2), toxin_order))]

gene_lookup <- read.csv('../../Figure_and_Scripts/z_data_files/gene_annotated_lookup.csv', header=T)
gene_lookup <- gene_lookup %>% 
  filter(Gene_ID %in% unique(corrMatres_filtered$column2)) %>% 
  arrange(match(Gene_ID, unique(corrMatres_filtered$column2))) %>% 
  dplyr::select(-c(TF,venom_family,Blair_candidate)) %>% 
  pivot_longer(cols = UPR_related:pioneer, 
               names_to = "Characteristic",
               values_to = "Yes_No") %>% 
  mutate(Yes_No = ifelse(Yes_No == 0, NA, 1))

remove_genes <- gene_lookup %>% 
  group_by(Gene_ID) %>%
  filter(all(is.na(Yes_No))) %>%
  pull(Gene_ID) %>% 
  unique()


p_corrMat <- corrMatres_filtered %>% 
  filter(!column2 %in% remove_genes) %>% # remove genes with no functional annotations from Blair
  ggplot() +
  theme_bw() +
  geom_tile(aes(x=row2, y=column2, fill=cor)) +
  theme(axis.text.x = element_text(angle=70, hjust=1, size=9),
        legend.text = element_text(size=10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size=10, angle = 90),
        legend.title.align = 0.5,
        legend.direction = "vertical",
        plot.margin = margin(5,0,5,0)) +
  scale_x_discrete(labels=toxin_order[sort(match(unique(corrMatres_filtered$row2), toxin_order))],
                   limits=toxin_order[sort(match(unique(corrMatres_filtered$row2), toxin_order))]) +
  scale_y_discrete(labels=c("CREB3", "FOS", "DDIT3", "ATF2", "ATF4", "BHLHA15", "CREB3L1", "CREG1", "CREM", "ERF", "KLF11", "NR4A1", "XBP1"),
                   limits=c("CREB3", "FOS", "DDIT3", "ATF2", "ATF4", "BHLHA15", "CREB3L1", "CREG1", "CREM", "ERF", "KLF11", "NR4A1", "XBP1")) +
  scale_fill_gradientn(colors = c("#084594", "#FFFFFF", "#990000"),
                       values = scales::rescale(c(-1, 0, 1), to = c(0, 1)),
                       space = "Lab", na.value = "transparent") +
  guides(fill = guide_colorbar(title = "Pearson's rho", barwidth = 1, title.position = "right")) +
  ggtitle('Toxin-TF expression correlation heatmap (all samples)')

p_balloon3 <- gene_lookup %>% 
  ggballoonplot(x = 'Characteristic', y = 'Gene_ID', size = 'Yes_No', fill = 'Characteristic', shape = 23, size.range = 4) +
  guides(size = 'none', fill = 'none') +
  scale_y_discrete(labels=c("CREB3", "FOS", "DDIT3", "ATF2", "ATF4", "BHLHA15", "CREB3L1", "CREG1", "CREM", "ERF", "KLF11", "NR4A1", "XBP1"),
                   limits=c("CREB3", "FOS", "DDIT3", "ATF2", "ATF4", "BHLHA15", "CREB3L1", "CREG1", "CREM", "ERF", "KLF11", "NR4A1", "XBP1")) +
  scale_x_discrete(labels=c('ASP related', 'ERK related', 'Pioneer', 'SP1 complex', 'UPR related')) +
  theme(plot.margin = margin(0,0,0,25),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_manual(values=brewer.pal(5, 'Dark2'))

p3 <- plot_grid(p_balloon3, p_corrMat, ncol = 2, align = 'h', rel_widths = c(0.2,1))

#### Plot Together ####
plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(1,0.7,0.8))
