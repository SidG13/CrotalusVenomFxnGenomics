# Written by BWP and continued by SSG #

# This code is used to take raw output following from a STAR alignment and 
# FeatureCounts counting to produce normalized RNA count matrices. Some basic 
# differential expression is also performed. 

library(DESeq2)
library(pheatmap)
library(viridis)
library(IHW)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(tidyverse)


# Reading in and preparing files ------------------------------------------
setwd('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/z_data_files')

name_conversion_table <- read.table('Cvv_GTF_to_converted_names_05.22.23.txt',header = T,stringsAsFactors = F)

sample_info <- read_tsv('rnaSampleInfo_06.29.21.txt')

rawCounts <- read.table('cvv_12SnakeVenomRNAseq_rawCounts_09.16.22.txt',stringsAsFactors = F,skip=1,header=T)
rawCounts <- rawCounts[,-c(2:7)]

# Generate simpler raw count table
rawCounts$Geneid <- name_conversion_table$converted_id_no_dups[match(rawCounts$Geneid, name_conversion_table$gtf_gene)] # Convert GTF names to simpler names

row.names(rawCounts) <- rawCounts$Geneid
rawCounts <- rawCounts[,c(-1)]
# rawCounts.simple <- rawCounts[rowSums( rawCounts != 0 ) > 30,]
colnames(rawCounts) <- colnames(rawCounts) %>% str_remove_all('..STAR_mapped.|Aligned.sortedByCoord.out.bam')
rawCounts <- rawCounts %>%  select(-contains('ODPE'),-contains('TDPE'),-contains('Un'),-contains('Stom'),-contains('Panc'),-contains('Skin'))

names(rawCounts)
names(rawCounts) %in% sample_info$Sample_ID 
head(rawCounts)


# DEseq2 Analyses ---------------------------------------------------------

##########################################
######## Tissue Type + Population ########
##########################################

### Note: want to test how expression varies across tissues across populations. By pasting tissue and population factors, we get all combinations:
paste0(sample_info$Ven_NonVen, sample_info$Population)

colData_TissueAndPop <- (DataFrame(batch= factor(sample_info$Batch), 
                                   condition= factor(paste0(sample_info$Ven_NonVen, sample_info$Population))))

dds.pop <- DESeqDataSetFromMatrix(countData = rawCounts,
                                  colData = colData_TissueAndPop,
                                  design = formula(~batch + condition))


keep <- rowSums(counts(dds.pop)) >= 30 # basic prefiltering, keep genes with >30 counts
dds.pop <- dds.pop[keep,]

dds.pop <- DESeq(dds.pop)

dds.pop <- dds.pop[which(mcols(dds.pop)$betaConv),]

resultsNames(dds.pop)

normcounts <- counts(dds.pop,normalized=TRUE)

### Vst PCA
vst.pop <- vst(dds.pop, blind=FALSE)

pcaData.pop <- plotPCA(vst.pop, intgroup=c("condition","batch"), returnData=TRUE)
percentVar.pop <- round(100 * attr(pcaData.pop, "percentVar"))

ggplot(pcaData.pop, aes(PC1, PC2, color=condition, label=name, shape=batch)) +
  geom_point(size=3) +
  geom_text_repel() +
  xlab(paste0("PC1: ",percentVar.pop[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.pop[2],"% variance")) +
  coord_fixed() + theme_minimal()

### Visually inspect PCA for outliers that do not cluster with similar tissues and remove
### Here, Skn_11, Skn_12, Lvr_13, and Pnc_13, Pnc_2 appear to be outliers
### Additionally, Pancreas appears most similar to venom gland samples while skin is different from all other tissues...
### Going to therefore use pancreas for comparisons to venom tissues

rawCounts.simple.noOutliers <- rawCounts %>% select(-Skn_11,-Skn_12,-Lvr_13,-Pnc_13,-Pnc_2)

sample_info.noOutliers <- sample_info %>% 
  filter(str_detect(Sample_ID,'Skn_11|Skn_12|Lvr_13|Pnc_13|Pnc_2',negate = T)) %>% 
  mutate(TissueType2 = ifelse(Tissue == 'Pnc','Pancreas',ifelse(Ven_NonVen == 'Ven','Venom','Ignore'))) %>% 
  mutate(Lineage2 = ifelse(Lineage == 'viridis', 'viridis', 'nonviridis'))

sum(names(rawCounts.simple.noOutliers) != sample_info.noOutliers$Sample_ID) #check to make sure columns in count table still match sample_info rows

# Re-run DEseq2 looking at venom population contrasts first
tfs <- read.csv('gene_annotated_lookup.csv', header = T)
tfs <- tfs %>% filter(TF == 1) %>% select(Gene_ID) %>% pull() # get list of all TFs
# Blair's 72 candidates TF motifs from 2022 GR paper
blair_72_tfs_jaspar <- c("Arid3a", "Arnt", "ATF2", "ATF4", "ATF6", "BARX2", "Bhlha15", "CLOCK", "CREB3", "CREB3L1", "Creb3l2", "CREM", "CTCF", "Dlx4", "EHF", "ELF5", "ELK4", "ERF", "FIGLA", "FOS", "FOXC2", "FOXI1", "FOXO3", "FOXO4", "GLIS2", "GRHL1", "GRHL2", "HES6", "JUN", "KLF11", "KLF13", "KLF16", "MEIS1", "NFATC1", "NFE2L1", "NFIA", "NFIB", "NFIX", "NR4A1", "NR4A2", "OVOL1", "PITX2", "PLAG1", "POU6F1", "PPARD", "RARA", "RORA", "RORC", "RREB1", "SOX10", "SOX9", "SP1", "SPDEF", "SREBF1", "SREBF2", "TBX3", "TFAP2A", "TFAP4", "TP73", "XBP1", "ZBTB26", "ZBTB33", "ZBTB6", "ZNF341", "ZNF410", "ZNF652", "ZNF740", "FOSB", "Tcfcp2l1", "Nr1h3::Rxra", "Ddit3::Cebpa", "Irx2")

colData_TissueAndPop.noOut <- (DataFrame(batch= factor(sample_info.noOutliers$Batch), 
                                         condition= factor(paste0(sample_info.noOutliers$TissueType2, sample_info.noOutliers$Population))))

dds.pop <- DESeqDataSetFromMatrix(countData = rawCounts.simple.noOutliers, 
                                  colData = colData_TissueAndPop.noOut,
                                  design = formula(~batch + condition))

keep <- rowSums(counts(dds.pop)) >= 30 # basic prefiltering, keep genes with >30 counts
dds.pop <- dds.pop[keep,]
dds.pop <- DESeq(dds.pop)
dds.pop <- dds.pop[which(mcols(dds.pop)$betaConv),]
resultsNames(dds.pop)

normcounts.noOut <- as.data.frame(counts(dds.pop,normalized=TRUE))
#write.table(normcounts.noOut,'./RNAseq_NormalizedCounts_noOutliers_06.29.23.txt', quote =F, col.names = T, row.names = T, sep = '\t')

vsd.pop.noOut <- as.data.frame(assay(vst(dds.pop, blind=FALSE)))
#write.table(vsd.pop.noOut,'./RNAseq_VSTNormalizeCounts_noOutliers_05.23.23.txt', row.names = T, col.names = T, quote = F, sep = '\t')
vsd.pop.noOut[grepl('Venom', rownames(vsd.pop.noOut)), grepl('VG', colnames(vsd.pop.noOut))]
pheatmap(vsd.pop.noOut[grepl('Venom', rownames(vsd.pop.noOut)), grepl('LVG|RVG', colnames(vsd.pop.noOut))], cluster_cols = F)

### Venom North vs. Venom Mid
deRes.NvsM <- as.data.frame(results(dds.pop, contrast=c('condition','VenomNorth','VenomMid')))
ihwRes.NvsM <- ihw(pvalue ~ baseMean,  data = deRes.NvsM, alpha = 0.05)
rejections(ihwRes.NvsM)
deRes.NvsM$IHW_pvalue <- ihwRes.NvsM@df$adj_pvalue
deRes.NvsM <- deRes.NvsM[order(deRes.NvsM$IHW_pvalue),]
deRes.NvsM <- deRes.NvsM %>% mutate(signif = ifelse(IHW_pvalue < 0.05,'True','False'))

new_tfs1 <- deRes.NvsM %>% 
  filter(signif == 'True') %>% 
  filter(row.names(.) %in% tfs) %>% 
  filter(!row.names(.) %in% blair_72_tfs_jaspar) %>% 
  rownames_to_column(var = "RowNames") %>%
  pull(RowNames)

### Venom North vs. Venom South
deRes.NvsS <- as.data.frame(results(dds.pop, contrast=c('condition','VenomNorth','VenomSouth')))
ihwRes.NvsS <- ihw(pvalue ~ baseMean,  data = deRes.NvsS, alpha = 0.05)
rejections(ihwRes.NvsS)
deRes.NvsS$IHW_pvalue <- ihwRes.NvsS@df$adj_pvalue
deRes.NvsS <- deRes.NvsS[order(deRes.NvsS$IHW_pvalue),]
deRes.NvsS <- deRes.NvsS %>% mutate(signif = ifelse(IHW_pvalue < 0.05,'True','False'))

new_tfs2 <- deRes.NvsS %>% 
  filter(signif == 'True') %>% 
  filter(row.names(.) %in% tfs) %>% 
  filter(!row.names(.) %in% blair_72_tfs_jaspar) %>% 
  rownames_to_column(var = "RowNames") %>%
  pull(RowNames)


### Venom Mid vs. Venom South
deRes.MvsS <- as.data.frame(results(dds.pop, contrast=c('condition','VenomMid','VenomSouth')))
ihwRes.MvsS <- ihw(pvalue ~ baseMean,  data = deRes.MvsS, alpha = 0.05)
rejections(ihwRes.MvsS)
deRes.MvsS$IHW_pvalue <- ihwRes.MvsS@df$adj_pvalue
deRes.MvsS <- deRes.MvsS[order(deRes.MvsS$IHW_pvalue),]
deRes.MvsS <- deRes.MvsS %>% mutate(signif = ifelse(IHW_pvalue < 0.05,'True','False'))

new_tfs3 <- deRes.MvsS %>% 
  filter(signif == 'True') %>% 
  filter(row.names(.) %in% tfs) %>% 
  filter(!row.names(.) %in% blair_72_tfs_jaspar) %>% 
  rownames_to_column(var = "RowNames") %>%
  pull(RowNames)

### Rerun DESeq2 to analyze viridis vs non-viridis TF differences
colData_TissueAndPop.noOut <- (DataFrame(batch= factor(sample_info.noOutliers$Batch), 
                                         condition= factor(paste0(sample_info.noOutliers$Lineage2, sample_info.noOutliers$TissueType2))))

dds.pop <- DESeqDataSetFromMatrix(countData = rawCounts.simple.noOutliers, 
                                  colData = colData_TissueAndPop.noOut,
                                  design = formula(~batch + condition))

keep <- rowSums(counts(dds.pop)) >= 30 # basic prefiltering, keep genes with >30 counts
dds.pop <- dds.pop[keep,]
dds.pop <- DESeq(dds.pop)
dds.pop <- dds.pop[which(mcols(dds.pop)$betaConv),]
resultsNames(dds.pop)

### Viridis venom vs. Nonviridis venom
deRes.VvsNV <- as.data.frame(results(dds.pop, contrast=c('condition','viridisVenom','nonviridisVenom')))
ihwRes.VvsNV <- ihw(pvalue ~ baseMean,  data = deRes.VvsNV, alpha = 0.05)
rejections(ihwRes.VvsNV)
deRes.VvsNV$IHW_pvalue <- ihwRes.VvsNV@df$adj_pvalue
deRes.VvsNV <- deRes.VvsNV[order(deRes.VvsNV$IHW_pvalue),]
deRes.VvsNV <- deRes.VvsNV %>% mutate(signif = ifelse(IHW_pvalue < 0.05,'True','False'))

new_tfs4 <- deRes.VvsNV %>% 
  filter(signif == 'True') %>% 
  filter(row.names(.) %in% tfs) %>% 
  filter(!row.names(.) %in% blair_72_tfs_jaspar) %>% 
  rownames_to_column(var = "RowNames") %>%
  pull(RowNames)


# all new TFs to add to JASPAR matrix for TFBS scanning
union(c(new_tfs1, new_tfs2), c(new_tfs3, new_tfs4))