### Remake the grant figure of SVMP array with ATAC-peaks ###

library(BSgenome.Cviridis.custom.CroVir.noSeqNamesMods)
library(tidyverse)
library(rtracklayer)
library(ggcoverage)
library(scales)
library(gggenes)
library(cowplot)
library(viridis)
library(chromVAR)
library(Biostrings)
library(patchwork)
library(scales)
library(ggforce)
library(pheatmap)
library(seqinr)

#### SVMP Array plot (from Blair 2022 Genome Res.) ####
setwd('~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_VenomGeneRegulation_NEW_Aug2021/')

pri_venom_genes <- read_tsv('./data/venom_annotation/PriorityVenomGenes_08.02.21.txt',col_names = F)
exp <- read_tsv('./analysis/1_gene_expression/norm_counts/AllGenes_AvgMedianTPM_1DPE_08.02.21.tsv') %>% 
  filter(txid %in% pri_venom_genes$X6) %>%
  left_join(pri_venom_genes,by=c('txid'='X6')) %>% 
  mutate(gene = ifelse(str_detect(X7,'ADAM28',negate = T),str_replace_all(X7,'_',' '),X7))


all_info <- read_tsv('./data/annotation/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019.GeneEntriesOnly.gff',col_names = F) %>% 
  filter(str_detect(X9,'trnascan',negate = T)) %>% 
  mutate(tx_id = str_split_fixed(X9,';',4)[,3]) %>% 
  mutate(tx_id = str_remove_all(tx_id,'Crovir_Transcript_ID=')) %>% 
  filter(tx_id %in% pri_venom_genes$X6) %>% 
  left_join(pri_venom_genes,by=c('tx_id' = 'X6')) %>% 
  select(molecule = 1, gene = 16, start = 4, end = 5, strand = 7,tx_id) %>% 
  mutate(strand = ifelse(strand == '+','forward','reverse')) %>% 
  mutate(direction = ifelse(strand == 'forward',1,-1)) %>% 
  mutate(gene = ifelse(str_detect(gene,'ADAM28',negate = T),str_replace_all(gene,'_',' '),gene)) %>% 
  left_join(exp) %>% 
  mutate(gene = ifelse(str_detect(gene,'ADAM28'),paste('NVP: ',gene,sep = ''),gene)) %>% 
  mutate(prom_start = ifelse(strand=='forward',start,end))




SVMP.info <- all_info %>% 
  filter(str_detect(gene,'SVMP|ADAM'))



###
### ******* This is where you set the x-axis limits. I did 20k bases up/downstream of regions. 
###
#SVMP.reg.start <- min(SVMP.info$start)-20000
#SVMP.reg.start <- 13500000
SVMP.reg.start = floor(13901005 / 1000 ) * 1000
#SVMP.reg.start <- 8500000 # SVSP

#SVMP.reg.end <- max(SVMP.info$end)+20000
#SVMP.reg.end <- 15000000
SVMP.reg.end = ceiling(14424729 / 1000 ) * 1000 
#SVMP.reg.end <- 9000000 # SVSP




SVMP.reg.length <- paste(c(round((SVMP.reg.end-SVMP.reg.start)/1000,digits = 2),'kb'),collapse = ' ')

# Read in vPERs and super-enhancers

svmp.vpers <- 
  read_tsv('./analysis/6_ABC_Enhancers/ABC_output/_reformat/SVMP_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed',col_names = F) %>% 
  #read_tsv('./analysis/6_ABC_Enhancers/ABC_output/_reformat/SVSP_EnhancerPredictionsFull_VenomGenes_simple_newID_08.18.21.bed', col_names = F) %>% 
  select(molecule=1,start=2,end=3,id=4) %>%  
  mutate(gene = str_split_fixed(id, '_',2)[,2]) %>% 
  mutate(gene = str_split(gene,'\\.')) %>% 
  unnest(gene) %>% 
  mutate(gene = str_replace(gene,'\\_',' ')) %>% 
  left_join(SVMP.info,by='gene') %>% 
  mutate(type=' vPERs') %>% 
  select(molecule=1,start=2,end=3,id,gene,gene.start=7,gene.end=8,type,prom_start)

#### Get venom gene expression ####
VST_count_mat <- read.table('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/z_data_files/RNAseq_NormalizedCounts_noOutliers_06.29.23.txt', header=T)

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
  select(-contains('LVG'),-contains('RVG')) %>% 
  rownames_to_column(var = 'Gene') %>% 
  filter(grepl('SVMP|ADAM', Gene)) %>% 
  mutate(AverageExp = rowMeans(select(., starts_with("VG_")), na.rm = TRUE)) %>% 
  select(Gene, AverageExp)

VST_count_mat_filtered <- VST_count_mat_filtered %>% 
  mutate(Gene = gsub('Venom_SVMP', 'SVMP ', Gene)) %>% 
  mutate(Gene = gsub('Venom_', 'NVP: ', Gene))

# ADD TO SVMP INFO
SVMP.info <- SVMP.info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
SVMP.info <- SVMP.info %>% mutate(gene = gsub(' ','',gene))
#all_info <- all_info %>% left_join(VST_count_mat_filtered, by = c('gene'='Gene'))
#all_info <- all_info %>% mutate(gene = gsub(' ','',gene))

#### Plotting ####
p1 <- ggplot(SVMP.info,aes(xmin = start, xmax = end, y = str_remove(molecule,'scaffold\\-'), forward = direction, fill = log10(AverageExp+1))) +
  ggrepel::geom_text_repel(data = all_info %>% 
                             filter(str_detect(gene,'SVMP')) %>% # Change for SVSP-SVMP
                             mutate(start = (start + end)/2), 
                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = 0.5, size=3) +
  #  ggrepel::geom_text_repel(data = all_info %>% 
  #                             filter(str_detect(gene,'ADAM')) %>% 
  #                             mutate(start = (start + end)/2), 
  #                           aes(x = start, y = str_remove(molecule,'scaffold\\-'), label = gene), inherit.aes = F, nudge_y = -0.5,size=3,color='grey60') +
  geom_segment(aes(x=SVMP.reg.start,xend=SVMP.reg.end,y=str_remove(molecule,'scaffold\\-'),yend=str_remove(molecule,'scaffold\\-')),lwd=1,color='grey70') +
  geom_diagonal(inherit.aes = F,data=svmp.vpers,aes(x=prom_start,xend=start,y='mi1',yend=type,alpha = stat(index)),strength = -0.2,show.legend = F) + # Change for SVSP-SVMP
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(2, "mm"),show.legend = F) +
  geom_segment(inherit.aes = F,data=svmp.vpers,aes(x=SVMP.reg.start,xend=SVMP.reg.end,y=type,yend=type),lwd=1,color='grey70') +
  geom_point(inherit.aes = F, data=svmp.vpers, aes(x=(start+end)/2, y=type),size=2) +
  ylab('') +
  xlab('') +
  scale_fill_viridis_c(option = 'B') +
  scale_x_continuous(labels = comma,limits=c(SVMP.reg.start,SVMP.reg.end),expand=c(0,0)) +
  theme_classic(base_size = 14) +
  theme(axis.line.y = element_blank(),
        plot.title.position = 'plot',
        plot.title = element_text(color='black',face='bold',size = 14),
        axis.title.x=element_blank())


#### SVMP Array ATAC-seq with peaks ####
setwd('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/3_ATACseq_comparisons/6_GeneVignettes/SVMP_wholeArray')
peakfile <- paste('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/3_ATACseq_comparisons/1_chromVAR/merged_peaks/allSample_mergedPeaks.bed')
peaks <- getPeaks(peakfile, sort_peaks = FALSE)
peaks <- resize(peaks, width = 500, fix = "center")
peaks <- as.data.frame(peaks) %>% filter(seqnames == 'scaffold-mi1')

# Load data and plot
track.folder = '/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/3_ATACseq_comparisons/1_chromVAR/bigwigs/'

# Add BAM info
sample.meta = data.frame(SampleName = gsub('.bw', '' , grep('bw$', list.files(track.folder), value = T)))
sample.meta <- sample.meta %>% 
  mutate(Type = str_extract(SampleName, "^[^_]+")) %>% 
  mutate(Group = str_extract(SampleName, 'viridis|concolor|lutosus|cerberus')) %>% 
  filter(Group != 'lutosus') %>% # remove lutosus for now because its ATAC-seq quality is bad
  #filter(Type != 'CV1096') %>% # remove CV1096 for now because we don't have a genome
  arrange(match(Group, c('viridis', 'concolor', 'cerberus')))

# load regions of interest from bam files
chrom = "scaffold-mi1"
start_pos = floor(13901005 / 1000 ) * 1000
end_pos = ceiling(14424729 / 1000 ) * 1000 
track.df = LoadTrackFile(track.folder = track.folder, format = "bw", # change to bam
                         bamcoverage.path = '/Users/sidgopalan/miniconda3/bin/bamCoverage',
                         meta.info = sample.meta,
                         single.nuc = F, # change to T for BAM
                         region = paste0(chrom, ':', start_pos, '-', end_pos) # short hand for getting a region, change region to single.nuc.region for BAM,
)

score_variance <- track.df %>%
  group_by(seqnames, start, end) %>%
  summarise(score = var(score))
score_variance$Type <- 'Variance'
score_variance$Group <- 'Variance'

combined_df <- bind_rows(track.df, score_variance)
  

p2 <- ggplot() +
  theme_minimal() +
  geom_col(data = combined_df, aes(x = start, y = score, fill = Group)) +
  #annotate(geom = "rect",
  #         xmin = peaks$start,
  #         xmax = peaks$end,
  #         ymin = -Inf,
  #         ymax = +Inf,
  #z         alpha = 0.2) +
  facet_wrap(~factor(Type, levels = sample.meta$Type), ncol = 1, strip.position = "right", scales = 'free_y') +
  scale_x_continuous(limits = c(start_pos,end_pos), expand = c(0,0), labels = label_number(scale = 1e-6)) +
  scale_y_continuous(breaks = pretty_breaks(n=2)) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        plot.margin = margin(10,5,5,5),
        #strip.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 1, title = NULL)) +
  scale_fill_manual(values = c("viridis" = "black", # #2E8B58
                               "concolor" = "black", # #D9A528
                               "cerberus" = "black",
                               "Variance" = "blue"),
                    labels = c(expression(italic("C. viridis")),
                               expression(italic("C. concolor")),
                               expression(italic("C. cerberus")),
                               breaks = c("viridis", "concolor", "cerberus"))) +
  labs(x = paste0("pos. on ", chrom, ' (Mb)'), y = "ATAC-seq read density") +
  ggtitle("Accessibility at SVMP array")


plot_grid(p1, p2, nrow = 2, axis = 'lr', align = 'hv', rel_heights = c(0.2,1))

#### Make TFBS P/A plots ####
t <- read.table('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/2_TFBS_Presence_Absence/z_outputs/1_Ciiider_FP_tables/z_old/PER_ciiider_table_fp_Aug2023.txt', header = T)

all_combs <- t %>% 
  filter(gene == 'PER12_SVMP6') %>% 
  select(Sample.Name, TFBS_ID) %>% 
  tidyr::expand(Sample.Name, TFBS_ID)

t <- t %>% 
  filter(gene == 'PER12_SVMP6') %>% 
  select(Sample.Name, TFBS_ID) %>% 
  mutate(Presence = 1)

x <- full_join(all_combs, t) %>% mutate(Presence = if_else(is.na(Presence), 'Absent', 'Present'))

# Relevant TFBS (look at alignment and x):
# PER12_SVMP6:Arid3a:213-218 (SNP)
# PER12_SVMP6:FOS:209-219 (SNP)
# PER12_SVMP6:Bhlha15:239-246 (INDEL)
# PER12_SVMP6:FOXC2:232-243 (INDEL)
# PER12_SVMP6:GATA4:238-249 (INDEL)
# PER12_SVMP6:GATA6:238-250 (INDEL)
# PER12_SVMP6:Arid3a:237-242 (INDEL)


#x$TFBS_ID <- (sapply(strsplit(x$TFBS_ID,':'), function(x) x[2]))

# selected_TFBS <- x %>%
#   distinct(TFBS_ID) %>%
#   sample_n(20) %>%
#   pull(TFBS_ID) # Take a random sample of TFBS for visualization purposes

selected_TFBS <- c('PER12_SVMP6:Arid3a:213-218','PER12_SVMP6:FOS:209-219','PER12_SVMP6:Bhlha15:239-246', 'PER12_SVMP6:FOXC2:232-243', 'PER12_SVMP6:GATA4:238-249', 'PER12_SVMP6:GATA6:238-250', 'PER12_SVMP6:Arid3a:237-242')

# Filter 'x' to include only the selected TFBS_ID and all corresponding Sample.Name values
x_selected <- x %>%
  filter(TFBS_ID %in% selected_TFBS) %>% 
  mutate(Presence = ifelse(Sample.Name %in% c('CV0985_concolor_Other_F', 'CV1090_cerberus_Other_M'), 'Absent_variant', Presence))

ggplot(x_selected, aes(TFBS_ID, Sample.Name, fill=Presence)) + 
  geom_tile(color = 'black') +
  scale_fill_manual(values = c('grey70','grey41', 'lawngreen')) +
  scale_y_discrete(limits = rev(gsub('\\..*', '', sample.meta$SampleName))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
#### Make ATAC-seq variation graphs for SVMP genes ####
# ATACseq score var within specific enhancers or promoters of each gene? (box/whisker or violin). Could do within viridis and overall ones also separately (violins)
# p1 All promoters per gene
# p2 Viridis promoters per gene
# p3 All enhancers per gene
# p4 Viridis enhancers per gene

ATAC_mat <- read.table('~/Dropbox/CastoeLabFolder/projects/CVV_Gene_Regulation/_12SnakeVenomRegulation/data/atacseq/scaleFactorNormCounts.June2023.txt', header = T)
bed <- read.table('/Volumes/SeagatePortableDrive/12_snakeVenomVar/Figure_and_Scripts/3_ATACseq_comparisons/6_GeneVignettes/SVMP_wholeArray/SVMP_pr_enh_ATAC_scores/SVMP_CREs_Sept2023.bed', col.names = c('seqnames', 'start', 'end', 'name','type'))

ATAC_scores <- read.table('./SVMP_pr_enh_ATAC_scores/AvgCounts.SVMP_CRE_Sept2023.txt', header = T)
ATAC_scores <- ATAC_scores %>% rename_with(~sub("\\..*", "", .), everything())
ATAC_scores <- ATAC_scores %>% left_join(bed %>% select(start,end,name,type), by = c('start','end')) # Add bed info back

ATAC_scores <- ATAC_scores %>% 
  pivot_longer(cols = starts_with('CV'), values_to = 'ATAC_peak_score', names_to = 'SampleID') #%>% 
  #filter(grepl('viridis', SampleID)) # Keep only viridis


ATAC_scores %>% ungroup() %>% # enhancer/promoter # Change P1 to P4
ggplot(., aes(x = ATAC_peak_score, y = name)) +
  geom_boxplot(color="#0a0dc4", fill="#03dffc", alpha=0.2) + # or geom_boxplot
  theme_bw() + 
  ylab(NULL) +
  xlab("Accessibility\nscores") +
  guides(fill = 'none') +
  theme(axis.text.y = element_text(angle = 15, hjust = 1),
        plot.margin = margin(5,5,5,5)) + 
  scale_y_discrete(limits = c("Promoter_SVMP1", "Promoter_SVMP2", "Promoter_SVMP3", "Promoter_SVMP4", "Promoter_SVMP5", "Promoter_SVMP6", "Promoter_SVMP7", "Promoter_SVMP8", "Promoter_SVMP9", "Promoter_SVMP10", "Promoter_SVMP11",
                              "PER9_SVMP5", "PER10_SVMP5", "PER11_SVMP5", "PER12_SVMP6","PER13_SVMP7","PER14_SVMP8","PER15_SVMP10")) +
  ggtitle('All samples, CRE accessibility scores')



