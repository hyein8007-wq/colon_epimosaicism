# FIGURE 3 — APC promoter methylation and transcript usage
# Run from the GitHub repository root. No files are written.
################################################################################

# 1. Libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# 2. Input directory
input_dir <- "input/Fig3"
in3 <- function(x) file.path(input_dir,x)
read_tpm <- function() read_tsv(in3('HCs_woBlood_woHC19_exp_bc_TPM_250805.tsv'),show_col_types=FALSE) %>% tibble::column_to_rownames('gene_name')
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
# APC count inputs are in input/APC_exon_count/ (sibling of input/Fig3).
ss <- c('HC02','HC13','HC15','HC16')
count_files <- file.path('input','APC_exon_count',paste0(ss,'_apc_exon_spliced_count_v3.tsv'))
stopifnot(all(file.exists(count_files)))
apc_count <- bind_rows(lapply(count_files,function(f) read_tsv(f,show_col_types=FALSE))) %>% filter(!grepl('_P',bam)) %>% mutate(sample=sub('\\.mark\\.bam$','',bam),promoter1A_me=sample %in% c('HC02-02','HC02-17','HC13-01','HC15-08','HC16-09'),prop_1A=exon2/(exon1+exon2),sample_id=sub('-.*','',sample))
print(ggboxplot(filter(apc_count,exon1+exon2>10),x='promoter1A_me',y='prop_1A',add='jitter')+stat_compare_means()+theme_bw())
# Methylation heatmap requires input/Fig3/cpg_annotation.tsv (chrom,pos,region)
# and input/Fig3/clone_order.tsv (donor,sample_id,order).
stopifnot(file.exists(in3('clone_order.tsv')),file.exists(in3('cpg_annotation.tsv')))
orders <- read_tsv(in3('clone_order.tsv'),show_col_types=FALSE) %>% filter(donor %in% ss) %>% arrange(match(donor,ss),order) %>% pull(sample_id)
cpg <- read_tsv(in3('cpg_annotation.tsv'),show_col_types=FALSE) %>% filter(as.character(chrom)=='5') %>% select(pos,region) %>% distinct(pos,.keep_all=TRUE)
m <- read_tsv(in3('methyl_frac_all_chr5.tsv.gz'),show_col_types=FALSE) %>% filter(between(pos,112038115,112075277)) %>% left_join(cpg,by='pos') %>% arrange(desc(pos))
clones <- intersect(orders,names(m)); stopifnot(length(clones)>0)
ha <- rowAnnotation(CpG=m$region,col=list(CpG=c(Island='#054d35',Shore='#abd4c9',Shelf='#f5e889',OpenSea='#e8ebea')))
h <- Heatmap(as.matrix(m[,clones,drop=FALSE]),name='ME_pct',col=colorRamp2(c(-2,0,.4,.6,1),c('#cccccc','#1465b5','#f5f29f','#f5f29f','#b51914')),left_annotation=ha,column_split=sub('-.*','',clones),cluster_rows=FALSE,cluster_columns=FALSE,cluster_column_slices=FALSE,show_row_names=FALSE)
draw(h)
