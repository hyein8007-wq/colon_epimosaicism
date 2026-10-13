# FIGURE 3 — MGMT promoter methylation and expression
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
# Required: input/Fig3/cpg_annotation.tsv (chrom,pos,cpg_context) and clone_order.tsv (donor,sample_id,order).
stopifnot(file.exists(in3('clone_order.tsv')),file.exists(in3('cpg_annotation.tsv')))
orders <- read_tsv(in3('clone_order.tsv'),show_col_types=FALSE) %>% filter(donor=='HC50') %>% arrange(order) %>% pull(sample_id)
cpg <- read_tsv(in3('cpg_annotation.tsv'),show_col_types=FALSE) %>% filter(as.character(chrom)=='10') %>% select(pos,cpg_context) %>% distinct(pos,.keep_all=TRUE)
m <- read_tsv(in3('class_converted_chr10.tsv.gz'),show_col_types=FALSE)
names(m) <- gsub('_EM$|_E_merge$|_E$|-EM$|-E$','',names(m))
m <- m %>% filter(between(pos,131264599,131266401)) %>% left_join(cpg,by='pos') %>% arrange(pos)
clones <- intersect(orders,names(m)); stopifnot(length(clones)>0)
ha <- HeatmapAnnotation(CpG=m$cpg_context,col=list(CpG=c(Island='#360c2a',Shore='#986b8c',Shelf='#cfbecb',OpenSea='#e9e7e6')))
h <- Heatmap(t(as.matrix(m[,clones,drop=FALSE])),name='ME_pct',col=colorRamp2(c(-2,0,.5,1),c('#cccccc','#003e73','#f2d8aa','#a51f27')),top_annotation=ha,cluster_rows=FALSE,cluster_columns=FALSE,show_row_names=FALSE)
draw(h)
tpm <- read_tpm(); clones <- intersect(orders,colnames(tpm))
gg <- tibble(sample_id=clones,exp=as.numeric(tpm['MGMT',clones])) %>% mutate(state=ifelse(sample_id %in% c('D2-21','DS-3','DS-1','D3-28','PS-21','D2-9'),'Half-methyl','Unmethyl'))
print(ggboxplot(gg,x='state',y='exp',add='jitter')+stat_compare_means()+coord_cartesian(ylim=c(0,16)))
