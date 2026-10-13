# FIGURE 3 — HGNC gene-group enrichment
# Run from the GitHub repository root. No files are written.
################################################################################

# 1. Libraries
library(dplyr)
library(tidyr)
library(readr)
library(data.table)
library(ggplot2)

# 2. Input directory
input_dir <- "input/Fig3"
in3 <- function(x) file.path(input_dir,x)

# 3. Load and annotate early methylation regions
read_early <- function(){
  anno <- data.table::fread(in3('element_clone_raw_avg_promoteranno.tsv.gz')) %>% select(region_id,segment_id,promoter_gene_name,promoter_start,promoter_end,gene_promoter) %>% distinct()
  x <- data.table::fread(in3('clone_state_mode_segment_for_phylo.tsv.gz')) %>% left_join(anno)
  x$chrom <- as.character(x$chrom)
  gtf <- readr::read_tsv(in3('Homo_sapiens.GRCh37.87.gtf.gz'),comment='#',col_names=FALSE,show_col_types=FALSE) %>% filter(X3=='gene') %>% transmute(chrom=as.character(X1),gene_name=sub('.*gene_name "([^"]+)".*','\\1',X9),gene_biotype=ifelse(grepl('gene_biotype "',X9),sub('.*gene_biotype "([^"]+)".*','\\1',X9),ifelse(grepl('gene_source "',X9),sub('.*gene_source "([^"]+)".*','\\1',X9),NA_character_))) %>% mutate(priority=ifelse(gene_biotype=='pseudogene',1,2)) %>% arrange(chrom,gene_name,priority) %>% distinct(chrom,gene_name,.keep_all=TRUE) %>% select(-priority)
  x <- x %>% rename(gene_name=promoter_gene_name) %>% left_join(gtf,by=c('chrom','gene_name'))
  names(x) <- sub('_E_merge$|_E$','',names(x))
  clones <- read_tsv(in3('clone_info_260917.tsv'),show_col_types=FALSE) %>% filter(WGS=='Y',EM_seq=='Y',!grepl('HC14',sampleID)) %>% pull(sampleID)
  x %>% select(all_of(c(names(x)[!grepl('-',names(x))],intersect(clones,names(x)))))
}

# 4. Load RNA expression and sample information
read_tpm <- function() read_tsv(in3('HCs_woBlood_woHC19_exp_bc_TPM_250805.tsv'),show_col_types=FALSE) %>% tibble::column_to_rownames('gene_name')
read_samples <- function() read_tsv(in3('final_rna_sample_list.tsv'),show_col_types=FALSE) %>% filter(donor!='HC14') %>% pull(sample_id)

# Requires the original tested-gene universe (gene_name, gene_chrom, gene_biotype).
early_me <- read_early(); hgnc <- data.table::fread(in3('hgnc_total.txt'))
stopifnot(file.exists(in3('anno_with_meth_final.tsv')))
background <- read_tsv(in3('anno_with_meth_final.tsv'),show_col_types=FALSE)
gene_early_anno <- bind_rows(early_me %>% filter(gene_biotype!='pseudogene',!gene_name %in% c('Y_RNA','snoU13')) %>% distinct(gene_name,.keep_all=TRUE) %>% transmute(chrom,gene_name,group='Early'),background %>% filter(gene_biotype!='pseudogene') %>% distinct(gene_name,.keep_all=TRUE) %>% transmute(chrom=as.character(gene_chrom),gene_name,group='non-Early'))
hgnc_map <- bind_rows(hgnc %>% transmute(gene_name=symbol,symbol,match_type='Current'),hgnc %>% filter(!is.na(prev_symbol),prev_symbol!='') %>% separate_rows(prev_symbol,sep='\\|') %>% transmute(gene_name=trimws(prev_symbol),symbol,match_type='Previous')) %>% arrange(factor(match_type,levels=c('Current','Previous'))) %>% distinct(gene_name,.keep_all=TRUE)
gene_early <- gene_early_anno %>% left_join(hgnc_map,by='gene_name') %>% filter(!is.na(symbol)) %>% distinct(symbol,.keep_all=TRUE)
family_anno <- hgnc %>% select(symbol,gene_group) %>% filter(!is.na(gene_group),gene_group!='') %>% separate_rows(gene_group,sep='\\|') %>% mutate(gene_group=trimws(gene_group)) %>% distinct(symbol,gene_group)
universe_n <- n_distinct(gene_early$symbol); early_n <- n_distinct(gene_early$symbol[gene_early$group=='Early'])
family_result <- family_anno %>% inner_join(select(gene_early,symbol,group),by='symbol') %>% group_by(gene_group) %>% summarise(family_total=n_distinct(symbol),early_n_family=n_distinct(symbol[group=='Early']),early_genes=paste(sort(unique(symbol[group=='Early'])),collapse=', '),.groups='drop') %>% mutate(expected=early_n*family_total/universe_n,fold_enrichment=early_n_family/expected) %>% rowwise() %>% mutate(p_value=fisher.test(matrix(c(early_n_family,family_total-early_n_family,early_n-early_n_family,universe_n-family_total-early_n+early_n_family),nrow=2))$p.value) %>% ungroup() %>% mutate(FDR=p.adjust(p_value,'BH'))
group_order <- c('HOXL subclass homeoboxes','NKL subclass homeoboxes and pseudogenes','PRD class homeoboxes and pseudogenes','TALE class homeoboxes and pseudogenes','KRAB domain containing','Zinc fingers C2H2-type','5S ribosomal RNAs','Antisense RNAs','H2B histones')
set.seed(123)
examples <- family_anno %>% filter(gene_group %in% group_order) %>% inner_join(gene_early %>% filter(group=='Early') %>% distinct(symbol),by='symbol') %>% group_by(gene_group) %>% slice_sample(n=3) %>% summarise(example_genes=paste(sort(symbol),collapse=', '),.groups='drop')
plot_df <- family_result %>% filter(gene_group %in% group_order) %>% left_join(examples,by='gene_group') %>% mutate(gene_group=factor(gene_group,levels=rev(group_order)))
p <- ggplot(plot_df,aes(fold_enrichment,gene_group))+geom_col(width=.65)+geom_vline(xintercept=1,linetype='dashed')+geom_text(aes(label=paste0(early_n_family,'/',family_total)),hjust=-.15,size=3)+geom_text(aes(x=max(plot_df$fold_enrichment)*1.25,label=example_genes),hjust=0,size=3,fontface='italic')+scale_x_continuous(limits=c(0,max(plot_df$fold_enrichment)*2))+labs(x='EEC-associated gene enrichment (fold)',y=NULL)+theme_classic()
print(p)
