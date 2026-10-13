# FIGURE 3 — EEC promoter and gene proportions
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

# Requires the original promoter–gene mapping (promoter_id, gene_name).
early_me <- read_early()
reg <- read_tsv(in3('ensembl_regulatory_sorted.grch37.gff.gz'),col_names=FALSE,show_col_types=FALSE)
names(reg) <- c('chrom','build','feature','start','end','score','strand','phase','info')
reg <- reg %>% filter(feature=='promoter',chrom %in% as.character(1:22)) %>% mutate(reg_id=sub('^id=([^;]+).*','\\1',info)) %>% distinct(chrom,start,end,.keep_all=TRUE)
stopifnot(file.exists(in3('promoter_gene_annotation.tsv')))
p_g_anno <- read_tsv(in3('promoter_gene_annotation.tsv'),show_col_types=FALSE)
prom_candi <- intersect(unique(p_g_anno$promoter_id),reg$reg_id)
prom_early <- early_me %>% transmute(chrom,start=promoter_start+1,end=promoter_end,early=TRUE) %>% distinct()
prom <- reg %>% left_join(prom_early,by=c('chrom','start','end')) %>% filter(reg_id %in% prom_candi) %>% mutate(group=ifelse(!is.na(early),'EEE','Non-EEE'))
all_gene <- p_g_anno %>% filter(promoter_id %in% prom_candi) %>% pull(gene_name) %>% unique() %>% na.omit()
early_gene <- p_g_anno %>% filter(promoter_id %in% prom$reg_id[prom$group=='EEE']) %>% pull(gene_name) %>% unique() %>% na.omit()
df <- bind_rows(prom %>% count(group) %>% mutate(type='Promoters'),tibble(type='Genes',group=c('EEE','Non-EEE'),n=c(length(early_gene),length(setdiff(all_gene,early_gene))))) %>% group_by(type) %>% mutate(prop=n/sum(n)) %>% ungroup() %>% mutate(type=factor(type,levels=c('Promoters','Genes')))
p <- ggplot(df,aes(type,n,fill=group))+geom_col(width=.6)+geom_text(aes(label=paste0(format(n,big.mark=','),'\n(',round(100*prop,1),'%)')),position=position_stack(vjust=.5),size=4)+labs(x=NULL,y='Count',fill='EEC')+theme_classic()
print(p)
