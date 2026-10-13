# FIGURE 3 — Representative methylation–expression associations
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

early_me <- as.data.frame(read_early()); normalized_count <- read_tpm(); common_cells <- read_samples()
common_cells <- intersect(common_cells,intersect(names(early_me),colnames(normalized_count)))
donor <- ifelse(grepl('^HC[0-9]+-',common_cells),sub('-.*','',common_cells),'HC50')
res <- lapply(seq_len(nrow(early_me)),function(i){
  gene <- early_me$gene_name[i]; if(is.na(gene)||!gene %in% rownames(normalized_count)) return(NULL)
  df <- data.frame(cell=common_cells,donor=donor,state=as.numeric(unlist(early_me[i,common_cells])),expr=as.numeric(normalized_count[gene,common_cells])) %>% filter(complete.cases(.)) %>% group_by(donor) %>% mutate(expr_rank=(rank(expr,ties.method='average')-1)/pmax(n()-1,1)) %>% ungroup()
  good <- df %>% group_by(donor) %>% summarise(meth_ok=sum(table(state)>=2)>=2,expr_ok=sum(expr>=1)>=2,.groups='drop') %>% filter(meth_ok & expr_ok) %>% pull(donor)
  df <- filter(df,donor %in% good); if(nrow(df)<4||n_distinct(df$state)<2) return(NULL)
  fit <- summary(lm(expr_rank~state,data=df))$coefficients
  data.frame(row=i,gene_name=gene,beta=fit['state','Estimate'],p_value=fit['state','Pr(>|t|)'],n_cell=nrow(df),n_donor=n_distinct(df$donor))
}) %>% bind_rows()
best_df <- res %>% group_by(gene_name) %>% group_modify(~{neg <- filter(.x,beta<0); if(nrow(neg)) slice_min(neg,p_value,n=1,with_ties=FALSE) else slice_min(.x,p_value,n=1,with_ties=FALSE)}) %>% ungroup() %>% mutate(adj_p=p.adjust(p_value,'BH'))
volcano <- best_df %>% filter(!is.na(beta),!is.na(adj_p)) %>% mutate(status=case_when(adj_p<.05 & beta>0~'Positive',adj_p<.05 & beta<0~'Negative',TRUE~'NS'),neg_log10_fdr=-log10(pmax(adj_p,.Machine$double.xmin)))
p <- ggplot(volcano,aes(beta,neg_log10_fdr,color=status))+geom_point(size=.8,alpha=.5)+geom_hline(yintercept=-log10(.05),linetype='dashed')+geom_vline(xintercept=0,linetype='dashed',color='grey50')+scale_color_manual(values=c(Negative='#2166AC',NS='grey75',Positive='#B2182B'))+coord_cartesian(xlim=c(-1,1))+labs(x='Regression coefficient (β)',y='-log10(FDR)',color=NULL)+theme_bw()+theme(panel.grid=element_blank())

# 7. Plot representative genes (Fig. 3d–e)
donor_levels <- c('HC01','HC02','HC08','HC09','HC12','HC13','HC15','HC16','HC17','HC18','HC19','HC20','HC21','HC50')
donor_colors <- setNames(scales::hue_pal()(length(donor_levels)),donor_levels)
plots <- lapply(c('MAP3K7CL','HOXD8'),function(g){
  x <- best_df %>% filter(gene_name==g) %>% slice(1); stopifnot(nrow(x)==1)
  df <- data.frame(cell=common_cells,donor=donor,state=as.numeric(unlist(early_me[x$row,common_cells])),expr=as.numeric(normalized_count[g,common_cells])) %>% filter(complete.cases(.)) %>% group_by(donor) %>% mutate(expr_rank=(rank(expr,ties.method='average')-1)/pmax(n()-1,1)) %>% ungroup()
  good <- df %>% group_by(donor) %>% summarise(meth_ok=sum(table(state)>=2)>=2,expr_ok=sum(expr>=1)>=2,.groups='drop') %>% filter(meth_ok & expr_ok) %>% pull(donor)
  df <- filter(df,donor %in% good) %>% mutate(state=factor(state,levels=c(0,.5,1),labels=c('U/U','U/M','M/M')))
  ggplot(df,aes(state,expr,color=donor))+geom_jitter(width=.15,height=0,size=1.5,alpha=.8)+scale_color_manual(values=donor_colors)+labs(x='Epigenotype',y='TPM',color='Donor',title=paste0(g,'  β = ',round(x$beta,3),', P = ',format.pval(x$p_value,digits=2)))+theme_classic()
})
for(p in plots) print(p)
