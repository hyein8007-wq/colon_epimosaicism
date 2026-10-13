################################################################################
# Figure 4j. ASE status of genes with hypomethylated (U/U) promoters
# Run from the repository root. No output files are written.
################################################################################

# 1. Libraries
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(grid)

# 2. Input data and matched female RNA clones
input_dir <- "input/Fig4"
final_clones <- read_tsv(file.path(input_dir,"final_rna_sample_list.tsv"),show_col_types=FALSE) %>% filter(sex=="F") %>% pull(sample_id)
m <- read_tsv(file.path(input_dir,"chrX_cgi_total_summary_pattern.tsv"),show_col_types=FALSE) %>% filter(par=="nonPAR")
ase_res_final <- read_tsv(file.path(input_dir,"ase_res_female_chrX_final_260919_v2.tsv"),show_col_types=FALSE)
common_clones <- intersect(final_clones,names(m))
stopifnot(length(common_clones)>0,all(c("genomic_region","gene_name","pattern") %in% names(m)),all(c("gene_name","sample","final_group") %in% names(ase_res_final)))

# 3. Promoter CGI methylation: gene-by-clone U/U calls
m_long <- m %>% filter(genomic_region=="promoter",!is.na(gene_name),gene_name!="",gene_name!=".") %>%
  separate_rows(gene_name,sep=",") %>% mutate(gene_name=str_trim(gene_name)) %>%
  pivot_longer(all_of(common_clones),names_to="sample",values_to="methylation") %>% filter(!is.na(methylation))
uu_gene_clone <- m_long %>% group_by(gene_name,sample) %>% summarise(UU=any(methylation<0.2),.groups="drop") %>% filter(UU)

# 4. Match U/U gene-clone pairs to ASE status
ase_check <- ase_res_final %>% filter(sample %in% common_clones) %>% group_by(gene_name,sample) %>%
  summarise(n_group=n_distinct(final_group),final_group=first(final_group),.groups="drop")
stopifnot(all(ase_check$n_group==1))
ase_gene_clone <- ase_check %>% select(-n_group)
uu_ase <- uu_gene_clone %>% left_join(ase_gene_clone,by=c("gene_name","sample")) %>% mutate(final_group=replace_na(final_group,"w/oSNP"))
stopifnot(anyDuplicated(u u_gene_clone[c("gene_name","sample")])==0)

# 5. Assign one CGI pattern per gene (A > B > C > D)
pattern_info_gene <- m_long %>% filter(!is.na(pattern),pattern %in% c("A","B","C","D")) %>%
  distinct(gene_name,pattern) %>% mutate(priority=match(pattern,c("A","B","C","D"))) %>%
  arrange(gene_name,priority) %>% distinct(gene_name,.keep_all=TRUE) %>% select(gene_name,pattern)

# 6. Classify ASE per gene (BAE > MAE > LowDepth > w/oSNP)
gene_ase <- uu_ase %>% group_by(gene_name) %>%
  summarise(ASE=case_when(any(final_group=="BAE")~"BAE",any(final_group=="MAE")~"MAE",any(final_group=="LowDepth")~"LowDepth",TRUE~"w/oSNP"),n_UU=n(),.groups="drop") %>%
  left_join(pattern_info_gene,by="gene_name") %>% filter(!is.na(pattern)) %>%
  mutate(pattern=factor(pattern,levels=c("A","B","C","D")),ASE=factor(ASE,levels=c("w/oSNP","LowDepth","MAE","BAE")))
stopifnot(anyDuplicated(gene_ase$gene_name)==0)

# 7. Figure 4j: proportion of ASE classes within each CGI pattern
plot_j <- gene_ase %>% count(pattern,ASE,.drop=FALSE) %>% group_by(pattern) %>% mutate(prop=if_else(sum(n)>0,n/sum(n),0)) %>% ungroup()
p_j <- ggplot(plot_j,aes(x=pattern,y=prop,fill=ASE)) + geom_col(width=0.85) +
  scale_fill_manual(values=c("BAE"="#A99BC2","MAE"="#493277","LowDepth"="#C7C7C7","w/oSNP"="#969696"),breaks=c("BAE","MAE","LowDepth","w/oSNP")) +
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.25,0.5,0.75,1),labels=sprintf("%.2f",c(0,0.25,0.5,0.75,1)),expand=c(0,0)) +
  labs(x="Pattern",y="Proportion",fill="ASE",title=paste0("No. of genes = ",n_distinct(gene_ase$gene_name))) +
  theme_classic(base_size=16) +
  theme(axis.line=element_line(linewidth=0.8),axis.ticks=element_line(linewidth=0.8),axis.title=element_text(size=17),axis.text=element_text(size=14,color="black"),legend.position="top",legend.title=element_text(face="bold"),legend.key.height=unit(0.8,"lines"),legend.key.width=unit(0.8,"lines"),plot.title=element_text(size=14,hjust=0.5))
print(p_j)
