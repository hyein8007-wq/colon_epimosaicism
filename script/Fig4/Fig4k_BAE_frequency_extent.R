################################################################################
# Figure 4k. BAE frequency and extent in promoter U/U clones
# Run from the repository root. No output files are written.
################################################################################

# 1. Libraries
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)

# 2. Input data and matched female RNA clones
input_dir <- "input/Fig4"
final_clones <- read_tsv(file.path(input_dir,"final_rna_sample_list.tsv"),show_col_types=FALSE) %>% filter(sex=="F") %>% pull(sample_id)
m <- read_tsv(file.path(input_dir,"chrX_cgi_total_summary_pattern.tsv"),show_col_types=FALSE) %>% filter(par=="nonPAR")
ase_res_final <- read_tsv(file.path(input_dir,"ase_res_female_chrX_final_260919_v2.tsv"),show_col_types=FALSE)
common_clones <- intersect(final_clones,names(m))
stopifnot(length(common_clones)>0,all(c("genomic_region","gene_name","pattern") %in% names(m)),all(c("gene_name","sample","final_group","VAF") %in% names(ase_res_final)))

# 3. Identify promoter U/U gene-clone pairs
m_long <- m %>% filter(genomic_region=="promoter",!is.na(gene_name),gene_name!="",gene_name!=".") %>%
  separate_rows(gene_name,sep=",") %>% mutate(gene_name=str_trim(gene_name)) %>%
  pivot_longer(all_of(common_clones),names_to="sample",values_to="methylation") %>% filter(!is.na(methylation))
uu_gene_clone <- m_long %>% group_by(gene_name,sample) %>% summarise(UU=any(methylation<0.2),.groups="drop") %>% filter(UU)

# 4. Match ASE status in the same clones
ase_check <- ase_res_final %>% filter(sample %in% common_clones) %>% group_by(gene_name,sample) %>%
  summarise(n_group=n_distinct(final_group),final_group=first(final_group),.groups="drop")
stopifnot(all(ase_check$n_group==1),anyDuplicated(uu_gene_clone[c("gene_name","sample")])==0)
ase_gene_clone <- ase_check %>% select(-n_group)
uu_ase <- uu_gene_clone %>% left_join(ase_gene_clone,by=c("gene_name","sample")) %>% mutate(final_group=replace_na(final_group,"w/oSNP"))
stopifnot(anyDuplicated(uu_ase[c("gene_name","sample")])==0)

# 5. Assign gene-level CGI pattern (A > B > C > D)
pattern_info_gene <- m_long %>% filter(!is.na(pattern),pattern %in% c("A","B","C","D")) %>%
  distinct(gene_name,pattern) %>% mutate(priority=match(pattern,c("A","B","C","D"))) %>%
  arrange(gene_name,priority) %>% distinct(gene_name,.keep_all=TRUE) %>% select(gene_name,pattern)

# 6. BAE frequency among ASE-informative U/U clones
bae_freq <- uu_ase %>% group_by(gene_name) %>%
  summarise(n_UU=n(),n_informative=sum(final_group %in% c("BAE","MAE")),n_BAE=sum(final_group=="BAE"),
            n_MAE=sum(final_group=="MAE"),n_LowDepth=sum(final_group=="LowDepth"),n_wSNP=sum(final_group=="w/oSNP"),
            BAE_freq=ifelse(n_informative>0,n_BAE/n_informative,NA_real_),.groups="drop")

# 7. BAE extent: mean minor VAF across BAE SNPs, then across matched U/U clones
bae_clone_vaf <- ase_res_final %>% filter(sample %in% common_clones,final_group=="BAE",!is.na(VAF),VAF>=0,VAF<=1) %>%
  mutate(minor_VAF=pmin(VAF,1-VAF)) %>% group_by(gene_name,sample) %>%
  summarise(clone_BAE_extent=mean(minor_VAF,na.rm=TRUE),.groups="drop")
uu_bae_extent <- uu_ase %>% filter(final_group=="BAE") %>% select(gene_name,sample) %>%
  inner_join(bae_clone_vaf,by=c("gene_name","sample"))
bae_extent <- uu_bae_extent %>% group_by(gene_name) %>%
  summarise(BAE_extent=mean(clone_BAE_extent,na.rm=TRUE),n_BAE_extent_clones=n(),.groups="drop")

# 8. Prepare Figure 4k data
plot_k <- bae_freq %>% filter(n_informative>0) %>% left_join(bae_extent,by="gene_name") %>%
  left_join(pattern_info_gene,by="gene_name") %>%
  mutate(BAE_extent=ifelse(n_BAE==0,0,BAE_extent),n_BAE_extent_clones=ifelse(n_BAE==0,0,n_BAE_extent_clones),
         pattern=factor(pattern,levels=c("A","B","C","D")))
stopifnot(anyDuplicated(plot_k$gene_name)==0)
print(plot_k %>% summarise(n_gene=n(),n_with_BAE=sum(n_BAE>0),n_MAE_only=sum(n_BAE==0 & n_MAE>0)))
print(plot_k %>% filter(!is.na(pattern)) %>% group_by(pattern) %>%
        summarise(n_gene=n(),n_with_BAE=sum(n_BAE>0),n_without_BAE=sum(n_BAE==0),
                  median_BAE_freq=median(BAE_freq,na.rm=TRUE),median_BAE_extent=median(BAE_extent,na.rm=TRUE),.groups="drop"))

# 9. Figure 4k: BAE frequency versus extent
p_k <- plot_k %>% filter(!is.na(BAE_freq),!is.na(BAE_extent),!is.na(pattern)) %>%
  mutate(draw_order=match(as.character(pattern),c("D","C","B","A"))) %>% arrange(draw_order) %>%
  ggplot(aes(x=BAE_freq,y=BAE_extent,color=pattern)) +
  stat_density_2d(aes(color=NULL),color="grey75",linewidth=0.6,bins=5) +
  geom_point(data=~filter(.x,BAE_freq!=0 | BAE_extent!=0) %>% arrange(draw_order),size=3,alpha=0.8) +
  geom_jitter(data=~filter(.x,BAE_freq==0 & BAE_extent==0) %>% arrange(draw_order),size=3,alpha=0.8,width=0.015,height=0.006) +
  scale_color_manual(values=c("A"="#b82860","B"="#e4a126","C"="#0095a2","D"="#253e51")) +
  scale_x_continuous(limits=c(-0.03,1),breaks=seq(0,1,0.25)) +
  scale_y_continuous(limits=c(-0.015,0.5),breaks=seq(0,0.5,0.1)) +
  labs(x="BAE frequency in U/U clones",y="BAE extent in U/U clones",color="Pattern") +
  theme_classic(base_size=16) +
  theme(axis.line=element_line(linewidth=0.8),axis.ticks=element_line(linewidth=0.8),
        axis.title=element_text(size=17),axis.text=element_text(size=14,color="black"),
        legend.position="right",legend.title=element_text(face="bold"))
print(p_k)
