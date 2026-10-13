# FIGURE 4 - SOURCE-PRESERVING EXTRACT (NOT YET STANDALONE)
# Source: original_XCI_methylation.R; original line ranges: 1381-1699
# This extract preserves the original calculations and figure commands.
# It depends on objects from earlier sections and original server inputs.
# DO NOT publish as independently reproducible until dependencies are resolved.


# ===== ORIGINAL LINES 1381-1699 =====

# 1. Non-PAR X genes
gtf <- readr::read_tsv("input/Homo_sapiens.GRCh37.87.gtf.gz",comment="#",col_names=F) %>%
  dplyr::filter(X3=="gene") %>% dplyr::rename(chr=X1,start=X4,end=X5) %>%
  tidyr::separate(X9,into=c("x1","gene_id","x2","gene_version","x3","gene_name","x4","gene_source"),sep="\"") %>%
  dplyr::select(chr,start,end,gene_name,gene_id) %>%
  dplyr::filter(chr=="X",!is.na(gene_name),gene_name!="",!(start<=2699520 & end>=60001),!(start<=155260560 & end>=154931044)) %>%
  dplyr::distinct(gene_name,.keep_all=TRUE) %>%
  dplyr::mutate(X_arm=ifelse((start+end)/2<60600000,"Xp","Xq"))
nonpar_genes <- gtf$gene_name

# 2. Previous XCI annotation
gtex <- readr::read_tsv("input/Nature_GTEx_Tukiainene_2017.txt") %>%
  dplyr::rename(gene_name=`Gene name`,XCI_status=`Combined XCI status`) %>%
  dplyr::select(gene_name,XCI_status) %>% dplyr::distinct()

# 3. RNA ASE summary
rna_summary <- ase_gene_clone %>% dplyr::filter(gene_name %in% nonpar_genes) %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise(ASE=dplyr::case_when(any(final_group=="BAE")~"BAE",any(final_group=="MAE")~"MAE",any(final_group=="LowDepth")~"LowDepth",TRUE~"w/oSNP"),n_BAE_RNA=sum(final_group=="BAE"),n_MAE_RNA=sum(final_group=="MAE"),n_LowDepth_RNA=sum(final_group=="LowDepth"),n_informative_RNA=n_BAE_RNA+n_MAE_RNA,BAE_freq_RNA=ifelse(n_informative_RNA>0,n_BAE_RNA/n_informative_RNA,NA_real_),.groups="drop")

rna_extent <- ase_res_final %>% dplyr::filter(gene_name %in% nonpar_genes,sample %in% common_clones,final_group=="BAE",!is.na(VAF),VAF>=0,VAF<=1) %>%
  dplyr::mutate(minor_VAF=pmin(VAF,1-VAF)) %>% dplyr::group_by(gene_name,sample) %>%
  dplyr::summarise(clone_extent=mean(minor_VAF,na.rm=TRUE),.groups="drop") %>% dplyr::group_by(gene_name) %>%
  dplyr::summarise(BAE_extent_RNA=mean(clone_extent,na.rm=TRUE),.groups="drop")

rna_summary <- rna_summary %>% dplyr::left_join(rna_extent,by="gene_name") %>%
  dplyr::mutate(BAE_extent_RNA=dplyr::case_when(n_informative_RNA==0~NA_real_,n_BAE_RNA==0~0,TRUE~BAE_extent_RNA))

# 4. Promoter U/U summary
uu_summary <- uu_gene_clone %>% dplyr::filter(gene_name %in% nonpar_genes) %>%
  dplyr::group_by(gene_name) %>% dplyr::summarise(n_UU=dplyr::n_distinct(sample),UU_freq=n_UU/278,.groups="drop")

# 5. Methylation pattern
pattern_summary <- pattern_info_gene %>% dplyr::filter(gene_name %in% nonpar_genes) %>%
  dplyr::select(gene_name,Pattern=pattern)

# 6. ASE among U/U clones
uu_metrics <- bae_freq %>% dplyr::filter(gene_name %in% nonpar_genes) %>%
  dplyr::select(gene_name,n_BAE,n_MAE,n_informative,BAE_freq) %>%
  dplyr::rename(n_BAE_UU=n_BAE,n_MAE_UU=n_MAE,n_informative_UU=n_informative,BAE_freq_UU=BAE_freq) %>%
  dplyr::left_join(bae_extent %>% dplyr::filter(gene_name %in% nonpar_genes) %>% dplyr::select(gene_name,BAE_extent_UU=BAE_extent),by="gene_name") %>%
  dplyr::mutate(BAE_extent_UU=dplyr::case_when(n_informative_UU==0~NA_real_,n_BAE_UU==0~0,TRUE~BAE_extent_UU))

# 7. Non-PAR gene sets
rna_bae <- ase_gene_clone %>% dplyr::filter(gene_name %in% nonpar_genes,final_group=="BAE") %>% dplyr::distinct(gene_name) %>% dplyr::pull(gene_name)
meth_uu <- uu_gene_clone %>% dplyr::filter(gene_name %in% nonpar_genes) %>% dplyr::distinct(gene_name) %>% dplyr::pull(gene_name)
reported_escape <- gtex %>% dplyr::filter(gene_name %in% nonpar_genes,XCI_status %in% c("escape","variable")) %>% dplyr::distinct(gene_name) %>% dplyr::pull(gene_name)
all_genes <- Reduce(union,list(rna_bae,meth_uu,reported_escape))

# 8. Sheet 1: gene-level integration
supp8_summary <- gtf %>% dplyr::filter(gene_name %in% all_genes) %>%
  dplyr::mutate(X_arm=ifelse((start+end)/2<60600000,"Xp","Xq"),RNA_BAE=gene_name %in% rna_bae,Promoter_UU=gene_name %in% meth_uu,Reported_escape_variable=gene_name %in% reported_escape,
                Venn_group=dplyr::case_when(RNA_BAE & Promoter_UU & Reported_escape_variable~"All three",RNA_BAE & Promoter_UU~"RNA BAE + Promoter U/U",RNA_BAE & Reported_escape_variable~"RNA BAE + Reported",Promoter_UU & Reported_escape_variable~"Promoter U/U + Reported",RNA_BAE~"RNA BAE only",Promoter_UU~"Promoter U/U only",Reported_escape_variable~"Reported only")) %>%
  dplyr::left_join(rna_summary,by="gene_name") %>% dplyr::left_join(uu_summary,by="gene_name") %>%
  dplyr::left_join(pattern_summary,by="gene_name") %>% dplyr::left_join(uu_metrics,by="gene_name") %>%
  dplyr::left_join(gtex,by="gene_name") %>% dplyr::arrange(X_arm,start) %>%
  dplyr::select(gene_name,gene_id,chr,start,end,X_arm,RNA_BAE,Promoter_UU,Reported_escape_variable,Venn_group,ASE,n_BAE_RNA,n_MAE_RNA,n_LowDepth_RNA,n_informative_RNA,BAE_freq_RNA,BAE_extent_RNA,Pattern,n_UU,UU_freq,n_BAE_UU,n_MAE_UU,n_informative_UU,BAE_freq_UU,BAE_extent_UU,XCI_status)

# 9. Sheet 2: SNP-level ASE results
supp8_ase_raw <- ase_res_final

# 10. Excel
wb8 <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb8,"BAE_UU_Report_integration")
openxlsx::writeData(wb8,"BAE_UU_Report_integration",supp8_summary,withFilter=TRUE)
openxlsx::freezePane(wb8,"BAE_UU_Report_integration",firstRow=TRUE)

openxlsx::addWorksheet(wb8,"ASE_variant_results")
openxlsx::writeData(wb8,"ASE_variant_results",supp8_ase_raw,withFilter=TRUE)
openxlsx::freezePane(wb8,"ASE_variant_results",firstRow=TRUE)

openxlsx::saveWorkbook(wb8,file.path(outdir,"supplementary_table_10_chrX_gene_summary.xlsx"),overwrite=TRUE)

######################################################################################





























## RNA BAE vs methylation U/U ===================================================

# 동일 matched-clone universe에서 한 clone이라도 BAE인 gene
rna_bae <- ase_gene_clone %>% filter(final_group=="BAE") %>% distinct(gene_name) %>% pull(gene_name)

# 동일 matched-clone universe에서 한 clone이라도 promoter U/U인 gene
meth_uu <- uu_gene_clone %>% distinct(gene_name) %>% pull(gene_name)

# overlap 확인
cat("RNA BAE:",length(rna_bae),"\n","Methylation U/U:",length(meth_uu),"\n","Both:",length(intersect(rna_bae,meth_uu)),"\n","RNA only:",length(setdiff(rna_bae,meth_uu)),"\n","Methylation only:",length(setdiff(meth_uu,rna_bae)),"\n")

# Venn diagram
library(ggVennDiagram)
venn_list <- list("RNA BAE"=rna_bae,"Promoter U/U"=meth_uu)
p_venn <- ggVennDiagram(venn_list,label_alpha=0,edge_size=0.8,label="count")+scale_fill_gradient(low="white",high="grey70")+theme_void(base_size=16)+theme(legend.position="none")
p_venn

pdf("output/extended_venn_260919.pdf")
p_venn
dev.off()


## U/U–BAE clone-level matching =================================================

# matching = U/U clone 중 informative ASE가 있는 clone에서 BAE가 관찰된 비율
matching_all <- uu_ase %>% group_by(gene_name) %>% summarise(n_UU=n(),n_BAE=sum(final_group=="BAE"),n_MAE=sum(final_group=="MAE"),n_LowDepth=sum(final_group=="LowDepth"),n_wSNP=sum(final_group=="w/oSNP"),n_informative=n_BAE+n_MAE,matching=ifelse(n_informative>0,n_BAE/n_informative,NA_real_),.groups="drop")
matching <- matching_all %>% left_join(pattern_info_gene,by="gene_name")
matching %>% count(gene_name) %>% filter(n>1) # 반드시 0 rows

# RNA와 methylation 양쪽에서 검출된 gene의 matching 정도
plot_match <- matching %>% filter(gene_name %in% intersect(rna_bae,meth_uu),!is.na(matching),!is.na(pattern)) %>% mutate(pattern=factor(pattern,levels=c("A","B","C","D")))

p_match <- ggplot(plot_match,aes(x=matching,y=pattern))+
  geom_jitter(aes(size=n_informative),height=0.15,width=0,alpha=0.7)+
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25))+
  labs(x="U/U–BAE matching frequency",y="Methylation pattern",size="No. informative\nU/U clones")+
  theme_classic(base_size=16)
p_match

pdf("output/extended_venn_matching_260919.pdf")
p_match
dev.off()


## GTEx XCI annotation
gtex <- read_tsv("input/Nature_GTEx_Tukiainene_2017.txt") %>% dplyr::rename(gene_name=`Gene name`,XCI_status=`Combined XCI status`) %>% dplyr::select(gene_name,XCI_status) %>% distinct()

# gene-level methylation pattern(A > B > C > D)과 GTEx XCI status matching
gtex_pattern <- pattern_info_gene %>% inner_join(gtex,by="gene_name") %>% filter(XCI_status %in% c("inactive","variable","escape")) %>% mutate(pattern=factor(pattern,levels=c("A","B","C","D")),XCI_status=factor(XCI_status,levels=c("inactive","variable","escape")))

# pattern별 GTEx XCI status proportion
plot_gtex <- gtex_pattern %>% count(pattern,XCI_status,.drop=F) %>% group_by(pattern) %>% mutate(prop=n/sum(n)) %>% ungroup()

# 실제 gene 수와 proportion 확인
gtex_pattern %>% count(pattern,XCI_status) %>% group_by(pattern) %>% mutate(prop=n/sum(n))

# Pattern별 GTEx XCI status
p_gtex <- ggplot(plot_gtex,aes(x=pattern,y=prop,fill=XCI_status))+
  geom_col(width=0.75)+
  scale_y_continuous(limits=c(0,1),labels=scales::percent_format(),expand=c(0,0))+
  labs(x="Promoter methylation pattern",y="Proportion of genes",fill="Reported XCI status")+
  theme_classic(base_size=16)+
  theme(axis.line=element_line(linewidth=0.8),axis.ticks=element_line(linewidth=0.8),axis.title=element_text(size=17),axis.text=element_text(size=14,color="black"),legend.position="right",legend.title=element_text(face="bold"))


pdf("output/extended_GTEx_260919.pdf")
p_gtex
dev.off()


## Supplementary Table 8 =================================================================

# published XCI status
gtex <- read_tsv("input/Nature_GTEx_Tukiainene_2017.txt") %>% dplyr::rename(gene_name=`Gene name`,XCI_status=`Combined XCI status`) %>% dplyr::select(gene_name,XCI_status) %>% distinct()

# RNA BAE 또는 methylation U/U 중 하나 이상에서 검출된 모든 gene
all_detected_genes <- union(rna_bae,meth_uu)

# 최종 supplementary table; BAE_freq는 matching과 동일하며 BAE_extent는 BAE clone의 mean minor VAF
supp_table <- tibble(gene_name=all_detected_genes) %>%
  mutate(RNA_BAE=gene_name %in% rna_bae,Methylation_UU=gene_name %in% meth_uu,detection=case_when(RNA_BAE & Methylation_UU~"Both",RNA_BAE~"RNA only",Methylation_UU~"Methylation only")) %>%
  left_join(pattern_info_gene,by="gene_name") %>%
  left_join(bae_freq %>% dplyr::select(gene_name,n_UU,n_BAE,n_MAE,n_LowDepth,n_wSNP,n_informative,BAE_freq),by="gene_name") %>%
  left_join(bae_extent %>% dplyr::select(gene_name,BAE_extent,n_BAE_extent_clones),by="gene_name") %>%
  left_join(gtex,by="gene_name") %>%
  dplyr::select(gene_name,detection,RNA_BAE,Methylation_UU,pattern,n_UU,n_BAE,n_MAE,n_LowDepth,n_wSNP,n_informative,BAE_freq,BAE_extent,n_BAE_extent_clones,XCI_status) %>%
  arrange(factor(detection,levels=c("Both","RNA only","Methylation only")),factor(pattern,levels=c("A","B","C","D")),gene_name)

# 반드시 gene당 1 row
supp_table %>% count(gene_name) %>% filter(n>1)
supp_table

# Build supplementary table
outdir <- "SUPPLEMENTARY_TABLE_EXTERNAL"

# DEG results
deg_fm <- read_tsv("input/chrX_deg_250827.tsv") %>%
  dplyr::rename(log2FC_FM=log2FoldChange) %>%
  dplyr::select(gene_name,baseMean,log2FC_FM,lfcSE,stat,pvalue,padj)

# Sheet 1: BAE-U/U methylation integration + DEG
supp8_summary <- supp_table %>%
  left_join(deg_fm,by="gene_name") %>%
  dplyr::select(gene_name,detection,RNA_BAE,Methylation_UU,pattern,n_UU,n_BAE,n_MAE,n_LowDepth,n_wSNP,n_informative,BAE_freq,BAE_extent,n_BAE_extent_clones,XCI_status,baseMean,log2FC_FM,lfcSE,stat,pvalue,padj) %>%
  dplyr::rename("Reported_status" = "XCI_status")

# Sheet 2: SNP-level ASE results
supp8_ase_raw <- ase_res_final

# Export
wb8 <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb8,"BAE_UU_integration")
openxlsx::writeData(wb8,"BAE_UU_integration",supp8_summary,withFilter=TRUE)
openxlsx::freezePane(wb8,"BAE_UU_integration",firstRow=TRUE)

openxlsx::addWorksheet(wb8,"ASE_variant_results")
openxlsx::writeData(wb8,"ASE_variant_results",supp8_ase_raw,withFilter=TRUE)
openxlsx::freezePane(wb8,"ASE_variant_results",firstRow=TRUE)

openxlsx::saveWorkbook(wb8,file.path(outdir,"supplementary_table_8.xlsx"),overwrite=TRUE)


################################################################################

# Xp vs Xq Escape

################################################################################

gg <- gg %>%
  dplyr::mutate(
    ESC_freq = ifelse(ASE_final == "MAE", NA, ESC_freq),
    ESC_extent = ifelse(ASE_final == "MAE", NA, ESC_extent)
  )
gg <- left_join(gg, cgi_total_summary_pattern %>% dplyr::select(promoter_id, cgi_id, unmethyl_freq) %>% distinct(promoter_id, .keep_all = T))
gg_tmp <- gg %>% dplyr::filter(position != "Centromere")

# Pattern별 Xp, Xq 비율
plot_df <- gg_tmp %>%
  dplyr::filter(pattern != "None") %>%
  dplyr::count(pattern, position) %>%
  dplyr::group_by(pattern) %>%
  dplyr::mutate(freq = n / sum(n))

ggplot(plot_df,
       aes(x = pattern, y = freq, fill = position)) +
  geom_col(width = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(
    Xp = "#4E79A7",
    Xq = "#E15759"
  )) +
  theme_classic() +
  labs(x = "Pattern", y = "Composition")

# Pattern 고려하지 않은버전
plot_df <- gg_tmp %>%
  tidyr::pivot_longer(
    cols=c(ESC_extent,ESC_freq,unmethyl_freq),
    names_to="Feature",
    values_to="Value"
  )

ggplot(plot_df,
       aes(position,Value,fill=position))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.15,size=1,alpha=0.6)+
  facet_wrap(~Feature,scales="free_y")+
  theme_classic()+
  theme(legend.position="none")

bar_df <- gg_tmp %>%
  dplyr::filter(pattern!="None") %>%
  dplyr::count(pattern, position) %>%
  dplyr::group_by(pattern) %>%
  dplyr::mutate(Value=n/sum(n),
                Feature="Pattern")

box_df <- gg_tmp %>%
  tidyr::pivot_longer(
    cols=c(ESC_extent,ESC_freq,unmethyl_freq),
    names_to="Feature",
    values_to="Value"
  )

p1 <- ggplot(box_df %>% dplyr::filter(Feature=="ESC_extent"),
             aes(position,Value,fill=position))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.15,size=1,alpha=0.6)+
  theme_classic()+theme(legend.position="none")

p2 <- ggplot(box_df %>% dplyr::filter(Feature=="ESC_freq"),
             aes(position,Value,fill=position))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.15,size=1,alpha=0.6)+
  theme_classic()+theme(legend.position="none")

p3 <- ggplot(box_df %>% dplyr::filter(Feature=="unmethyl_freq"),
             aes(position,Value,fill=position))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.15,size=1,alpha=0.6)+
  theme_classic()+theme(legend.position="none")

p4 <- ggplot(bar_df,
             aes(pattern,Value,fill=position))+
  geom_col(width=0.7)+
  scale_y_continuous(labels=scales::percent)+
  theme_classic()

pdf("output/fig4_escape_Xpq_260808.pdf", width = 12, height = 25)
(p1|p2|p3|p4)
dev.off()



################################################################################
