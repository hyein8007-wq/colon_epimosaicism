################################################################################
# Extended Data. Integrated analysis of XCI escape and promoter hypomethylation
################################################################################

# 1. Libraries
library(dplyr)
library(readr)
library(ggplot2)
library(ggVennDiagram)
library(ComplexHeatmap)
library(circlize)
library(grid)

# 2. Load reference annotations
input_dir <- "input/ExtendedData"

gtf <- read_tsv(file.path(input_dir,"Homo_sapiens.GRCh37.87.gtf.gz"),comment="#",col_names=FALSE,show_col_types=FALSE) %>%
  filter(X3=="gene",X1=="X") %>%
  transmute(chr=X1,pos1=X4,pos2=X5,
            gene_id=sub('.*gene_id "([^"]+)".*','\\1',X9),
            gene_name=sub('.*gene_name "([^"]+)".*','\\1',X9)) %>%
  filter(grepl('gene_name "',gene_name,fixed=TRUE)==FALSE,
         !(pos1<=2699520 & pos2>=60001),
         !(pos1<=155260560 & pos2>=154931044)) %>%
  distinct(gene_name,.keep_all=TRUE)

gtex <- read_tsv(file.path(input_dir,"Nature_GTEx_Tukiainene_2017.txt"),show_col_types=FALSE) %>%
  rename(gene_name=`Gene name`,XCI_status=`Combined XCI status`) %>%
  select(gene_name,XCI_status) %>% distinct()

nonpar_genes <- gtf$gene_name

# 3. Identify genes with RNA BAE, promoter U/U or reported XCI escape
rna_bae <- ase_gene_clone %>%
  filter(final_group=="BAE",gene_name %in% nonpar_genes) %>%
  distinct(gene_name) %>% pull(gene_name)

meth_uu <- uu_gene_clone %>%
  filter(gene_name %in% nonpar_genes) %>%
  distinct(gene_name) %>% pull(gene_name)

reported_escape <- gtex %>%
  filter(XCI_status %in% c("escape","variable"),gene_name %in% nonpar_genes) %>%
  distinct(gene_name) %>% pull(gene_name)

venn_list <- list("RNA BAE"=rna_bae,"Promoter U/U"=meth_uu,
                  "Reported escape/variable"=reported_escape)
all_genes <- Reduce(union,venn_list)

cat("RNA BAE:",length(rna_bae),"\n")
cat("Promoter U/U:",length(meth_uu),"\n")
cat("Reported escape/variable:",length(reported_escape),"\n")

# 4. Venn diagram
p_venn <- ggVennDiagram(venn_list,label_alpha=0,edge_size=0.8,label="count") +
  scale_fill_gradient(low="white",high="grey70") +
  theme_void(base_size=16) + theme(legend.position="none")

print(p_venn)

# 5. Summarize RNA ASE
rna_summary <- ase_gene_clone %>%
  group_by(gene_name) %>%
  summarise(ASE=case_when(any(final_group=="BAE")~"BAE",
                          any(final_group=="MAE")~"MAE",
                          any(final_group=="LowDepth")~"LowDepth",TRUE~"w/oSNP"),
            n_BAE=sum(final_group=="BAE"),n_MAE=sum(final_group=="MAE"),
            n_informative=n_BAE+n_MAE,
            BAE_freq=ifelse(n_informative>0,n_BAE/n_informative,NA_real_),.groups="drop")

rna_extent <- ase_res_final %>%
  filter(sample %in% common_clones,final_group=="BAE",!is.na(VAF),between(VAF,0,1)) %>%
  mutate(minor_VAF=pmin(VAF,1-VAF)) %>%
  group_by(gene_name,sample) %>%
  summarise(clone_extent=mean(minor_VAF),.groups="drop") %>%
  group_by(gene_name) %>%
  summarise(BAE_extent=mean(clone_extent),.groups="drop")

rna_summary <- rna_summary %>%
  left_join(rna_extent,by="gene_name") %>%
  mutate(BAE_extent=case_when(n_informative==0~NA_real_,n_BAE==0~0,TRUE~BAE_extent))

# 6. Summarize promoter U/U methylation and matched RNA ASE
uu_freq <- uu_gene_clone %>%
  count(gene_name,name="n_UU") %>%
  mutate(UU_freq=n_UU/278)

uu_metrics <- bae_freq %>%
  left_join(bae_extent,by="gene_name") %>%
  mutate(BAE_extent=case_when(n_informative==0~NA_real_,n_BAE==0~0,TRUE~BAE_extent))

# 7. Order genes by X-chromosome position
gene_info <- gtf %>%
  filter(gene_name %in% all_genes) %>%
  mutate(gene_mid=(pos1+pos2)/2,arm=ifelse(gene_mid<60600000,"Xp","Xq")) %>%
  arrange(pos1)

gene_order_p <- gene_info %>% filter(arm=="Xp") %>% pull(gene_name)
gene_order_q <- gene_info %>% filter(arm=="Xq") %>% pull(gene_name)

cat("Venn union:",length(all_genes),"\n")
cat("Heatmap genes:",nrow(gene_info),"\n")
cat("Xp:",length(gene_order_p),"\n")
cat("Xq:",length(gene_order_q),"\n")

# 8. Heatmap colors
ase_col <- c("BAE"="#A99BC2","MAE"="#493277","LowDepth"="#C7C7C7","w/oSNP"="#F2F2F2")
pattern_col <- c("A"="#b82860","B"="#e4a126","C"="#284054","D"="#54704c")
xci_col <- c("inactive"="#253e51","variable"="#e4a126","escape"="#b82860")

# 9. Heatmap functions
get_values <- function(df,genes,column) df[[column]][match(genes,df$gene_name)]

make_bar_heatmap <- function(values,label,max_value=1){
  mat <- matrix(values,nrow=1,dimnames=list(label,names(values)))
  Heatmap(mat,name=label,col=colorRamp2(c(0,max_value),c("white","white")),
          na_col="white",cluster_rows=FALSE,cluster_columns=FALSE,
          show_column_names=FALSE,show_heatmap_legend=FALSE,height=unit(10,"mm"),
          cell_fun=function(j,i,x,y,w,h,fill){
            grid.rect(x=x,y=unit(0,"npc"),width=w,height=unit(1,"npc"),
                      just=c("center","bottom"),gp=gpar(fill="white",col="grey90",lwd=0.2))
            if(!is.na(values[j])) grid.rect(x=x,y=unit(0,"npc"),width=w*0.75,
                                            height=unit(values[j]/max_value,"npc"),just=c("center","bottom"),
                                            gp=gpar(fill="grey35",col=NA))
          })
}

make_xci_heatmap <- function(genes,arm){
  rna_ase <- get_values(rna_summary,genes,"ASE")
  rna_freq <- get_values(rna_summary,genes,"BAE_freq")
  rna_extent_v <- get_values(rna_summary,genes,"BAE_extent")
  meth_freq <- get_values(uu_freq,genes,"UU_freq")
  meth_freq[is.na(meth_freq)] <- 0
  meth_bae_freq <- get_values(uu_metrics,genes,"BAE_freq")
  meth_bae_extent <- get_values(uu_metrics,genes,"BAE_extent")
  pattern <- get_values(pattern_info_gene,genes,"pattern")
  pattern[!genes %in% meth_uu] <- NA
  reported <- get_values(gtex,genes,"XCI_status")
  reported[!reported %in% c("inactive","variable","escape")] <- NA
  
  categorical_heatmap <- function(values,label,colors,show_names=FALSE){
    mat <- matrix(values,nrow=1,dimnames=list(label,genes))
    Heatmap(mat,name=paste0(label,"_",arm),col=colors,na_col="white",
            cluster_rows=FALSE,cluster_columns=FALSE,
            show_column_names=show_names,column_names_rot=90,
            column_names_gp=gpar(fontsize=5),column_names_centered=TRUE,
            height=unit(4,"mm"))
  }
  
  names(rna_freq) <- names(rna_extent_v) <- names(meth_freq) <- genes
  names(meth_bae_freq) <- names(meth_bae_extent) <- genes
  
  ht1 <- categorical_heatmap(rna_ase,"ASE",ase_col)
  ht2 <- make_bar_heatmap(rna_freq,"RNA BAE frequency")
  ht3 <- make_bar_heatmap(rna_extent_v,"RNA BAE extent",0.5)
  ht4 <- categorical_heatmap(pattern,"Pattern",pattern_col)
  ht5 <- make_bar_heatmap(meth_freq,"U/U frequency")
  ht6 <- make_bar_heatmap(meth_bae_freq,"BAE frequency (U/U)")
  ht7 <- make_bar_heatmap(meth_bae_extent,"BAE extent (U/U)",0.5)
  ht8 <- categorical_heatmap(reported,"Reported XCI",xci_col,TRUE)
  
  ht1 %v% ht2 %v% ht3 %v% ht4 %v% ht5 %v% ht6 %v% ht7 %v% ht8
}

# 10. Visualize Xp and Xq
ht_p <- make_xci_heatmap(gene_order_p,"Xp")
ht_q <- make_xci_heatmap(gene_order_q,"Xq")

grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1)))

pushViewport(viewport(layout.pos.row=1))
draw(ht_p,merge_legends=TRUE,heatmap_legend_side="right",newpage=FALSE)
popViewport()

pushViewport(viewport(layout.pos.row=2))
draw(ht_q,merge_legends=TRUE,heatmap_legend_side="right",newpage=FALSE)
popViewport(2)