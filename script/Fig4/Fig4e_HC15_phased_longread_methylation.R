################################################################################
# Figure 4e. Haplotype-specific CpG island methylation in HC15
################################################################################

# 1. Libraries
library(dplyr)
library(readr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# 2. Input data
input_dir <- "input/Fig4/Phased_HC15"

ss <- c("HC15-17","HC15-10","HC15-22","HC15-05","HC15-16","HC15-09",
        "HC15-21","HC15-14","HC15-20","HC15-04","HC15-23","HC15-19",
        "HC15-13","HC15-01","HC15-12","HC15-11","HC15-08","HC15-06",
        "HC15-18","HC15-15","HC15-02","HC15-07","HC15-03")

cgi_promoter_total <- read_tsv("input/Fig4/cgi_promoter_total.tsv",show_col_types=FALSE)

# 3. Read haplotype-specific methylation
read_phased <- function(s){
  hp1 <- read_tsv(file.path(input_dir,paste0(s,"_hap1.bismark.cov.cpgi.txt")),
                  col_names=c("chrom","pos","pos2","pct","me_count","unme_count"),show_col_types=FALSE) %>%
    filter(me_count+unme_count>0) %>%
    select(chrom,pos,pos2,pct) %>%
    rename(!!paste0(s,"_hp1_pct"):=pct)
  
  hp2 <- read_tsv(file.path(input_dir,paste0(s,"_hap2.bismark.cov.cpgi.txt")),
                  col_names=c("chrom","pos","pos2","pct","me_count","unme_count"),show_col_types=FALSE) %>%
    filter(me_count+unme_count>0) %>%
    select(chrom,pos,pos2,pct) %>%
    rename(!!paste0(s,"_hp2_pct"):=pct)
  
  full_join(hp1,hp2,by=c("chrom","pos","pos2")) %>% drop_na()
}

# 4. Merge HC15 clones
hp_list <- lapply(ss,read_phased)
hp <- Reduce(function(x,y) full_join(x,y,by=c("chrom","pos","pos2")),hp_list)

hp <- hp %>%
  left_join(cgi_promoter_total %>% select(chrom,pos,pos2,gene_name,promoter_id),
            by=c("chrom","pos","pos2"))

# 5. Prepare heatmap matrix
ht <- hp %>%
  filter(!is.na(promoter_id),promoter_id!=".") %>%
  arrange(pos) %>%
  select(contains("_pct"),promoter_id,gene_name)

ht[is.na(ht)] <- -200
hp2_cols <- grep("_hp2_pct$",colnames(ht),value=TRUE)
ht[hp2_cols] <- lapply(ht[hp2_cols],function(x) -abs(x))

meth_cols <- grep("_pct$",colnames(ht),value=TRUE)
clone_group <- sub("_hp[12]_pct$","",meth_cols)

# 6. Heatmap
hm <- Heatmap(t(as.matrix(ht[,meth_cols,drop=FALSE])),
              col=colorRamp2(c(-200,-100,0,100),c("#dddddd","#07406b","#ffffff","#de3174")),
              column_split=factor(ht$gene_name,levels=unique(ht$gene_name)),
              row_split=factor(clone_group,levels=unique(clone_group)),
              cluster_rows=FALSE,cluster_columns=FALSE,
              column_title_gp=gpar(fontsize=8),column_title_rot=90,row_title_rot=0,
              heatmap_legend_param=list(title="Methylation(%)"),
              show_row_names=TRUE)

# 7. Visualization
draw(hm)