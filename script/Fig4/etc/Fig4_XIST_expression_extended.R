# Archived source excerpt. Use Fig4_XIST_expression_revised.R instead.
# FIGURE 4 - SOURCE-PRESERVING EXTRACT (NOT YET STANDALONE)
# Source: original_XCI_methylation.R; original line ranges: 582-701
# This extract preserves the original calculations and figure commands.
# It depends on objects from earlier sections and original server inputs.
# DO NOT publish as independently reproducible until dependencies are resolved.


# ===== ORIGINAL LINES 582-701 =====
###############################################################################

# + XIST expression & CGI methylation level

###############################################################################


## XIST expression ============================================================

library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ggpubr)
library(tidyverse)
library(data.table)

# 1. Collect raw RSEM XIST TPM
bam_dir <- "/home/users/team_projects/hcolon_heterogeneity/01_LINE1/02_RNA/02_bam"
rna_clones <- read_tsv("input/final_rna_sample_list.tsv") %>% .$sample_id
male_donors <- c("HC01","HC02","HC09","HC10","HC13")
xist_id <- "ENSG00000229807"

hc_dirs <- list.dirs(bam_dir,recursive=FALSE,full.names=TRUE)
hc_dirs <- hc_dirs[grepl("^HC",basename(hc_dirs))]
files <- unlist(lapply(hc_dirs,function(d) list.files(file.path(d,"rsem"),pattern="\\.genes\\.results$",full.names=TRUE)))
files <- files[which(!grepl("HC19|HC22", files))] %>% .[which(!grepl("_P", .))] %>% .[which(!grepl("-B", .))]

xist_df <- lapply(files,function(f){
  x <- fread(f,select=c("gene_id","TPM"))
  x <- x[gene_id==xist_id]
  if(nrow(x)==0) return(NULL)
  data.frame(sample=sub("\\.genes\\.results$","",basename(f)),XIST_TPM=x$TPM)
}) %>% bind_rows() %>%
  filter(sample %in% rna_clones) %>%
  mutate(donor=ifelse(grepl("^HC",sample),sub("-.*","",sample),"HC50"),
         sex=ifelse(donor %in% male_donors,"Male","Female"),
         sex=factor(sex,levels=c("Male","Female")),
         log2_XIST=log2(XIST_TPM+0.01))

# 2. QC
cat("RSEM samples:",nrow(xist_df),"\n")
cat("Final RNA clones:",length(rna_clones),"\n")
cat("Missing:",length(setdiff(rna_clones,xist_df$sample)),"\n")
print(setdiff(rna_clones,xist_df$sample))
print(table(xist_df$donor,xist_df$sex))

# 3. Male vs Female
g_xist <-ggplot(xist_df,aes(x=sex,y=log2_XIST,fill=sex)) +
  geom_boxplot(width=0.5,outlier.shape=NA,alpha=0.7) +
  geom_jitter(width=0.15,size=1.5,alpha=0.6) +
  ggpubr::stat_compare_means(method="wilcox.test",label="p.format") +
  labs(x=NULL,y=expression(log[2]*"(XIST TPM + 0.01)")) +
  theme_classic() +
  theme(legend.position="none")

# ggplot(xist_df,aes(x=sample,y=log2_XIST,fill=sex)) +
#  geom_boxplot(width=0.5,outlier.shape=NA,alpha=0.7) +
#  geom_jitter(width=0.15,size=1.5,alpha=0.6) +
#  facet_grid(~donor) +
#  #ggpubr::stat_compare_means(method="wilcox.test",label="p.format") +
#  labs(x=NULL,y=expression(log[2]*"(XIST TPM + 1)")) +
#  theme_classic() +
#  theme(legend.position="none")


## CpG island methylation ======================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Non-PAR X CpGs
epi_x <- cpg_x %>%
  dplyr::filter(par!="PAR",male_WT!=0.5,cgi_id %in% cgi_total_summary_pattern$cgi_id) %>%
  dplyr::select(contains("-"))

male_donors <- c("HC01","HC02","HC09","HC10","HC13","HC19")

# 2. Proportion calculation
epi_prop <- epi_x %>%
  dplyr::mutate(cpg_id=dplyr::row_number()) %>%
  tidyr::pivot_longer(-cpg_id,names_to="sample",values_to="epigenotype") %>%
  dplyr::filter(!is.na(epigenotype),epigenotype %in% c(0,0.5,1)) %>%
  dplyr::mutate(donor=ifelse(grepl("^HC",sample),sub("-.*","",sample),"HC50"),
                sex=ifelse(donor %in% male_donors,"Male","Female"),
                sex=factor(sex,levels=c("Male","Female")),
                epigenotype=factor(epigenotype,levels=c(0,0.5,1),labels=c("0","0.5","1"))) %>%
  dplyr::count(sample,donor,sex,epigenotype,name="n") %>%
  dplyr::group_by(sample,donor,sex) %>%
  dplyr::mutate(prop=n/sum(n)) %>%
  dplyr::ungroup()

# 3. Sample order
sample_order <- epi_prop %>% dplyr::distinct(sample,donor,sex) %>% dplyr::arrange(sex,donor,sample) %>% dplyr::pull(sample)
donor_order <- epi_prop %>% dplyr::distinct(donor,sex) %>% dplyr::arrange(sex,donor) %>% dplyr::pull(donor)
epi_prop <- epi_prop %>% dplyr::mutate(sample=factor(sample,levels=sample_order),donor=factor(donor,levels=donor_order))

# 4. Figure
epi_prop$epigenotype <- factor(epi_prop$epigenotype, levels = c(1, 0.5, 0)) 
epi_prop$donor <- factor(epi_prop$donor, levels = c("HC01", "HC02", "HC09", "HC10", "HC13","HC19", "HC08", "HC12", "HC14", "HC15", "HC16", "HC17", "HC18","HC20", "HC21", "HC50")) 
g_x_me <- ggplot(epi_prop,aes(sample,prop,fill=epigenotype)) +
  geom_col(width=1) +
  scale_fill_manual(values = c(
    "0" = "#003e73",
    "0.5" = "#f2d8aa",
    "1" = "#a51f27"
  )) +
  facet_grid(~donor,scales="free_x",space="free_x") +
  scale_y_continuous(labels=scales::percent_format()) +
  labs(x=NULL,y="Proportion of CpGs",fill="Epigenotype") +
  theme_classic() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        strip.background=element_blank(),strip.text=element_text(face="bold"),
        panel.spacing.x=unit(0.15,"lines"),legend.position="top")

pdf("output/Extended_fig_XISTexp_CGIepigenotype.pdf")
print(g_xist)
print(g_x_me)
dev.off()
