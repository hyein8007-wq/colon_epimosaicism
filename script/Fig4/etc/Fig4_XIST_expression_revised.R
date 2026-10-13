################################################################################
# XIST expression: uncorrected RSEM TPM, 374 clones (2026-10-13)
# Run from github/ repository root. Input file is linked under input/.
################################################################################
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Load uncorrected TPM (do NOT use exp_bc_TPM for this analysis)
expr <- fread("input/exp_RSEM_TPM_374clones_261013.tsv")
stopifnot("gene_name" %in% names(expr),sum(expr$gene_name=="XIST",na.rm=TRUE)==1)
rna_info <- fread("input/final_rna_sample_list.tsv")
stopifnot("sample_id" %in% names(rna_info))
rna_clones <- unique(rna_info$sample_id)

# 2. Extract clone-level XIST TPM
xist <- expr[gene_name=="XIST",] %>% select(-gene_name) %>%
  pivot_longer(everything(),names_to="sample",values_to="XIST_TPM") %>%
  filter(sample %in% rna_clones,!grepl("^HC19|^HC22|_P|-B",sample)) %>%
  mutate(donor=ifelse(grepl("^HC",sample),sub("-.*","",sample),"HC50"),
         sex=factor(ifelse(donor %in% c("HC01","HC02","HC09","HC10","HC13"),"Male","Female"),levels=c("Male","Female")),
         log2_XIST=log2(as.numeric(XIST_TPM)+0.01))
stopifnot(!anyDuplicated(xist$sample),!anyNA(xist$XIST_TPM))
print(xist %>% group_by(sex) %>% summarise(n=n(),median_TPM=median(XIST_TPM),.groups="drop"))
print(wilcox.test(log2_XIST~sex,data=xist))

# 3. Plot
dir.create("output",showWarnings=FALSE)
p <- ggplot(xist,aes(sex,log2_XIST,fill=sex)) +
  geom_boxplot(width=0.5,outlier.shape=NA,alpha=0.7) +
  geom_jitter(width=0.15,size=1.5,alpha=0.6) +
  labs(x=NULL,y="log2(XIST RSEM TPM + 0.01)") +
  theme_classic() + theme(legend.position="none")
ggsave("output/XIST_RSEM_TPM.pdf",p,width=4,height=4)
