################################################################################
# Figure 4c. Genomic distribution of X-chromosome epigenotypes
################################################################################

# 1. Libraries
library(dplyr)
library(readr)
library(ggplot2)

# 2. Load data
cpg_x <- read_tsv("input/Fig4/class_converted_MF_chrX_260917.tsv.gz",show_col_types=FALSE)

# 3. Classify epigenotypes
cpg_x <- cpg_x %>%
  filter(male_WT!=0.5) %>%
  mutate(type=factor(paste0(male_WT,"_",female_WT),levels=c("0_0","0_0.5","0_1","1_0","1_0.5","1_1")))

# 4. Summarize genomic annotations
df_region <- cpg_x %>%
  count(type,genomic_region,name="Freq") %>%
  group_by(type) %>% mutate(prop=Freq/sum(Freq)) %>% ungroup() %>%
  mutate(genomic_region=factor(genomic_region,levels=c("intergenic","gene_body","promoter")))

df_context <- cpg_x %>%
  count(type,cpg_context,name="Freq") %>%
  group_by(type) %>% mutate(prop=Freq/sum(Freq)) %>% ungroup() %>%
  mutate(cpg_context=factor(cpg_context,levels=c("Island","Shore","Shelf","OpenSea")))

# 5. Visualization
g1 <- ggplot(df_region,aes(type,prop,fill=genomic_region)) +
  geom_col() +
  scale_fill_manual(values=c("promoter"="#543a91","gene_body"="#ab98d9","intergenic"="#e0e0e0")) +
  labs(x="Type",y="Proportion",fill="Genomic region") +
  coord_flip() + theme_classic()

g2 <- ggplot(df_context,aes(type,prop,fill=cpg_context)) +
  geom_col() +
  scale_fill_manual(values=c("Island"="#360c2a","Shore"="#986b8c","Shelf"="#cfbecb","OpenSea"="#e9e7e6")) +
  labs(x="Type",y="Proportion",fill="CpG context") +
  coord_flip() + theme_classic()

print(g1)
print(g2)