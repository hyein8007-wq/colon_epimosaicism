################################################################################
# Figure 4b. X-chromosome epigenotype transitions between males and females
################################################################################

# 1. Libraries
library(dplyr)
library(readr)
library(ggplot2)
library(ggalluvial)

# 2. Load data
cpg_x <- read_tsv("input/Fig4/class_converted_MF_chrX_260917.tsv.gz",show_col_types=FALSE)

# 3. Classify PAR and non-PAR regions
cpg_x <- cpg_x %>%
  mutate(par=ifelse(between(pos,60001,2699520) | between(pos,154931044,155260560),"PAR","nonPAR")) %>%
  filter(!(par=="nonPAR" & male_WT==0.5))

# 4. Summarize epigenotype transitions
df_transition <- cpg_x %>%
  count(par,male_WT,female_WT,name="Freq") %>%
  mutate(male_WT=factor(male_WT,levels=c(1,0.5,0)),female_WT=factor(female_WT,levels=c(1,0.5,0))) %>%
  group_by(par) %>% mutate(prop=Freq/sum(Freq),label=sprintf("%.1f%%",prop*100)) %>% ungroup()

# 5. Sankey plots
epi_colors <- c("0"="#003e73","0.5"="#f2d8aa","1"="#a51f27")

plot_sankey <- function(region){
  ggplot(filter(df_transition,par==region),aes(axis1=male_WT,axis2=female_WT,y=Freq)) +
    geom_alluvium(aes(fill=female_WT),width=1/6) +
    geom_stratum(width=1/5,fill="grey85",color="black") +
    geom_text(stat="alluvium",aes(label=label),size=3,color="black") +
    scale_fill_manual(values=epi_colors) +
    scale_x_discrete(limits=c("Male WT","Female WT"),expand=c(0.1,0.1)) +
    labs(x=NULL,y="Count",fill="Female WT") +
    theme_minimal() +
    theme(panel.grid=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank())
}

p_par <- plot_sankey("PAR")
p_nonpar <- plot_sankey("nonPAR")

print(p_par)
print(p_nonpar)