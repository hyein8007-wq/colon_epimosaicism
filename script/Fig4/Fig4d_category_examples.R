################################################################################
# Figure 4. CpG island methylation patterns across female clones
################################################################################

# 1. Libraries
library(dplyr)
library(readr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

# 2. Load data
cpg_x <- read_tsv("input/Fig4/class_converted_MF_chrX_260917.tsv.gz",show_col_types=FALSE)
cgi_pattern <- read_tsv("input/Fig4/unmethyl_x_cgis_pattern_manual_v3_260918.tsv",show_col_types=FALSE)
clone_info <- read_tsv("input/Fig4/clone_info_260917.tsv",show_col_types=FALSE)

# 3. Select U/U CpG islands and female clones
sample_f <- c("HC08","HC12","HC14","HC15","HC16","HC17","HC18","HC20","HC21","HC50")
cgis <- cpg_x %>% filter(unmethyl_freq>0) %>% pull(cgi_id) %>% unique()
cgis <- cgis[!is.na(cgis)]

female_clones <- clone_info %>%
  filter(WGS=="Y",EM_seq=="Y",sampleID %in% colnames(cpg_x)) %>%
  filter(sub("-.*","",sampleID) %in% sample_f | grepl("^D[0-9]+-|^PS-",sampleID)) %>%
  pull(sampleID)

# 4. Prepare CpG methylation matrix
hm_input <- cpg_x %>%
  filter(cgi_id %in% cgis) %>%
  arrange(pos) %>%
  select(cgi_id,pos,all_of(female_clones))

meth_mat <- as.matrix(hm_input[,female_clones,drop=FALSE])
storage.mode(meth_mat) <- "numeric"
meth_mat[is.na(meth_mat)] <- -1

# 5. CpG island methylation heatmaps
meth_colors <- colorRamp2(c(-1,0,0.5,1),c("#cccccc","#003e73","#f2d8aa","#a51f27"))

plot_cgi <- function(id){
  idx <- which(hm_input$cgi_id==id)
  if(!length(idx)) return(NULL)
  
  mat <- t(meth_mat[idx,,drop=FALSE])
  donor <- sub("-.*","",rownames(mat))
  donor[!donor %in% sample_f] <- "HC50"
  
  Heatmap(mat,name="Methylation",col=meth_colors,
          cluster_rows=FALSE,cluster_columns=FALSE,
          show_row_names=FALSE,show_column_names=FALSE,
          row_split=factor(donor,levels=sample_f),
          column_title=id,column_title_gp=gpar(fontsize=12),
          row_title_rot=0,row_title_gp=gpar(fontsize=8),
          use_raster=TRUE,raster_quality=4,
          width=unit(200,"mm"),height=unit(90,"mm"))
}

# 6. Visualize a representative CpG island
cgi_id <- cgis[1]
draw(plot_cgi(cgi_id))