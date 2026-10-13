################################################################################
# Figure 4d. X-chromosome CpG island methylation landscape
################################################################################

# 1. Libraries and input
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(data.table)
library(ape)
library(ggtree)
library(ComplexHeatmap)
library(circlize)
library(grid)

input_dir <- "input/Fig4"
cpg_x <- read_tsv(file.path(input_dir,"class_converted_MF_chrX_260917.tsv.gz"),show_col_types=FALSE)
cgi_summary <- read_tsv(file.path(input_dir,"chrX_cgi_total_summary_pattern.tsv"),show_col_types=FALSE)
ase <- read_tsv(file.path(input_dir,"ase_res_female_chrX_final_260919_v2.tsv"),show_col_types=FALSE)
rna_samples <- read_tsv(file.path(input_dir,"final_rna_sample_list.tsv"),show_col_types=FALSE)
anno_info <- read_tsv(file.path(input_dir,"hm_anno_info.tsv"),show_col_types=FALSE)
allele_df <- read_tsv(file.path(input_dir,"allele_info_260918.tsv"),show_col_types=FALSE) %>% rename(sample_id=clone_id)
clone_info <- read_tsv(file.path(input_dir,"clone_info_260917.tsv"),show_col_types=FALSE)

# 2. Clone ordering from the original analysis
sample_m <- c("HC01","HC02","HC09","HC10","HC13","HC19")
sample_f <- c("HC08","HC12","HC14","HC15","HC16","HC17","HC18","HC20","HC21","HC50")
clone_ids <- clone_info %>% filter(WGS=="Y",EM_seq=="Y") %>% pull(sampleID)
orders <- list("HC01" = c("HC01-21", "HC01-15", "HC01-12", "HC01-09", "HC01-13", "HC01-02", "HC01-18", "HC01-20", "HC01-11", "HC01-06",
                          "HC01-17", "HC01-19", "HC01-04", "HC01-05", "HC01-10", "HC01-07", "HC01-03", "HC01-08", "HC01-16", "HC01-14"),
               "HC02" = c("HC02-04", "HC02-01", "HC02-08", "HC02-14", "HC02-02", "HC02-12", "HC02-11", "HC02-17", "HC02-10", "HC02-07", 
                          "HC02-21", "HC02-15", "HC02-22", "HC02-03", "HC02-16", "HC02-13", "HC02-20", "HC02-05", "HC02-19", "HC02-09", "HC02-06", "HC02-18"),
               "HC09" = c("HC09-10", "HC09-08", "HC09-02", "HC09-01", "HC09-20", "HC09-17", "HC09-16", "HC09-11", "HC09-13", "HC09-15", 
                          "HC09-19", "HC09-04", "HC09-06", "HC09-21", "HC09-05", "HC09-03", "HC09-14", "HC09-07", "HC09-18", "HC09-09", "HC09-12"),
               "HC10" = c("HC10-19", "HC10-05", "HC10-07", "HC10-11", "HC10-10", "HC10-17", "HC10-14", "HC10-15", "HC10-06", "HC10-16" , 
                          "HC10-03", "HC10-09", "HC10-18", "HC10-04", "HC10-13", "HC10-02", "HC10-08", "HC10-20", "HC10-01", "HC10-12"),
               "HC13" = c("HC13-19", "HC13-17", "HC13-14", "HC13-21", "HC13-09", "HC13-07", "HC13-05", "HC13-12", "HC13-08", "HC13-01", 
                          "HC13-16", "HC13-03", "HC13-22", "HC13-20", "HC13-02", "HC13-11", "HC13-18", "HC13-06", "HC13-13", "HC13-10", "HC13-15", "HC13-04"),
               "HC08" = c("HC08-15", "HC08-08", "HC08-21", "HC08-18", "HC08-20", "HC08-04", "HC08-17", "HC08-13", "HC08-23", "HC08-12", 
                          "HC08-09", "HC08-14", "HC08-07", "HC08-03", "HC08-19", "HC08-10", "HC08-22", "HC08-01", "HC08-11", "HC08-06", "HC08-16","HC08-05"),
               "HC12" = c("HC12-18", "HC12-04", "HC12-11", "HC12-23","HC12-07", "HC12-08", "HC12-22", "HC12-09", "HC12-17", "HC12-10", "HC12-05", #07앞 "HC12-23"/ 05앞 10/20 앞 24 
                          "HC12-21", "HC12-14", "HC12-15", "HC12-03", "HC12-19", "HC12-12", "HC12-16", "HC12-24", "HC12-20"),
               "HC14" = c("HC14-14", "HC14-13", "HC14-05", "HC14-01", "HC14-03", "HC14-04", "HC14-11", "HC14-12", "HC14-17", "HC14-16"),
               "HC15" = c("HC15-17", "HC15-10", "HC15-22", "HC15-05", "HC15-16", "HC15-09", "HC15-21", "HC15-14", "HC15-20", "HC15-04", 
                          "HC15-23", "HC15-19", "HC15-13", "HC15-01", "HC15-12", "HC15-11", "HC15-08", "HC15-06", "HC15-18", "HC15-15", "HC15-02", "HC15-07", "HC15-03"),
               "HC16" = c("HC16-08", "HC16-01", "HC16-19", "HC16-05", "HC16-17", "HC16-18", "HC16-12", "HC16-14", "HC16-04", "HC16-10", "HC16-11",
                          "HC16-09", "HC16-15", "HC16-13", "HC16-02", "HC16-07", "HC16-06", "HC16-16", "HC16-03"),
               "HC17" = c("HC17-16", "HC17-07", "HC17-12", "HC17-03", "HC17-17", "HC17-02", "HC17-01", "HC17-18", "HC17-08", "HC17-15", 
                          "HC17-05", "HC17-19", "HC17-13", "HC17-04"), #231031: HC17-17샘플 없는것 확인해서 다시 집어넣음. #241107: HC17-02 앞에 있던 HC17-17 제거. overall escaping
               "HC18" = c("HC18-20", "HC18-13", "HC18-05", "HC18-09", "HC18-17", "HC18-10", "HC18-07", "HC18-16", "HC18-19", #"HC18-01", 
                          "HC18-11", "HC18-15", "HC18-12", "HC18-18", "HC18-03", "HC18-14"),
               "HC19" = clone_ids[grepl("HC19", clone_ids)],
               "HC20" = c("HC20-01", "HC20-03", "HC20-06", "HC20-16", "HC20-05", "HC20-04", "HC20-14", "HC20-10", "HC20-07", "HC20-18", 
                          "HC20-17", "HC20-08", "HC20-22", "HC20-15", "HC20-09", "HC20-02", "HC20-13", "HC20-20", "HC20-19", "HC20-21", "HC20-12", "HC20-11"),
               "HC21" = c("HC21-07", "HC21-02", "HC21-11", "HC21-18", "HC21-14", "HC21-01", "HC21-21", "HC21-19", "HC21-06", "HC21-20", 
                          "HC21-15", "HC21-10", "HC21-23", "HC21-16", "HC21-05", "HC21-22", "HC21-12", "HC21-17", "HC21-09", "HC21-13", "HC21-04", "HC21-08", "HC21-03"),
               "HC50" = read.newick("input/Fig4/HC50_colon_final_250703.length.nwk") %>% ggtree() %>% .$data %>% dplyr::filter(isTip) %>% arrange(y) %>% pull(label)
)
male_clones <- unlist(orders[sample_m],use.names=FALSE)
female_clones <- unlist(orders[sample_f],use.names=FALSE)
stopifnot(!anyDuplicated(c(male_clones,female_clones)),all(c(male_clones,female_clones) %in% names(cpg_x)))

# 3. Reconstruct CpG island IDs and restrict to non-PAR X chromosome
cgi_bed <- read_tsv(file.path(input_dir,"cpgIslandExt.sort.nochr.range.bed"),col_names=FALSE,show_col_types=FALSE) %>%
  transmute(chr=as.character(X1),start=X2,end=X3) %>% filter(chr=="X") %>% mutate(cgi_id=paste0("cgi_",row_number()))
if(!"cgi_id" %in% names(cpg_x)){
  cpg_dt <- as.data.table(cpg_x) %>% copy()
  cpg_dt[,`:=`(chr="X",pos_start=pos,pos_end=pos)]
  bed_dt <- as.data.table(cgi_bed)
  setkey(bed_dt,chr,start,end); setkey(cpg_dt,chr,pos_start,pos_end)
  overlap <- foverlaps(cpg_dt,bed_dt,by.x=c("chr","pos_start","pos_end"),by.y=c("chr","start","end"),type="within",nomatch=NA)
  stopifnot(nrow(overlap)==nrow(cpg_dt))
  cpg_dt[,cgi_id:=overlap$cgi_id]
  cpg_dt[,c("chr","pos_start","pos_end"):=NULL]
  cpg_x <- as_tibble(cpg_dt)
}
cpg_x <- cpg_x %>%
  mutate(par=ifelse(between(pos,60001,2699520)|between(pos,154931044,155260560),"PAR","nonPAR"),
         position=case_when(gene_name=="XIST"~"XIST",pos<58100000~"Xp",pos>63000000~"Xq",TRUE~"Centromere")) %>%
  filter(par=="nonPAR",male_WT!=0.5,cpg_context=="Island") %>% arrange(pos)
stopifnot(!anyDuplicated(cpg_x$pos),all(cpg_x$cgi_id[!is.na(cpg_x$cgi_id)] %in% cgi_summary$cgi_id))

# 4. Match female RNA and methylation clones
final_clones <- rna_samples %>% filter(sex=="F") %>% pull(sample_id)
common_clones <- intersect(final_clones,intersect(female_clones,names(cgi_summary)))
stopifnot(length(common_clones)>0)
m_long <- cgi_summary %>%
  filter(par=="nonPAR",genomic_region=="promoter",!is.na(gene_name),!gene_name %in% c("",".")) %>%
  separate_rows(gene_name,sep=",") %>% mutate(gene_name=trimws(gene_name)) %>%
  pivot_longer(all_of(common_clones),names_to="sample",values_to="methylation") %>% filter(!is.na(methylation))
uu_gene_clone <- m_long %>% group_by(gene_name,sample) %>% summarise(UU=any(methylation<0.2),.groups="drop") %>% filter(UU)
ase_gene_clone <- ase %>% filter(sample %in% common_clones) %>% group_by(gene_name,sample) %>%
  summarise(n_group=n_distinct(final_group),final_group=first(final_group),.groups="drop")
stopifnot(all(ase_gene_clone$n_group==1))
uu_ase <- uu_gene_clone %>% left_join(select(ase_gene_clone,-n_group),by=c("gene_name","sample")) %>%
  mutate(final_group=replace_na(final_group,"w/oSNP"))

# 5. Calculate BAE frequency and extent among U/U clones
bae_freq <- uu_ase %>% group_by(gene_name) %>%
  summarise(n_UU=n(),n_informative=sum(final_group %in% c("BAE","MAE")),n_BAE=sum(final_group=="BAE"),
            n_LowDepth=sum(final_group=="LowDepth"),BAE_freq=ifelse(n_informative>0,n_BAE/n_informative,NA_real_),.groups="drop")
bae_extent <- ase %>% filter(sample %in% common_clones,final_group=="BAE",!is.na(VAF),between(VAF,0,1)) %>%
  mutate(minor_VAF=pmin(VAF,1-VAF)) %>% group_by(gene_name,sample) %>%
  summarise(clone_extent=mean(minor_VAF),.groups="drop") %>%
  inner_join(uu_ase %>% filter(final_group=="BAE") %>% select(gene_name,sample),by=c("gene_name","sample")) %>%
  group_by(gene_name) %>% summarise(BAE_extent=mean(clone_extent),.groups="drop")
bae_info <- bae_freq %>% left_join(bae_extent,by="gene_name") %>%
  mutate(BAE_extent=ifelse(n_BAE==0 & n_informative>0,0,BAE_extent))

# 6. Prepare CGI-level annotations using the original annotation table
required_anno <- c("gene_name","baseMean","log2FC","lfcSE","stat","pvalue","padj","ASE","ASE_rescue","meanExp")
stopifnot(all(required_anno %in% names(anno_info)),!anyDuplicated(anno_info$gene_name))
hm_anno <- cgi_summary %>%
  separate_rows(gene_name,sep=",") %>% mutate(gene_name=trimws(gene_name)) %>%
  left_join(anno_info %>% select(all_of(required_anno)),by="gene_name") %>%
  left_join(bae_info,by="gene_name") %>%
  arrange(is.na(ASE),is.na(ASE_rescue),is.na(BAE_freq),desc(meanExp)) %>%
  distinct(cgi_id,.keep_all=TRUE) %>%
  mutate(ASE_final=ifelse(is.na(ASE),"w/oSNP",ifelse(!is.na(ASE_rescue),ASE_rescue,ASE)),
         BAE_freq_plot=case_when(!is.na(BAE_freq)~BAE_freq,
                                 n_UU>0 & n_informative==0 & n_LowDepth>0~-1,
                                 n_UU>0 & n_informative==0~-2,TRUE~-2),
         BAE_extent_plot=case_when(!is.na(BAE_extent)~BAE_extent,
                                   n_UU>0 & n_informative==0 & n_LowDepth>0~-1,
                                   n_UU>0 & n_informative==0~-2,TRUE~-2))

# 7. Construct aligned methylation matrices
cpg_x <- cpg_x %>% left_join(hm_anno %>% select(cgi_id,pattern,unmethyl_freq,genomic_region,ASE_final,BAE_freq_plot,BAE_extent_plot,log2FC),by="cgi_id",suffix=c("","_anno"))
# Use the CGI-level annotations, rather than any pre-existing CpG-level fields.
for(nm in c("pattern","unmethyl_freq","genomic_region")){
  alt <- paste0(nm,"_anno")
  if(alt %in% names(cpg_x)) { cpg_x[[nm]] <- cpg_x[[alt]]; cpg_x[[alt]] <- NULL }
}
male_mat <- t(as.matrix(cpg_x[,male_clones,drop=FALSE]))
female_mat <- t(as.matrix(cpg_x[,female_clones,drop=FALSE]))
storage.mode(male_mat) <- "numeric"; storage.mode(female_mat) <- "numeric"
male_mat[is.na(male_mat)] <- -1; female_mat[is.na(female_mat)] <- -1
stopifnot(identical(colnames(male_mat),colnames(female_mat)))

# 8. Heatmap annotations and colors
meth_col <- colorRamp2(c(-1,0,0.5,1),c("#cccccc","#003e73","#f2d8aa","#a51f27"))
position <- ifelse(cpg_x$cgi_id %in% c("cgi_108","cgi_250","cgi_564","cgi_844"),"Example",cpg_x$position)
ha <- HeatmapAnnotation(
  Position=position,Region=cpg_x$genomic_region,Unmethyl_freq=cpg_x$unmethyl_freq,
  Pattern=cpg_x$pattern,ASE=cpg_x$ASE_final,BAE_freq=cpg_x$BAE_freq_plot,
  BAE_extent=cpg_x$BAE_extent_plot,log2FC_FM=cpg_x$log2FC,
  col=list(Position=c("Xp"="#e4c9f0","Xq"="#f0c9dd","Centromere"="#bfbfbf","PAR"="#581082","XIST"="#e06666","Example"="#000000"),
           Region=c("promoter"="#006674","gene_body"="#87c4c3","intergenic"="#e0e0e0"),
           Unmethyl_freq=colorRamp2(c(0,1),c("#ffffff","#003e73")),
           Pattern=c("A"="#b82860","B"="#e4a126","C"="#284054","D"="#54704c","None"="#e0e0e0"),
           ASE=c("MAE"="#472a7a","BAE"="#bea4ed","LowDepth"="#e0e0e0","w/oSNP"="#adadad"),
           BAE_freq=colorRamp2(c(-2,-1,0,1),c("#adadad","#e0e0e0","#ffffff","#c42565")),
           BAE_extent=colorRamp2(c(-2,-1,0,0.5),c("#adadad","#e0e0e0","#ffffff","#c42565")),
           log2FC_FM=colorRamp2(c(-4,-1.5,0,1.5),c("#dbdbdb","#002b8f","#ffffff","#b30202"))),
  annotation_height=unit(rep(2,8),"mm"))
allele <- allele_df$allele[match(female_clones,allele_df$sample_id)]
ha2 <- rowAnnotation(XCI=allele,col=list(XCI=c("A"="#db1860","B"="#053f68")))

# 9. Draw landscape heatmap
hm1 <- Heatmap(male_mat,name="Male methylation",cluster_rows=FALSE,cluster_columns=FALSE,
               show_row_names=FALSE,show_column_names=FALSE,row_split=sub("-.*","",male_clones),
               col=meth_col,top_annotation=ha,use_raster=TRUE,raster_quality=4,
               row_title_gp=gpar(fontsize=8),width=unit(200,"mm"),height=unit(40,"mm"))
female_donor <- sub("-.*","",female_clones)
female_donor[!female_donor %in% sample_f] <- "HC50"
hm2 <- Heatmap(female_mat,name="Female methylation",cluster_rows=FALSE,cluster_columns=FALSE,
               show_row_names=FALSE,show_column_names=FALSE,row_split=factor(female_donor,levels=sample_f),
               col=meth_col,left_annotation=ha2,use_raster=TRUE,raster_quality=4,
               row_title_gp=gpar(fontsize=8),width=unit(200,"mm"),height=unit(90,"mm"))
draw(hm1 %v% hm2)
