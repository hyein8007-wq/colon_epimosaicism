################################################################################
# Figure 2m. Genomic feature enrichment of early methylation regions
# Run from the repository root. No output files are written.
################################################################################

# 1. Libraries
library(data.table)
library(dplyr)
library(tidyr)
library(readr)
library(GenomicRanges)
library(IRanges)
library(rtracklayer)
library(ggplot2)

# 2. Input data
input_dir <- "input/Fig2"
reference_dir <- "reference/Fig2"
early <- fread(file.path(input_dir,"region_only_cov_wide.tsv.gz"),select=c("region_id","chrom","start","end","n_cpg","n_patients_covered"))
background <- fread(file.path(input_dir,"fig2_enrichemnt_background_261002.tsv"))
regulatory <- read_tsv(file.path(reference_dir,"ensembl_regulatory_sorted.grch37.gff.gz"),col_names=FALSE,show_col_types=FALSE)
colnames(regulatory) <- c("chr","regulatory_build","func.anno","start","end","X6","X7","X9","info")
regulatory <- regulatory %>% separate(info,into=c("reg_id","."),sep=";",extra="drop") %>% mutate(reg_id=gsub("id=","",reg_id)) %>% filter(chr %in% as.character(1:22)) %>% distinct(chr,start,end,func.anno,.keep_all=TRUE)
promoter <- as.data.table(regulatory %>% filter(func.anno=="promoter") %>% transmute(chrom=paste0("chr",chr),start,end))
enhancer <- as.data.table(regulatory %>% filter(func.anno=="enhancer") %>% transmute(chrom=paste0("chr",chr),start,end))
ctcf <- fread(file.path(reference_dir,"ENCFF663OEV.hg19.bed"),header=FALSE)[,.(chrom=V1,start=V2,end=V3)]
icr_raw <- fread(file.path(reference_dir,"icr_hg19.tsv"))
icr <- icr_raw[,.(chrom=CHROM,start=START,end=END)]
early[,`:=`(chrom=paste0("chr",sub("^chr","",chrom)),width=end-start+1)]
background[,chrom:=paste0("chr",sub("^chr","",chrom))]
stopifnot(all(c("region_id","perm_id","chrom","start","end","width","n_cpg") %in% names(background)))
stopifnot(all(file.exists(file.path(reference_dir,c("E075-H3K4me3.imputed.pval.signal.bigwig","E075-H3K27ac.imputed.pval.signal.bigwig","E075-H3K27me3.imputed.pval.signal.bigwig","E075-DNase.imputed.pval.signal.bigwig","wgEncodeUwRepliSeqImr90WaveSignalRep1.bigWig")))))

# 3. Matched background generation (optional; use existing background by default)
# Original matching: chromosome and number of CpGs, excluding overlaps with early regions.
make_background <- function(early,bg_cpg,n_perm=1000){
  early <- copy(as.data.table(early)); bg_cpg <- copy(as.data.table(bg_cpg))
  setkey(early,chrom,start,end)
  bg_list <- vector("list",nrow(early))
  for(i in seq_len(nrow(early))){
    if(i%%100==0) message(i," / ",nrow(early))
    chr <- early$chrom[i]; ncpg <- early$n_cpg[i]; tmp <- bg_cpg[chrom==chr][order(cpos)]
    n_window <- nrow(tmp)-ncpg+1
    if(is.na(n_window) || n_window<=0) next
    candidate <- data.table(chrom=chr,start=tmp$cpos[seq_len(n_window)],end=tmp$cpos[seq.int(ncpg,n_window+ncpg-1)])
    candidate[,`:=`(width=end-start+1,cid=.I)]; setkey(candidate,chrom,start,end)
    ov <- foverlaps(candidate,early,nomatch=0L)
    if(nrow(ov)>0) candidate <- candidate[!cid %in% unique(ov$cid)]
    candidate[,cid:=NULL]
    if(nrow(candidate)==0) next
    candidate <- candidate[sample(.N,n_perm,replace=(.N<n_perm))]
    candidate[,`:=`(region_id=early$region_id[i],perm_id=seq_len(.N),n_cpg=ncpg)]
    setcolorder(candidate,c("region_id","perm_id","chrom","start","end","width","n_cpg"))
    bg_list[[i]] <- candidate
  }
  rbindlist(bg_list)
}

# To regenerate the background rather than use the supplied file, uncomment:
# bg_cpg <- fread(file.path(input_dir,"CG_cpos.tsv"))
# setnames(bg_cpg,"chr","chrom")
# bg_cpg[,chrom:=paste0("chr",sub("^chr","",chrom))]
# set.seed(1234)
# background <- make_background(early,bg_cpg,n_perm=1000)

# 4. Background quality control
qc <- merge(background,early[,.(region_id,early_width=width,early_n_cpg=n_cpg)],by="region_id")
print(table(qc$n_cpg==qc$early_n_cpg))
print(summary(qc$width/qc$early_width))

# 5. Genomic annotation and BigWig signal functions
calc_annotation <- function(regions,feature,window=10000,id_col="region_id"){
  regions <- as.data.table(regions); feature <- as.data.table(feature)
  gr_region <- GRanges(seqnames=regions$chrom,ranges=IRanges(start=pmax(1L,regions$start-window),end=regions$end+window),region=regions[[id_col]])
  gr_feature <- GRanges(seqnames=feature$chrom,ranges=IRanges(start=feature$start,end=feature$end))
  hit <- findOverlaps(gr_region,gr_feature)
  if(length(hit)==0) return(data.table(region=regions[[id_col]],value=0L))
  out <- data.table(region=S4Vectors::mcols(gr_region)$region[S4Vectors::queryHits(hit)])[,.(value=.N),by=region]
  out <- merge(unique(data.table(region=regions[[id_col]])),out,by="region",all.x=TRUE)
  out[is.na(value),value:=0L]; setorder(out,region); out
}
calc_bigwig_signal <- function(regions,gr_signal,window=10000,id_col="region_id"){
  regions <- as.data.table(regions)
  gr_query <- GRanges(seqnames=regions$chrom,ranges=IRanges(start=pmax(1L,regions$start-window),end=regions$end+window),region=regions[[id_col]])
  hit <- findOverlaps(gr_query,gr_signal)
  if(length(hit)==0) return(data.table(region=regions[[id_col]],value=0))
  q <- S4Vectors::queryHits(hit); s <- S4Vectors::subjectHits(hit)
  overlap_bp <- pmin(end(gr_query)[q],end(gr_signal)[s])-pmax(start(gr_query)[q],start(gr_signal)[s])+1L
  dt <- data.table(region=S4Vectors::mcols(gr_query)$region[q],score=S4Vectors::mcols(gr_signal)$score[s],bp=overlap_bp)
  out <- dt[,.(value=sum(score*bp)/sum(bp)),by=region]
  out <- merge(data.table(region=regions[[id_col]]),out,by="region",all.x=TRUE)
  out[is.na(value),value:=0]; out
}

# 6. Permutation-based fold enrichment
calc_feature <- function(early,background,feature,feature_name,feature_type=c("annotation","bigwig"),window=10000){
  feature_type <- match.arg(feature_type)
  perm_ids <- sort(unique(background$perm_id)); exp_mean <- numeric(length(perm_ids))
  if(feature_type=="annotation"){
    observed <- calc_annotation(early,feature,window,"region_id"); obs_mean <- mean(observed$value,na.rm=TRUE)
    for(i in seq_along(perm_ids)){
      if(i%%100==0) message(feature_name,": ",i," / ",length(perm_ids))
      bg <- background[perm_id==perm_ids[i]]
      exp_mean[i] <- mean(calc_annotation(bg,feature,window,"region_id")$value,na.rm=TRUE)
    }
  }else{
    message(feature_name,": importing BigWig")
    gr_signal <- rtracklayer::import(feature)
    observed <- calc_bigwig_signal(early,gr_signal,window,"region_id"); obs_mean <- mean(observed$value,na.rm=TRUE)
    for(i in seq_along(perm_ids)){
      if(i%%100==0) message(feature_name,": ",i," / ",length(perm_ids))
      bg <- background[perm_id==perm_ids[i]]
      exp_mean[i] <- mean(calc_bigwig_signal(bg,gr_signal,window,"region_id")$value,na.rm=TRUE)
      if(i%%100==0) gc()
    }
    rm(gr_signal); gc()
  }
  expected_mean <- mean(exp_mean); fold <- obs_mean/expected_mean
  p_upper <- (sum(exp_mean>=obs_mean)+1)/(length(exp_mean)+1)
  p_lower <- (sum(exp_mean<=obs_mean)+1)/(length(exp_mean)+1)
  fold_perm <- obs_mean/exp_mean
  data.table(feature=feature_name,observed=obs_mean,expected=expected_mean,fold=fold,log2FE=log2(fold),z_score=(obs_mean-expected_mean)/sd(exp_mean),p_value=min(1,2*min(p_upper,p_lower)),ci_low=quantile(fold_perm,.025,na.rm=TRUE),ci_high=quantile(fold_perm,.975,na.rm=TRUE),perm_mean=list(exp_mean),perm_fold=list(fold_perm))
}

# 7. Run enrichment tests (original analysis uses the first 100 permutations)
background100 <- background[perm_id<=100]
res <- rbindlist(list(
  calc_feature(early,background100,promoter,"Promoter","annotation",10000),
  calc_feature(early,background100,enhancer,"Enhancer","annotation",10000),
  calc_feature(early,background100,ctcf,"CTCF","annotation",10000),
  calc_feature(early,background100,icr,"ICR","annotation",10000),
  calc_feature(early,background100,file.path(reference_dir,"E075-H3K4me3.imputed.pval.signal.bigwig"),"H3K4me3","bigwig",10000),
  calc_feature(early,background100,file.path(reference_dir,"E075-H3K27ac.imputed.pval.signal.bigwig"),"H3K27ac","bigwig",10000),
  calc_feature(early,background100,file.path(reference_dir,"E075-H3K27me3.imputed.pval.signal.bigwig"),"H3K27me3","bigwig",10000),
  calc_feature(early,background100,file.path(reference_dir,"E075-DNase.imputed.pval.signal.bigwig"),"DHS","bigwig",10000),
  calc_feature(early,background100,file.path(reference_dir,"wgEncodeUwRepliSeqImr90WaveSignalRep1.bigWig"),"RepliSeq_IMR90","bigwig",100000)
))
print(res[,.(feature,observed,expected,fold,log2FE,z_score,p_value)])

# 8. Figure 2m
plot_df <- copy(res)
plot_df[,`:=`(ci_low_log2=log2(ci_low),ci_high_log2=log2(ci_high))]
plot_df[,feature:=reorder(feature,log2FE)]
g1 <- ggplot(plot_df,aes(x=log2FE,y=feature))+
  geom_vline(xintercept=0,linetype=2,colour="grey55")+
  geom_col(width=.55,fill="grey45")+
  geom_errorbarh(aes(xmin=ci_low_log2,xmax=ci_high_log2),height=.18,linewidth=.8)+
  geom_text(aes(label=sprintf("%.2f",fold),hjust=ifelse(log2FE>=0,-.15,1.15)),size=3.8)+
  labs(x="log2(Observed / Expected)",y=NULL)+theme_bw(base_size=13)+
  theme(axis.text.y=element_text(size=11),axis.title.x=element_text(size=13),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),panel.grid.major.x=element_line(colour="grey92"),legend.position="none")
g2 <- ggplot(plot_df,aes(x=log2FE,y=feature,label=sprintf("%.2f×",fold)))+
  geom_vline(xintercept=0,linetype=2,colour="grey65",linewidth=.6)+
  geom_errorbarh(aes(xmin=ci_low_log2,xmax=ci_high_log2),height=0,linewidth=.8,colour="grey35")+
  geom_point(size=3.2,shape=21,fill="white",stroke=1)+
  geom_text(nudge_x=.06,hjust=0,size=3.7)+
  scale_x_continuous(breaks=c(-.5,0,.5,1,1.5,2,2.5),expand=expansion(mult=c(.03,.12)))+
  labs(x=expression(log[2]*"(Observed / Expected)"),y=NULL)+theme_classic(base_size=13)+
  theme(axis.text.y=element_text(size=11,colour="black"),axis.text.x=element_text(colour="black"),axis.title.x=element_text(size=13,margin=margin(t=8)),axis.line.y=element_blank(),axis.ticks.y=element_blank())
print(g1)
print(g2)

# 9. Supplementary table (in-memory only)
supp_enrichment <- as.data.frame(res) %>%
  mutate(Category=case_when(feature %in% c("Promoter","Enhancer","CTCF","ICR")~"Genomic annotation",feature %in% c("H3K4me3","H3K27ac","H3K27me3","DHS")~"Epigenomic signal",grepl("^RepliSeq",feature)~"Replication timing"),
         Data_source=case_when(feature %in% c("Promoter","Enhancer")~"Ensembl",feature=="CTCF"~"ENCODE",feature=="ICR"~"Human Imprintome",feature %in% c("H3K4me3","H3K27ac","H3K27me3","DHS")~"NIH Roadmap Epigenomics",grepl("^RepliSeq",feature)~"ENCODE/UW"),
         Reference_ID=case_when(feature %in% c("Promoter","Enhancer")~"Ensembl Regulatory Build (GRCh37)",feature=="CTCF"~"ENCFF663OEV",feature=="ICR"~"Human Imprintome",feature %in% c("H3K4me3","H3K27ac","H3K27me3","DHS")~"E075",feature=="RepliSeq_IMR90"~"IMR90"),
         Biosample=case_when(feature %in% c("H3K4me3","H3K27ac","H3K27me3","DHS")~"Colonic mucosa",feature=="CTCF"~"Colon mucosa",feature=="RepliSeq_IMR90"~"IMR90 fibroblast",TRUE~NA_character_),
         Tissue_origin=case_when(feature %in% c("H3K4me3","H3K27ac","H3K27me3","DHS","CTCF")~"Colon",feature=="RepliSeq_IMR90"~"Fibroblast",TRUE~NA_character_),
         Assay=case_when(feature %in% c("Promoter","Enhancer","ICR")~"Genomic annotation",feature %in% c("CTCF","H3K4me3","H3K27ac","H3K27me3")~"ChIP-seq",feature=="DHS"~"DNase-seq",grepl("^RepliSeq",feature)~"Repli-seq"),
         Assembly="GRCh37",Window_bp=ifelse(grepl("^RepliSeq",feature),100000,10000)) %>%
  transmute(Feature=feature,Category,Data_source,Reference_ID,Biosample,Tissue_origin,Assay,Assembly,Window_bp,Observed=observed,Expected=expected,Fold_enrichment=fold,Log2_fold_enrichment=log2FE,Z_score=z_score,Empirical_P_value=p_value,CI_low=ci_low,CI_high=ci_high)
print(supp_enrichment)
