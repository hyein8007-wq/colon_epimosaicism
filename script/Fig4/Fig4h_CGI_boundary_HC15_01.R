################################################################################
# Figure 4. Hypomethylated CpG island boundary (HC15-01)
# Run from the repository root; no files are written.
################################################################################

# 1. Libraries and input
library(dplyr)
library(readr)
library(data.table)
library(ggplot2)
library(patchwork)

input_dir <- "input/Fig4"
cpg_x <- read_tsv(file.path(input_dir,"class_converted_MF_chrX_260917.tsv.gz"),show_col_types=FALSE)
cgi_total_summary_pattern <- read_tsv(file.path(input_dir,"chrX_cgi_total_summary_pattern.tsv"),show_col_types=FALSE)
cgi_bed <- read_tsv(file.path(input_dir,"cpgIslandExt.sort.nochr.range.bed"),col_names=FALSE,show_col_types=FALSE) %>%
  transmute(chr=as.character(X1),start=as.integer(X2),end=as.integer(X3)) %>%
  filter(chr=="X") %>% mutate(cgi_id=paste0("cgi_",row_number()))

# 2. Assign CpGs to islands, preserving the original BED-based CGI identifiers
if(!"cgi_id" %in% names(cpg_x)){
  cpg_dt <- as.data.table(cpg_x); cpg_dt[,`:=`(row_id=.I,chr="X",pos_start=as.integer(pos),pos_end=as.integer(pos))]
  bed_dt <- as.data.table(cgi_bed)
  setkey(bed_dt,chr,start,end); setkey(cpg_dt,chr,pos_start,pos_end)
  overlaps <- foverlaps(cpg_dt,bed_dt,by.x=c("chr","pos_start","pos_end"),by.y=c("chr","start","end"),type="within",nomatch=NA)
  stopifnot(!anyDuplicated(overlaps$row_id))
  cpg_x$cgi_id <- overlaps$cgi_id[match(seq_len(nrow(cpg_x)),overlaps$row_id)]
}
if(!"par" %in% names(cpg_x)) cpg_x <- cpg_x %>% mutate(par=if_else(between(pos,60001,2699520)|between(pos,154931044,155260560),"PAR","nonPAR"))

# 3. Select islands hypomethylated in HC15-01
stopifnot(all(c("cgi_id","HC15-01") %in% names(cgi_total_summary_pattern)),all(c("pos","HC15-01","cpg_context","cgi_id","par") %in% names(cpg_x)))
target_cgi <- cgi_total_summary_pattern %>% filter(!is.na(`HC15-01`),`HC15-01`<0.2) %>% pull(cgi_id)
dt <- cpg_x %>% filter(par!="PAR") %>% select(pos,`HC15-01`,cpg_context,cgi_id) %>% as.data.table()
setorder(dt,pos)
dt[,cpg_idx:=.I]
dt[,is_island:=!is.na(cgi_id)]
dt[,beta_state:=fifelse(is.na(`HC15-01`),NA_character_,fifelse(`HC15-01`==0,"0",fifelse(`HC15-01`==0.5,"0.5","1")))]
island_bounds <- dt[cgi_id %in% target_cgi,.(start_idx=min(cpg_idx),end_idx=max(cpg_idx)),by=cgi_id]
stopifnot(nrow(island_bounds)>0)

# 4. Sample up to 1,000 non-island CpGs on each side; stop at the next island
res_list <- lapply(seq_len(nrow(island_bounds)),function(i){
  start <- island_bounds$start_idx[i]; end <- island_bounds$end_idx[i]
  left_idx <- integer(); if(start>1) for(j in seq(start-1,1,-1)){if(dt$is_island[j]) break; left_idx <- c(left_idx,j); if(length(left_idx)>=1000) break}
  right_idx <- integer(); if(end<nrow(dt)) for(j in seq(end+1,nrow(dt))){if(dt$is_island[j]) break; right_idx <- c(right_idx,j); if(length(right_idx)>=1000) break}
  left_dt <- copy(dt[left_idx]); left_dt[,rel_idx:=cpg_idx-start]
  right_dt <- copy(dt[right_idx]); right_dt[,rel_idx:=cpg_idx-end]
  rbind(left_dt,right_dt,fill=TRUE)
})
plot_dt <- rbindlist(res_list,fill=TRUE)
stopifnot(nrow(plot_dt)>0)

# 5. Count methylation states and CpG contexts by relative index
meth_dt <- plot_dt[!is.na(beta_state),.(N=.N),by=.(rel_idx,beta_state)]
meth_dt <- merge(CJ(rel_idx=sort(unique(meth_dt$rel_idx)),beta_state=c("0","0.5","1")),meth_dt,all.x=TRUE)
meth_dt[is.na(N),N:=0]
context_dt <- plot_dt[,.(N=.N),by=.(rel_idx,cpg_context)]
context_dt <- merge(CJ(rel_idx=sort(unique(context_dt$rel_idx)),cpg_context=sort(unique(plot_dt$cpg_context))),context_dt,all.x=TRUE)
context_dt[is.na(N),N:=0]

# 6. Visualization (no PDF or other output files)
p1 <- ggplot(meth_dt,aes(rel_idx,N,fill=beta_state,group=beta_state)) +
  geom_area(position="fill") + scale_fill_manual(values=c("0"="#2166ac","0.5"="#fddbc7","1"="#b2182b")) +
  theme_classic() + labs(x=NULL,y="Epigenotype")
p2 <- ggplot(context_dt,aes(rel_idx,N,fill=cpg_context,group=cpg_context)) +
  geom_area(position="fill") + theme_classic() + labs(x="CpG index from island boundary",y="CpG context")
print(p1/p2)
