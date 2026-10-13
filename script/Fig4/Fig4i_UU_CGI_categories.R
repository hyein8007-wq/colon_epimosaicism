################################################################################
# Figure 4. Hypomethylated CpG islands on chromosome X
# Run from the GitHub repository root; no output files are generated.
################################################################################

# 1. Libraries and inputs
library(dplyr)
library(readr)
library(data.table)

input_dir <- "input/Fig4"
cpg_x <- read_tsv(file.path(input_dir,"class_converted_MF_chrX_260917.tsv.gz"),show_col_types=FALSE)
clone_info <- read_tsv(file.path(input_dir,"clone_info_260917.tsv"),show_col_types=FALSE)
cpg_bed <- read_tsv(file.path(input_dir,"cpgIslandExt.sort.nochr.range.bed"),col_names=FALSE,show_col_types=FALSE) %>%
  transmute(chr=as.character(X1),start=as.integer(X2),end=as.integer(X3)) %>% filter(chr=="X") %>% mutate(cgi_id=paste0("cgi_",row_number()))
manual_pattern <- read_tsv(file.path(input_dir,"unmethyl_x_cgis_pattern_manual_v3_260918.tsv"),show_col_types=FALSE) %>% select(cgi_id,pattern) %>% distinct(cgi_id,.keep_all=TRUE)

# 2. Identify male and female clone columns (order does not affect means)
sample_m <- c("HC01","HC02","HC09","HC10","HC13","HC19")
sample_f <- c("HC08","HC12","HC14","HC15","HC16","HC17","HC18","HC20","HC21","HC50")
clone_ids <- clone_info %>% filter(WGS=="Y",EM_seq=="Y") %>% pull(sampleID)
donor_of <- function(x) sub("-.*$","",x)
male_cols <- intersect(names(cpg_x),clone_ids[donor_of(clone_ids) %in% sample_m])
female_cols <- intersect(names(cpg_x),clone_ids[donor_of(clone_ids) %in% sample_f])
sample_cols <- c(male_cols,female_cols)
stopifnot(length(male_cols)>0,length(female_cols)>0,!anyDuplicated(sample_cols))

# 3. Assign CpG island IDs using genomic overlap
# Preserve the original CpG row order explicitly; BED coordinates are 0-based half-open.
# The original code used start/end directly; convert BED start to 1-based inclusive.
stopifnot(all(c("pos","cpg_context","male_WT","par","promoter_id","promoter_chrom","promoter_start","promoter_end","gene_name","genomic_region") %in% names(cpg_x)))
cpg_dt <- as.data.table(cpg_x)
cpg_dt[,row_id:=.I]
cpg_dt[,`:=`(chr="X",pos_start=as.integer(pos),pos_end=as.integer(pos))]
bed_dt <- as.data.table(cpg_bed)
bed_dt[,start:=start+1L]
setkey(bed_dt,chr,start,end)
matched <- foverlaps(cpg_dt,bed_dt,by.x=c("chr","pos_start","pos_end"),by.y=c("chr","start","end"),type="within",nomatch=NA)
stopifnot(!anyDuplicated(matched$row_id))
setorder(matched,row_id)
cpg_dt[,cgi_id:=matched$cgi_id]
cpg_dt[,c("row_id","chr","pos_start","pos_end"):=NULL]
cpg_x <- as_tibble(cpg_dt)

# 4. Summarize island-level methylation across clones
cgi_total_summary <- cpg_x %>% filter(par=="nonPAR",cpg_context=="Island",male_WT!=0.5,!is.na(cgi_id)) %>%
  group_by(cgi_id) %>% summarise(n_cpg=n(),promoter_id=first(promoter_id),promoter_chrom=first(promoter_chrom),
                                 promoter_start=min(promoter_start),promoter_end=max(promoter_end),gene_name=first(gene_name),
                                 genomic_region=first(genomic_region),par=first(par),across(all_of(sample_cols),~mean(.x,na.rm=TRUE)),.groups="drop") %>%
  mutate(unmethyl_freq=rowMeans(across(all_of(female_cols),~.x<0.2),na.rm=TRUE)) %>%
  left_join(cpg_bed %>% select(cgi_id,start,end),by="cgi_id") %>%
  rename(n_cpg_cgi=n_cpg,cgi_start=start,cgi_end=end)

# 5. Attach manual A-D patterns and restrict to islands with >=2 CpGs
cgi_total_summary_pattern <- cgi_total_summary %>% left_join(manual_pattern,by="cgi_id") %>% filter(n_cpg_cgi>1)
cpg_x <- cpg_x %>% left_join(cgi_total_summary %>% select(cgi_id,n_cpg_cgi,unmethyl_freq),by="cgi_id")

# 6. Inspect results (no files saved)
print(table(cgi_total_summary_pattern$pattern,useNA="ifany"))
print(cgi_total_summary_pattern %>% select(cgi_id,cgi_start,cgi_end,n_cpg_cgi,par,genomic_region,promoter_id,gene_name,pattern,unmethyl_freq))
