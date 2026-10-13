################################################################################
# Figure 4f. XCI timing inferred from shared mutations and blood VAF
# Run from the repository root. No files are written.
################################################################################

# 1. Libraries and input paths
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(ape)
library(phangorn)
library(data.table)
library(patchwork)

tree_dir <- "input/Fig4/tree_F"
vaf_dir <- "input/Fig4/blood_vaf_F"
allele_file <- "input/Fig4/allele_info_260918.tsv"

# 2. Read XCI allele assignments and select trees
allele_info <- read_tsv(allele_file,show_col_types=FALSE)
stopifnot(all(c("clone_id","allele") %in% names(allele_info)))
df_allele <- allele_info %>% transmute(sample_id=as.character(clone_id),allele=as.character(allele)) %>% distinct(sample_id,.keep_all=TRUE)

files <- list.files(tree_dir,pattern="edit\\.nwk$",full.names=FALSE)
files <- setdiff(files,c("HC50_241017.length.edit.nwk","HC50.length.edit.nwk","HC01.length.edit.nwk","HC02.length.edit.nwk","HC13.length.edit.nwk"))
hc50_file <- "HC50_colon_final_260919_woConta.length.nwk"
if(file.exists(file.path(tree_dir,hc50_file))) files <- c(files,hc50_file)
if(!length(files)) stop("No eligible tree files found in ",tree_dir)

# 3. Count shared mutations along the ancestral path of each clone pair
shared_from_tree <- function(f){
  tree <- ape::read.tree(file.path(tree_dir,f))
  if(is.null(tree$edge.length)) stop("Missing branch lengths: ",f)
  sample <- substr(f,1,4)
  tips <- tree$tip.label
  if(length(tips)<2) return(tibble())
  paths <- ape::nodepath(tree)
  names(paths) <- as.character(seq_along(paths))
  pairs <- combn(seq_along(tips),2,simplify=FALSE)
  bind_rows(lapply(pairs,function(pair){
    a <- paths[[as.character(pair[1])]]; b <- paths[[as.character(pair[2])]]
    k <- 0L
    while(k<min(length(a),length(b)) && a[k+1L]==b[k+1L]) k <- k+1L
    shared_nodes <- head(a,k)
    shared_edges <- which(tree$edge[,2] %in% tail(shared_nodes,-1))
    tibble(sample=sample,clone1=sort(tips[pair])[1],clone2=sort(tips[pair])[2],shared_length=sum(tree$edge.length[shared_edges]))
  }))
}

df_shared_length_final <- bind_rows(lapply(files,shared_from_tree)) %>% mutate(pairID=paste(clone1,clone2,sep=","))
df_shared_length_allele_final <- df_shared_length_final %>%
  left_join(df_allele %>% rename(clone1=sample_id,allele1=allele),by="clone1") %>%
  left_join(df_allele %>% rename(clone2=sample_id,allele2=allele),by="clone2") %>%
  mutate(shared_allele=case_when(is.na(allele1)|is.na(allele2)~NA_character_,allele1==allele2~"same",TRUE~"different"))

# 4. Shared XCI allele plots (original plotting definitions)
p_bar <- df_shared_length_allele_final %>% filter(shared_length<=30) %>%
  ggplot(aes(shared_length,fill=shared_allele)) + geom_bar(position="fill") +
  scale_fill_brewer(palette="Blues") + geom_vline(xintercept=18,color="red") +
  labs(x="Shared somatic mutations",y="Proportion",fill="Shared allele") + theme_classic()

p_density <- df_shared_length_allele_final %>% filter(shared_length<30,sample!="HC50") %>%
  ggplot(aes(shared_length,fill=shared_allele)) +
  geom_density(position="fill",bw=2,adjust=0.8) + scale_fill_brewer(palette="Blues") +
  geom_vline(xintercept=18,color="red") + coord_cartesian(xlim=c(0,30)) + theme_classic()

# 5. Infer blood VAF for each clone pair
# As in the source: direct = nth-highest shared blood VAF; otherwise inferred from descendants.
get_blood_vaf <- function(s){
  message("Processing ",s)
  df <- df_shared_length_allele_final %>% filter(sample==s)
  tree_file <- file.path(tree_dir,paste0(s,".length.edit.nwk"))
  vaf_file <- file.path(vaf_dir,paste0(s,".variant.merge.mark.edit.vaf.tsv"))
  if(!file.exists(tree_file) || !file.exists(vaf_file)){
    warning("Skipping ",s,": missing standard tree or blood VAF input")
    return(NULL)
  }
  tree <- ape::read.tree(tree_file)
  var <- data.table::fread(vaf_file)
  stopifnot(all(c("shared_sampleID","VAF") %in% names(var)))
  sl <- strsplit(ifelse(is.na(var$shared_sampleID),"",var$shared_sampleID),",",fixed=TRUE)
  mem <- vapply(tree$tip.label,function(tip) vapply(sl,function(x) tip %in% x,logical(1)),logical(nrow(var)))
  if(is.null(dim(mem))) mem <- matrix(mem,ncol=length(tree$tip.label),dimnames=list(NULL,tree$tip.label))
  colnames(mem) <- tree$tip.label
  mm <- ape::mrca(tree,full=FALSE)
  if(!all(c(df$clone1,df$clone2) %in% rownames(mm))) stop("Tree tip mismatch: ",s)
  df$node <- mapply(function(x,y) mm[x,y],df$clone1,df$clone2)
  desc <- setNames(lapply(unique(df$node),function(n) tree$tip.label[phangorn::Descendants(tree,n,type="tips")[[1]]]),unique(df$node))
  res <- lapply(seq_len(nrow(df)),function(i){
    c1 <- df$clone1[i]; c2 <- df$clone2[i]; n <- df$shared_length[i]
    if(n==0) return(c(direct=0.5,inferred=NA_real_))
    shared <- mem[,c1] & mem[,c2]
    v <- sort(var$VAF[shared & !is.na(var$VAF) & var$VAF>0],decreasing=TRUE)
    direct <- if(n>=1 && n==floor(n) && length(v)>=n) v[n] else NA_real_
    if(!is.na(direct)) return(c(direct=direct,inferred=NA_real_))
    d <- desc[[as.character(df$node[i])]]
    idx <- rowSums(mem[,d,drop=FALSE])>0 & !shared & !is.na(var$VAF) & var$VAF>0
    c(direct=NA_real_,inferred=if(any(idx)) max(var$VAF[idx]) else NA_real_)
  })
  res <- do.call(rbind,res)
  df %>% mutate(direct_vaf=res[,"direct"],inferred_vaf=res[,"inferred"],
                blood_vaf=case_when(shared_length==0~0.5,!is.na(direct_vaf)~direct_vaf,!is.na(inferred_vaf)~inferred_vaf,TRUE~0),
                vaf_status=case_when(shared_length==0~"root",!is.na(direct_vaf)~"direct",!is.na(inferred_vaf)~"inferred",TRUE~"undetected"))
}

df_blood_vaf <- bind_rows(lapply(unique(df_shared_length_allele_final$sample),get_blood_vaf))
if(!nrow(df_blood_vaf)) stop("No blood VAF data matched the tree inputs")
donor_vaf <- df_blood_vaf %>% filter(shared_length<=30) %>% group_by(sample,shared_length) %>% summarise(vaf=mean(blood_vaf),.groups="drop")
vaf_summary <- donor_vaf %>% group_by(shared_length) %>% summarise(mean_vaf=mean(vaf),se=sd(vaf)/sqrt(n()),n=n(),.groups="drop")

# 6. Blood VAF plots and combined figure
p_vaf <- ggplot(vaf_summary,aes(shared_length,mean_vaf)) +
  geom_errorbar(aes(ymin=mean_vaf-se,ymax=mean_vaf+se),width=0.15,na.rm=TRUE) +
  geom_point(size=2.5) + geom_smooth(method="loess",se=FALSE,span=0.8,linewidth=1) +
  geom_vline(xintercept=18,linetype=2) + coord_cartesian(ylim=c(0,0.5)) +
  labs(x="Shared somatic mutations",y="Mean blood VAF") + theme_classic()

p_donor_vaf <- ggplot(donor_vaf,aes(shared_length,vaf,group=sample)) +
  geom_line(alpha=0.2) + geom_point(alpha=0.3,size=1) +
  geom_errorbar(data=vaf_summary,aes(x=shared_length,ymin=mean_vaf-se,ymax=mean_vaf+se),width=0.15,inherit.aes=FALSE,na.rm=TRUE) +
  geom_point(data=vaf_summary,aes(shared_length,mean_vaf),size=2.5,inherit.aes=FALSE) +
  geom_smooth(data=vaf_summary,aes(shared_length,mean_vaf),method="loess",se=FALSE,span=0.8,linewidth=1,inherit.aes=FALSE) +
  geom_vline(xintercept=18,linetype=2) + coord_cartesian(ylim=c(0,0.5)) +
  labs(x="Shared somatic mutations",y="Blood VAF") + theme_classic()

p_combined <- (p_bar + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())) / p_vaf + plot_layout(heights=c(1,1))
print(p_combined)
