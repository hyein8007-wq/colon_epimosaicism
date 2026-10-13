################################################################################
# Figure 5g. Application of EEC-guided phylogenetic refinement
# Run from the GitHub repository root. No files are written.
################################################################################

# 1. Libraries
library(data.table)
library(dplyr)
library(ape)
library(phangorn)

# 2. Input data
input_dir <- "input/Fig5/application"
tree_files <- list.files(input_dir,pattern="^HC[0-9]+\\.length\\.edit\\.nwk$",full.names=TRUE)
eec_file <- file.path(input_dir,"clone_state_mode_segment_for_phylo.tsv.gz")
coverage_file <- file.path(input_dir,"region_element_cov_wide.tsv.gz")
stopifnot(length(tree_files)>0,all(file.exists(c(eec_file,coverage_file))))

# The EEC methylation-state matrix is required in addition to the trees and coverage table.
early_seg_total <- fread(eec_file)
seg_candi <- fread(coverage_file)
stopifnot(all(c("chrom","start","end","region_id","segment_id") %in% names(early_seg_total)))
stopifnot(all(c("chrom","start","end","n_patients_covered") %in% names(seg_candi)))

# 3. Candidate topology scoring and phylogenetic refinement
# Original functions: identify polytomies, enumerate candidates, assess supporting EECs,
# calculate Fitch parsimony, prediction confidence and refined topologies.

split_top <- function(x){
  depth <- 0
  start <- 1
  out <- c()
  
  for(i in seq_len(nchar(x))){
    ch <- substr(x,i,i)
    if(ch=="(") depth <- depth+1
    if(ch==")") depth <- depth-1
    if(ch=="," & depth==0){
      out <- c(out,substr(x,start,i-1))
      start <- i+1
    }
  }
  
  out <- c(out,substr(x,start,nchar(x)))
  trimws(out)
  
}


get_leaves <- function(x){
  leaves <- regmatches(x, gregexpr("[^,():]+(?=:)", x, perl = TRUE))[[1]]
  unique(leaves)
}


find_multifurcation <- function(x){
  res <- list()
  recurse <- function(s){
    if(!grepl("^\\(",s))
      return()
    
    depth <- 0
    end <- NA
    
    for(i in seq_len(nchar(s))){
      ch <- substr(s,i,i)
      if(ch=="(") depth <- depth+1
      if(ch==")"){
        depth <- depth-1
        if(depth==0){
          end <- i
          break
        }
      }
    }
    inside <- substr(s,2,end-1)
    child <- split_top(inside)
    if(length(child)>=3){
      before <- lapply(child,get_leaves)
      res[[length(res)+1]] <<- before
    }
    for(cc in child){
      recurse(cc)
    }
  }
  recurse(x)
  res
}



make_candidates <- function(before){
  
  n <- length(before)
  stopifnot(n >= 3)
  out <- list()
  idx <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      left <- sort(c(before[[i]], before[[j]]))
      right <- sort(unlist(before[-c(i,j)]))
      out[[idx]] <- list(
        left = left,
        right = right,
        before = before,
        after_string = paste0(
          "(",
          paste(left, collapse=","),
          ") | (",
          paste(right, collapse=","),
          ")"
        )
      )
      idx <- idx+1
    }
  }
  out
}




find_multifurcation_node <- function(tree, before){
  
  edge <- tree$edge
  tips <- tree$tip.label
  Ntip <- length(tips)
  
  get_descendants <- function(node){
    child <- edge[edge[,1]==node,2]
    out <- c()
    for(ch in child){
      if(ch <= Ntip){
        out <- c(out, tips[ch])
      }else{
        out <- c(out, get_descendants(ch))
      }
    }
    unique(out)
  }
  
  target <- sort(unique(unlist(before)))
  
  internal_nodes <- unique(edge[,1])
  for(node in internal_nodes){
    
    child <- edge[edge[,1]==node,2]
    
    if(length(child) != length(before))
      next
    
    child_sets <- lapply(child,function(x){
      if(x <= Ntip){tips[x]}
      else{get_descendants(x)}
    })
    
    child_sets <- lapply(child_sets, sort)
    
    before2 <- lapply(before, sort)
    
    matched <- rep(FALSE,length(before2))
    
    for(cs in child_sets){
      idx <- which(
        sapply(before2,function(x)
          identical(x,cs)))
      if(length(idx)==1)
        matched[idx] <- TRUE
    }
    
    if(all(matched)){
      return(node)
    }
    
  }
  
  stop("Cannot find corresponding multifurcation node.")
  
}



get_descendants <- function(tree, node){
  
  edge <- tree$edge
  tips <- tree$tip.label
  Ntip <- length(tips)
  
  recurse <- function(x){
    
    child <- edge[edge[,1] == x, 2]
    out <- c()
    for(ch in child){
      if(ch <= Ntip){
        out <- c(out, tips[ch])
      }else{
        out <- c(out, recurse(ch))
      }
    }
    out
  }
  
  sort(unique(recurse(node)))
  
}



make_reduced_tree <- function(tree,
                              candidate){
  reduced_tree <- ape::drop.tip(tree, candidate$left)
  reduced_tree
}


score_candidate <- function(early_seg,
                            candidate,
                            clone_cols,
                            min_frac = 0.8,
                            outside_frac = 0.8){
  
  left <- candidate$left
  
  outside <- setdiff(clone_cols,left)
  
  support <- apply(early_seg[, clone_cols],1,function(x){
    
    
    left_state <- na.omit(x[left])
    
    if(length(left_state) < ceiling(length(left) * min_frac))
      return(FALSE)
    
    if(length(unique(left_state)) != 1)
      return(FALSE)
    
    left_value <- as.character(unique(left_state))
    
    
    outside_state <- na.omit(x[outside])
    
    if(length(outside_state) < ceiling(length(outside) * min_frac))
      return(FALSE)
    tb <- table(outside_state)
    
    if(sum(tb == max(tb)) > 1)
      return(FALSE)
    major_frac <- max(tb) / sum(tb)
    
    if(major_frac < outside_frac)
      return(FALSE)
    major_value <- names(tb)[which.max(tb)]
    
    if(left_value == major_value)
      return(FALSE)
    TRUE
  }
  )
  list(
    support_n = sum(support),
    support_region = early_seg[support, ]
  )
}



calculate_parsimony <- function(region,
                                reduced_tree,
                                min_frac = 0.8){
  
  state <- unlist(region[reduced_tree$tip.label])
  names(state) <- reduced_tree$tip.label
  
  if(sum(!is.na(state)) < ceiling(length(state)*min_frac))
    return(NA)
  
  keep <- !is.na(state)
  tree_tmp <- reduced_tree
  
  if(any(!keep)){
    tree_tmp <- ape::drop.tip(tree_tmp, names(state)[!keep])
    state <- state[keep]
  }
  
  dat <- matrix(as.character(state), ncol=1)
  rownames(dat) <- names(state)
  dat <- phangorn::phyDat(dat, type="USER", levels=c("0","0.5","1"))
  
  score <- phangorn::parsimony(tree=tree_tmp, data=dat, method="fitch", site="pscore")
  as.numeric(score)
  
}



get_node_shared_mutation <- function(tree, node){
  edge <- tree$edge
  edge.length <- tree$edge.length
  root <- setdiff(edge[,1], edge[,2])
  total <- 0
  cur <- node
  repeat{
    idx <- which(edge[,2] == cur)
    if(length(idx)==0) break
    total <- total + edge.length[idx]
    parent <- edge[idx,1]
    if(parent == root) break
    cur <- parent
  }
  total
}

evaluate_candidate <- function(early_seg,
                               candidate,
                               clone_cols,
                               tree,
                               min_frac = 0.8,
                               outside_frac = 0.8){
  reduced_tree <- make_reduced_tree(tree = tree, candidate = candidate)
  res <- score_candidate(
    early_seg = early_seg,
    candidate = candidate,
    clone_cols = clone_cols,
    min_frac = min_frac,
    outside_frac = outside_frac
  )
  support_region <- res$support_region
  parsimony <- rep(NA_real_, nrow(support_region))
  if(nrow(support_region) > 0){
    for(i in seq_len(nrow(support_region))){
      parsimony[i] <- calculate_parsimony(
        region = support_region[i, ],
        reduced_tree = reduced_tree,
        min_frac = min_frac
      )
    }
  }
  support_region$parsimony <- parsimony
  score_valid <- parsimony[!is.na(parsimony)]
  summary_tbl <- table(factor(ifelse(score_valid >= 3,"3+",as.character(score_valid)),
                              levels = c("0","1","2","3+")))
  list(
    candidate = candidate,
    reduced_tree = reduced_tree,
    support_n = nrow(support_region),
    summary = summary_tbl,
    support_region = support_region
  )
}


analyze_all_multifurcations <- function(tree,
                                        nwk,
                                        early_seg,
                                        clone_cols,
                                        min_frac=0.8,
                                        outside_frac=0.8){
  before_all <- find_multifurcation(nwk)
  out <- vector("list", length(before_all))
  for(case_n in seq_along(before_all)){
    message("Analyzing multifurcation ", case_n, " / ", length(before_all))
    before_case <- before_all[[case_n]]
    node <- find_multifurcation_node(tree=tree, before=before_case)
    shared_mutation <- get_node_shared_mutation(tree=tree, node=node)
    cand <- make_candidates(before_case)
    result <- vector("list", length(cand))
    for(i in seq_along(cand)){
      result[[i]] <- evaluate_candidate(
        early_seg=early_seg,
        candidate=cand[[i]],
        clone_cols=clone_cols,
        tree=tree,
        min_frac=min_frac,
        outside_frac=outside_frac
      )
    }
    summary_df <- bind_rows(lapply(seq_along(result), function(i){
      sm <- result[[i]]$summary
      data.frame(
        Candidate=paste0("Candidate ",i),
        Total=result[[i]]$support_n,
        P0=as.numeric(sm["0"]),
        P1=as.numeric(sm["1"]),
        P2=as.numeric(sm["2"]),
        P3plus=as.numeric(sm["3+"]),
        stringsAsFactors=FALSE
      )
    }))
    all_scores <- unlist(lapply(result, function(x) x$support_region$parsimony))
    all_scores <- all_scores[!is.na(all_scores)]
    candidate_total <- sapply(result, function(x) x$support_n)
    if(length(all_scores)==0){
      min_parsimony <- NA_real_
      min_count <- rep(0,length(result))
      best_candidate <- NA_integer_
    }else{
      min_parsimony <- min(all_scores)
      min_count <- sapply(result, function(x)
        sum(x$support_region$parsimony==min_parsimony, na.rm=TRUE))
      best_candidate <- which.max(min_count)
    }
    ord <- order(min_count, decreasing=TRUE)
    top_candidate <- if(length(ord)>=1) ord[1] else NA_integer_
    second_candidate <- if(length(ord)>=2) ord[2] else NA_integer_
    top_count <- if(!is.na(top_candidate)) min_count[top_candidate] else 0
    second_count <- if(!is.na(second_candidate)) min_count[second_candidate] else 0
    confidence <- if(top_count + second_count == 0) NA_real_ else top_count/(top_count + second_count)
    out[[case_n]] <- list(
      case_n=case_n,
      node=node,
      before=before_case,
      candidates=cand,
      result=result,
      summary=summary_df,
      shared_mutation=shared_mutation,
      min_parsimony=min_parsimony,
      min_count=min_count,
      candidate_total=candidate_total,
      best_candidate=best_candidate,
      top_candidate=top_candidate,
      second_candidate=second_candidate,
      top_count=top_count,
      second_count=second_count,
      confidence=confidence
    )
  }
  out
}


resolve_one_multifurcation <- function(
  tree,
  before_case,
  selected_tips
){
  
  
  target_node <- find_multifurcation_node(tree = tree, before = before_case)
  edge <- tree$edge
  tips <- tree$tip.label
  Ntip <- length(tips)
  
  
  get_node_tips <- function(node){
    
    if(node <= Ntip){
      return(tips[node])
    }
    children <- edge[edge[, 1] == node, 2]
    sort(unique(unlist(lapply(children,get_node_tips))))
  }
  
  
  target_children <- edge[edge[, 1] == target_node,2]
  selected_children <- target_children[
    vapply(
      target_children,
      function(ch){
        child_tips <- get_node_tips(ch)
        all(child_tips %in% selected_tips)
      },
      logical(1)
    )
  ]
  
  
  if(length(selected_children) < 2){
    stop(
      "Could not identify at least two direct child lineages ",
      "for the selected candidate."
    )
  }
  
  
  if(length(selected_children) == length(target_children)){
    stop(
      "Selected candidate contains all direct children ",
      "of the multifurcation."
    )
  }
  
  
  build_newick <- function(node){
    
    if(node <= Ntip){
      return(tips[node])
    }
    
    
    children <- edge[edge[, 1] == node, 2]
    
    
    if(node == target_node){
      
      other_children <- setdiff(children, selected_children)
      
      selected_strings <- vapply(selected_children, build_newick, character(1))
      
      selected_clade <- paste0("(",paste(selected_strings,collapse = ","),")")
      
      other_strings <- vapply(other_children, build_newick, character(1))
      
      
      return(paste0("(",paste(c(selected_clade,other_strings),collapse = ","),")"))
    }
    
    
    child_strings <- vapply(children, build_newick, character(1))
    paste0("(", paste( child_strings, collapse = "," ), ")")
    
  }
  
  
  root_node <- setdiff(unique(edge[, 1]), unique(edge[, 2]))
  
  if(length(root_node) != 1){
    stop("Could not uniquely identify the root node.")
  }
  
  
  new_nwk <- paste0(build_newick(root_node), ";")
  
  
  refined_tree <- ape::read.tree(text = new_nwk)
  refined_tree
}



make_refined_tree <- function(tree,
                              analysis,
                              min_shared_mutation=5,
                              min_confidence=NULL){
  refined_tree <- tree
  for(k in seq_along(analysis)){
    x <- analysis[[k]]
    if(is.null(x$best_candidate) || length(x$best_candidate)==0 || is.na(x$best_candidate)){
      message("M", k, ": no candidate selected -> skipped")
      next
    }
    if(is.na(x$shared_mutation) || x$shared_mutation <= min_shared_mutation){
      message("M", k, ": shared mutation <= ", min_shared_mutation, " -> skipped")
      next
    }
    if(!is.null(min_confidence) && (is.na(x$confidence) || x$confidence < min_confidence)){
      message("M", k, ": confidence < ", min_confidence, " -> skipped")
      next
    }
    selected_candidate <- x$candidates[[x$best_candidate]]
    message("M", k, ": Candidate ", x$best_candidate,
            " selected | shared=", x$shared_mutation,
            " | confidence=", round(x$confidence,3))
    refined_tree <- resolve_one_multifurcation(
      tree=refined_tree,
      before_case=x$before,
      selected_tips=selected_candidate$left
    )
  }
  refined_tree
}

find_refined_node <- function(
  tree,
  selected_tips
){
  
  selected_tips <- intersect(
    selected_tips,
    tree$tip.label
  )
  
  if(length(selected_tips) < 2){
    return(NA_integer_)
  }
  
  ape::getMRCA(
    tree,
    selected_tips
  )
}


# 4. Apply refinement to all donors
ss <- sub("\\.length\\.edit\\.nwk$","",basename(tree_files))
all_analysis <- setNames(vector("list",length(ss)),ss)
all_trees <- setNames(vector("list",length(ss)),ss)
result_list <- vector("list",length(ss))
refined_trees <- setNames(vector("list",length(ss)),ss)

for(i in seq_along(ss)){
  s <- ss[i]
  message("Analyzing ",s," (",i,"/",length(ss),")")
  nwk <- paste(readLines(tree_files[i],warn=FALSE),collapse="")
  tree <- ape::read.tree(text=nwk)
  clone_cols <- tree$tip.label
  stopifnot(all(clone_cols %in% names(early_seg_total)),s %in% names(seg_candi))

  # Exclude regions observed exclusively in the donor being analyzed.
  seg_filtout <- seg_candi %>% filter(n_patients_covered==1,!is.na(.data[[s]]),.data[[s]]!=0) %>% select(chrom,start,end) %>% distinct() %>% mutate(filt=TRUE)
  early_seg <- early_seg_total %>% select(chrom,start,end,region_id,segment_id,all_of(clone_cols)) %>% anti_join(seg_filtout,by=c("chrom","start","end"))

  analysis <- analyze_all_multifurcations(tree=tree,nwk=nwk,early_seg=early_seg,clone_cols=clone_cols,min_frac=0.8,outside_frac=0.8)
  all_analysis[[s]] <- analysis
  all_trees[[s]] <- tree
  refined_trees[[s]] <- make_refined_tree(tree=tree,analysis=analysis,min_shared_mutation=5,min_confidence=0.7)

  result_list[[i]] <- bind_rows(lapply(analysis,function(x){
    selected <- !is.na(x$best_candidate) && !is.na(x$shared_mutation) && x$shared_mutation>5 && !is.na(x$confidence) && x$confidence>=0.7
    data.frame(Donor=s,Polytomy=paste0("M",x$case_n),Shared_somatic_mutations=x$shared_mutation,
               No_candidate_topologies=length(x$candidates),Selected_topology=if(is.na(x$best_candidate)) NA_character_ else x$candidates[[x$best_candidate]]$after_string,
               Top_EEC_support=x$top_count,Second_EEC_support=x$second_count,Prediction_confidence=x$confidence,
               Minimum_parsimony=x$min_parsimony,Refined=ifelse(selected,"Yes","No"))
  }))
}

# 5. Inspect application results
result_df <- bind_rows(result_list)
print(result_df)
print(result_df %>% filter(!is.na(Minimum_parsimony)) %>% summarise(n=n()))
print(result_df %>% filter(Shared_somatic_mutations>5,Prediction_confidence>=0.7) %>% count(Donor,name="n_refined"))

# 6. Inspect HC15 application (Figure 5g example)
hc15_result <- result_df %>% filter(Donor=="HC15")
print(hc15_result)
if("HC15" %in% names(all_trees)){
  old_par <- par(no.readonly=TRUE)
  par(mfrow=c(1,2),mar=c(2,2,3,2))
  ape::plot.phylo(all_trees[["HC15"]],use.edge.length=FALSE,main="HC15: original",cex=0.7)
  ape::plot.phylo(refined_trees[["HC15"]],use.edge.length=FALSE,main="HC15: EEC-refined",cex=0.7)
  par(old_par)
}
