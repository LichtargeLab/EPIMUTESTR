#
# Get the significant genes from 10 times EPIMUTESTR replications
#
rm(list = ls())
library(dplyr)

# Function for cbind data frames with different rows
cbindPad <- function(...){
  args <- list(...)
  n <- sapply(args,nrow)
  mx <- max(n)
  pad <- function(x, mx){
    if (nrow(x) < mx){
      nms <- colnames(x)
      padTemp <- matrix(NA, mx - nrow(x), ncol(x))
      colnames(padTemp) <- nms
      if (ncol(x)==0) {
        return(padTemp)
      } else {
        return(rbind(x,padTemp))
      }
    }
    else{
      return(x)
    }
  }
  rs <- lapply(args,pad,mx)
  return(do.call(cbind,rs))
}

setwd("path_to_the_direcotry")
cosmic_file <- read.table("Census_allMon Mar 29 18_06_49 2021.tsv", sep = "\t", header = T, colClasses = "character")
cosmic_gene <- cosmic_file$Gene.Symbol

cancers_list <- c("ACC",  "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", 
                  "KIRC", "KIRP", "LAML", "LGG",  "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", 
                  "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "PAN")

setwd("path_to_the_direcotry")

cancers_SigGenes <- data.frame()
all_genes <- NULL
for (ca in cancers_list){
  Genes_reps <- NULL
  for (re in 1:10){
    # My genes
    load(paste0(ca, "_scores_", re, ".RData"))
    my_gene_rank <- get(paste0(ca, "_scores_", re))
    geneRank_fltr <- sort(my_gene_rank[which(my_gene_rank>0)], T)
    
    # Back-Transformation (square root)
    trans_fun <- function(x){x^(1/2)}
    
    geneRank_fltr_trans <- trans_fun(geneRank_fltr)
    
    # calculate the P value for each gene
    Pvals <- NULL
    for (p in geneRank_fltr_trans){
      Pvals <- c(Pvals, pnorm(p, mean(geneRank_fltr_trans), sd(geneRank_fltr_trans), lower.tail = F))
    }
    
    Qvals <- p.adjust(Pvals, method = "BH")
    g <- names(geneRank_fltr)[which(Qvals < 0.1)]
    p <- Pvals[which(Qvals < 0.1)]
    q <- Qvals[which(Qvals < 0.1)]
    
    Genes_reps <- rbind(Genes_reps, data.frame(gene=g, pvalue=p, qvalue=q))
  }
  
  # filter genes with frequency of at least 2
  Genes_reps$qvalue <- as.numeric(as.character(Genes_reps$qvalue))
  Genes_reps_fltr <- Genes_reps %>%
    group_by(gene) %>%
    mutate(freq = n()) %>% 
    filter(freq >= 5) %>%
    filter(qvalue == min(qvalue)) %>%
    ungroup()
  Genes_reps_fltr <- as.data.frame(Genes_reps_fltr)
  # sort based on average qvalue
  Genes_reps_sorted <- Genes_reps_fltr[order(Genes_reps_fltr$qvalue), ]
  
  cancers_SigGenes <- cbindPad(cancers_SigGenes, Genes_reps_sorted)
  
  if (ca != "PAN"){
    all_genes <- rbind(all_genes, data.frame(gene=as.character(Genes_reps_sorted$gene), qvalue=Genes_reps_sorted$qvalue))
  } else {
    pan_genes <- data.frame(gene=as.character(Genes_reps_sorted$gene), qvalue=Genes_reps_sorted$qvalue)
  }
}

# PanCancer genes standalone
length(setdiff(unique(pan_genes$gene), unique(all_genes$gene)))

idx <- seq(1, 136, by = 4)
names(cancers_SigGenes)[idx] <- cancers_list
rownames(cancers_SigGenes) <- NULL


colnames(all_genes) <- c("gene", "qvalue")
all_genes$qvalue <- as.numeric(all_genes$qvalue)
all_genes_fltr <- all_genes %>%
  group_by(gene) %>%
  filter(qvalue == min(qvalue)) %>%
  ungroup()
all_genes_sorted <- unique(all_genes_fltr[order(all_genes_fltr$qvalue), ])

final_genes <- c(as.character(all_genes_sorted$gene), setdiff(pan_genes$gene, all_genes_sorted$gene))

setwd("path_to_the_directory")
write.table(cancers_SigGenes, file = "EPIMUTESTR_SigGenes.tsv", quote = F, row.names = F, sep = "\t", na = "")

write.table(data.frame(Gene=final_genes), file = "EPIMUTESTR_FinalGenes_list.tsv", quote = F, row.names = F, sep = "\t", na = "")


