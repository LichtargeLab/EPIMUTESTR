#
# Criteria for success 
#

##########################################################################################
#### ------- Create 428 imputed “gold standard” putative cancer driver genes -------- ####
library(readxl)
setwd("path_to_the_results_direcotry")
# Read EPIMUTESTR final gene list
EPI_finalgenes <- read.table("EPIMUTESTR_FinalGenes_list.tsv", header = TRUE, sep = "\t", colClasses = "character")
EPI_finalgenes <- EPI_finalgenes$Gene

# Read 10 sources of cancer driver genes
gs_table <- read_excel("gs_table.xlsx")
gs_genes <- as.vector(as.matrix(gs_table))
gs_genes <- na.omit(gs_genes)
gs_common_genes <- names(which(table(gs_genes)>=2)) # genes that are co-occur at least twice

EPIMUTESTR_novel_genes <- setdiff(EPI_finalgenes, gs_common_genes)
write.table(data.frame(Gene=EPIMUTESTR_novel_genes), "EPIMUTESTR_NovelGenes.tsv", row.names = F, quote = F, sep = "\t")

# create a table to show the gs agreement
gs_overlap_mat <- matrix(NA, length(gs_common_genes), 12)
rownames(gs_overlap_mat) <- gs_common_genes
colnames(gs_overlap_mat) <- c("dNdScv", "Wheeler_2018", "TUSON", "MUSIC", "MutSigCV", "MutSig2CV", 
                              "Bailey_2018", "2020+", "2020", "COSMIC", "agreement", "EPIMUTESTR")


for (i in 1:length(gs_common_genes)){
  for (j in 1:10){
    gs_overlap_mat[i, j] <- ifelse(is.element(gs_common_genes[i], na.omit(as.data.frame(gs_table)[, j])), 1, NA)
  }
  gs_overlap_mat[i, "agreement"] <- sum(as.numeric(gs_overlap_mat[i, ]), na.rm = TRUE)
  gs_overlap_mat[i, "EPIMUTESTR"] <- ifelse(is.element(gs_common_genes[i], EPI_finalgenes), 1, NA)
}


gs_overlap_mat[is.na(gs_overlap_mat)] <- ""
write.table(gs_overlap_mat, "GS_overlap_matrix.tsv", quote = FALSE, sep = "\t")

########################################################################################
#### ------- Enrich all 33 cancer types plus pan-cancer using DOSE library -------- ####
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("clusterProfiler")
# browseVignettes("clusterProfiler")

library(clusterProfiler)

# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

# BiocManager::install("DOSE")
# browseVignettes("DOSE")

library(DOSE)

library(rJava)
library(xlsx)

cancers_list <- c("ACC",  "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", 
                  "KIRC", "KIRP", "LAML", "LGG",  "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", 
                  "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "PAN")

setwd("path_to_the_results_direcotry")

all_genelists <- read.table("EPIMUTESTR_SigGenes.tsv", header = T, sep = "\t", colClasses = "character",na.strings=c("",NA)) 
EPIMUTESTR_final_genes <- read.table("EPIMUTESTR_FinalGenes_list.tsv", header = T, sep = "\t", colClasses = "character",na.strings=c("",NA)) 
EPIMUTESTR_novel_genes <- read.table("EPIMUTESTR_NovelGenes.tsv", header = T, sep = "\t", colClasses = "character")


ncg_genes_list <- NULL
novel_genes_list <- NULL
for (ca in cancers_list){
  ca_genelist <- na.omit(all_genelists[, ca])
  # In order to use enrichDGN library we need ENTREZID from clusterProfiler library
  # Get Genes Entrez ID
  GeneEID <- bitr(ca_genelist, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  # Use network of cancer genes to enrich cancer genes 
  ncg <- enrichNCG(GeneEID$ENTREZID, readable = TRUE)
  ncg_result <- ncg@result[, c("geneID", "Count", "Description", "pvalue", "qvalue")]
  
  ncg_genes <- unique(unlist(strsplit(ncg_result$geneID, split = "/")))
  ncg_genes_list <- c(ncg_genes_list, ncg_genes)
  
  novel_genes <- setdiff(ca_genelist, ncg_genes)
  novel_genes_list <- c(novel_genes_list, novel_genes)
  
  write.csv(ncg_result, paste0(ca, "_enrich_NCG.csv"), col.names = T, row.names = F)
  write.csv(ncg_result, paste0(ca, "_enrich_NCG.csv"))
  
}

write.table(unique(ncg_genes_list), "EPIMUTESTR_genes_in_NCG.tsv", row.names = F, col.names = F, append = TRUE, sep = "\n")
write.table(unique(novel_genes_list), "EPIMUTESTR_novel_genes_NCG.tsv", row.names = F, col.names = F, append = TRUE, sep = "\n")


###################################################
#### ------- PubMed literature search -------- ####
setwd("path_to_the_results_directory")
novel_genes <- read.table("EPIMUTESTR_NovelGenes.tsv", header = TRUE, colClasses = "character", sep = "\t")
EPIMUTESTR_final_genes <- read.table("EPIMUTESTR_FinalGenes_list.tsv", header = T, sep = "\t", colClasses = "character",na.strings=c("",NA)) 
entire_genes <- read.csv("whole_genome.csv", header = FALSE, colClasses = "character")

# First, PubMed_search_all.py has to run
setwd("path_to_the_pubmed_search_query_output_direcotry")
pubmed_search <- NULL
for (gene in EPIMUTESTR_final_genes$Gene){
  query <- read.table(paste0(gene, ".tsv"), header = TRUE)
  assoc <- ifelse(is.numeric(query$No.literature_w_Cancer), query$No.literature_w_Cancer, 0)
  greater10 <- ifelse(assoc>=10, 1, 0)
  between1_10 <- ifelse((assoc<10 & assoc>=1), 1, 0)
  pubmed_search <- rbind(pubmed_search, c(gene, assoc, greater10, between1_10))
}

colnames(pubmed_search) <- c("Gene", "PubMed_query", "In_PubMed>10", "In_PubMed_1_10")

setwd("path_to_the_results_direcotry")
write.table(pubmed_search, "pubmed_FinalGenes.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

setwd("path_to_the_pubmed_search_query_output_direcotry")
pubmed_search <- NULL
for (gene in entire_genes$V1){
  query <- read.table(paste0(gene, ".tsv"), header = TRUE)
  assoc <- ifelse(is.numeric(query$No.literature_w_Cancer), query$No.literature_w_Cancer, 0)
  greater10 <- ifelse(assoc>=10, 1, 0)
  between1_10 <- ifelse((assoc<10 & assoc>=1), 1, 0)
  pubmed_search <- rbind(pubmed_search, c(gene, assoc, greater10, between1_10))
}

colnames(pubmed_search) <- c("Gene", "PubMed_query", "In_PubMed>10", "In_PubMed_1_10")

setwd("path_to_the_results_direcotry")
write.table(pubmed_search, "pubmed_wholegenome.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

####################################################################
#### ------- Robustness/Power analysis by downsampling -------- ####
library(parallel)
library(foreach)
library(doParallel)

# Assign cores
no_cores <- detectCores() - 2

# Initiate cluster
cl <- makeCluster(no_cores, outfile="")
registerDoParallel(cl)

setwd("path_to_the_data_direcotry")
# Percentage of missing values for each gene
pMiss <- function(x){sum(is.na(x))/length(x)*100}

# Read transcripts that are invalid
trans <- read.table("invalid_trans", header = FALSE, sep = "\t")
invalid_trans <- as.vector(as.character(trans$V1))

# Read all transcript in one file
valid_trans <- read.table("all_nm_ids.txt", header = FALSE, colClasses = "character")
all_trans <- as.vector(as.character(valid_trans$V1))

# Giving a tumor name to begin the analysis (e.g. LUAD)
setwd("path_to_the_results_direcotry")
ca_maf <- read.table("LUAD_maf_cleaned.tsv", header = TRUE,  sep = "\t", colClasses = "character")
signature_file <- read.csv("LUAD_signature.csv", header = FALSE, colClasses = "character")
ca_mat <- get(load("LUAD_max_feat_mat.RData"))

subsam_prob <- seq(0.05, 1, 0.05)
foreach (prob=1:length(subsam_prob)) %dopar%{
  epi_genes_list <- list()
  for (i in 1:10){
    replicate_genes <- NULL
    for (re in 1:10){
      # Create control matrix
      ctrl_mat <- matrix(NA, nrow(ca_mat), ncol(ca_mat))
      rownames(ctrl_mat) <- rownames(ca_mat)
      colnames(ctrl_mat) <- colnames(ca_mat)
      
      
      rand_muts <- NULL
      genes <- colnames(ca_mat)
      for (gene in genes) {
        Transcript <- as.character(unlist(strsplit(ca_maf$NM[which(ca_maf$GENE == gene)], ";"))[1])
        
        setwd("path_to_AllNucelotideChanges_file_direcotry")
        if (length(Transcript) == 0 || is.element(Transcript, invalid_trans) || !is.element(Transcript, all_trans)){
          Trans_gene <- read.table(paste0(sample(all_trans, 1)), header = FALSE, sep = "\t", colClasses = "character")
        } else {
          Trans_gene <- read.table(paste0(Transcript, ".random"), header = FALSE, sep = "\t", colClasses = "character")
        }
        
        # Accounting for mutation signature
        # There is no FS and Indel mutations in the random files, so we go with all
        Trans_gene$RefAlt <- gsub('[[:digit:]]+', '', Trans_gene$V5)
        Trans_lines <- NULL
        for (allele in unique(Trans_gene$RefAlt)){
          mut_signature <- signature_file[which(signature_file$V1==allele), "V2"]
          mut_idx <- which(Trans_gene$RefAlt==allele)
          Trans_lines <- c(Trans_lines, sample(mut_idx, as.numeric(mut_signature)*nrow(Trans_gene), replace = T))
        }
        Trans_gene <- Trans_gene[Trans_lines, ]
        
        pool <- Trans_gene$V7[!Trans_gene$V7 %in% c("silent", "no_action", "STOP_LOSS")]
        while(length(pool)==0){
          Trans_gene <- read.table(paste0(sample(all_trans, 1)), header = FALSE, sep = "\t", colClasses = "character")
          pool <- Trans_gene$V7[!Trans_gene$V7 %in% c("silent", "no_action", "STOP_LOSS")]
        }
        
        num_muts <- length(which(!is.na(ca_mat[, gene])))
        if (length(pool) >= num_muts) {
          rand_muts <- as.numeric(as.character(sample(pool, num_muts)))
        } else {
          rand_muts <- as.numeric(as.character(sample(pool, num_muts, replace = TRUE)))
        }
        muts_idx <- which(!is.na(ca_mat[, gene]))
        ctrl_mat[muts_idx, gene] <- rand_muts
      }
      
      # Replace the control Sample IDs
      rownames(ctrl_mat) <- paste("Sample_", 1:nrow(ctrl_mat), sep = "")
      
      # Select the common genes
      common_genes <- intersect(colnames(ca_mat), colnames(ctrl_mat))
      case_fltr <- ca_mat[, common_genes]
      control_fltr <- ctrl_mat[, common_genes]
      
      # subsampling
      random_subsample <- sample(nrow(case_fltr), as.numeric(subsam_prob[prob])*nrow(case_fltr))
      case_fltr_subsampled <- case_fltr[random_subsample, common_genes]
      control_fltr_subsampled <- control_fltr[random_subsample, common_genes]
      
      # combine case and control
      cancers_combined <- rbind(case_fltr_subsampled, control_fltr_subsampled)
      
      # Create phenotype
      case_class <- rep(1, dim(case_fltr_subsampled)[1])
      control_class <- rep(0, dim(control_fltr_subsampled)[1])
      case_control_class <- c(case_class, control_class)
      
      # Relief feature selection
      fold_data <- data.frame(cancers_combined, class = as.factor(case_control_class))
      # Percentage of missing values for each gene
      pMiss <- function(x){sum(is.na(x))/length(x)*100}
      
      # removing columns with missing value
      fold_data_fltr <- fold_data[ ,which(apply(fold_data, 2, pMiss) != 100)]
      fold_data_fltr <- fold_data_fltr[which(apply(t(fold_data_fltr), 2, pMiss) != 100), ]
      # InfGain, MDL, Gini, Accuracy, DKMcost, MDLsmp
      genes_score <- CORElearn::attrEval("class", fold_data_fltr, estimator = "ReliefFequalK", 
                                         kNearestEqual = floor((dim(fold_data_fltr)[1]-1)*0.154))
      
      # find significant genes
      geneRank_fltr <- sort(genes_score[which(genes_score>0)], T)
      
      # Back-Transformation (square root)
      trans_fun <- function(x){x^(1/2)}
      
      geneRank_fltr_trans <- trans_fun(geneRank_fltr)
      
      # calculate the p-value for each gene
      Pvals <- NULL
      for (p in geneRank_fltr_trans){
        Pvals <- c(Pvals, pnorm(p, mean(geneRank_fltr_trans), sd(geneRank_fltr_trans), lower.tail = F))
      }
      
      Qvals <- p.adjust(Pvals, method = "BH")
      genes <- names(geneRank_fltr)[which(Qvals < 0.1)]
      
      replicate_genes <- c(replicate_genes, genes)
      
    }
    epi_genes_list[[i]] <- names(which(table(replicate_genes)>=5))
    
  }
  # Save scores
  setwd("path_to_the_results_direcotry")
  save(epi_genes_list, file = paste0("LUAD_EPI_genes_list_",subsam_prob[prob],".RData"))
}

# Finish
stopCluster(cl)


##############################################
##### dNdScv (Martincorena et al., 2017) #####
# -------------------------------------------#
# http://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/dNdScv.html
# library(devtools); install_github("im3sanger/dndscv")
library("dndscv")

nt = c("A","C","G","T")
subsam_prob <- seq(0.05, 1, 0.05)
dnds_gs_overlap <- list()
dnds_genes_lists <- list()
dnds_gs_pval <- NULL
setwd("path_to_the_data_direcotry")
mutations_file <- read.table("BLCA_mutations.tsv", header = T, sep = "\t", colClasses = "character")
load("BLCA_max_feat_mat.RData")
for (prob in subsam_prob) {
  dnds_genes_list <- list()
  for (k in 1:10){
    random_subsample <- sample(rownames(BLCA_max_feat_mat), round(prob*nrow(BLCA_max_feat_mat)))
    mutations_file_fltr <- mutations_file[which(mutations_file$Tumor_Sample_Barcode%in%random_subsample), ]
    mutations_file_fltr <- mutations_file_fltr[-which(!(mutations_file_fltr$Reference_Allele %in% nt)|
                                                        !(mutations_file_fltr$Tumor_Seq_Allele2 %in% nt)), ]
    
    mutations_file_fltr <- unique(mutations_file_fltr)
    dndsout = dndscv(mutations_file_fltr)#, refdb="output_refcds.rda")#, cv=NULL)
    dnds_res <- strsplit(dndsout$sel_cv[dndsout$sel_cv$qallsubs_cv<0.1, c("gene_name")], split =":")
    dnds_genes <- NULL
    for (i in 1:length(dnds_res)){
      dnds_genes <- c(dnds_genes, dnds_res[[i]][1])
    }
    dnds_genes_list[[k]] <- dnds_genes
  }
  # Save dnds genes
  setwd("path_to_the_results_direcotry")
  save(dnds_genes_list, file = paste0("LUAD_dnds_genes_list_",prob,".RData"))
}










