#
# Replicate analysis 10 times with different synthetic control
#
# Saeid Parvandeh June 2020
# -------------------------

library(parallel)
library(foreach)
library(doParallel)
library(Rcpp)

# Assign cores
no_cores <- detectCores() - 2

# Initiate cluster
cl <- makeCluster(no_cores, outfile="")
registerDoParallel(cl)

cancers_list <- c("ACC",  "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH",
                  "KIRC", "KIRP", "LAML", "LGG",  "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG",
                  "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "PAN")

foreach (ca=1:length(cancers_list)) %dopar% {
  # read cancer maf file
  setwd("path_to_the_data_directory")
  ca_maf <- read.table(paste0(cancers_list[ca], "_maf_cleaned.tsv"), header = TRUE, fill = TRUE, colClasses = "character")
  # Extract genes and samples id
  all_genes <- unique(ca_maf$GENE)
  all_samples <- unique(ca_maf$Tumor_Sample_Barcode)

  # Removing samples may cause problem
  samples_fltr <- NULL
  for (s in all_samples){
    num_sample_var <- length(which(ca_maf$Tumor_Sample_Barcode==s))
    if(num_sample_var > 1000){
      samples_fltr <- c(samples_fltr, s)
    }
  }

  # removing samples with high number of mutations with cutoff point of 1000
  if(!is.null(samples_fltr)){
    ca_maf_fltr <- ca_maf[!ca_maf$Tumor_Sample_Barcode %in% samples_fltr, ]
  } else {
    ca_maf_fltr <- ca_maf
  }
  dim(ca_maf_fltr)

  # update the sample ids
  all_samples <- unique(ca_maf_fltr$Tumor_Sample_Barcode)

  # Fill the feature matrix with nonsynonymous and stopgain EA scores
  #### using Rcpp ####
  Rcpp::sourceCpp("Create_design_matrix.cpp")
  ca_maf_fltr$ACTION <- as.numeric(ca_maf_fltr$ACTION)
  cancer_ft_mat <- DesignMatrix(ca_maf_fltr, all_genes, all_samples)
  rownames(cancer_ft_mat) <- all_samples
  colnames(cancer_ft_mat) <- all_genes

  # Percentage of missing values for each gene
  pMiss <- function(x){sum(is.na(x))/length(x)*100}

  # removing columns with missing value
  cancer_na.out <- cancer_ft_mat[ ,which(apply(cancer_ft_mat, 2, pMiss) != 100)]
  cancer_na.out <- cancer_na.out[which(apply(t(cancer_na.out), 2, pMiss) != 100), ]

  # Create control matrix
  control_ft_mat <- matrix(NA, nrow(cancer_na.out), ncol(cancer_na.out))
  rownames(control_ft_mat) <- rownames(cancer_na.out)
  colnames(control_ft_mat) <- colnames(cancer_na.out)

  # rename to associated cancer -- Max EA score
  assign(paste0(cancers_list[ca], "_max_feat_mat"), cancer_na.out)
  # save the matrix
  setwd("path_to_the_results_direcotry")
  save(list = noquote(paste0(cancers_list[ca], "_max_feat_mat")),
       file = paste0(cancers_list[ca], "_max_feat_mat.RData"))
  
  # Read transcripts that are invalid
  trans <- read.table("invalid_trans", header = FALSE, sep = "\t")
  invalid_trans <- as.vector(as.character(trans$V1))
  
  # Read all transcript in one file
  valid_trans <- read.table("all_nm_ids.txt", header = FALSE, colClasses = "character")
  all_trans <- as.vector(as.character(valid_trans$V1))

  # Create random control from all single nocleotide changes by incorporating molecular signature
  signature_file <- read.csv(paste0(cancers_list[ca], "_signature.csv"), header = F, colClasses = "character")
  rand_muts <- NULL
  genes <- colnames(cancer_na.out)
  for (gene in genes) {
    Transcript <- as.character(unlist(strsplit(ca_maf_fltr$NM[which(ca_maf_fltr$GENE == gene)], ";"))[1])
    
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

    num_muts <- length(which(!is.na(cancer_na.out[, gene])))
    if (length(pool) >= num_muts) {
      rand_muts <- as.numeric(as.character(sample(pool, num_muts)))
    } else {
      rand_muts <- as.numeric(as.character(sample(pool, num_muts, replace = TRUE)))
    }
    muts_idx <- which(!is.na(cancer_na.out[, gene]))
    control_ft_mat[muts_idx, gene] <- rand_muts
  }

  # Replace the control Sample IDs
  rownames(control_ft_mat) <- paste("Sample_", 1:nrow(control_ft_mat), sep = "")

  # rename to associated control -- Control EA score
  assign(paste0(cancers_list[ca], "_ctrl_feat_mat"), control_ft_mat)
  # save the matrix
  setwd("path_to_the_results_directory")
  save(list = noquote(paste0(cancers_list[ca], "_ctrl_feat_mat")),
       file = paste0(cancers_list[ca], "_ctrl_feat_mat.RData"))
  
  # Select the common genes
  common_genes <- intersect(colnames(cancer_na.out), colnames(control_ft_mat))
  case_fltr <- cancer_na.out[, common_genes]
  control_fltr <- control_ft_mat[, common_genes]
  
  # combine case and control
  cancers_combined <- rbind(case_fltr, control_fltr)
  
  # Create phenotype
  case_class <- rep(1, dim(case_fltr)[1])
  control_class <- rep(0, dim(control_fltr)[1])
  case_control_class <- c(case_class, control_class)
  
  # Relief feature selection
  fold_data <- data.frame(cancers_combined, class = as.factor(case_control_class))
  # Save scores
  assign(paste0(cancers_list[ca], "_scores_1"), genes_score)
  save(list = noquote(paste0(cancers_list[ca], "_scores_1")),
       file = paste0(cancers_list[ca], "_scores_1.RData"))
}

# Finish
stopCluster(cl)



