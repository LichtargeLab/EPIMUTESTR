# 
# Preprocessing the TCGA cancer MC3 
# author: Saeid Parvandeh - June 2020
# 
library(parallel)
library(foreach)
library(doParallel)

# Assign cores
no_cores <- detectCores() - 2

# Initiate cluster
cl <- makeCluster(no_cores)
registerDoParallel(cl)

cancers_list <- c("ACC",  "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", 
                  "KIRC", "KIRP", "LAML", "LGG",  "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", 
                  "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "PAN")

foreach (ca=1:length(cancers_list)) %dopar% {
  # read cancer MC3 maf file
  setwd("path_to_the_data_direcotry")
  ca_maf <- read.delim(paste0(cancers_list[ca], "_EA.maf"), fill = T, header = T, colClasses = "character")

  ca_maf_fltr <- ca_maf[, c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode", "Chromosome",
                                    "Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele1",
                                    "Tumor_Seq_Allele2","dbSNP_RS", "GENE", "NM", "SUB", "ACTION")]

  ca_maf_fltr <- ca_maf_fltr[-which(ca_maf_fltr$ACTION == ""), ]
  ca_maf_fltr <- ca_maf_fltr[-which(ca_maf_fltr$ACTION == "-"), ]
  ca_maf_fltr$ACTION[which(ca_maf_fltr$ACTION %in% c("no_STOP", "STOP_gain", "STOP_loss", "START_loss", "STOP"))] <- "100"

  # Addressing the iso-form issue by taking a transcript with maximum EA score
  for (l in 1:nrow(ca_maf_fltr)){
    EA <- as.numeric(unlist(strsplit(ca_maf_fltr$ACTION[l], ";")))
    if (any(EA, na.rm = TRUE)){
      max_nm <- which.max(EA)
      iso_ea <- EA
      iso_sub <- unlist(strsplit(ca_maf_fltr$SUB[l], ";"))
      iso_nm <- unlist(strsplit(ca_maf_fltr$NM[l], ";"))
      if (length(iso_ea) == 1 & length(iso_sub) == 1 & length(iso_nm) == 1){
        ca_maf_fltr$ACTION[l] <- iso_ea
        ca_maf_fltr$SUB[l] <- iso_sub
        ca_maf_fltr$NM[l] <- iso_nm
      } else {
        ca_maf_fltr$ACTION[l] <- iso_ea[max_nm]
        ca_maf_fltr$SUB[l] <- ifelse(length(iso_sub)==length(iso_ea), iso_sub[max_nm], iso_sub)
        ca_maf_fltr$NM[l] <- ifelse(length(iso_nm)==length(iso_ea), iso_nm[max_nm], iso_nm)
      }
    }
  }

  # write all mutations
  write.table(ca_maf_fltr, paste0(cancers_list[ca], "_mutations.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
  ca_maf_fltr <- read.table(paste0(cancers_list[ca], "_mutations.tsv"), header = TRUE, colClasses = "character")
  
  # Calculate mutations signature
  cancer_signature <- ca_maf_fltr[!ca_maf_fltr$Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins"), ]
  cancer_signature$RefAlt <- paste(cancer_signature$Reference_Allele,cancer_signature$Tumor_Seq_Allele2, sep = "")
  allele_table <- table(cancer_signature$RefAlt)
  allele_table <- allele_table[!grepl("-", names(allele_table))]
  mut_signature <- NULL
  for (allele in names(allele_table)){
    mut_signature <- rbind(mut_signature, c(allele, allele_table[allele]/nrow(cancer_signature)))
  }
  write.table(mut_signature, paste0(cancers_list[ca], "_signature.csv"), quote = F, row.names = F, col.names = F, sep = ",")
  
  ca_maf_fltr <- ca_maf_fltr[ca_maf_fltr$Variant_Classification %in% c("Missense_Mutation",
                                                                       "Nonsense_Mutation",
                                                                       "Frame_Shift_Del",
                                                                       "Frame_Shift_Ins",
                                                                       "Nonstop_Mutation",
                                                                       "Splice_Site",
                                                                       "Translation_Start_Site"), ]
  
  # filter samples based on Bailey's paper
  # Exclude hypermutated samples - greater than 1,000 mutations and interqurtile 
  setwd("path_to_the_data_direcotry")
  hypermutated_samples <- read.table("hypermutated_samples.txt", header = TRUE, colClasses = "character")
  HM_samples <- hypermutated_samples[which(hypermutated_samples$CODE==cancers_list[ca]), "Tumor_Sample_Barcode"]
  if (length(HM_samples)!=0){
    ca_maf_fltr <- ca_maf_fltr[!ca_maf_fltr$Tumor_Sample_Barcode%in%HM_samples, ]}
  # Exclude pathology review samples
  pathology_review_samples <- read.delim("tumor_samples_excluded_pathology_review.txt", header = TRUE, colClasses = "character")
  PR_samples <- pathology_review_samples[which(pathology_review_samples$acronym==cancers_list[ca]), "sample_barcode"]
  ca_maf_fltr$Tumor_Sample_12digits <- substr(ca_maf_fltr$Tumor_Sample_Barcode, 1, 16)
  if (length(PR_samples)!=0) {
    ca_maf_fltr <- ca_maf_fltr[!ca_maf_fltr$Tumor_Sample_12digits%in%PR_samples, ]}
  

  ca_maf_fltr <- ca_maf_fltr[!is.na(as.numeric(ca_maf_fltr$ACTION)), ]
  setwd("path_to_the_results_direcotry")
  write.table(ca_maf_fltr, paste0(cancers_list[ca], "_maf_cleaned.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")
}

# Finish
stopCluster(cl)
