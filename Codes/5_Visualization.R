#
# Seris offplots for EPIMUTESTR results 
# Author: Saeid Parvandeh Sep, 2021

###########################################################################
####### ---- Evaluate the cancer significant genes with Depmap ----- ######
# First, DepMap analysis has to run to generate AllFigureData.csv
setwd("path_to_the_results_direcotry")
depmap_data <- read.csv("AllFigureData.csv")
cce_OG  <- depmap_data$Ceres.Score[which(depmap_data$Mut..Group=="Moderate EA"&depmap_data$GeneSet=="CCE")]
cce_other  <- depmap_data$Ceres.Score[which(depmap_data$Mut..Group=="Other"&depmap_data$GeneSet=="CCE")]
ce_OG   <- depmap_data$Ceres.Score[which(depmap_data$Mut..Group=="Moderate EA"&depmap_data$GeneSet=="CE")]
ce_other   <- depmap_data$Ceres.Score[which(depmap_data$Mut..Group=="Other"&depmap_data$GeneSet=="CE")]
ne_OG   <- depmap_data$Ceres.Score[which(depmap_data$Mut..Group=="Moderate EA"&depmap_data$GeneSet=="NE")]
ne_other   <- depmap_data$Ceres.Score[which(depmap_data$Mut..Group=="Other"&depmap_data$GeneSet=="NE")]
cos_OG  <- depmap_data$Ceres.Score[which(depmap_data$Mut..Group=="Moderate EA"&depmap_data$GeneSet=="COSMIC GS")]
cos_other  <- depmap_data$Ceres.Score[which(depmap_data$Mut..Group=="Other"&depmap_data$GeneSet=="COSMIC GS")]
rand_OG <- depmap_data$Ceres.Score[which(depmap_data$Mut..Group=="Moderate EA"&depmap_data$GeneSet=="Random Genes")]
rand_other <- depmap_data$Ceres.Score[which(depmap_data$Mut..Group=="Other"&depmap_data$GeneSet=="Random Genes")]
EPI_indv_OG <- depmap_data$Ceres.Score[which(depmap_data$Mut..Group=="Moderate EA"&depmap_data$GeneSet=="ML Individual Cancer (Pooled) Oncogene")]
EPI_indv_others <- depmap_data$Ceres.Score[which(depmap_data$Mut..Group=="Other"&depmap_data$GeneSet=="ML Individual Cancer (Pooled) TumorSupressors")]
EPI_pan_OG <- depmap_data$Ceres.Score[which(depmap_data$Mut..Group=="Moderate EA"&depmap_data$GeneSet=="ML PAN (All) Cancer Oncogenes")]
EPI_pan_other <-depmap_data$Ceres.Score[which(depmap_data$Mut..Group=="Other"&depmap_data$GeneSet=="ML Individual Cancer (Pooled) TumorSupressors")]

Ceres_Score <- c(matrix(ne_OG),
                 matrix(rand_OG),
                 matrix(ce_OG),
                 matrix(cos_OG), 
                 matrix(EPI_indv_OG), 
                 matrix(EPI_pan_OG))

Gene_Set <- c(rep("NE", length(c(ne_OG))), rep("Random", length(c(rand_OG))), 
              rep("CE", length(c(ce_OG))), rep("COSMIC", length(c(cos_OG))),
              rep("EPI. Indiv.", length(c(EPI_indv_OG))),
              rep("EPI. PAN", length(c(EPI_pan_OG))))

Category <- factor(c(rep("OncoGenes", length(ne_OG)), 
                     rep("OncoGenes", length(rand_OG)), 
                     rep("OncoGenes", length(ce_OG)), 
                     rep("OncoGenes", length(cos_OG)), 
                     rep("OncoGenes", length(EPI_indv_OG)), 
                     rep("OncoGenes", length(EPI_pan_OG))), 
                   levels = c("OncoGenes")) 

boxplot_df <- data.frame(Gene_Set, Category, Ceres_Score)


library(ggplot2)
library(ggpubr)

# Perorm pairwise comparisons
compare_means(Ceres_Score ~ Gene_Set,  data = boxplot_df,method = "t.test")
# Visualize: Specify the comparisons you want
my_comparisons <- list( c("Random", "CE"), c("COSMIC", "EPI. Indiv."), c("Random", "COSMIC"),  
                        c("Random", "EPI. Indiv."), c("Random", "EPI. PAN") )
ggboxplot(boxplot_df, x = "Gene_Set", y = "Ceres_Score",
          color = "Gene_Set", palette = "jco")+ 
  geom_hline(yintercept=-0.5, linetype='dashed', col = 'red', lwd=1)+
  labs(title="DepMap - Large scale screens oncogenes comparison") +
  theme(legend.position = "none") +
  scale_y_continuous(breaks=seq(-2.5, 3, 0.5)) +
  stat_compare_means(method = "t.test", method.args = list(alternative="greater"), comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 4)     # Add global p-value

######################################################################################
####### ---- Plot the robustness comparison between EPIMUTESTR and dNdScv ----- ######
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


ca <- "COAD" # BLCA, COAD, HNSC, LUAD, UCEC
setwd("path_to_the_data_direcotry")
known_genes_file <- read.table("gs_genes.tsv", header = F, colClasses = "character")
known_genes <- known_genes_file$V1
subsam_prob <- seq(0.05, 1, 0.05)
overlap_mean <- NULL
overlap_sd <- NULL
subsam_num <- NULL
overlap_genes_df <- data.frame()
ca_mat <- get(load(paste0(ca, "_max_feat_mat.RData")))
setwd("path_to_the_results_direcotry")
for (prob in subsam_prob){
  robust_genes <- get(load(paste0(ca, "_EPI_genes_list_", prob, ".RData")))
  subsam_num <- c(subsam_num,  round(prob*(nrow(ca_mat)*2)))
  overlap_num <- NULL
  overlap_genes <- NULL
  for (i in 1:10){
    overlap_num <- c(overlap_num, length(intersect(robust_genes[[i]], known_genes)))
    overlap_genes <- c(overlap_genes, robust_genes[[i]])
  }
  overlap_genes_df <- cbindPad(overlap_genes_df, data.frame(sort(table(overlap_genes), T)))
  overlap_mean <- c(overlap_mean, mean(overlap_num))
  overlap_sd <- c(overlap_sd, sd(overlap_num))
}
# write robustness genes into file
names(overlap_genes_df)[seq(1, 40, by = 2)] <- subsam_num
# write.csv(overlap_genes_df, paste0(ca, "_robustness.csv"), quote = F, row.names = F, na = "")
plot(subsam_num, overlap_mean, xaxt="n", col = "blue", type="o", pch=19, lty=1, cex.lab = 1.3, cex.main = 2,
     ylim=range(c(overlap_mean-overlap_sd), (overlap_mean+overlap_sd)), cex.axis = 1.3,
     main = paste0(ca),xlab = "No. of samples", ylab = "")#"No. significant genes")
# hack: we draw arrows but with very special “arrowheads”
arrows(subsam_num, overlap_mean-overlap_sd, subsam_num, overlap_mean+overlap_sd,
       length=0.05, angle=90, code=3)
axis(1, at = subsam_num, las=2, cex.axis = 1.2)

subsam_prob <- seq(0.05, 1, 0.05)
overlap_mean <- NULL
overlap_sd <- NULL
subsam_num <- NULL
for (prob in subsam_prob){
  robust_genes <- get(load(paste0("dnds_genes_list_", prob, ".RData")))
  subsam_num <- c(subsam_num,  round(prob*(nrow(ca_mat)*2)))
  overlap_num <- NULL
  for (i in 1:10){
    overlap_num <- c(overlap_num, length(intersect(robust_genes[[i]], known_genes)))
  }
  overlap_mean <- c(overlap_mean, mean(overlap_num))
  overlap_sd <- c(overlap_sd, sd(overlap_num))
}
points(subsam_num, overlap_mean, col="green", pch=17)
lines(subsam_num, overlap_mean, col="green",lty=2)
arrows(subsam_num, overlap_mean-overlap_sd, subsam_num, overlap_mean+overlap_sd,
       length=0.05, angle=90, code=3)
legend(100,25,legend=c("EPIMUTESTR","dNdScv"), col=c("blue","green"), pch=c(19,17),lty=c(1,2), ncol=1, cex = 1.3)


#####################################################################
###### ---- Plot the genes association in NCG and PubMed ----- ######
setwd("path_to_the_results_direcotry")
# EPIMUTESTR results
final_genelist <- read.table("EPIMUTESTR_FinalGenes_list.tsv", header = TRUE, colClasses = "character")
EPI_finalgenes <- final_genelist$Gene
EPI_novel_genes <- read.table("EPIMUTESTR_NovelGenes.tsv", sep = "\t", header = T, colClasses = "character")
# NCG results
ncg_genes <- read.table("EPIMUTESTR_genes_in_NCG.tsv", sep = "\t", colClasses = "character")
ncg_genes <- unique(ncg_genes$V1)
# PubMed query results
pubmed_entire_genes <- read.table("pubmed_wholegenome.tsv", sep = "\t", header = TRUE, colClasses = "character")
pubmed_final_genes <- read.table("pubmed_FinalGenes.tsv", sep = "\t", header = TRUE, colClasses = "character")
pubmed_genes <- pubmed_final_genes$Gene[which(pubmed_final_genes$In_PubMed.10==1)] # | pubmed_novel_genes$In_PubMed.10==1

# Percent stacked barchart
library(ggplot2)
Category <- c(rep(c("EPI. genes-NCG", "All genes-NCG"), 2), rep(c("EPI. genes-PubMed", "All genes-PubMed"), 3))
Groups <- c(rep("In NCG", 2), rep("Not in NCG", 2), rep("In PubMed > 10", 2), rep("In PubMed 1-10", 2), rep("Not in PubMed", 2))
Values <- c(length(intersect(EPI_finalgenes, ncg_genes)), 2372, 
            length(EPI_finalgenes) - length(intersect(EPI_finalgenes, ncg_genes)),
            18696 - 2372, 
            length(which(pubmed_final_genes$In_PubMed.10==1)), length(which(pubmed_entire_genes$In_PubMed.10==1)),
            length(which(pubmed_final_genes$In_PubMed_1_10==1)), length(which(pubmed_entire_genes$In_PubMed_1_10==1)),
            length(EPI_finalgenes) - sum(length(which(pubmed_final_genes$In_PubMed.10==1)),length(which(pubmed_final_genes$In_PubMed_1_10==1))),
            length(pubmed_entire_genes$Gene) - sum(length(which(pubmed_entire_genes$In_PubMed.10==1)), length(which(pubmed_final_genes$In_PubMed_1_10==1))))
barchat_data <- data.frame(Category, Groups, Values)
barchat_data$Category <- factor(barchat_data$Category, 
                                levels=c("EPI. genes-NCG", "All genes-NCG", "EPI. genes-PubMed", "All genes-PubMed"))
a<-barchat_data$Values[barchat_data$Category%in%"EPI. genes-NCG"]
b<-barchat_data$Values[barchat_data$Category%in%"All genes-NCG"]
c<-barchat_data$Values[barchat_data$Category%in%"EPI. genes-PubMed"]
d<-barchat_data$Values[barchat_data$Category%in%"All genes-PubMed"]
barchat_data$Fraction <- c(a[1]/sum(a), b[1]/sum(b), a[2]/sum(a), b[2]/sum(b), c[1]/sum(c), d[1]/sum(d), c[2]/sum(c), d[2]/sum(d), c[3]/sum(c), d[3]/sum(d))

ggplot(barchat_data, aes(fill=Groups, y=Fraction, x=Category)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity") +
  scale_fill_manual(values = c("grey30", "grey20", "grey30", "grey60", "grey60")) +
  geom_text(aes(y=Fraction, label=paste0(Values)), position = position_stack(vjust = 0.5, reverse = TRUE), color="white", size=6)+
  annotate("text", x = 1.5, y = 1.05, label = "p < 3.8e-116") +
  annotate("text", x = 3.5, y = 1.05, label = "p < 1.9e-43") +
  theme_bw() +  theme(legend.title = element_blank(), 
                      axis.text=element_text(size=12), 
                      legend.text=element_text(size=12),
                      axis.title.x = element_text(size = 16),
                      axis.title.y = element_text(size = 16)) 


########################################################
###### ---- Pie chart for genes association ----- ######

EPI_support_mat <- data.frame(matrix(NA, nrow(EPI_novel_genes), 5))
rownames(EPI_support_mat) <- EPI_novel_genes$Gene
colnames(EPI_support_mat) <- c("PubMed", "NCG", "two_supports", "one_support","novel")


supports <- c("PubMed", "NCG")
for (i in 1:nrow(EPI_novel_genes)){
  for (support in supports){
    if (support == "PubMed") {
      EPI_support_mat[i, "PubMed"] <- ifelse(is.element(EPI_novel_genes$Gene[i], pubmed_genes), 1, NA)
    } else if (support == "NCG") {
      EPI_support_mat[i, "NCG"] <- ifelse(is.element(EPI_novel_genes$Gene[i], ncg_genes), 1, NA)
    }
  }
}

support_vec <- rowSums(EPI_support_mat, na.rm = TRUE)

for (i in 1:nrow(EPI_novel_genes)){
  if(support_vec[i] == 0) {
    EPI_support_mat[i, "novel"] <- 1
  } else if(support_vec[i] == 1) {
    EPI_support_mat[i, "one_support"] <- 1
  } else if(support_vec[i] == 2) {
    EPI_support_mat[i, "two_supports"] <- 1
  }
}

EPI_overlap_matrix <- read.table("EPIMUTESTR_overlap_matrix.tsv", sep = "\t", header = T, colClasses = "character",
                                 row.names = 1)
EPI_overlap_matrix[EPI_overlap_matrix==""] <- NA
pie_df <- data.frame(FR=c(length(na.omit(EPI_overlap_matrix$all_agreement)), 
                          length(na.omit(EPI_overlap_matrix$nine_agreement)),
                          length(na.omit(EPI_overlap_matrix$eight_agreemnet)),
                          length(na.omit(EPI_overlap_matrix$seven_agreemnet)),
                          length(na.omit(EPI_overlap_matrix$six_agreemnet)),
                          length(na.omit(EPI_overlap_matrix$five_agreemnet)),
                          length(na.omit(EPI_overlap_matrix$four_agreement)),
                          length(na.omit(EPI_overlap_matrix$three_agreement)),
                          length(na.omit(EPI_overlap_matrix$two_agreement)),
                          length(na.omit(EPI_support_mat$two_supports)),
                          length(na.omit(EPI_support_mat$one_support)),
                          length(na.omit(EPI_support_mat$novel))))

# Margin: down, left, up, right 
par(mar=c(4,0,4,3))
with(pie_df,pie(FR, radius=1.0, clockwise = TRUE, col = colorRampPalette(c("black", "white"))(12), cex=1.5,
                labels=paste0(FR, "\n", (round(FR/nrow(EPI_overlap_matrix), digits = 2))*100, "%")))

##############################################
###### ---- Novel genes filtering ----- ######
# calculate the missense/truncating ratio for final gene list
setwd("path_to_the_results_direcotry")
all_genelists <- read.table("EPIMUTESTR_SigGenes.tsv", header = T, sep = "\t", colClasses = "character",na.strings=c("",NA)) 
final_genelist <- read.table("EPIMUTESTR_FinalGenes_list.tsv", header = T, sep = "\t", colClasses = "character",na.strings=c("",NA)) 
final_genelist <- final_genelist$Gene

# sort genes by lowest qvalue
cancers_list <- c("ACC",  "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", 
                  "KIRC", "KIRP", "LAML", "LGG",  "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", 
                  "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "PAN")
integrated_insights <- NULL
for (ca in cancers_list){
  ca_genelist <- all_genelists[, c(which(colnames(all_genelists)==ca), which(colnames(all_genelists)==ca)+2)]
  if(any(ca_genelist[,ca]%in%final_genelist)){
    ca_genelist_fltr <- ca_genelist[which(ca_genelist[,ca]%in%final_genelist), ]
    colnames(ca_genelist_fltr) <- c("gene", "qvalue")
    ca_genelist_fltr$tumor <- rep(ca, nrow(ca_genelist_fltr))
    integrated_insights <- rbind(integrated_insights, ca_genelist_fltr)
  }
}

# Truncating mutations
# Missense_Mutation (No), Nonsense_Mutation (Yes), Nonstop_Mutation (Yes), 
# Frame_Shift_Ins (Yes), Frame_Shift_Del (Yes), Translation_Start_Site (Yes), 
# Intron (No), 3'UTR (No), 3'Flank (No), 5'Flank (No), 5'UTR (No), RNA (No), 
# Silent (No), Splice_Region (No), Splice_Site (Yes)
setwd("path_to_the_data_direcotry")
pan_maf <- read.delim("PAN_EA.maf", header = TRUE, colClasses = "character")
pan_maf_fltr <- pan_maf[-which(pan_maf$Variant_Classification=="Silent"), ]
whole_genome <- unique(pan_maf_fltr$GENE)

gof_lof_ratio <- NULL
for (gene in final_genelist){
  pan_maf_gene_fltr <- pan_maf_fltr[which(pan_maf_fltr$GENE==gene), ]
  missense <- nrow(pan_maf_gene_fltr[which(pan_maf_gene_fltr$Variant_Classification%in%
                                             c("Missense_Mutation")), ])
  truncated <- nrow(pan_maf_gene_fltr[which(pan_maf_gene_fltr$Variant_Classification%in%
                                              c("Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Ins", 
                                                "Frame_Shift_Del", "Translation_Start_Site", "Splice_Site")), ])
  gof_lof_ratio <- rbind(gof_lof_ratio, data.frame(gene=gene,mis.trun_ratio=missense/truncated))
}

integrated_insights$mis.truc <- gof_lof_ratio$mis.trun_ratio[match(integrated_insights$gene, gof_lof_ratio$gene)]
# This is suplementary Table 2
write.table(integrated_insights, "integrated_insights.tsv", col.names = T, sep = "\t", quote = FALSE, row.names = F)

# Add information to novel genes
novel_genes <- read.table("EPIMUTESTR_NovelGenes.tsv", sep = "\t", header = T, colClasses = "character")
integrated_insights <- read.table("integrated_insights.tsv", sep = "\t", header = T, colClasses = "character")
integrated_insights_fltr <- integrated_insights[which(integrated_insights$gene%in%novel_genes$Gene), ]
setwd("path_to_the_data_direcotry")
novel_genes_info <- NULL
for(gene in integrated_insights_fltr$gene){
  tumors <- strsplit(integrated_insights_fltr[which(integrated_insights_fltr$gene==gene), "tumor"], ",")
  for (ca in tumors[[1]]) {
    cancer_file <- read.table(paste0(ca, "_maf_cleaned.tsv"), sep = "\t", header = T, colClasses = "character")
    cancer_file$GENE[which(cancer_file$GENE%in%c("FPGT-TNNI3K;TNNI3K;FPGT-TNNI3K",
                                                 "TNNI3K;FPGT-TNNI3K",
                                                 "FPGT-TNNI3K;FPGT;FPGT-TNNI3K;FPGT;FPGT"))] <- "FPGT"
    cancer_file_fltr <- cancer_file[which(cancer_file$GENE==gene), c("ACTION", "Variant_Classification", "Tumor_Sample_Barcode",
                                                                     "Chromosome","Start_Position","End_Position","Reference_Allele",
                                                                     "Tumor_Seq_Allele2","Tumor_Seq_Allele1","dbSNP_RS", "NM", "SUB")]
    
    novel_genes_info <- rbind(novel_genes_info, c(gene, ca, 
                                                  length(cancer_file_fltr$ACTION),
                                                  paste(cancer_file_fltr$Variant_Classification, collapse = ","),
                                                  paste(cancer_file_fltr$Tumor_Sample_Barcode, collapse = ","),
                                                  paste(cancer_file_fltr$NM, collapse = ","),
                                                  paste(cancer_file_fltr$Chromosome, collapse = ","),
                                                  paste(cancer_file_fltr$Start_Position, collapse = ","),
                                                  paste(cancer_file_fltr$End_Position, collapse = ","),
                                                  paste(cancer_file_fltr$Reference_Allele, collapse = ","),
                                                  paste(cancer_file_fltr$Tumor_Seq_Allele2, collapse = ","),
                                                  paste(cancer_file_fltr$Tumor_Seq_Allele1, collapse = ","),
                                                  paste(cancer_file_fltr$dbSNP_RS, collapse = ","),
                                                  paste(cancer_file_fltr$SUB, collapse = ","),
                                                  paste(cancer_file_fltr$ACTION, collapse = ",")))
  }
}
colnames(novel_genes_info) <- c("gene", "type", "freq", "mutation_type", "sample_id", "NM_id", 
                                "Chromosome","Start_Position","End_Position","Reference_Allele",
                                "Tumor_Seq_Allele2","Tumor_Seq_Allele1", "dbSNP_RS", "SUB_id", "EA")

# Plot the original and filtered EA scores distributions for novel genes
all_eas <- NULL
all_eas_genes <- NULL
all_less30 <- NULL
less30_genes <- NULL
all_stop95 <- NULL
stop95_genes <- NULL
all_eas_fltr <- NULL
all_more30 <- NULL
more30_genes <- NULL
for (i in 1:nrow(novel_genes_info)){
  item <- novel_genes_info$EA[i]
  eas <- as.numeric(strsplit(item, ",")[[1]])
  all_eas <- c(all_eas, eas)
  all_eas_genes <- c(all_eas_genes , novel_genes_info$gene[i])
  if(all(eas<30)){
    all_less30 <- c(all_less30, eas)
    less30_genes <- c(less30_genes , novel_genes_info$gene[i])
  } else if ((sum(eas==100)/length(eas))>0.95){
    all_stop95 <- c(all_stop95, eas)
    stop95_genes <- c(stop95_genes , novel_genes_info$gene[i])
  } else {
    all_eas_fltr <- c(all_eas_fltr, eas)
  }
  if ((sum(eas>30)/length(eas))>0.80){
    all_more30 <- c(all_more30, eas)
    more30_genes <- c(more30_genes , novel_genes_info$gene[i])
  } else next
}
integrated_insights_fltr <- integrated_insights[which(integrated_insights$gene%in%more30_genes), ]
all_more30[all_more30==100] <- 101
all_eas[all_eas==100] <- 101
par(mfrow = c(1,2))
hist(all_eas, breaks = seq(0, 110, 10),
     ylab = "Frequency", 
     xlab = "EA", 
     main = paste0(length(unique(all_eas_genes)), " novel genes"), 
     xlim = c(0, 110))#, ylim = c(0,1000))
hist(all_more30, breaks = seq(0, 110, 10),
     ylab = "Frequency", 
     xlab = "EA", 
     main = paste0(length(unique(more30_genes)), " filtered novel genes"), 
     xlim = c(0, 110))#, ylim = c(0,1000))
ks.test(all_eas, all_more30, alternative = "greater")$p.value

