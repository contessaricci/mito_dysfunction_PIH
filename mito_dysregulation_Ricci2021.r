##############################################################################
###                                                                        ###
### AUTHOR: CONTESSA A. RICCI, PhD                                         ###
### MANUSCRIPT: DYSREGULATION OF MITOCHONDRIA-MEDIATED MATERNAL-FETAL      ###
###             INTERACTIONS IN HYPERTENSIVE DISORDERS OF PREGNANCY        ###
###             DOI: (TBD)                                                 ###
### STUDY PURPOSE: Reanalysis of longitudinal maternal RNAseq data from    ###
###                peripheral blood plasma and endpoint fetal RNAseq data  ###
###                from placenta at delivery. Goal is to understand        ###
###                consequences of mitochondrial gene dysregulation by     ###
###                examining expression patterns of genes the dysregulated ###
###                mitochondrial genes are known to interact with          ###
### DATA: GEO accession GSE154377 (maternal longitudinal data)             ###
###       GEO accession GSE114691 (fetal placental)                        ###
###       Accessible via NCBI GEO datasets                                 ###
###                                                                        ###
### NOTES: Script assumes all data has been compiled and is in correct     ###
###        format. This script utilizes some files that have been compiled ###
###        via Python.                                                     ###
###                                                                        ###
##############################################################################

###############################################
###                                         ###
###       EFFECT SIZE MUTATIONAL LOAD       ###
###                                         ###
###############################################
# RUN EFFECT SIZE USING COHEN'S D
# OUTLIER REMOVAL AND T-TEST CARRIED OUT IN GRAPH PRISM

install.packages("effsize")
library(effsize)
mutational_loads <- read.csv("Mutational_load_outliersrm.csv")
cohen.d(mutational_loads$Mutations, mutational_loads$Treatment)
#d estimate: -0.5997953 (medium)
#95 percent confidence interval:
#  lower      upper 
#-1.4936515  0.2940609 

############################################################################
###                                                                      ###
### MATERNAL LONGITUDINAL ANALYSIS ###### MATERNAL LONGITUDINAL ANALYSIS ###
###                                                                      ###
############################################################################


##########################################################
### COVARIATE EXAMINATION ###### COVARIATE EXAMINATION ###
##########################################################
# CORRELATION MATRIX FOR MATERNAL PATIENT CHARACTERISTICS TO DETERMINE
#     AUTOCORRELATION AMONG COVARIATES (**FOR DESEQ2 MODEL LATER**)

covars <- read.csv("maternal_covariates.csv", row.names = 1)

library("car")
qqPlot(covars$Age)
qqPlot(covars$BMI)
qqPlot(covars$Baby_weight_g)
qqPlot(covars$Baby_length_cm) #not normal
qqPlot(covars$Placental_weight_g)
qqPlot(covars$Delivery)

names(covars)[names(covars) == "Baby_weight_g"] <- "Infant wt"
names(covars)[names(covars) == "Baby_length_cm"] <- "Infant lgth"
names(covars)[names(covars) == "Baby_head_circumference_cm"] <- "Head circ"
names(covars)[names(covars) == "Placental_weight_g"] <- "Placenta wt"

library(corrplot)
library(Hmisc)
corr_matrix <- rcorr(as.matrix(covars[,c(3,4,6,7,8,9,10)]), type = "spearman")
corrplot(corr_matrix$r, type="upper", order="hclust", 
         p.mat = corr_matrix$P, sig.level = 0.05, insig = "blank")
#infant weight correlated with placental weight, head circumference, infant length
#age correlated with BMI

# COLLAPSE CORRELATED VARIABLES INTO SUMMARY STATISTICS
#     CREATE UTERO EIGEN VARIABLE
install.packages("factoextra")
library(factoextra)
pca_utero <- prcomp(covars[,c(6:9)], scale = TRUE)
fviz_eig(pca_utero) #PC 1 explains >60% variation
loadings <- as.data.frame(pca_utero$x)
covars$utero_eigen <- loadings$PC1

#     CREATE AGE-BMI EIGEN VARIABLE
pca_ageBMI <- prcomp(covars[,c(3:4)], scale = TRUE)
fviz_eig(pca_ageBMI) #PC 1 explains ~70% variation
loadings <- as.data.frame(pca_ageBMI$x)
covars$ageBMI_eigen <- loadings$PC1

# GROUNDTRUTH THAT SUMMARY STATISTICS REMOVE AUTOCORRELATION
corr_matrix <- rcorr(as.matrix(covars[,c(10:12)]), type = "spearman")
corrplot(corr_matrix$r, type="upper", order="hclust", 
         p.mat = corr_matrix$P, sig.level = 0.05, insig = "blank")
#no correlation present

# WRITE NEW FILE
write.csv(covars, "maternal_covariates_autocorr_rm.csv", row.names = TRUE)


##################################################
### OUTLIER DETECTION ###### OUTLIER DETECTION ###
##################################################
install.packages("ggfortify")
library(ggfortify)
counts_raw <- read.csv("maternal_gestation_raw_counts.csv", row.names="Transcripts")

# SEPARATE GESTATION AND DELIVERY
# REMOVE SAMPLE NT1 (DOES NOT HAVE 1ST TRIMESTER COLLECTION + TCGSA REQUIRES BALANCED, FULL DATA)
gest_counts_raw <- counts_raw[,-c(1,2,3,7,11,15,19,23,27,31,3,35,39,43,47,51,55,59,63,67)]
del_counts_raw <- counts_raw[,c(7,11,15,19,23,27,31,35,39,43,47,51,55,59,63,67)]

# OUTLIER DETECTION VIA MEDIAN ABSOLUTE DISTANCE (MAD) -- OUTLIERS INDICATED IN RED
library(magrittr)
library(ggplot2)

#   OUTLIER DETECTION GESTATION TIMEPOINTS
pca_subjs <- prcomp(t(gest_counts_raw))
PCs <- pca_subjs$x
qplot(PCs[, 1], PCs[, 2])
apply(PCs, 2, function(x) which( abs(x - mean(x)) > (6 * sd(x)) )) %>%
  Reduce(union, .)
ind.out <- apply(PCs, 2, function(x) which( (abs(x - median(x)) / sd(x)) > 6 )) %>%
  Reduce(union, .) %>%
  print()
col <- rep("black", nrow(PCs)); col[ind.out] <- "red"
qplot(PCs[, 1], PCs[, 2], color = I(col), size = I(2))
#no outliers detected

#   OUTLIER DETECTION DELIVERY TIMEPOINT
pca_subjs <- prcomp(t(del_counts_raw))
PCs <- pca_subjs$x
qplot(PCs[, 1], PCs[, 2])
apply(PCs, 2, function(x) which( abs(x - mean(x)) > (6 * sd(x)) )) %>%
  Reduce(union, .)
ind.out <- apply(PCs, 2, function(x) which( (abs(x - median(x)) / sd(x)) > 6 )) %>%
  Reduce(union, .) %>%
  print()
col <- rep("black", nrow(PCs)); col[ind.out] <- "red"
qplot(PCs[, 1], PCs[, 2], color = I(col), size = I(2))
#no outliers detected


######################################################################
### EXPLORE STRUCTURING OF DATA ###### EXPLORE STRUCTURING OF DATA ###
### PURPOSE: determine which covariates are most influential       ###
###             at each timepoint on gene expression. These will   ###
###             be used as covariates for DEseq models             ###
######################################################################

metafile <- read.csv("maternal_DESeq_metafile.csv", row.names = 1)

# SEPARATE GESTATION AND DELIVERY
dim_df <- as.data.frame(dim(metafile))
covars_gest <- subset(metafile[c(4:dim_df[1,]),], Collection != "Delivery")
covars_del <- subset(metafile[c(4:dim_df[1,]),], Collection == "Delivery")

# SPECIFY SAMPLES (COLUMNS) FOR EACH CONTRAST
names(gest_counts_raw)
T1 <- gest_counts_raw[,c(1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46)]
T2 <- gest_counts_raw[,c(2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47)]
T3 <- gest_counts_raw[,c(3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48)]
DEL <- del_counts_raw

# CONFIRM PC1 EXPLAINS ADEQUATE AMOUNT OF VARIATION
#   TRIMESTER 1
pca_T1 <- prcomp(T1, scale = TRUE)
fviz_eig(pca_T1) #PC 1 explains > 80% variation for T1
#   TRIMESTER 2
pca_T2 <- prcomp(T2, scale = TRUE)
fviz_eig(pca_T2) #PC 1 explains ~ 80% variation for T2
#   TRIMESTER 3
pca_T3 <- prcomp(T3, scale = TRUE)
fviz_eig(pca_T3) #PC 1 explains ~ 80% variation for T3
#   DELIVERY
pca_DEL <- prcomp(DEL, scale = TRUE)
fviz_eig(pca_DEL) #PC 1 explains > 80% variation for delivery

# BIND PC1 LOADINGS FROM GENE EXPRESSION PCA TO EACH COVARIATE TIMEPOINT DATAFRAME
#   TRIMESTER 1
covars_T1 <- subset(metafile[-c(1:3),], Collection == "Trimester 1")
cbind(rownames(covars_T1), rownames(pca_T1$rotation)) #confirm samples are in correct order
covars_T1$loadings_T1 <- as.data.frame(pca_T1$rotation)$PC1
#   TRIMESTER 2
covars_T2 <- subset(metafile[-c(1:3),], Collection == "Trimester 2")
cbind(rownames(covars_T2), rownames(pca_T2$rotation)) #confirm samples are in correct order
covars_T2$loadings_T2 <- as.data.frame(pca_T2$rotation)$PC1
#   TRIMESTER 3
covars_T3 <- subset(metafile[-c(1:3),], Collection == "Trimester 3")
cbind(rownames(covars_T3), rownames(pca_T3$rotation)) #confirm samples are in correct order
covars_T3$loadings_T3 <- as.data.frame(pca_T3$rotation)$PC1
#   DELIVERY
covars_DEL <- subset(metafile[-c(1:3),], Collection == "Delivery")
cbind(rownames(covars_DEL), rownames(pca_DEL$rotation)) #confirm samples are in correct order
covars_DEL$loadings_DEL <- as.data.frame(pca_DEL$rotation)$PC1


# CONDUCT LINEAR MODELS TO DETERMINE INFLUENCE OF COVARIATES ON GENE EXPRESSION PROFILES
#   (note: replace interaction term with covariate of interest)
#   (note: asteriks [*] next to covariate indicates inclusion in DESeq model)
#   TRIMESTER 1
summary(lm(loadings_T1 ~ ageBMI_eigen * Condition, data = covars_T1))
# Delivery_mode *
# Fetal_sex
# T1_collect
# T2_collect
# T3_collect
# Delivery_wks *
# utero_eigen *
# ageBMI_eigen

#   TRIMESTER 2
summary(lm(loadings_T2 ~ utero_eigen * Condition, data = covars_T2))
# Delivery_mode
# Fetal_sex
# T1_collect *
# T2_collect
# T3_collect
# Delivery_wks
# utero_eigen *
# ageBMI_eigen

#   TRIMESTER 3
summary(lm(loadings_T3 ~ ageBMI_eigen * Condition, data = covars_T3))
# Delivery_mode
# Fetal_sex
# T1_collect
# T2_collect
# T3_collect *
# Delivery_wks
# utero_eigen
# ageBMI_eigen *

#   DELIVERY
summary(lm(loadings_DEL ~ T3_collect * Condition, data = covars_DEL))
# Delivery_mode
# Fetal_sex
# T1_collect
# T2_collect
# T3_collect *
# Delivery_wks *
# utero_eigen
# ageBMI_eigen


##################################################################
###     CONDUCT DESEQ ANALYSES ### CONDUCT DESEQ ANALYSES      ###
### Contrasts of interest:                                     ###
###    Trimester 1 NT v PE                                     ###
###    Trimester 2 NT v PE                                     ###
###    Trimester 3 NT v PE                                     ###
###    Delivery NT v PE                                        ###
###                                                            ###
### PURPOSE: Detect mitochondrial genes in DEGs (mtDEGs)       ###
###                                                            ###
### NOTES: will pool all mtDEGs for TcGSA gene sets            ###
##################################################################

# THERE ARE FOUR RAW COUNTS FILES WITH CORRESPONDING METAFILES:
#     - T1    covars_T1
#     - T2    covars_T2
#     - T3    covars_T3
#     - DEL   covars_DEL

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

# ENSURE THAT ROW NAMES OF METAFILES MATCH COLUMN NAMES OF COUNTFILES
#   TRIMESTER 1
cbind(rownames(covars_T1), colnames(T1)) #make sure sample names are in same order
rownames(covars_T1) <- colnames(T1)
#   TRIMESTER 2
cbind(rownames(covars_T2), colnames(T2)) #make sure sample names are in same order
rownames(covars_T2) <- colnames(T2)
#   TRIMESTER 3
cbind(rownames(covars_T3), colnames(T3)) #make sure sample names are in same order
rownames(covars_T3) <- colnames(T3)
#   DELIVERY
cbind(rownames(covars_DEL), colnames(DEL)) #make sure sample names are in same order
rownames(covars_DEL) <- colnames(DEL)


# RUN DESEQ
#   (note: countData is count file [T1, T2, T3, DEL];
#          colData is corresponding metafile [covars_T1, covars_T2, covars_T3, covars_DEL])
#   (note: unremark covariates to include them in the DESeq model)
dds <- DESeqDataSetFromMatrix(countData = T2, colData = covars_T2, design = ~Condition +
                                #Delivery_mode +
                                #Fetal_sex +
                                T1_collect +
                                #T2_collect +
                                #T3_collect +
                                #Delivery_wks +
                                utero_eigen +
                                ageBMI_eigen
                                )
dds <- DESeq(dds) #run DESeq

# WRITE DESEQ FILES - WILL USE:
#                         - ALL CONTRAST RESULTS FOR EACH GESTATION TIME POINT [IPA]
#                         - ONLY SIGNIFICANT CONTRAST RESULTS FOR EACH TIME POINT (GESTATION TIME
#                           POINTS AND AT DELIVERY) [MITO GENE AND MITO INTERACTION GENE SEARCHING + 
#                           GO FUNCTIONAL ENRICHMENTS]

#conduct contrasts and write files
Contr_results <- results(dds, contrast=c("Condition", "Normotensive", "PE/HTN"))
write.csv(Contr_results, "T1_Contr.csv") #full contrasts results; rename as necessary
sigContr <- subset(Contr_results, Contr_results$padj <= 0.05)
write.csv(del_sigContr, "T1_sigContr.csv") #only significant DEGs; rename as necessary



##########################################################################
###                                                                    ###
###                   NORMALIZE ALL GESTATION BY CPM                   ###
###                           * MATERNAL ONLY *                        ###
###                                                                    ###
###                                                                    ###
### PURPOSE:                                                           ###
###     Total count normalization across entire gestation for TcGSA    ###
###                                                                    ###
########################################################################## 

# ALL GESTATION COUNTS: gest_counts_raw

# CONFIRM LIBRARY SIZES NEED TO BE CONTROLLED FOR
lib_sizes = as.data.frame(seq(1,length(gest_counts_raw)))
colnames(lib_sizes) <- "col_num"
lib_sizes$libsize = 0
for(i in seq(1,length(gest_counts_raw), by = 1)){
  lib_sizes[i,2] = sum(gest_counts_raw[,i])
}
# plot distribution of library sizes
ggplot(data = lib_sizes, aes(x=col_num, y=libsize)) + 
  geom_bar(stat="identity")

# CORRECT FOR LIBRARY SIZE BY CPM
gest_CPM <- gest_counts_raw
for(i in seq(1,length(gest_CPM), by = 1)){
  gest_CPM[,i] = gest_CPM[,i]/(sum(gest_CPM[,i])/1000000)
}

# CONFIRM LIBRARY SIZES ARE CONSISTENT ACROSS SAMPLES
lib_sizes = as.data.frame(seq(1,length(gest_CPM)))
colnames(lib_sizes) <- "col_num"
lib_sizes$libsize = 0
for(i in seq(1,length(gest_CPM), by = 1)){
  lib_sizes[i,2] = sum(gest_CPM[,i])
}
# plot distribution of normalized library sizes
ggplot(data = lib_sizes, aes(x=col_num, y=libsize)) + 
  geom_bar(stat="identity")

# WRITE NORMALIZED COUNTS FILE
write.csv(gest_CPM, "maternal_gestation_normalized_counts_TcGSA.csv")


#################################################################
###                                                           ###
### SEARCH DEGs FROM MATERNAL (GESTATION AND AT DELIVERY) AND ### 
###   FETAL (AT DELIVERY) SAMPLES FOR mtDEGs AND RESPECTIVE   ###
###                   INTERACTION GENES.                      ###
###                                                           ###
###  USE PYTHON SCRIPTS:                                      ###
###                 - find_mtDEGs.py                          ###
###                 - find_mtDEG-INTXs.py                     ###
###                                                           ###
### PURPOSE: -  Find mtDEGs and respective mtDEG interaction  ### 
###             genes (mtDEG-INTXs) present for all samples   ###
###             and timepoints of interest                    ###
###          -  Compile GMT files (actual and null model) for ###
###             TcGSA on maternal gestation                   ###
###          -  Write mtDEG-INTXs expression files for        ###
###             maternal and fetal at delivery for gene       ###
###             ontology enrichment and plotting purposes     ###
#################################################################


#################################################################
###                                                           ###
###   RUN TCGSA ON NORMALIZED GESTATION COUNTS AND GMT FILE   ###
###               MADE USING PYTHON SCRIPTS                   ###
###                                                           ###
#################################################################

# note: outliers have already been identified and removed
# note: NT1 samples are removed due to incomplete sampling (not sampled during 1st trimester)
BiocManager::install("multtest")
library(GSA)
require(TcGSA)


# READ IN GMT FILE
gmt_file <- GSA.read.gmt("gestation_GMT.txt")
#gmt_file <- GSA.read.gmt("null_GMT.txt") # (note:  unremark and read in for null model testing)
head(gmt_file) #make sure was imported in the correct format

# READ IN NORMALIZED EXPRESSION DATA
gestation_expr <- read.csv("maternal_gestation_normalized_counts_TcGSA.csv", row.names = 1)

# READ IN EXPERIMENTAL DESIGN DATA
exp_design <- read.csv("TcGSA_experimental_design.csv")

# RUN TcGSA USING DEFAULT PARAMETERS
tcgsa_result <- TcGSA.LR(expr = gestation_expr, 
                         gmt = gmt_file, 
                         design = exp_design,
                         subject_name = "Subject_no", 
                         time_name = "Collection")

summary(tcgsa_result) #view results summary
TcGSA::signifLRT.TcGSA(tcgsa_result)$mixedLRTadjRes #view significant gene sets
model_check <- as.data.frame(TcGSA::multtest.TcGSA(tcgsa_result)) #check model quality

#cluster significant TcGSA gene sets
clust <- TcGSA::clustTrend(tcgs = tcgsa_result, 
                           expr = tcgsa_result$Estimations,
                           Subject_ID = exp_design$Subject_no,
                           TimePoint = exp_design$Collection,
                           group.var = exp_design$Condition,
                           group_of_interest="PE/HTN",
                           ref="Normotensive")

#visualize TcGSA results
#   (note: gene set names on y-axis and collection names on x-axis were modified in illustrator
#          due to inability to customize using built-in visualization function for TcGSA package)
plot(x = tcgsa_result, expr = tcgsa_result$Estimations,
     Subject_ID = exp_design$Subject_no,
     TimePoint = exp_design$Collection,
     group_of_interest = "Condition",
     clust_trends = clust,
     legend.breaks = seq(from = -2,to = 2, by = 0.01), 
     time_unit = "Collection ",
     subtitle = "Normotensive v PE/HTN", cex.label.row = 1, cex.label.col = 0.7, cex.main = 0.7,
     heatmap.width = 0.5, dendrogram.size = 0.1, heatKey.size = 0.5,margins = c(4,6))

#######################################################
###                                                 ###
###     MAKE GO PLOTS USING DATAFRAMES MADE BY      ###
###           make_GOplots_dataframes.py            ###
###                                                 ###
### PURPOSE: REDUCE TERMS IN GO ENRICHMENTS FOR     ###
###          EFFECTIVE VISUALIZATION                ###
###                                                 ###
#######################################################

install.packages("GOplot")
library(GOplot)

maternal_GOdata <- read.csv("maternal_delivery_GOplot_data.csv")
fetal_GOdata <- read.csv("fetal_delivery_GOplot_data.csv")
head(maternal_GOdata)
head(fetal_GOdata)

# calculate z-score with directionality
maternal_GOdata$zscore <- ifelse(maternal_GOdata$logFC < 0, qnorm(maternal_GOdata$adj_pval, lower.tail = T),qnorm(maternal_GOdata$adj_pval, lower.tail = F))
fetal_GOdata$zscore <- ifelse(fetal_GOdata$logFC == "Down", qnorm(fetal_GOdata$adj_pval, lower.tail = T),qnorm(fetal_GOdata$adj_pval, lower.tail = F))
#write.csv(maternal_GOdata,"maternal_delivery_GOplot_data_zscore.csv", row.names = F) #uncomment to write dataframe with zscores
#write.csv(fetal_GOdata,"fetal_delivery_GOplot_data_zscore.csv", row.names = F) #uncomment to write dataframe with zscores

# collapse terms based on 99% similarity of shared genes
maternalBP_reduced <- reduce_overlap(subset(maternal_GOdata[,c(1,3:9)], category == "BP"), overlap = 0.99)
maternalCC_reduced <- reduce_overlap(subset(maternal_GOdata[,c(1,3:9)], category == "CC"), overlap = 0.99)
fetalBP_reduced <- reduce_overlap(subset(fetal_GOdata[,c(1,3:9)], category == "BP"), overlap = 0.99)
fetalCC_reduced <- reduce_overlap(subset(fetal_GOdata[,c(1,3:9)], category == "CC"), overlap = 0.99)

# make dataframes for plotting
maternalBP_reduced$Individual <- "Maternal"
maternalCC_reduced$Individual <- "Maternal"
fetalBP_reduced$Individual <- "Fetal"
fetalCC_reduced$Individual <- "Fetal"

both_BP <- rbind(maternalBP_reduced, fetalBP_reduced)
both_CC <- rbind(maternalCC_reduced, fetalCC_reduced)

BP_plot <- both_BP %>%
  ggplot( aes(x=Individual, y=term, size = count, color=zscore)) +
  scale_colour_gradient(low = "cyan4", high = "chocolate1") + 
  theme_minimal()+theme(axis.title.x=element_blank(),
                        axis.title.y=element_blank())+
  geom_point(alpha=0.8)+
  scale_size(range = c(0.1, 10), name="count") + ylab("")

CC_plot <- both_CC %>%
  ggplot( aes(x=Individual, y=term, size = count, color=zscore)) +
  scale_colour_gradient(low = "cyan4", high = "chocolate1") + 
  theme_minimal()+theme(axis.title.x=element_blank(),
                        axis.title.y=element_blank())+
  geom_point(alpha=0.8)+
  scale_size(range = c(0.1, 10), name="count") + ylab("")
