# Import ----
library(stringr)
library(Biobase)
library(BiocGenerics)
library(ggplot2)
library(limma)
library(oligo)
library(tools)
library(oligo)
library(variancePartition)
library(GEOquery)
# Initial variables ----
# Get genetic information from database
GSE124691 <- getGEO('GSE124691')
# File stored at: 
#   /var/folders/1b/dyj_c0r54wncq4rq60scgf9m0000gn/T//RtmpiB711P/GPL19057.soft

# extract GSE Object
GSEOBJ <- GSE124691$GSE124691_series_matrix.txt.gz

#extract data about gene assay and subjects
assayexprs <- exprs(GSEOBJ)
feats1 <- GSEOBJ@featureData
pheno1 <-GSEOBJ@phenoData@data


# # Check which attributes are prime for this experiment
findRelevantQualitativeVars <- function(phenotypeData, linebreak = '\n __________________ \n', print_unique = TRUE){
  for (ii in seq_along(attributes(phenotypeData)$names)) {
    n <- attributes(phenotypeData)$names[ii];
    m_unique <- unique(phenotypeData[ii])
    m_all <- phenotypeData[ii]
    totalrows <- nrow(phenotypeData[ii])
    
    if (print_unique){
      mp <- m_unique
      is_unique = 'Unique'
    }else{
      mp <- m_all
      is_unique = ''
    }
    
    if (nrow(m_unique)>1){
      cat('ATTRIBUTE NAME:',n,'\n')
      if (nrow(mp)<50){
        print(paste0(is_unique,' Attribute Options: \n'))
        print(mp)
      }else{
        print(paste0('Has ', nrow(m_unique),'/', nrow(m_all),' unique rows. Too many to print'))
      }
      cat('Number of unique attribute options out of total', nrow(m_unique),'/',nrow(m_all),'\n')
      print(linebreak)
    }
  }
}
# Example usage: findRelevantQualitativeVars(pheno1)
findRelevantQualitativeVars(pheno1)


# Create list of matrices -------------------------------------------------------------------
gsm_mat_list <- vector("list", 7)
gsm_id_list <- vector("list", 7)
dfinit = FALSE

i <- 1
for (gsmID in (pheno1$geo_accession)){
  gsm <- getGEO(gsmID, destdir = "GSE124691/GSM")
  
  supfile2ex <- gsm@header$supplementary_file_2
  fname <- tail(unlist(strsplit(supfile2ex,"/")),1)
  fname <- str_sub(fname, end=-4)
  # print(fname)
  supmatrix <- readMM(paste0('GSE124691/GSE124691_RAW/',fname))
  gsm_mat_list[[i]] <- supmatrix
  gsm_id_list[[i]] <- gsmID
  i <- i+1
}

# Exclude non-relevant sample ----
# sample 3 was excluded but happens later in the ode
pheno_concise <- pheno1[c('title','source_name_ch1','characteristics_ch1.1',
                          'characteristics_ch1.3', 'description','description.1',
                          'tissue:ch1','time point:ch1')]
included_GSM = c(2,3,4,6,7)
pheno_concise_EXC <- pheno_concise[included_GSM]
gsm_mat_list_EXC <- gsm_mat_list[included_GSM]
gsm_id_list_EXC <- gsm_id_list[included_GSM]


# Create full matrix ------
# make one large matrix and vector of column names (GSM___)
for (i in seq_along((gsm_mat_list_EXC))){
  ID <- gsm_id_list_EXC[[i]]
  mat <- gsm_mat_list_EXC[[i]]
  
  # list of ID's as wide as matrix
  IDL <- (rep(ID,dim(mat)[[2]]))
  print(i)
  if (i==1){
    #init list of IDs and large matrix
    gsm_mat_combo_EXC <- mat
    gsm_id_combo_EXC <- IDL
  }else{
    
    gsm_mat_combo_EXC <- cbind2(gsm_mat_combo_EXC, mat)
    gsm_id_combo_EXC <- c(gsm_id_combo_EXC,IDL)
  }
}
gsm_id_combo_EXC <- factor(gsm_id_combo_EXC)


# get genes from file ----
genesfname <- 'GSE124691/GSE124691_Genes.tsv'
genenames <- read.table(file = genesfname, sep = '\t', header = FALSE)

# Create scaled phenotype/annotation matrix
# GSM_pheno_EXC = samples x 4
# gsm_mat_combo_EXC = genes x samples
# gsm_id_combo_EXC = vector of length(samples)

first_run=TRUE
for (ID in gsm_ID_EXC) {
  if (first_run){
    first_run=FALSE
    ind = which(pheno_concise2_EXC$indv == ID)
    GSM_pheno_EXC <- pheno_concise2_EXC[ind,]
  } else {
    ind = which(pheno_concise2_EXC$indv == ID)
    GSM_pheno_EXC <- rbind(GSM_pheno_EXC,pheno_concise2_EXC[ind,])
  }
}


# Create dense matrix for use in genomics software ----
# GSM
GSM_MAT_EXC <- as.matrix(gsm_mat_combo_EXC) #Matrix of counts
GSM_ID_EXC <- as.character(gsm_id_combo_EXC) #Sample ID's
colnames(GSM_MAT_EXC) <- GSM_ID_EXC
rownames(GSM_MAT_EXC) <- genenames$V2


# Using Limma-Voom -----
# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

# Create DGEList object
d0 <- DGEList(as.matrix(GSM_MAT_EXC))

# Calculate downstream normalization factors
d0 <- calcNormFactors(d0)

# Filter low-expressed genes 
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left

# remove GSM3543445 because it is not a replicant
LCb2 <- which(GSM_pheno_EXC$indv='GSM3543445')
GSM_pheno_EXC2 <- GSM_pheno_EXC[-LCb2,]
dim(GSM_pheno_EXC2)
d <- d[,-LCb2]
dim(d)

# Correct replicant/batch number
Batch3 <- which(GSM_pheno_EXC2$batch_num==3)
GSM_pheno_EXC2$batch_num[Batch3] <- 2


# RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR -----
# https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#differential-expression-analysis

# Create DGEList object
# see left script
dim(d)

# Add Group (cell type) and lane (batch/rep) to samples
gp <- as.factor(GSM_pheno_EXC2$source_name_ch1)

levels(gp) <- c("LCMV_GP66+","TILs_CD44+")
# Add group (cell tissue type) to count matrix d
d$samples$group <- gp
# Add lane (batch) to count matrix d
d$samples$lane <- as.factor(GSM_pheno_EXC2$batch_num)

# get genes to match 
d$genes <- genenames$V2[-drop,]

# 5. Data pre-processing -----
# Transformations from the raw-scale ----
cpm <- cpm(d)
lcpm <- cpm(d, log=TRUE)

# 5.3 Normalising gene expression distributions ----
d <- calcNormFactors(d, method = "TMM")



# Unsupervised clustering of samples (PCA)----
# Computationally expensive skip graphing if possible
group <- gp
lane <- as.factor(GSM_pheno_EXC2$batch_num)

lcpm <- cpm(d, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)

plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")

plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")

# Differential expression analysis ----
# 6.1 Creating a design matrix and contrasts
levels(group) <- c("LCMV_GP66","TILs_CD44")
glint<-interaction(group,lane,sep='_L')
head(glint)


design <- model.matrix(~0+group+glint)
colnames(design) <- gsub("group", "", colnames(design))
colnames(design) <- gsub("glint", "", colnames(design))
head(design)

contr.matrix <- makeContrasts(
  LCMVvsTIL = LCMV_GP66 - TILs_CD44, 
  # LCMV__L1vsL2 = LCMV_GP66_L1 - LCMV_GP66_L2,
  # TIL__L1vsL2 = TILs_CD44_L1 - TILs_CD44_L2,
  levels = colnames(design))
contr.matrix

# 6.2 Removing heteroscedascity from count data ----
par(mfrow=c(1,2))
v <- voom(d, design, plot=TRUE)
head(v)


vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

# Examining the number of DE genes ----
print('Based on p-value alone')
summary(decideTests(efit))

# The treat method (McCarthy and Smyth 2009) can be used to calculate 
# p-values from empirical Bayes moderated t-statistics with a minimum 
# log-FC requirement
LogFCreq = 2
tfit <- treat(vfit, lfc=LogFCreq)
dt <- decideTests(tfit)
print('Based on p-value with minimum Log-FC that is significantly greater than 1')
summary(dt)

# 6.5 Examining individual DE genes from top to bottom ----
  # topTreat arranges genes from smallest to largest adjusted p-value
LCMV.vs.TIL <- topTreat(tfit, coef=1, n=Inf)
head(LCMV.vs.TIL)
# plotMD(tfit, column=1, status=dt, main=colnames(tfit)[1]) 
       # xlim=c(-8,13))

# interactive plot
# glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
#          side.main="ID", counts=lcpm, groups=group, launch=FALSE)

# make use of significant genes
sig_genes <- rownames(dt)[which(dt!=0)]
length(sig_genes)

# index significant genes
indices <- match(sig_genes, d$genes)
head(sig_genes)
head(d$genes[indices])

# plot best gene ----
nullhyptest <- function(pval, a=.05){
  if (pval<a){
    hyp <- 'rejected'
  }else{
    hyp <- 'accepted'
  }
  return(hyp)
}

# 5 looks good, ~6, 7, 9, 12
for (geneindex in indices){
  most_significant_gene <- geneindex
  datalg <- lcpm[most_significant_gene,]
  
  gene1 <- data.frame('logcounts' = datalg, 'CellnBatch' = glint, 'Cell_Type' = group)
  
  ## Statistical test
  # Distribution: non-normal
  # Input var: Nominal
  # Outcome var: Quantitative non-normal/Discrete
  # Best test: Independent 2 group Whitney-Man U test
  
  #Organize by batch/rep and cell tissue
  X_by_glint <- split(gene1, gene1$CellnBatch)
  #LCMV rep
  LCMV_rep1 <- X_by_glint$LCMV_GP66_L1$logcounts
  LCMV_rep2 <- X_by_glint$LCMV_GP66_L2$logcounts
  #TIL rep
  TIL_rep1 <- X_by_glint$TILs_CD44_L1$logcounts
  TIL_rep2 <- X_by_glint$TILs_CD44_L2$logcounts
  
  # independent 2-group Mann-Whitney U Test 
  w_LCMV <- wilcox.test(LCMV_rep1,LCMV_rep2) # For Numeric Data
  w_TIL <- wilcox.test(TIL_rep1,TIL_rep2) # For Numeric Data
  
  hyp_LCMV <- paste0("LCMV batch null hypothesis ",nullhyptest(w_LCMV$p.value)," with p-value=",round(w_LCMV$p.value,5))
  hyp_TIL <- paste0("TIL batch null hypothesis ",nullhyptest(w_TIL$p.value)," with p-value=",round(w_TIL$p.value,5))
  
  plot_title <- paste0("Boxplot of gene \"",d$genes[most_significant_gene],"\" expression in \nLCMV_GP66+ and TIL CD4+ cells\n")
  plot_caption <- paste0("Hypothesis testing across batch within cell tissue type\n",hyp_LCMV,'\n',hyp_TIL)
  
  p <- ggplot(gene1, aes(x=CellnBatch, y=logcounts, fill=Cell_Type)) + 
    # geom_violin(trim=TRUE) +
    geom_boxplot(width=0.3)+
    labs(title=plot_title,caption = plot_caption, x="Cell and Batch", y = "Log Transformed Counts Per Million") +
    scale_x_discrete(limits=c("LCMV_GP66_L1","LCMV_GP66_L2", "TILs_CD44_L1", "TILs_CD44_L2"))
  
  # Save File
  genename_ = d$genes[most_significant_gene]
  path_ = paste0("./LogFC_2_plots/",genename_,"_boxplot.png")
  # png(file=path_, width=600, height=350)
  p
  # dev.off()
  ggsave(
    path_,
    p,
    width = 6,
    height = 3.5,
    dpi = 1200
  )
}




