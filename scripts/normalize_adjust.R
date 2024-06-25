library(stringr)
library(data.table)
library(readr)
library(CePa)
library(stringr)
library(plyr)
library(dplyr)
library(limma)
library(matrixStats)
library(tidyverse)
library(edgeR)

`%!in%` = Negate(`%in%`)

###################
## upload kallisto counts
###################

txi.gene = readRDS('gtex_kallisto_genes.rds')
counts_all = txi.gene$counts

###################
## upload meta data
###################

keep_samples = readRDS('gtex_combined_meta.rds')

gam_list = read.csv('gametologs_in_genome.csv')

######################################
## normalize, filter, and adjust expression data
######################################

for (i in 1:length(tissue.list)){
  
  print(paste("Now analyzing:",tissue.list[i]))
  
  # get current data sets 
  m_now = subset(keep_samples, SMTSD == tissue.list[i])
  
  # only keep samples with all data available
  m_now = m_now[complete.cases(m_now[,c('SEX','AGE','SMTSISCH','SMRIN','SMNTRNRT')]),]

  if(length(m_now$SAMPID) == 0) {
    print('not enough samples') } else {
    
    # get count matrix for tissue
    counts = counts_all[,which(colnames(counts_all) %in% m_now$SAMPID)]
    
    ##################
    ## males & females
    ##################
    
    # calculate TMM normalization factors
    d0 <- DGEList(counts)
    d0 <- calcNormFactors(d0, method = 'TMM')

    # remove low expresssed genes (mean CPM < 5)
    cutoff = 5
    mean_cpm = rowMeans(cpm(d0))
    drop = names(which(mean_cpm < cutoff))
    d = d0[which(rownames(d0) %!in% drop),]
    rm(d0)
    exp_genes = rownames(d)
    
    # design covariate matrix
    indx = match(colnames(d$counts),m_now$SAMPID) 
    m_now = m_now[indx,]
    m_now$scaleAGE = scale(m_now$AGE)
    m_now$scaleSMTSISCH = scale(m_now$SMTSISCH)
    m_now$scaleSMRIN = scale(m_now$SMRIN)
    m_now$scaleSMNTRNRT = scale(m_now$SMNTRNRT)
    design = model.matrix(~ SEX + scaleAGE + scaleSMTSISCH + scaleSMRIN + scaleSMNTRNRT, data=m_now, na.action=na.pass)

    # normalize data 
    y = voom(d, design, plot = T)
    exp = y$E
    
    # E	= numeric matrix of normalized expression values on the log2 scale
    # weights	= numeric matrix of inverse variance weights
    # design = design matrix
    # lib.size = numeric vector of total normalized library sizes
    
    # prune samples based on mean (outside 3 SDs)
    samp_cor <- cor(exp, use="everything", method="spearman")
    rownames(samp_cor) <- colnames(exp)
    colnames(samp_cor) <- colnames(exp)
    samp_cor_means <- apply(samp_cor, 1, mean)
    samp_cor_z <- scale(samp_cor_means)
    colnames(samp_cor_z) <- "means"
    samples_pruned <- rownames_to_column(as.data.frame(samp_cor_z)) %>% filter(means >= -3) %>% filter(means <= 3)
    m_now = subset(m_now, SAMPID %in% samples_pruned$rowname)
    samp = m_now$SAMPID
    
    # save gene list
    saveRDS(exp_genes, paste(tissue.list[i],"_genes_ALL.rds",sep=""))
    
    # save sample list
    saveRDS(samp, paste(tissue.list[i],"_samples_ALL.rds",sep=""))
    
    # prune and export voom for sex-biased expression
    design = model.matrix(~ SEX + scaleAGE + scaleSMTSISCH + scaleSMRIN + scaleSMNTRNRT, data=m_now, na.action=na.pass)
    y = voom(d[,which(colnames(d) %in% m_now$SAMPID)], design, plot = T)
    print(table(colnames(y$E) == m_now$SAMPID))
    
    saveRDS(y, file=paste(tissue.list[i],"_normalized_exp_ALL.rds",sep=""))
    
    ################
    ## males only
    ################
    
    # prune calculate residual values for co-expression within each sex
    m = subset(m_now, SEX == 1) 
    
    # calculate TMM normalization factors
    d0 <- DGEList(counts[,which(colnames(counts) %in% m$SAMPID)])
    d0 <- calcNormFactors(d0, method = 'TMM')
    
    # remove low expresssed genes (mean CPM < 5)
    cutoff = 5
    mean_cpm = rowMeans(cpm(d0))
    drop = names(which(mean_cpm < cutoff))
    d = d0[which(rownames(d0) %!in% drop),]

    # add back gams (mean CPM >= 1)
    missing_gams = gam_list$ensembl_id[which(gam_list$ensembl_id %!in% rownames(d))]
    missing_gamsE = mean_cpm[missing_gams]
    add_gams = names(missing_gamsE[which(missing_gamsE >= 1)])
    add_gamsE = d0[add_gams,]
    d = rbind(d, add_gamsE)
    rm(d0)
    exp_genesM = rownames(d)
    
    # design covariate matrix
    indx = match(colnames(d$counts),m$SAMPID) 
    m = m[indx,]
    m$scaleAGE = scale(m$AGE)
    m$scaleSMTSISCH = scale(m$SMTSISCH)
    m$scaleSMRIN = scale(m$SMRIN)
    m$scaleSMNTRNRT = scale(m$SMNTRNRT)
    
    design = model.matrix(~ scaleAGE + scaleSMTSISCH + scaleSMRIN + scaleSMNTRNRT, data=m, na.action=na.pass)
    ym = voom(d, design, plot = T)
    tryCatch(print(table(colnames(ym$E) == m$SAMPID)), error=function(err) NA)
    
    # adjust expression levels for age + technical effects
    exp_now_adj_males = tryCatch(removeBatchEffect(ym, covariates=design), error=function(err) NA)
    
    # export residual expression
    saveRDS(exp_now_adj_males, file = paste(tissue.list[i],'_adjusted_exp_MALES.rds',sep=""))
    
    ################
    ## females only
    ################
    
    # prune calculate residual values for co-expression within each sex
    f = subset(m_now, SEX == 2)
    
    # calculate TMM normalization factors
    d0 <- DGEList(counts[,which(colnames(counts) %in% f$SAMPID)])
    d0 <- calcNormFactors(d0, method = 'TMM')
    
    # remove low expresssed genes (mean CPM < 5)
    cutoff = 5
    mean_cpm = rowMeans(cpm(d0))
    drop = names(which(mean_cpm < cutoff))
    d = d0[which(rownames(d0) %!in% drop),]
    
    # add back gams (mean CPM >= 1)
    missing_gams = gam_list$ensembl_id[which(gam_list$ensembl_id %!in% rownames(d))]
    missing_gamsE = mean_cpm[missing_gams]
    add_gams = names(missing_gamsE[which(missing_gamsE >= 1)])
    add_gamsE = d0[add_gams,]
    d = rbind(d, add_gamsE)
    rm(d0)
    exp_genesF = rownames(d)
    
    # design covariate matrix
    indx = match(colnames(d$counts),f$SAMPID) 
    f = f[indx,]
    f$scaleAGE = scale(f$AGE)
    f$scaleSMTSISCH = scale(f$SMTSISCH)
    f$scaleSMRIN = scale(f$SMRIN)
    f$scaleSMNTRNRT = scale(f$SMNTRNRT)
    
                                 design = model.matrix(~ scaleAGE + scaleSMTSISCH + scaleSMRIN + scaleSMNTRNRT, data=f, na.action=na.pass)
    yf = tryCatch(voom(d[,which(colnames(d) %in% f$SAMPID)], design, plot = T), error=function(err) NA)
    tryCatch(print(table(colnames(yf$E) == f$SAMPID)), error=function(err) NA)
    
    # adjust expression levels for age + technical effects
    exp_now_adj_females = tryCatch(removeBatchEffect(yf, covariates=design), error=function(err) NA)
    
    # export residual expression
    saveRDS(exp_now_adj_females, file = paste(tissue.list[i],'_adjusted_exp_FEMALES.rds',sep=""))
    
    
    }}

