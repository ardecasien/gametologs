library(limma)
library(matrixStats)
library(tidyverse)
library(dcanr)
library(data.table)
library(readr)
library(CePa)
library(stringr)
library(plyr)
library(dplyr)
library(Hmisc)
library(spqn)

`%!in%` = Negate(`%in%`)

###################
## upload meta data
###################

keep_samples = readRDS('gtex_combined_meta.rds')
table(keep_samples$SMTSD, keep_samples$SEX)

keep_samples %>% group_by(SEX) %>% summarise(mean=mean(AGE))
s = data.frame(keep_samples %>% group_by(SMTSD, SEX) %>% summarise(n=n()))
mean(subset(s, SEX == 1)$n, na.rm = T)
mean(subset(s, SEX == 2)$n, na.rm = T)

gam_list = read.csv('gametologs_in_genome.csv')

#######################
# get co-expression 
# get coupled co-expression
# IN MALES
#######################

m = subset(pheno, SEX == 1)

for (i in 1:length(tissue.list)){
  
  print(paste("Now analyzing:",tissue.list[i]))
  
  # get current data set
  m_now = subset(keep_samples, SMTSD == tissue.list[i] & SAMPID2 %in% m$SAMPID2)
  
  # limit to samples with all data available
  m_now = m_now[complete.cases(m_now[,c('AGE','SMTSISCH','SMRIN','SMNTRNRT')]),]

  if(length(m_now$SAMPID) == 0) {
    print('not enough samples') } else {
    
    # load adjusted expression
    exp_now = readRDS(file = paste(tissue.list[i],'_adjusted_exp_MALES.rds',sep=""))
    exp_now = exp_now[,which(colnames(exp_now) %in% m_now$SAMPID)]
    
    print("estimating co-expression")
    
    # get coexpression matrix
    coexp_matrix = cor(t(exp_now), use= 'everything', method = 'spearman')
    # coexp_matrix_sig = rcorr(t(exp_now), type = 'spearman')
    
    # save  
    saveRDS(coexp_matrix, file = paste(tissue.list[i],'_male_coexp.rds',sep=""))

    print("estimating normalized co-expression")
    
    # get normalized coexpression matrix
    avg_exp = rowMeans(exp_now)
    
    # examine the mean-correlation relationship
    # higher co-expression variance for highly exp genes
    pdf(file = paste(tissue.list[i],'_male_mean_correlation_raw.pdf',sep=""))
    plot_signal_condition_exp(coexp_matrix, avg_exp, signal=0)
    plot_signal_condition_exp(coexp_matrix, avg_exp, signal=0.001)
    IQR_spqn_list <- get_IQR_condition_exp(coexp_matrix, avg_exp)
    plot_IQR_condition_exp(IQR_spqn_list)
    #par(mfrow = c(3,3))
    #for(j in c(1:8,10)){
    #  qqplot_condition_exp(coexp_matrix, avg_exp, j, j)}
    IQR_unlist <- unlist(lapply(1:10, function(ii) IQR_spqn_list$IQR_cor_mat[ii, ii:10]))
    plot(rep(IQR_spqn_list$grp_mean, times = 1:10),
         IQR_unlist,
         xlab="min(average(log2CPM))", ylab="IQR", cex.lab=1.5, cex.axis=1.2, col="blue")
    dev.off()
    
    ## normalize co-expression matrix
    ## ngrp = # bins in each row/column to be used to partition the correlation matrix
    ## size_grp = size of the outer bins to be used to approximate the distribution of the inner bins, in order to smooth the normalization
    ## the product of size_grp and ngrp must be equal or larger than than the row/column # of cor_mat
    ## ref_grp = location of the reference bin on the diagonal, whose distribution will be used as target distribution in the normalization
    dim(coexp_matrix)
    cor_m_spqn = normalize_correlation(coexp_matrix, ave_exp=avg_exp, 
                                       ngrp=20, size_grp=1000, ref_grp=18)
    rownames(cor_m_spqn) = colnames(cor_m_spqn) = rownames(coexp_matrix)
    
    # examine new the mean-correlation relationship
    
    pdf(file = paste(tissue.list[i],'_male_mean_correlation_norm.pdf',sep=""))
    plot_signal_condition_exp(cor_m_spqn, avg_exp, signal=0)
    plot_signal_condition_exp(cor_m_spqn, avg_exp, signal=0.001)
    IQR_spqn_list <- get_IQR_condition_exp(cor_m_spqn, avg_exp)
    plot_IQR_condition_exp(IQR_spqn_list)
    #par(mfrow = c(3,3))
    #for(j in c(1:8,10)){
    #  qqplot_condition_exp(cor_m_spqn, avg_exp, j, j)}
    IQR_unlist <- unlist(lapply(1:10, function(ii) IQR_spqn_list$IQR_cor_mat[ii, ii:10]))
    plot(rep(IQR_spqn_list$grp_mean, times = 1:10),
         IQR_unlist,
         xlab="min(average(log2CPM))", ylab="IQR", cex.lab=1.5, cex.axis=1.2, col="blue")
    dev.off()
    
    # save
    saveRDS(cor_m_spqn, file = paste(tissue.list[i],'_male_coexp_norm.rds',sep=""))
    
    print("estimating coupled co-expression")
    
    # coupled coexpression matrices for gametologues
    diag(coexp_matrix) = NA
    coexp_matrix_gams = coexp_matrix[,which(colnames(coexp_matrix) %in% gam_list$ensembl_id)]
    coup_coexp_matrix = cor.pairs(coexp_matrix_gams, cor.method = 'spearman')
    
    diag(cor_m_spqn) = NA
    coexp_matrix_gams_norm = cor_m_spqn[,which(colnames(cor_m_spqn) %in% gam_list$ensembl_id)]
    coup_coexp_matrix_norm = cor.pairs(coexp_matrix_gams_norm, cor.method = 'spearman')
    
    # save  
    saveRDS(coup_coexp_matrix, file = paste(tissue.list[i],'_male_coup_coexp.rds',sep=""))
    saveRDS(coup_coexp_matrix_norm, file = paste(tissue.list[i],'_male_coup_coexp_norm.rds',sep=""))
    
    # remove outputs
    rm(coexp_matrix)
    rm(cor_m_spqn)
    rm(coup_coexp_matrix)
    rm(coup_coexp_matrix_norm)
    
    }}

#######################
# get co-expression 
# IN FEMALES
#######################

f = subset(pheno, SEX == 2)

for (i in 1:length(tissue.list)){
  
  print(paste("Now analyzing:",tissue.list[i]))
  
  # get current data set
  m_now = subset(keep_samples, SMTSD == tissue.list[i] & SAMPID2 %in% f$SAMPID2)
  
  # limit to samples with all technical data available
  m_now = m_now[complete.cases(m_now[,c('AGE','SMTSISCH','SMRIN','SMNTRNRT')]),]
  
  if(length(m_now$SAMPID) == 0) {
    print('not enough samples') } else {
      
      # load adjusted expression
      exp_now = readRDS(file = paste(tissue.list[i],'_adjusted_exp_FEMALES.rds',sep=""))
      exp_now = exp_now[,which(colnames(exp_now) %in% m_now$SAMPID)]
      
      # get coexpression matrix
      coexp_matrix = cor(t(exp_now), use= 'everything', method = 'spearman')
      
      # save coexpression 
      saveRDS(coexp_matrix, file = paste(tissue.list[i],'_female_coexp.rds',sep=""))
      
      # remove outputs
      rm(coexp_matrix)

    }}

    

  
