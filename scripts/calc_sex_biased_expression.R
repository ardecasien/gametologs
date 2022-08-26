`%!in%` = Negate(`%in%`)

library(limma)
library(edgeR)

###################
## upload meta data
###################

keep_samples = readRDS('gtex_combined_meta.rds')

##################
## calculate sex-biased expression
##################

for (i in 1:length(tissue.list)){
  
  print(paste("Now analyzing:",tissue.list[i]))
  
  # get current data sets 
  # limit to samples with all technical data available
  m_now = subset(keep_samples, SMTSD == tissue.list[i])
  m_now = m_now[complete.cases(m_now[,c('SEX','AGE','SMTSISCH','SMRIN','SMNTRNRT')]),]
  m_now$SEX = as.factor(m_now$SEX)
  
  if(length(m_now$SAMPID) <= 1) {
    print('not enough samples') } else {
    
    # load normalized expression
    y = readRDS(file=paste(tissue.list[i],"_normalized_exp_ALL.rds",sep=""))
    
    m_now = m_now[which(m_now$SAMPID %in% colnames(y)),]
    idx = match(colnames(y),m_now$SAMPID)
    m_now = m_now[idx,]

    # design model
    design = model.matrix(~ 0 + SEX + AGE + SMTSISCH + SMRIN + SMNTRNRT, data=m_now, na.action=na.pass)
    
    # run model
    fit = lmFit(y, design)
    contr = makeContrasts(SEX1 - SEX2, levels = colnames(coef(fit)))
    tmp = contrasts.fit(fit, contr)
    tmp = eBayes(tmp)
    
    # save results for tissue
    saveRDS(tmp, file = paste(tissue.list[i],'_sex_effects.rds',sep=""))
    
    }}

#################
## mashr
#################

library(limma)
library(edgeR)
library(mashr)

Bhat = data.frame()
Shat = data.frame()

for (i in 1:length(tissue.list)){
  
  print(paste("Now analyzing:",tissue.list[i]))
  
  sexnow = readRDS(paste(tissue.list[i],'_sex_effects.rds',sep=""))
  
  beta = sexnow$coefficients
  se.coef = sqrt(sexnow$s2.post) * sexnow$stdev.unscaled
  
  Bhat = merge(Bhat, beta, by = 'row.names', all = T)
  rownames(Bhat) = Bhat$Row.names
  Bhat = Bhat[,-1, drop=FALSE]
  
  Shat = merge(Shat, se.coef, by = 'row.names', all = T)
  rownames(Shat) = Shat$Row.names
  Shat = Shat[,-1, drop=FALSE]
  
}

colnames(Bhat) = colnames(Shat) = short.listnow

saveRDS(Bhat, file = 'Bhat.rds')
saveRDS(Shat, file = 'Shat.rds')

Bhat[is.na(Bhat)] = 0
Shat[is.na(Shat)] = 100

Bhat = as.matrix(Bhat)
Shat = as.matrix(Shat)

# Create the mashr data object
mash.data = mash_set_data(Bhat,Shat)

# Compute canonical covariance matrices
U.c = cov_canonical(mash.data)

m.1by1 = mash_1by1(mash.data)
strong.subset = get_significant_results(m.1by1, thresh = 0.01)
length(strong.subset) 

# The code below computes covariance matrix on a null dataset 
# to learn correlation structure among null tests

# Get random subset (random choose half of all genes)
set.seed(88)
random.subset = sample(1:nrow(Bhat),ceiling(nrow(Bhat)/2))

# Set temporary objects in order to estimate null correlation structure
temp = mash_set_data(Bhat[random.subset,],Shat[random.subset,])
temp.U.c = cov_canonical(temp)
Vhat = estimate_null_correlation_simple(temp)

mash.random = mash_set_data(Bhat[random.subset,],Shat[random.subset,],V=Vhat)
mash.strong = mash_set_data(Bhat[strong.subset,],Shat[strong.subset,], V=Vhat)

# Perform PCA and extreme deconvolution to obtain data-driven covariances
U.pca = cov_pca(mash.strong,5)
U.ed = cov_ed(mash.strong, U.pca)

# Fit mash model
U.c = cov_canonical(mash.random)
m.r = mash(mash.random, Ulist = c(U.ed,U.c), outputlevel = 1)
m = mash(mash.data, g=get_fitted_g(m.r), fixg=TRUE)

saveRDS(m, 'mashr_sex.rds')

##############
## plot mash results
##############

m = readRDS(file = 'mashr_sex.rds')

library(mashr)
library(ggplot2)

mash.beta = get_pm(m)
mash.lfsr = get_lfsr(m)

write.csv(mash.beta, file = 'mash.beta.csv')
write.csv(mash.lfsr, file = 'mash.lfsr.csv')

print(length(get_significant_results(m, thresh = 0.05))) # 10463

bard = data.frame(get_estimated_pi(m))
bard$val = rownames(bard)
bard$val = factor(bard$val, levels = bard$val)
ggplot(bard, aes(x=val, y=get_estimated_pi.m.)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  ylab("Estimated mixture proportions") +
  xlab("") +
  theme(axis.text.x = element_text(angle=45, hjust=1,size =7))

print(get_pairwise_sharing(m))

pl = cor(mash.beta)
corrplot::corrplot.mixed(pl)
      
  
