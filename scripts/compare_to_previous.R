setwd("~/Desktop/Papers/GTEX/v8_remapped")

library(ggplot2)
library(matrixStats)
library(reshape2)

########################
## TPM
########################

## load current data TPM

txi.gene = readRDS('gtex_kallisto_genes.rds')
tpm_all = txi.gene$abundance

keep_samples = readRDS('gtex_combined_meta.rds')
gam_list = read.csv('gametologs_in_genome.csv')
gams = unique(gam_list[,c('ensembl_id','common_name')])

m = subset(keep_samples, SEX == 1)
mtpm = tpm_all[which(rownames(tpm_all) %in% gam_list$ensembl_id),which(colnames(tpm_all) %in% m$SAMPID)]

out = data.frame()

for(i in 1:length(tissue.list)){
  
  print(tissue.list[i])
  samples_now = subset(m, SMTSD == tissue.list[i])
  tpmnow = mtpm[,which(colnames(mtpm) %in% samples_now$SAMPID)]
  medtpm = data.frame(rowMedians(tpmnow))
  colnames(medtpm) = 'median_tpm'
  medtpm$ensembl_id = rownames(tpmnow)
  medtpm$region = tissue.list[i]
  out = rbind(out, medtpm)

}

out = merge(out, gams, by = 'ensembl_id')
View(out)
out$key = paste(out$region, out$common_name, sep = "")

## load Godfrey data adjusted TPM

god_tpm = read.csv('godfrey_median_tpm.csv')
god_tpm = melt(god_tpm)
colnames(god_tpm) = c('region','common_name','median_tpm')
god_tpm$key = paste(god_tpm$region, god_tpm$common_name, sep = "")

## compare TPM data

comp = merge(out, god_tpm, by = 'key')

ggplot(comp, aes(x = median_tpm.x, y = median_tpm.y)) +
  facet_wrap(~common_name.x, scales = 'free') +
  geom_point() +
  xlab("Current Study TPM") +
  ylab("Godfrey et al. 2020 TPM") +
  geom_smooth(method = 'lm') +
  theme_classic()

ggplot(comp, aes(x = median_tpm.x, y = median_tpm.y)) +
  facet_wrap(~region.x, scales = 'free') +
  geom_point() +
  xlab("Current Study TPM") +
  ylab("Godfrey et al. 2020 TPM") +
  geom_smooth(method = 'lm') +
  theme_classic()

stat = data.frame()
tis = levels(as.factor(comp$region.x))
for(i in 1:length(tis)){
  now = subset(comp, region.x == tis[i])
  cc = cor.test(now$median_tpm.x, now$median_tpm.y)
  stat[i,1] = cc$estimate
  stat[i,2] = cc$p.value
  stat[i,3] = tis[i]
}
colnames(stat) = c('rho','p','tissue')
View(stat)
write.csv(stat, file = 'godfrey_comparison_tpm_regions.csv')

stat = data.frame()
tis = levels(as.factor(comp$common_name.x))
for(i in 1:length(tis)){
  now = subset(comp, common_name.x == tis[i])
  cc = cor.test(now$median_tpm.x, now$median_tpm.y)
  stat[i,1] = cc$estimate
  stat[i,2] = cc$p.value
  stat[i,3] = tis[i]
}
colnames(stat) = c('rho','p','tissue')
View(stat)
write.csv(stat, file = 'godfrey_comparison_tpm_gams.csv')

########################
## coexpression
########################

library(reshape2)

# load current coexpression data

coexp_plot = readRDS('coexp_plot_43tissues_norm.rds')
coexp_plot2 = melt(coexp_plot, id = 'region')
coexp_plot2$value = as.numeric(coexp_plot2$value)
coexp_plot2$key = paste(coexp_plot2$region, coexp_plot2$variable, sep = "")

# load godfrey coexpression data

god = read.csv('godfrey_gam_coexp.csv')
short.list.god = c('Adipose(Sub)','Adipose(Vis)','Adrenal','Artery(Aor)','Artery(Cor)','Artery(Tib)','Brain(Amy)','Brain(Cblm1)','Brain(Cort)','Brain(Hip)','Brain(Hyp)','Brain(Caud)','Brain(SBN)','Mammary','Colon(Sig)','Colon(Trans)','Esoph(Muc)','Esoph(Musc)','Heart(Atr)','Heart(Ven)','Liver','Lung','Nerve(Tib)','Pancreas','Pituitary','Prostate','Salivary','Muscle','Skin(NoSun)','Ileum','Spinal(C1)','Spleen','Stomach','Testes','Thyroid')
god$tissue = as.factor(god$tissue)
levels(god$tissue) = short.list.god
god$pair = gsub("/"," & ",god$pair)
god$key = paste(god$tissue, god$pair, sep = "")

god_pairs = data.frame(god %>% group_by(pair) %>% summarise(mean = mean(Spearman.r)))
god_pairs = data.frame(god %>% group_by(pair) %>% summarise(mean = median(Spearman.r)))
View(god_pairs)

# compare coexpression values

comp = merge(coexp_plot2, god, by = 'key')

ggplot(comp, aes(x=value, y=Spearman.r)) +
  geom_point() +
  xlim(c(-0.5,1)) +
  ylim(c(-0.5,1)) +
  geom_smooth(method = 'lm') +
  geom_abline(slope=1, intercept = 0, linetype='dashed') +
  xlab('Current') +
  ylab('Godfrey et al. 2020') +
  theme_classic()

cor.test(comp$value, comp$Spearman.r, method = 'pearson')
cor.test(comp$value, comp$Spearman.r, method = 'spearman')

ggplot(comp, aes(x=value, y=Spearman.r)) +
  geom_point() +
  facet_wrap(~region) +
  xlim(c(-0.5,1)) +
  ylim(c(-0.5,1)) +
  geom_smooth(method = 'lm') +
  geom_abline(slope=1, intercept = 0, linetype='dashed') +
  xlab('Current') +
  ylab('Godfrey et al. 2020') +
  theme_classic()

ggplot(comp, aes(x=value, y=Spearman.r)) +
  geom_point() +
  facet_wrap(~pair) +
  xlim(c(-0.5,1)) +
  ylim(c(-0.5,1)) +
  geom_smooth(method = 'lm') +
  geom_abline(slope=1, intercept = 0, linetype='dashed') +
  xlab('Current') +
  ylab('Godfrey et al. 2020') +
  theme_classic()

comp$region = as.factor(comp$region)
tissues_now = levels(comp$region)
out = data.frame()
for(i in 1:length(tissues_now)){
  now = subset(comp, region == tissues_now[i])
  cc = cor.test(now$value, now$Spearman.r, method = 'pearson')
  out[i,1] = tissues_now[i]
  out[i,2] = cc$estimate
  out[i,3] = cc$p.value
  }
colnames(out) = c("Tissue","Correlation","P-value")
View(out)

table(out$Correlation > 0)
table(out$`P-value` < 0.05)
min(out$Correlation)
max(out$Correlation)
mean(out$Correlation)
median(out$Correlation)

write.csv(out, file = 'godfrey_comparison_correlations.csv')

comp$pair = as.factor(comp$pair)
pair_now = levels(comp$pair)
out = data.frame()
for(i in 1:length(pair_now)){
  now = subset(comp, pair == pair_now[i])
  cc = tryCatch(cor.test(now$value, now$Spearman.r, method = 'pearson'), error=function(err) NA)
  out[i,1] = pair_now[i]
  out[i,2] = tryCatch(cc$estimate, error=function(err) NA)
  out[i,3] = tryCatch(cc$p.value, error=function(err) NA)
}
colnames(out) = c("Tissue","Correlation","P-value")
View(out)

table(out$Correlation > 0)
table(out$`P-value` < 0.05)
min(out$Correlation, na.rm = T)
max(out$Correlation, na.rm = T)
mean(out$Correlation, na.rm = T)
median(out$Correlation, na.rm = T)

write.csv(out, file = 'godfrey_comparison_correlations_pairs.csv')

comp$dif = comp$value - comp$Spearman.r

ggplot(comp, aes(x = pair, y = dif)) +
  geom_boxplot() +
  theme_classic() +
  xlab("") +
  ylab("Co-expression Difference (Current - Godfrey)") + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(comp, aes(x = region, y = dif)) +
  geom_boxplot() +
  theme_classic() +
  xlab("") +
  ylab("Co-expression Difference (Current - Godfrey)") + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

new = data.frame(region = coexp_plot2$region, pair = coexp_plot2$variable, coexp = coexp_plot2$value, study = 'Current')
godf = data.frame(region = god$tissue, pair = god$pair, coexp = god$Spearman.r, study = 'Godfrey')

comp2 = rbind(new, godf)
comp2 = subset(comp2, region %in% short.list.god)

ggplot(comp2, aes(x = study, y = coexp)) +
  geom_boxplot() +
  facet_wrap(~region) +
  xlab('') +
  ylab('Co-expression') +
  theme_classic()

ggplot(comp2, aes(x = study, y = coexp)) +
  geom_boxplot() +
  facet_wrap(~pair) +
  xlab('') +
  ylab('Co-expression') +
  theme_classic()

library(car)

comp2$study = as.factor(comp2$study)

comp2$region = as.factor(comp2$region)
tissues_now = levels(comp2$region)
out = data.frame()
for(i in 1:length(tissues_now)){
  print(tissues_now[i])
  datanow = subset(comp2, region == tissues_now[i])
  tt = t.test(coexp ~ study, data = datanow)
  ll = leveneTest(coexp ~ study, data = datanow)
  bb = bartlett.test(coexp ~ study, data = datanow)
  out[i,1] = tissues_now[i]
  out[i,2] = tt$p.value
  out[i,3] = tt$estimate[1]
  out[i,4] = tt$estimate[2]
  out[i,5] = ll$`Pr(>F)`[1]
  out[i,6] = bb$p.value
  out[i,7] = var(subset(datanow, study == 'Current')$coexp, na.rm = T)
  out[i,8] = var(subset(datanow, study == 'Godfrey')$coexp, na.rm = T)
  out[]
}
colnames(out) = c("Tissue","T-test P","Mean Current","Mean Godfrey",
                  "Levene P","Bartlett P","Variance Current","Variance Godfrey")
View(out)
table(out$`Mean Current` > out$`Mean Godfrey`)
table(out$`Variance Current` > out$`Variance Godfrey`)
write.csv(out, 'godfrey_coexpression_comparison.csv')

comp2$pair = as.factor(comp2$pair)
pair_now = levels(comp2$pair)
out = data.frame()
for(i in 1:length(pair_now)){
  print(pair_now[i])
  datanow = subset(comp2, pair == pair_now[i])
  tt = t.test(coexp ~ study, data = datanow)
  ll = leveneTest(coexp ~ study, data = datanow)
  out[i,1] = pair_now[i]
  out[i,2] = tt$p.value
  out[i,3] = tt$estimate[1]
  out[i,4] = tt$estimate[2]
  out[i,5] = ll$`Pr(>F)`[1]
  out[i,6] = var(subset(datanow, study == 'Current')$coexp, na.rm = T)
  out[i,7] = var(subset(datanow, study == 'Godfrey')$coexp, na.rm = T)
  out[]
}
colnames(out) = c("Tissue","T-test P","Mean Current","Mean Godfrey",
                  "Levene P","Variance Current","Variance Godfrey")
View(out)
table(out$`Mean Current` > out$`Mean Godfrey`)
table(out$`Variance Current` > out$`Variance Godfrey`)

write.csv(out, 'godfrey_coexpression_comparison_pairs.csv')


