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
library(DescTools)
library(egg)
library(stringr)
library(moments)
library(splitstackshape)
library(pvclust)
library(tidytext)
library(ggrepel)
library(reshape2)
library(egg)

#### load ####

`%!in%` = Negate(`%in%`)

gam_list = read.csv('gametologs_in_genome.csv')
gam_list$pair = as.factor(gam_list$pair)
gams = data.frame()
for(i in 1:length(levels(gam_list$pair))){
  gams[i,1] = levels(gam_list$pair)[i]
  xgene = subset(gam_list, pair == levels(gam_list$pair)[i] & Gametolog == 'X')$common_name
  ygene = subset(gam_list, pair == levels(gam_list$pair)[i] & Gametolog == 'Y')$common_name
  gams[i,2] = paste(xgene, ygene, sep = " & ")
}
colnames(gams) = c('pair','gene')
gams

# sample attributes (N=22951)
attrib = read.delim('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
head(attrib)
keep_samples = subset(attrib, SMAFRZE == "RNASEQ")

#### end ####

#### combine files #### 

out_xy = data.frame()
out_xy_sub = data.frame()
out_xy_mean = data.frame()
out_scdced = data.frame()
resamp_xy = data.frame()

out_xy_auto = data.frame()
out_xy_mean_auto = data.frame()
out_scdced_auto = data.frame()
resamp_xy_auto = data.frame()

for (i in 1:length(tissue.list)) {
  
  print(tissue.list[i])
  
  out_xy_now = readRDS(paste(nospace.list[i],'_xy_coupling_res.rds',sep = ""))
  out_xy = rbind(out_xy, out_xy_now)

  out_xy_now = readRDS(paste(nospace.list[i],'_xy_coupling_res_subsamp.rds',sep = ""))
  out_xy_now$gene = rownames(out_xy_now)
  rownames(out_xy_now) = NULL
  out_xy_sub = rbind(out_xy_sub, out_xy_now)

  out_xy_now_mean = readRDS(paste(nospace.list[i],'_xy_coupling_mean_res.rds',sep = ""))
  out_xy_now_mean$tissue = nospace.list[i]
  out_xy_mean = rbind(out_xy_mean, out_xy_now_mean)

  out_scdced_now = readRDS(paste(nospace.list[i],'_sexchr_dep_ced_res.rds',sep = ""))
  out_scdced = rbind(out_scdced, out_scdced_now)

  out_resamp = readRDS(paste(nospace.list[i],'_xy_resamp.rds',sep = ""))
  resamp_xy = rbind(resamp_xy, out_resamp)
  
  out_xy_now = readRDS(paste(nospace.list[i],'_xy_coupling_res_auto.rds',sep = ""))
  out_xy_auto = rbind(out_xy_auto, out_xy_now)
  
  out_xy_now_mean = readRDS(paste(nospace.list[i],'_xy_coupling_mean_res_auto.rds',sep = ""))
  out_xy_now_mean$tissue = nospace.list[i]
  out_xy_mean_auto = rbind(out_xy_mean_auto, out_xy_now_mean)
  
  out_scdced_now = readRDS(paste(nospace.list[i],'_sexchr_dep_ced_res_auto.rds',sep = ""))
  out_scdced_auto = rbind(out_scdced_auto, out_scdced_now)
  
  out_resamp = readRDS(paste(nospace.list[i],'_xy_resamp_auto.rds',sep = ""))
  resamp_xy_auto = rbind(resamp_xy_auto, out_resamp)
  

  }

out_scdced = merge(out_scdced, gams, by = 'pair')
out_xy = merge(out_xy, gams, by = 'pair')
out_xy_sub = merge(out_xy_sub, gams, by = 'pair')
out_xy_sub$gene.x = out_xy$gene.x
resamp_xy = merge(resamp_xy, gams, by = 'pair')

out_scdced_auto = merge(out_scdced_auto, gams, by = 'pair')
out_xy_auto = merge(out_xy_auto, gams, by = 'pair')
resamp_xy_auto = merge(resamp_xy_auto, gams, by = 'pair')

saveRDS(out_xy, file = 'out_xy.rds')
saveRDS(out_xy_sub, file = 'out_xy_sub.rds') 
saveRDS(out_xy_mean, file = 'out_xy_mean.rds')
saveRDS(out_scdced, file = 'out_scdced.rds')
saveRDS(resamp_xy, file = 'resamp_xy.rds')

saveRDS(out_xy_auto, file = 'out_xy_auto.rds')
saveRDS(out_xy_mean_auto, file = 'out_xy_mean_auto.rds')
saveRDS(out_scdced_auto, file = 'out_scdced_auto.rds')
saveRDS(resamp_xy_auto, file = 'resamp_xy_auto.rds')

#### end #### 

#### load and output to Tables S7 and S9 #### 

out_xy = readRDS(file = 'out_xy.rds') # per gene per tissue/pair # signed diff coupling & p-value
out_xy_sub = readRDS(file = 'out_xy_sub.rds') # above, subsampled to N=66 males per region
out_xy_mean = readRDS(file = 'out_xy_mean.rds') # per gene per tissue # signed diff coupling & p-value
out_scdced = readRDS(file = 'out_scdced.rds') # per tissue/pair # signed/abs diff coupling & p-value
resamp_xy = readRDS(file = 'resamp_xy.rds') # per tissue/pair # resample signed/abs diff coupling 

out_xy = readRDS(file = 'out_xy_auto.rds') # per gene per tissue/pair # signed diff coupling & p-value
out_xy_mean = readRDS(file = 'out_xy_mean_auto.rds') # per gene per tissue # signed diff coupling & p-value
out_scdced = readRDS(file = 'out_scdced_auto.rds') # per tissue/pair # signed/abs diff coupling & p-value
resamp_xy = readRDS(file = 'resamp_xy_auto.rds') # per tissue/pair # resample signed/abs diff coupling 

# Table S7

out_xy = out_xy[complete.cases(out_xy),]
length(unique(paste(out_xy$region, out_xy$gene.y)))
out = dcast(out_xy, gene.x ~ gene.y + region, value.var = 'diff')
write.csv(out, file = 'out_xy.csv') # Supplementary Table 7

# Table S9

out_xy_sub = out_xy_sub[complete.cases(out_xy_sub),]
length(unique(paste(out_xy_sub$region, out_xy_sub$gene.y)))
out = dcast(out_xy_sub, gene.x ~ gene.y + region, value.var = 'diff')
write.csv(out, file = 'out_xy_sub.csv') # Supplementary Table 9

#### end #### 

#### estimate CI and p-values ####

resamp_xy = resamp_xy[complete.cases(resamp_xy),]
resamp_xy$meanit[which(resamp_xy$meanit == Inf | resamp_xy$meanit == -Inf)] = NA
resamp_xy = resamp_xy[complete.cases(resamp_xy),]
cival = data.frame()
for(i in 1:length(nospace.list)){
  for(j in 1:length(gams$gene)){
    now = subset(resamp_xy, tissue == nospace.list[i] & gene == gams$gene[j])
    m = mean(now$meanit, na.rm = T)
    sd = sd(now$meanit, na.rm = T)
    qCI = quantile(now$meanit, probs = c(0.025, 0.975))
    lowCI = m - 1.96*sd
    highCI = m + 1.96*sd
    outnow = data.frame(tissue = nospace.list[i], gene = gams$gene[j],
                        mean = m,
                        lowCI = lowCI, highCI = highCI)
    cival = rbind(cival, outnow)
  }}
cival = cival[complete.cases(cival),]
cival$sig = ifelse(cival$lowCI < 0 & cival$highCI < 0, 'sig', 'NS')
cival$sig = ifelse(cival$lowCI > 0 & cival$highCI > 0, 'sig', cival$sig)
table(cival$sig)
cival$se = (cival$highCI - cival$lowCI)/(2 * 1.96)
cival$z = abs(cival$mean/cival$se)
cival$p = exp(-0.717*cival$z - 0.416*(cival$z^2))
cival$padj = p.adjust(cival$p)
table(cival$padj<0.05)

cival_mean = data.frame()
for(i in 1:length(nospace.list)){
  now = subset(resamp_xy, tissue == nospace.list[i])
  m = mean(now$meanit, na.rm = T)
  sd = sd(now$meanit, na.rm = T)
  qCI = quantile(now$meanit, probs = c(0.025, 0.975))
  lowCI = m - 1.96*sd
  highCI = m + 1.96*sd
  outnow = data.frame(tissue = nospace.list[i], 
                      mean = m,
                      lowCI = lowCI, highCI = highCI)
  cival_mean = rbind(cival_mean, outnow)
}
cival_mean = cival_mean[complete.cases(cival_mean),]
cival_mean$sig = ifelse(cival_mean$lowCI < 0 & cival_mean$highCI < 0, 'sig', 'NS')
cival_mean$sig = ifelse(cival_mean$lowCI > 0 & cival_mean$highCI > 0, 'sig', cival_mean$sig)
table(cival_mean$sig)
cival_mean$se = (cival_mean$highCI - cival_mean$lowCI)/(2 * 1.96)
cival_mean$z = abs(cival_mean$mean/cival_mean$se)
cival_mean$p = exp(-0.717*cival_mean$z - 0.416*(cival_mean$z^2))
cival_mean$padj = p.adjust(cival_mean$p)
table(cival_mean$padj<0.05)

resamp_mf = readRDS('resamp_X.rds')
resamp_mf = merge(resamp_mf, subset(gam_list, pair != 14)[,c('ensembl_id','pair')], by = 'ensembl_id')
resamp_xy = readRDS('resamp_xy.rds')

both = merge(resamp_mf, resamp_xy, by = c('pair','tissue','it'))
both = both[complete.cases(both),]

cival2 = data.frame()
for(i in 1:length(nospace.list)){
  print(nospace.list[i])
  for(j in 1:length(gams$gene)){
    now = subset(both, tissue == nospace.list[i] & gene == gams$gene[j])
    
    m_mf = mean(now$diff, na.rm = T)
    sd_mf = sd(now$diff, na.rm = T)
    qCI_mf = quantile(now$diff, probs = c(0.025, 0.975))
    lowCI_mf = m_mf - 1.96*sd_mf
    highCI_mf = m_mf + 1.96*sd_mf
    
    m_xy = mean(now$meanabsit, na.rm = T)
    sd_xy = sd(now$meanabsit, na.rm = T)
    qCI_xy = quantile(now$meanabsit, probs = c(0.025, 0.975))
    lowCI_xy = m_xy - 1.96*sd_xy
    highCI_xy = m_xy + 1.96*sd_xy
    
    outnow = data.frame(tissue = nospace.list[i], gene = gams$gene[j],
                        mean_mf = m_mf, lowCI_mf = lowCI_mf, highCI_mf = highCI_mf,
                        mean_xy = m_xy, lowCI_xy = lowCI_xy, highCI_xy = highCI_xy)
    cival2 = rbind(cival2, outnow)
  }}

cival2$CI_test = ifelse(cival2$lowCI_mf > cival2$highCI_xy, 'mf > xy', 'none')
cival2$CI_test = ifelse(cival2$lowCI_xy > cival2$highCI_mf, 'xy > mf', cival2$CI_test)
table(cival2$CI_test)

#### end #### 

#### Figure 3A 4B S3 | Table S2 ####

out_scdced = readRDS('out_scdced.rds')
out_scdced = out_scdced[complete.cases(out_scdced),]
out_scdced$abs_diff[which(out_scdced$abs_diff == Inf)] = NA
out_scdced$diff[which(out_scdced$diff == Inf | out_scdced$diff == -Inf)] = NA
colnames(cival)[1] = colnames(cival2)[1] = 'region'
out_scdced = merge(out_scdced, cival, by = c('region','gene'))
out_scdced = merge(out_scdced, cival2, by = c('region','gene'))

# Table S2
csv = out_scdced
csv$gene = str_remove_all(csv$gene, pattern = fixed(" "))
write.csv(csv, file = 'sexchrdep_ced.csv')

# summarize
table(out_scdced$diff > 0)
table(out_scdced$region, out_scdced$diff > 0)
table(out_scdced$gene, out_scdced$diff > 0)

summ = data.frame()
for(i in 1:length(nospace.list)){
  now = subset(out_scdced, region == nospace.list[i])
  m = mean(now$diff, na.rm = T)
  ma = mean(now$abs_diff, na.rm = T)
  summ[i,1] = nospace.list[i]
  summ[i,2] = m
  summ[i,3] = ma
}
View(summ)

summ = data.frame()
for(i in 1:length(gams$gene)){
  now = subset(out_scdced, gene == gams$gene[i])
  m = mean(now$diff, na.rm = T)
  ma = max(now$abs_diff, na.rm = T)
  summ[i,1] = gams$gene[i]
  summ[i,2] = m
  summ[i,3] = ma
}
View(summ)

## Figure 3A S3

pl = out_scdced
pl$region = as.factor(pl$region)
levels(pl$region) = short.listnow
length(unique(paste(pl$region, pl$gene)))

mod = aov(abs_diff ~ region + gene, data = pl)
summary(mod)

pl$code = ifelse(pl$CI_test == 'xy > mf', "*", NA)
table(pl$code)
pl$code = ifelse(is.nan(pl$abs_diff), NA, pl$code)
pl$code = ifelse(pl$abs_diff == Inf, NA, pl$code)
pl$code2 = ifelse(pl$region %in% c('Testes','Prostate','Mammary'), "º", NA)

# Figure S3
cl = dcast(gene ~ region, value.var = 'abs_diff', data = pl)
rownames(cl) = cl$gene
cl = cl[,-1]
cl2 = pvclust(cl)
plot(cl2, cex = 0.5)
ord = cl2$hclust$order
names(ord) = levels(pl$region)
pl$region = gplots::reorder.factor(pl$region, new.order = ord)

clusMember = data.frame(cutree(cl2$hclust, h=0.6))
colnames(clusMember) = 'cluster'
clusMember$region = rownames(clusMember)
table(clusMember$cluster)
View(clusMember)

pl = merge(pl, clusMember, by = 'region')
pl$cluster = as.factor(pl$cluster) 
levels(pl$cluster) = c('E','A','C','F','D','B')
clus.ord = c(2,6,3,5,1,4) 
names(clus.ord) = levels(pl$cluster)
pl$cluster = reorder(pl$cluster, new.order = clus.ord)

pl$gene = as.factor(pl$gene)
levels(pl$gene)
levels(pl$gene) = c('DDX3X/Y','EIF1AX/Y','KDM5C/D','UTX/Y','NLGN4X/Y','PCDH11X/Y','PRKX/Y',
                    'RPS4X/Y1','SOX3/SRY','TBL1X/Y','TGIF2LX/Y','TMSB4X/Y','TXLNG/Y','USP9X/Y','ZFX/Y')

ggplot(pl, aes(y = region, x = gene, fill = abs_diff)) +
  geom_tile() +
  scale_fill_gradient2(high=tissue.colors[7],mid="white",low=tissue.colors[14],
                       na.value="gray95", midpoint=0, name = bquote(aCEFD[MX_MY])) +
  facet_grid(cluster ~ ., scale = 'free', space = 'free', drop = TRUE) +
  theme_article() +
  geom_text(label = pl$code, size = 6, vjust = 0.8) +
  geom_text(label = pl$code2, size = 2, vjust = 0.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_blank(),
        strip.text.y.right = element_text(angle = 0),
        legend.position = 'top',
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = 'grey95'))  

# Figure 4B

pl = out_scdced
pl$region = as.factor(pl$region)
levels(pl$region) = short.listnow
pl$diff = ifelse(pl$diff == Inf, NA, pl$diff)
pl$diff = ifelse(pl$diff == -Inf, NA, pl$diff)
levels(pl$region)
length(unique(paste(pl$region, pl$gene)))

# mod = aov(diff ~ region + gene, data = subset(pl, region %in% keep_regions))
mod = aov(diff ~ region + gene, data = pl)
summary(mod)
TukeyHSD(mod)

pl$code = ifelse(pl$padj < 0.05, "*", NA)
pl$code2 = ifelse(pl$p < 0.05, "sig", NA)
table(pl$code)
table(pl$code, pl$diff > 0)
table(pl$p<0.05, pl$diff > 0)

# Figure S3
cl = dcast(gene ~ region, value.var = 'diff', data = pl)
rownames(cl) = cl$gene
cl = cl[,-1]
cl2 = pvclust(cl)
plot(cl2, cex = 0.5)
ord = cl2$hclust$order
names(ord) = levels(pl$region)
pl$region = gplots::reorder.factor(pl$region, new.order = ord)

clusMember = data.frame(cutree(cl2$hclust, h=0.7))
colnames(clusMember) = 'cluster'
clusMember$region = rownames(clusMember)
View(clusMember)

pl = merge(pl, clusMember, by = 'region')
pl$cluster = as.factor(pl$cluster) 
levels(pl$cluster) = c('D','A','C','B')
clus.ord = c(2,4,3,1) 
names(clus.ord) = levels(pl$cluster)
pl$cluster = reorder(pl$cluster, new.order = clus.ord)

pl2 = complete(pl, region, gene)

av = pl2 %>% group_by(region) %>% summarise(diff = mean(diff, na.rm = T))
av$gene = 'mean'
av$pair = NA
av = av[,c(1,3,4,2)]
av[,c(5:25)] = NA
colnames(av) = colnames(pl)
pl2 = rbind(pl2, av)

pl2 = merge(pl2, clusMember, by = 'region')
pl2$cluster.y = as.factor(pl2$cluster.y) 
levels(pl2$cluster.y) = c('D','A','C','B')
clus.ord = c(2,4,3,1) 
names(clus.ord) = levels(pl2$cluster.y)
pl2$cluster.y = reorder(pl2$cluster.y, new.order = clus.ord)
pl2$code3 = ifelse(pl2$gene == 'mean','mean','')

pl2$gene = as.factor(pl2$gene)
levels(pl2$gene)
levels(pl2$gene) = c('DDX3X/Y','EIF1AX/Y','KDM5C/D','UTX/Y','mean',
                    'NLGN4X/Y','PCDH11X/Y','PRKX/Y',
                    'RPS4X/Y1','SOX3/SRY','TBL1X/Y','TGIF2LX/Y',
                    'TMSB4X/Y','TXLNG/Y','USP9X/Y','ZFX/Y')

ggplot(pl2, aes(y = region, x = gene, fill = diff)) +
  geom_tile() +
  geom_tile(data = pl2[which(pl2$code2 == 'sig'),], aes(y = region, x = gene), fill = "transparent", color = "black") +
  scale_fill_gradient2(high=tissue.colors[7],mid="white",low=tissue.colors[14],
                       na.value="grey95", midpoint=0, name = bquote(sCEFD[MX_MY])) +
  theme_article() +
  facet_grid(cluster.y ~ code3, scale = 'free', space = 'free', drop = TRUE) +
  geom_text(label = pl2$code, size = 6, vjust = 0.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text.x = element_blank(),
        strip.text.y.right = element_text(angle = 0),
        legend.position = 'top',
        legend.text = element_text(size = 6),
        axis.title = element_blank())  

#### end #### 

#### Figure 4A ####

pl = out_scdced
pl$region = as.factor(pl$region)
levels(pl$region) = short.listnow
pl$gene = as.factor(pl$gene)
pl$gene = droplevels(pl$gene)

mod = lm(abs_diff ~ diff, data = subset(pl, diff > 0))
su = summary(mod)
i1 = su$coefficients[1,1]
s1 = su$coefficients[2,1]
mod = lm(abs_diff ~ diff, data = subset(pl, diff < 0))
su = summary(mod)
i2 = su$coefficients[1,1]
s2 = su$coefficients[2,1]

pl2 = pl
pl2$code = ifelse(pl2$diff < 0, 'y', 'x')
levels(pl2$gene)
levels(pl2$gene) = c('DDX3X/Y','EIF1AX/Y','KDM5C/D','UTX/Y',
                     'NLGN4X/Y','PCDH11X/Y','PRKX/Y','RPS4X/Y1',
                     'SOX3/SRY','TBL1X/Y','TGIF2LX/Y','TMSB4X/Y',
                     'TXLNG/Y','USP9X/Y','ZFX/Y')
pl2$label = paste(pl2$gene, pl2$region, sep = " in ")
ggplot(pl2, aes(x = diff, y = abs_diff)) +
  geom_point(size = NA) + 
  geom_smooth(method = 'lm') +
  ylab(bquote(aCEFD[MX_MY])) +
  xlab(bquote(sCEFD[MX_MY])) +
  scale_x_continuous(breaks = c(-0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08)) +
  facet_wrap(~code, scales = 'free_y', ncol = 1) +
  geom_point(aes(x = diff, y = abs_diff, color = region, shape = gene), size = 3) +
  scale_color_manual(values = tissue.colors) +
  scale_shape_manual(values=c(8:25)) +
  theme_article() +
  coord_flip() +
  geom_label_repel(data=subset(pl2, region == 'Kidney(Cor)' & gene == 'TMSB4X/Y' |
                                 region == 'Kidney(Cor)' & gene == 'UTX/Y' |
                                 region == 'Testes' & gene == 'ZFX/Y' |
                                 region == 'Prostate' & gene == 'TBL1X/Y'), 
                   aes(x = diff, y = abs_diff, label = label),
                   size = 2,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme(axis.text = element_text(size = 14),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.box = 'vertical',
        legend.key.size = unit(0, 'lines'),
        strip.text.x = element_blank(),
        axis.title = element_text(size = 14))

#### end #### 

#### Figure 3C ####

sim = read.csv('similarity_measures.csv')
colnames(sim)[2] = 'gene'

out_scdced = readRDS('out_scdced.rds')
out_scdced$abs_diff[which(out_scdced$abs_diff == Inf)] = NA
out_scdced = out_scdced[complete.cases(out_scdced),]
table(out_scdced$gene)

comp = out_scdced %>% group_by(gene) %>% summarise(med = median(abs_diff, na.rm = T))

comp$med[which(comp$med == -Inf)] = NA
combo = merge(sim, comp, by = 'gene')
combo = subset(combo, gene != 'SOX3 & SRY')

cor.test(combo$med, 1-combo$prot_sim, method = 'pearson')
cor.test(combo$med, 1-combo$prot_sim, method = 'spearman')
cor.test(combo$med, 1-combo$DNA_sim, method = 'pearson')
cor.test(combo$med, 1-combo$DNA_sim, method = 'spearman')
cor.test(combo$med, 1-combo$prop_match_nongapped, method = 'pearson')
cor.test(combo$med, 1-combo$prop_match_nongapped, method = 'spearman')
cor.test(combo$med, 1-combo$prop_match_overall, method = 'pearson')
cor.test(combo$med, 1-combo$prop_match_overall, method = 'spearman')

combo$prom.div1 = 1 - combo$prop_match_nongapped
combo$prom.div2 = 1 - combo$prop_match_overall
combo$prot.div = 1 - combo$prot_sim
combo$dna.div = 1 - combo$DNA_sim

mod.prom1 = lm(med ~ prom.div1, data = combo)
mod.prom2 = lm(med ~ prom.div2, data = combo)
mod.prot = lm(med ~ prot.div, data = combo)
mod.dna = lm(med ~ dna.div, data = combo)

AIC(mod.prom1, mod.prom2, mod.prot, mod.dna)
BIC(mod.prom1, mod.prom2, mod.prot, mod.dna)

combo$gene = as.factor(combo$gene)
levels(combo$gene)
levels(combo$gene) = c('DDX3X/Y','EIF1AX/Y','KDM5C/D','UTX/Y',
                       'NLGN4X/Y','PCDH11X/Y','PRKX/Y','RPS4X/Y1',
                       'TBL1X/Y','TGIF2LX/Y','TMSB4X/Y',
                       'TXLNG/Y','USP9X/Y','ZFX/Y')

ggplot(combo, aes(x = 1-prot_sim, y = med, label = gene)) +
  geom_point(size = 4) +
  geom_smooth(method = 'lm') +
  xlab('Protein Divergence') +
  ylab(bquote("Mean" ~ aCEFD[MX_MY])) +
  theme_article() +
  annotate(geom="text", x=0.05, y=0.35, label="r = 0.46 \n p = 0.10", size = 5) +
  geom_label_repel(size = 4,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 18))

ggplot(combo, aes(x = 1-prop_match_nongapped, y = med, label = gene)) +
  geom_point(size = 4) +
  geom_smooth(method = 'lm') +
  xlab('Promoter Divergence') +
  ylab(bquote("Mean" ~ aCEFD[MX_MY])) +
  theme_article() +
  annotate(geom="text", x=0.06, y=0.25, label="r = 0.59 \n p = 0.03", size = 5) +
  geom_label_repel(size = 4,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 18))

#### end ####

#### Figure S7 ####

pl = out_scdced
pl = pl[complete.cases(pl),]
pl$region = as.factor(pl$region)
levels(pl$region) = nospace.list
pl$gene = as.factor(pl$gene)
pl$gene = droplevels(pl$gene)

# get gam expression (adjusted expression)

rel_exp = data.frame()
for (i in 1:length(tissue.list)){
  exp_now = readRDS(file = paste('/Users/decasienar/Desktop/Papers/GTEX/v8_remapped/',tissue.list[i],'_adjusted_exp_MALES.rds',sep=""))
    for(j in 1:length(levels(gam_list$pair))){
      x_gams = exp_now[which(rownames(exp_now) == subset(gam_list, pair == levels(gam_list$pair)[j] & Gametolog == 'X')$ensembl_id),]
      y_gams = exp_now[which(rownames(exp_now) == subset(gam_list, pair == levels(gam_list$pair)[j] & Gametolog == 'Y')$ensembl_id),]
      xname = subset(gam_list, pair == levels(gam_list$pair)[j] & Gametolog == 'X')$common_name
      yname = subset(gam_list, pair == levels(gam_list$pair)[j] & Gametolog == 'Y')$common_name
      rat = y_gams / x_gams
      rat[which(rat == Inf)] = x_gams[which(rat == Inf)] = y_gams[which(rat == Inf)] = NA
      rat[which(rat == -Inf)] = x_gams[which(rat == -Inf)] = y_gams[which(rat == -Inf)] = NA
      rat[which(rat == 0)] = x_gams[which(rat == 0)] = y_gams[which(rat == 0)] = NA
      xy_out = data.frame(tissue = nospace.list[i],
                          gene = paste(xname, "&", yname),
                          xy = mean(rat, na.rm = T),
                          sdxy = sd(rat, na.rm = T),
                          sdx = sd(x_gams, na.rm = T),
                          sdy = sd(y_gams, na.rm = T)
                          )
      rel_exp = rbind(rel_exp, xy_out)
    }}

# get gam expression (unadjusted tpm)

txi.gene = readRDS('gtex_kallisto_genes.rds')
tpm = txi.gene$abundance
rm(txi.gene)
keep_samples = readRDS('gtex_combined_meta.rds')
rel_exp = data.frame()

for (i in 1:length(tissue.list)){
  
  print(paste("Now analyzing:",tissue.list[i]))
  m_now = subset(keep_samples, SMTSD == tissue.list[i] & SEX == 1)
  m_now = m_now[complete.cases(m_now[,c('SEX','AGE','SMTSISCH','SMRIN','SMNTRNRT')]),]
  
  if(length(m_now$SAMPID) == 0) {
    print('not enough samples') } else {
      
      for(j in 1:length(levels(gam_list$pair))){
      x_tpm = tpm[which(rownames(tpm) == subset(gam_list, pair == levels(gam_list$pair)[j] & Gametolog == 'X')$ensembl_id),which(colnames(tpm) %in% m_now$SAMPID)]
      y_tpm = tpm[which(rownames(tpm) == subset(gam_list, pair == levels(gam_list$pair)[j] & Gametolog == 'Y')$ensembl_id),which(colnames(tpm) %in% m_now$SAMPID)]
      xname = subset(gam_list, pair == levels(gam_list$pair)[j] & Gametolog == 'X')$common_name
      yname = subset(gam_list, pair == levels(gam_list$pair)[j] & Gametolog == 'Y')$common_name
      rat = y_tpm / x_tpm
      rat[which(rat == Inf)] = x_tpm[which(rat == Inf)] = y_tpm[which(rat == Inf)] = NA
      rat[which(rat == -Inf)] = x_tpm[which(rat == -Inf)] = y_tpm[which(rat == -Inf)] = NA
      rat[which(rat == 0)] = x_tpm[which(rat == 0)] = y_tpm[which(rat == 0)] = NA
      xy_out = data.frame(tissue = nospace.list[i],
                          gene = paste(xname, "&", yname),
                          xy = median(rat, na.rm = T),
                          sdxy = sd(rat, na.rm = T),
                          sdx = sd(x_tpm, na.rm = T),
                          sdy = sd(y_tpm, na.rm = T),
                          meanx = mean(x_tpm, na.rm = T),
                          meany = mean(y_tpm, na.rm = T))
      rel_exp = rbind(rel_exp, xy_out)
}}}

rel_exp$key = paste(rel_exp$tissue, rel_exp$gene)
rel_exp$sd_diff = rel_exp$sdx - rel_exp$sdy
rel_exp$mean_diff = rel_exp$meanx - rel_exp$meany
rel_exp = subset(rel_exp, gene != 'AMELX & AMELY')
rel_exp$tissue = as.factor(rel_exp$tissue)
levels(rel_exp$tissue) = short.listnow

# get gam expression (unadjusted tpm - subtraction)

txi.gene = readRDS('gtex_kallisto_genes.rds')
tpm = txi.gene$abundance
rm(txi.gene)
keep_samples = readRDS('gtex_combined_meta.rds')
rel_exp = data.frame()

for (i in 1:length(tissue.list)){
  
  print(paste("Now analyzing:",tissue.list[i]))
  m_now = subset(keep_samples, SMTSD == tissue.list[i] & SEX == 1)
  m_now = m_now[complete.cases(m_now[,c('SEX','AGE','SMTSISCH','SMRIN','SMNTRNRT')]),]
  
  if(length(m_now$SAMPID) == 0) {
    print('not enough samples') } else {
      
      for(j in 1:length(levels(gam_list$pair))){
        x_tpm = tpm[which(rownames(tpm) == subset(gam_list, pair == levels(gam_list$pair)[j] & Gametolog == 'X')$ensembl_id),which(colnames(tpm) %in% m_now$SAMPID)]
        y_tpm = tpm[which(rownames(tpm) == subset(gam_list, pair == levels(gam_list$pair)[j] & Gametolog == 'Y')$ensembl_id),which(colnames(tpm) %in% m_now$SAMPID)]
        xname = subset(gam_list, pair == levels(gam_list$pair)[j] & Gametolog == 'X')$common_name
        yname = subset(gam_list, pair == levels(gam_list$pair)[j] & Gametolog == 'Y')$common_name
        rat = x_tpm - y_tpm
        rat[which(rat == Inf)] = x_tpm[which(rat == Inf)] = y_tpm[which(rat == Inf)] = NA
        rat[which(rat == -Inf)] = x_tpm[which(rat == -Inf)] = y_tpm[which(rat == -Inf)] = NA
        rat[which(rat == 0)] = x_tpm[which(rat == 0)] = y_tpm[which(rat == 0)] = NA
        xy_out = data.frame(tissue = nospace.list[i],
                            gene = paste(xname, "&", yname),
                            xy = rat,
                            indiv = m_now$SAMPID2)
        rel_exp = rbind(rel_exp, xy_out)
      }}}

rel_exp = subset(rel_exp, gene != 'AMELX & AMELY')
levels(rel_exp$tissue) = short.listnow
rel_exp2 = dcast(rel_exp, formula = indiv ~ gene + tissue, value.var = 'xy')
rel_exp2 = data.frame(rel_exp2)
rownames(rel_exp2) = rel_exp2$indiv
rel_exp2 = rel_exp2[,-1]

c1 = cor(rel_exp2, use = 'pairwise.complete.obs')

levels(out_xy$region) = short.listnow
out_xy = subset(out_xy, gene.y != 'AMELX & AMELY')
out_xy2 = dcast(out_xy, formula = gene.x ~ gene.y + region, value.var = 'diff')
out_xy2 = data.frame(out_xy2)
rownames(out_xy2) = out_xy2$gene.x
out_xy2 = out_xy2[,-1]
c2 = cor(out_xy2, use = 'pairwise.complete.obs')
c2 = c2[which(colnames(c2) %in% colnames(c1)), which(colnames(c2) %in% colnames(c1))]
idx = match(colnames(c2), colnames(c1))
c2 = c2[idx, idx]

c1p = c1[rowSums(is.na(c1)) != ncol(c1),rowSums(is.na(c1)) != ncol(c1)]
corrplot::corrplot(c1p, tl.cex = 0.2, order = 'hclust', addgrid.col = NA)
c1l = c1
c1l[lower.tri(c1l, diag = T)] = 1000
c1l = melt(c1l)
c1l = subset(c1l, value != 1000)
c1l = subset(c1l, value != 1)
c1l = subset(c1l, value != -1)

c2p = c2[rowSums(is.na(c2)) != ncol(c2),rowSums(is.na(c2)) != ncol(c2)]
corrplot::corrplot(c2p, tl.cex = 0.2, order = 'hclust', addgrid.col = NA)
c2l = c2
c2l[lower.tri(c2l,  diag = T)] = 1000
c2l = melt(c2l)
c2l = subset(c2l, value != 1000)

pl = merge(c1l, c2l, by = c('Var1', 'Var2'))
pl = separate(data = pl, col = Var1, into = c("pair1", "region1"), sep = "_")
pl = separate(data = pl, col = Var2, into = c("pair2", "region2"), sep = "_")
pl$pairmatch = ifelse(pl$pair1 == pl$pair2, 'same', 'diff')
pl$regionmatch = ifelse(pl$region1 == pl$region2, 'same', 'diff')

cor.test(pl$value.x, pl$value.y)
cor.test(subset(pl, pairmatch == 'same')$value.x, subset(pl, pairmatch == 'same')$value.y)
cor.test(subset(pl, regionmatch == 'same')$value.x, subset(pl, regionmatch == 'same')$value.y)

ggplot(pl, aes(x = value.x, y = value.y)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm', se = T) +
  theme_article()

ggplot(pl, aes(x = value.x, y = value.y, color = pairmatch)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm', se = T) +
  theme_article()

ggplot(subset(pl, pairmatch == 'same'), aes(x = value.x, y = value.y)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm', se = T) +
  theme_article()

ggplot(subset(pl, regionmatch == 'same'), aes(x = value.x, y = value.y)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm', se = T) +
  theme_article()

#rel_exp = subset(rel_exp, gene != 'TMSB4X & TMSB4Y')

ggplot(rel_exp, aes(x = gene, y = tissue, fill = log2(xy))) +
  geom_tile() +
  theme_article() +
  scale_fill_gradient2(high=tissue.colors[14],mid="white",low=tissue.colors[7],
                       na.value="grey95", midpoint=0, name = 'log2(Y/X expression') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = 'right')

ggplot(rel_exp, aes(x = gene, y = tissue, fill = sd_diff)) +
  geom_tile() +
  theme_article() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = 'right')

ggplot(rel_exp, aes(x = gene, y = tissue, fill = mean_diff)) +
  geom_tile() +
  theme_article() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = 'right')

cor.test(rel_exp$mean_diff, (rel_exp$sd_diff)^2, method = 'spearman')
ggplot(rel_exp, aes(x = mean_diff, y = sd_diff^2)) +
  geom_point(size = NA) +
  geom_smooth(method = 'lm', color = 'black') + 
  geom_point(aes(x = mean_diff, y = sd_diff^2, color = tissue, shape = gene), size = 3, alpha = 0.5) +
  geom_smooth(aes(x = mean_diff, y = sd_diff^2, color = tissue), method = 'lm', 
              se = FALSE, alpha = 1, size = 1) +
  scale_color_manual(values = tissue.colors) +
  #scale_shape_manual(values=c(8:25)) +
  scale_shape_manual(values=c(8:21,23:25)) +
  guides(shape=guide_legend(ncol=4,bycol=TRUE)) +
  geom_point(size = NA) +
  geom_smooth(method = 'lm', color = 'black',  size = 1.5) +
  xlab('X vs. Y \n mean expression') +
  ylab('X vs. Y \n expression variance') +
  theme_article() +
  theme(axis.text = element_text(size = 14),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.box = 'vertical',
        legend.key.size = unit(0, 'lines'),
        axis.title = element_text(size = 14))


pl$key = paste(pl$region, pl$gene)

combo = merge(rel_exp, pl, by = 'key')
combo = subset(combo, gene.x != 'TMSB4X & TMSB4Y')
#combo = subset(combo, key != 'SmallIntestineTerminalIleum PCDH11X & PCDH11Y')

combo$tissue = as.factor(combo$tissue)
levels(combo$tissue) = short.listnow

cor.test(combo$diff, log(combo$xy))

ggplot(combo, aes(x = log(xy), y = diff)) +
#ggplot(combo, aes(x = xy, y = diff)) +
  
  geom_point(size = NA) +
  geom_smooth(method = 'lm', color = 'black') + 
  geom_point(aes(x = log(xy), y = diff, color = tissue, shape = gene.x), size = 3, alpha = 0.5) +
  geom_smooth(aes(x = log(xy), y = diff, color = tissue), method = 'lm', 
              se = FALSE, alpha = 1, size = 1) +
  scale_color_manual(values = tissue.colors) +
  #scale_shape_manual(values=c(8:25)) +
  scale_shape_manual(values=c(8:21,23:25)) +  
  guides(shape=guide_legend(ncol=4,bycol=TRUE)) +
  geom_point(size = NA) +
  geom_smooth(method = 'lm', color = 'black',  size = 1.5) +
  xlab('Median Y/X TPM ratio (log)') +
  ylab(bquote(sCEFD[MX_MY])) +
  theme_article() +
  theme(axis.text = element_text(size = 14),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.box = 'vertical',
        legend.key.size = unit(0, 'lines'),
        axis.title = element_text(size = 14))

cor.test(combo$abs_diff, log(combo$xy))

ggplot(combo, aes(x = log(xy), y = abs_diff)) +
  #ggplot(combo, aes(x = xy, y = diff)) +
  
  geom_point(size = NA) +
  geom_smooth(method = 'lm', color = 'black') + 
  geom_point(aes(x = log(xy), y = abs_diff, color = tissue, shape = gene.x), size = 3, alpha = 0.5) +
  geom_smooth(aes(x = log(xy), y = abs_diff, color = tissue), method = 'lm', 
              se = FALSE, alpha = 1, size = 1) +
  scale_color_manual(values = tissue.colors) +
  #scale_shape_manual(values=c(8:25)) +
  scale_shape_manual(values=c(8:21,23:25)) +
  guides(shape=guide_legend(ncol=4,bycol=TRUE)) +
  geom_point(size = NA) +
  geom_smooth(method = 'lm', color = 'black',  size = 1.5) +
  xlab('Median Y/X TPM ratio (log)') +
  ylab(bquote(aCEFD[MX_MY])) +
  theme_article() +
  theme(axis.text = element_text(size = 14),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.box = 'vertical',
        legend.key.size = unit(0, 'lines'),
        axis.title = element_text(size = 14))

ggplot(combo, aes(x = log(xy), y = log(sdxy))) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm') +
  xlab('Median Y/X TPM ratio') +
  ylab('SD Y/X TPM ratio') +
  theme_article()

ggplot(combo, aes(x = log(xy), y = log(sdy))) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm') +
  xlab('Median Y/X TPM ratio') +
  ylab('SD Y expression') +
  theme_article()

ggplot(combo, aes(x = log(xy), y = log(sdx))) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm') +
  xlab('Median Y/X TPM ratio') +
  ylab('SD X expression') +
  theme_article()

cor.test(log(combo$sdy), y = combo$diff, method = 'spearman')
ggplot(combo, aes(x = log(sdy), y = diff)) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm') +
  xlab('SD Y expression') +
  ylab(bquote(sCEFD[MX_MY])) +
  theme_article()

cor.test(log(combo$sdx), y = combo$diff, method = 'spearman')
ggplot(combo, aes(x = log(sdx), y = diff)) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm') +
  xlab('SD X expression') +
  ylab(bquote(sCEFD[MX_MY])) +
  theme_article()

cor.test(log(combo$meany), y = combo$diff, method = 'spearman')
ggplot(combo, aes(x = log(meany), y = diff)) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm') +
  xlab('Mean Y expression') +
  ylab(bquote(sCEFD[MX_MY])) +
  theme_article()

cor.test(log(combo$meanx), y = combo$diff, method = 'spearman')
ggplot(combo, aes(x = log(meanx), y = diff)) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm') +
  xlab('Mean X expression') +
  ylab(bquote(sCEFD[MX_MY])) +
  theme_article()

cor.test(log(combo$sdx), y = log(combo$meanx), method = 'spearman')
ggplot(combo, aes(x = log(sdx), y = log(meanx))) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm') +
  xlab('SD Y expression') +
  ylab(bquote(sCEFD[MX_MY])) +
  theme_article()

ggplot(combo, aes(x = log(sdy), y = log(meany))) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm') +
  xlab('SD Y expression') +
  ylab(bquote(sCEFD[MX_MY])) +
  theme_article()

## higher Y vs X expression = higher X vs Y coupling
## higher Y vs X expression = more variable Y vs X expression = lower Y vs X coupling

#### end #### 

#### Figure S4 | compare all vs. autosomal only ####

out_scdced = readRDS(file = 'out_scdced.rds') 
out_scdced_auto = readRDS(file = 'out_scdced_auto.rds') 

pl = cbind(out_scdced[,c('pair','region','diff','abs_diff')], out_scdced_auto[,c('diff','abs_diff')])
colnames(pl) = c('pair','region','diff_all','abs_diff_all','diff_auto','abs_diff_auto')
pl$region = as.factor(pl$region)
levels(pl$region) = short.listnow

ggplot(pl, aes(x = diff_all, y = diff_auto, color = region)) +
  geom_point() +
  scale_color_manual(values = tissue.colors) +
  theme_article() +
  xlab('All genes') + ylab('Autosomal genes only') +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = 'lm', color = 'blue') +
  theme(legend.position = 'bottom',
        legend.title = element_blank())

ggplot(pl, aes(x = diff_all - diff_auto)) +
  geom_histogram(fill = tissue.colors[30], color = tissue.colors[40]) +
  theme_article() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlab('all genes - autosomal only') + ylab('') +
  theme(legend.position = 'bottom',
        legend.title = element_blank())

ggplot(pl, aes(x = abs_diff_all, y = abs_diff_auto, color = region)) +
  geom_point() +
  scale_color_manual(values = tissue.colors) +
  theme_article() +
  xlab('All genes') + ylab('Autosomal genes only') +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = 'lm', color = 'blue') +
  theme(legend.position = 'bottom',
        legend.title = element_blank())

ggplot(pl, aes(x = abs_diff_all - abs_diff_auto)) +
  geom_histogram(fill = tissue.colors[30], color = tissue.colors[40]) +
  theme_article() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlab('all genes - autosomal only') + ylab('') +
  theme(legend.position = 'bottom',
        legend.title = element_blank())

#### end #### 







