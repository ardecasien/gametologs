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
library(pvclust)
library(cluster)
library(ggpubr)

#### combine files ####

out_all = data.frame()
out_gam = data.frame()
mean_p = data.frame()
mean_p_strict = data.frame()
pergam_p = data.frame()
coexp_all = data.frame()
out_all_signed = data.frame()
out_gam_signed = data.frame()
resamp_XY = data.frame()

out_all_auto = data.frame()
out_gam_auto = data.frame()
mean_p_auto = data.frame()
mean_p_strict_auto = data.frame()
pergam_p_auto = data.frame()
coexp_all_auto = data.frame()
out_all_signed_auto = data.frame()
out_gam_signed_auto = data.frame()
resamp_XY_auto = data.frame()

out_all_X = data.frame()
out_gam_X = data.frame()
mean_p_X = data.frame()
pergam_p_X = data.frame()
coexp_all_X = data.frame()
resamp_X = data.frame()

out_all_X_auto = data.frame()
out_gam_X_auto = data.frame()
mean_p_X_auto = data.frame()
pergam_p_X_auto = data.frame()
coexp_all_X_auto = data.frame()
resamp_X_auto = data.frame()
  
for (i in 1:length(nospace.list)) {
  
  print(nospace.list[i])
  
  out_all_now = readRDS(paste(nospace.list[i],'_out_all_norm.rds',sep = ""))
  out_all = rbind(out_all, out_all_now)
  out_gam_now = readRDS(paste(nospace.list[i],'_out_gam_norm.rds',sep = ""))
  out_gam = rbind(out_gam, out_gam_now)
  mean_p_now = readRDS(paste(nospace.list[i],'_mean_p_norm.rds',sep = ""))
  mean_p = rbind(mean_p, mean_p_now)
  mean_p_strict_now = readRDS(paste(nospace.list[i],'_mean_p_strict_norm.rds',sep = ""))
  mean_p_strict = rbind(mean_p_strict, mean_p_strict_now)
  pergam_p_now = readRDS(paste(nospace.list[i],'_pergam_p_norm.rds',sep = ""))
  pergam_p = rbind(pergam_p, pergam_p_now)
  coexp_all_now = readRDS(paste(nospace.list[i],'_coexp_out_norm.rds',sep = ""))
  coexp_all = merge(coexp_all, coexp_all_now, by = 'row.names', all = T)
  rownames(coexp_all) = coexp_all$Row.names
  coexp_all = coexp_all[,-1]
  out_all_now = readRDS(paste(nospace.list[i],'_out_all_signed_norm.rds',sep = ""))
  out_all_signed = rbind(out_all_signed, out_all_now)
  out_gam_now = readRDS(paste(nospace.list[i],'_out_gam_signed_norm.rds',sep = ""))
  out_gam_signed = rbind(out_gam_signed, out_gam_now)
  resamp_now = readRDS(paste(nospace.list[i],'_out_gam_resamp_norm_XY.rds',sep=""))
  resamp_XY = rbind(resamp_XY, resamp_now)
  
  out_all_now = readRDS(paste(nospace.list[i],'_out_all_norm_auto.rds',sep = ""))
  out_all_auto = rbind(out_all_auto, out_all_now)
  out_gam_now = readRDS(paste(nospace.list[i],'_out_gam_norm_auto.rds',sep = ""))
  out_gam_auto = rbind(out_gam_auto, out_gam_now)
  mean_p_now = readRDS(paste(nospace.list[i],'_mean_p_norm_auto.rds',sep = ""))
  mean_p_auto = rbind(mean_p_auto, mean_p_now)
  mean_p_strict_now = readRDS(paste(nospace.list[i],'_mean_p_strict_norm_auto.rds',sep = ""))
  mean_p_strict_auto = rbind(mean_p_strict_auto, mean_p_strict_now)
  pergam_p_now = readRDS(paste(nospace.list[i],'_pergam_p_norm_auto.rds',sep = ""))
  pergam_p_auto = rbind(pergam_p_auto, pergam_p_now)
  coexp_all_now = readRDS(paste(nospace.list[i],'_coexp_out_norm_auto.rds',sep = ""))
  coexp_all = merge(coexp_all, coexp_all_now, by = 'row.names', all = T)
  rownames(coexp_all) = coexp_all$Row.names
  coexp_all = coexp_all[,-1]
  out_all_now = readRDS(paste(nospace.list[i],'_out_all_signed_norm_auto.rds',sep = ""))
  out_all_signed_auto = rbind(out_all_signed_auto, out_all_now)
  out_gam_now = readRDS(paste(nospace.list[i],'_out_gam_signed_norm_auto.rds',sep = ""))
  out_gam_signed_auto = rbind(out_gam_signed_auto, out_gam_now)
  resamp_now = readRDS(paste(nospace.list[i],'_out_gam_resamp_norm_XY_auto.rds',sep=""))
  resamp_XY_auto = rbind(resamp_XY_auto, resamp_now)
  
  out_all_now_X = readRDS(paste(nospace.list[i],'_out_all_norm_X.rds',sep = ""))
  out_all_X = rbind(out_all_X, out_all_now_X)
  out_gam_now_X = readRDS(paste(nospace.list[i],'_out_gam_norm_X.rds',sep = ""))
  out_gam_X = rbind(out_gam_X, out_gam_now_X)
  mean_p_now_X = readRDS(paste(nospace.list[i],'_mean_p_norm_X.rds',sep = ""))
  mean_p_X = rbind(mean_p_X, mean_p_now_X)
  pergam_p_now_X = readRDS(paste(nospace.list[i],'_pergam_p_norm_X.rds',sep = ""))
  pergam_p_X = rbind(pergam_p_X, pergam_p_now_X)
  coexp_all_now_X = readRDS(paste(nospace.list[i],'_coexp_out_norm_X.rds',sep = ""))
  coexp_all_X = merge(coexp_all_X, coexp_all_now_X, by = 'row.names', all = T)
  rownames(coexp_all_X) = coexp_all_X$Row.names
  coexp_all_X = coexp_all_X[,-1]
  resamp_now = readRDS(paste(nospace.list[i],'_out_gam_resamp_norm_X.rds',sep=""))
  resamp_X = rbind(resamp_X, resamp_now)
  
  out_all_now_X = readRDS(paste(nospace.list[i],'_out_all_norm_X_auto.rds',sep = ""))
  out_all_X_auto = rbind(out_all_X_auto, out_all_now_X)
  out_gam_now_X = readRDS(paste(nospace.list[i],'_out_gam_norm_X_auto.rds',sep = ""))
  out_gam_X_auto = rbind(out_gam_X_auto, out_gam_now_X)
  mean_p_now_X = readRDS(paste(nospace.list[i],'_mean_p_norm_X_auto.rds',sep = ""))
  mean_p_X_auto = rbind(mean_p_X_auto, mean_p_now_X)
  pergam_p_now_X = readRDS(paste(nospace.list[i],'_pergam_p_norm_X_auto.rds',sep = ""))
  pergam_p_X_auto = rbind(pergam_p_X_auto, pergam_p_now_X)
  coexp_all_now_X = readRDS(paste(nospace.list[i],'_coexp_out_norm_X.rds',sep = ""))
  coexp_all_X = merge(coexp_all_X, coexp_all_now_X, by = 'row.names', all = T)
  rownames(coexp_all_X) = coexp_all_X$Row.names
  coexp_all_X = coexp_all_X[,-1]
  resamp_now = readRDS(paste(nospace.list[i],'_out_gam_resamp_norm_X_auto.rds',sep=""))
  resamp_X_auto = rbind(resamp_X_auto, resamp_now)
  
}

mean_p$padj = p.adjust(mean_p$V1, method = 'BH')
mean_p_strict$padj = p.adjust(mean_p_strict$V1, method = 'BH')
mean_p_X$padj = p.adjust(mean_p_X$V1, method = 'BH')
pergam_p$padj = p.adjust(pergam_p$V3, method = 'BH')
pergam_p_X$padj = p.adjust(pergam_p_X$V3, method = 'BH')

mean_p_auto$padj = p.adjust(mean_p_auto$V1, method = 'BH')
mean_p_strict_auto$padj = p.adjust(mean_p_strict_auto$V1, method = 'BH')
mean_p_X_auto$padj = p.adjust(mean_p_X_auto$V1, method = 'BH')
pergam_p_auto$padj = p.adjust(pergam_p_auto$V3, method = 'BH')
pergam_p_X_auto$padj = p.adjust(pergam_p_X_auto$V3, method = 'BH')

saveRDS(out_all, file = 'out_all_norm.rds')
saveRDS(out_gam, file = 'out_gam_norm.rds')
saveRDS(mean_p, file = 'mean_p_norm.rds')
saveRDS(mean_p_strict, file = 'mean_p_strict_norm.rds')
saveRDS(pergam_p, file = 'pergam_p_norm.rds')
saveRDS(coexp_all, file = 'coexp_all.rds')
saveRDS(out_all_signed, file = 'out_all_signed_norm.rds')
saveRDS(out_gam_signed, file = 'out_gam_signed_norm.rds')
saveRDS(resamp_XY, file = 'resamp_XYsex.rds')

saveRDS(out_all_X, file = 'out_all_norm_X.rds')
saveRDS(out_gam_X, file = 'out_gam_norm_X.rds')
saveRDS(mean_p_X, file = 'mean_p_norm_X.rds')
saveRDS(pergam_p_X, file = 'pergam_p_norm_X.rds')
saveRDS(coexp_all_X, file = 'coexp_all_X.rds')
saveRDS(resamp_X, file = 'resamp_X.rds')

saveRDS(out_all_auto, file = 'out_all_norm_auto.rds')
saveRDS(out_gam_auto, file = 'out_gam_norm_auto.rds')
saveRDS(mean_p_auto, file = 'mean_p_norm_auto.rds')
saveRDS(mean_p_strict_auto, file = 'mean_p_strict_norm_auto.rds')
saveRDS(pergam_p_auto, file = 'pergam_p_norm_auto.rds')
saveRDS(coexp_all_auto, file = 'coexp_all_auto.rds')
saveRDS(out_all_signed_auto, file = 'out_all_signed_norm_auto.rds')
saveRDS(out_gam_signed_auto, file = 'out_gam_signed_norm_auto.rds')
saveRDS(resamp_XY_auto, file = 'resamp_XYsex_auto.rds')

saveRDS(out_all_X_auto, file = 'out_all_norm_X_auto.rds')
saveRDS(out_gam_X_auto, file = 'out_gam_norm_X_auto.rds')
saveRDS(mean_p_X_auto, file = 'mean_p_norm_X_auto.rds')
saveRDS(pergam_p_X_auto, file = 'pergam_p_norm_X_auto.rds')
saveRDS(coexp_all_X_auto, file = 'coexp_all_X_auto.rds')
saveRDS(resamp_X_auto, file = 'resamp_X_auto.rds')

#### end ####

#### load ####

out_all = readRDS('out_all_norm.rds')
out_gam = readRDS('out_gam_norm.rds')
#mean_p = readRDS('mean_p_norm.rds')
mean_p = readRDS('mean_p_strict_norm.rds')
pergam_p = readRDS('pergam_p_norm.rds')
coexp_all = readRDS('coexp_all.rds')
out_all_signed = readRDS('out_all_signed_norm.rds')
out_gam_signed = readRDS('out_gam_signed_norm.rds')
resamp_XY = readRDS('resamp_XYsex.rds')

out_all_X = readRDS('out_all_norm_X.rds')
out_gam_X = readRDS('out_gam_norm_X.rds')
mean_p_X = readRDS('mean_p_norm_X.rds')
pergam_p_X = readRDS('pergam_p_norm_X.rds')
coexp_all_X = readRDS('coexp_all_X.rds')
resamp_X = readRDS('resamp_X.rds')

# autosomal only

out_all = readRDS('out_all_norm_auto.rds')
out_gam = readRDS('out_gam_norm_auto.rds')
#mean_p = readRDS('mean_p_norm_auto.rds')
mean_p = readRDS('mean_p_strict_norm_auto.rds')
pergam_p = readRDS('pergam_p_norm_auto.rds')
coexp_all = readRDS('coexp_all_auto.rds')
out_all_signed = readRDS('out_all_signed_norm_auto.rds')
out_gam_signed = readRDS('out_gam_signed_norm_auto.rds')
resamp_XY = readRDS('resamp_XYsex_auto.rds')

out_all_X = readRDS('out_all_norm_X_auto.rds')
out_gam_X = readRDS('out_gam_norm_X_auto.rds')
mean_p_X = readRDS('mean_p_norm_X_auto.rds')
pergam_p_X = readRDS('pergam_p_norm_X_auto.rds')
coexp_all_X = readRDS('coexp_all_X_auto.rds')
resamp_X = readRDS('resamp_X_auto.rds')

#### end ####

#### Figures S3, 2C Table S2 ####

pl = data.frame()
all_norm = data.frame()
gams = unique(out_gam$gene)
for (i in 1:length(nospace.list)){
  print(nospace.list[i])
  gam = subset(out_gam, tissue == nospace.list[i])
  other = subset(out_all, tissue == nospace.list[i])
  all = rbind(gam, other)
  all$z = (all$diff - mean(all$diff)) / sd(all$diff)
  pl = rbind(pl, subset(all, gene %in% gams))
  all_norm = rbind(all_norm, all)
}

mod = aov(z ~ tissue + gene, data = pl)
mod = aov(diff ~ tissue + gene, data = pl)
summary(mod)
TukeyHSD(mod)

pl$tissue = as.factor(pl$tissue)
levels(pl$tissue)
colnames(pergam_p) = c('gene','tissue','p','padj')
pl = merge(pl, pergam_p, by = c('gene','tissue'))
pl$code = ifelse(pl$padj < 0.05, "*", NA)
pl$code2 = ifelse(pl$p < 0.05, "sig", NA)
levels(pl$tissue) = short.listnow
table(pl$code)
table(pl$code2)
table(is.na(pl$code))
table(is.na(pl$code2))

# Figure S3

cl = dcast(gene ~ tissue, value.var = 'z', data = pl)
rownames(cl) = cl$gene
cl = cl[,-1]
cl2 = pvclust(cl)
plot(cl2, cex = 0.5)
ord = cl2$hclust$order
names(ord) = levels(pl$tissue)
ord
levels(pl$tissue)
pl$tissue = reorder(pl$tissue, new.order = ord)
levels(pl$tissue)

clusMember = data.frame(cutree(cl2$hclust, h=0.6))
colnames(clusMember) = 'cluster'
clusMember$tissue = rownames(clusMember)
View(clusMember)

pl = merge(pl, clusMember, by = 'tissue')
pl$cluster = as.factor(pl$cluster) 
levels(pl$cluster)
levels(pl$cluster) = c('E','D','B','A','F','C','G')
clus.ord = c(4,3,6,2,1,5,7)
names(clus.ord) = levels(pl$cluster)
pl$cluster = reorder(pl$cluster, new.order = clus.ord)
levels(pl$cluster)

pl$gene = as.factor(pl$gene)
levels(pl$gene)
levels(pl$gene) = c('DDX3X/Y','EIF1AX/Y','KDM5C/D','UTX/Y','NLGN4X/Y','PCDH11X/Y','PRKX/Y',
                    'RPS4X/Y1','TBL1X/Y','TMSB4X/Y','TXLNG/Y','USP9X/Y','ZFX/Y')

pl2 = complete(pl, tissue, gene)
table(pl2$code2)

av = pl2 %>% group_by(tissue) %>% summarise(z = mean(z, na.rm = T))
av$gene = 'mean'
av$diff = NA
av$p = NA
av$padj = NA
av$code = NA
av$code2 = NA
av$cluster = NA
av = av[,c(1,3,4,2,5,6,7,8,9)]
pl2 = rbind(pl2, av)

pl2$code3 = ifelse(pl2$gene == 'mean','mean','')

pl2 = pl2[order(pl2$tissue),]
pl2 = pl2[complete.cases(pl2$z),]
for(i in 1:length(pl2$cluster)){
  if(is.na(pl2$cluster[i])) {pl2$cluster[i] = pl2$cluster[i-1]} else {pl2$cluster[i]}}

# Table S2

write.csv(pl, file = 'sex-dep-ced-per-tissue-and-gam.csv')

# Figure 2C 

ggplot(pl2, aes(y = tissue, x = gene, fill = z)) +
  geom_tile() +
  geom_tile(data = pl2[which(pl2$code2 == 'sig'),], aes(y = tissue, x = gene), fill = "transparent", color = "black") +
  scale_fill_gradient2(high=tissue.colors[7],
                       mid="white",
                       low=tissue.colors[14],
                       na.value="gray95", 
                       midpoint=0, name = bquote(atop("normalized",aCEFD[MXY_FXX]))) +
  facet_grid(cluster ~ code3, scale = 'free', space = 'free', drop = TRUE) +
  theme_article() +
  geom_text(label = pl2$code, size = 6, vjust = 0.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text.y.right = element_text(angle = 0),
        strip.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = 'top',
        panel.background = element_rect(fill = 'gray95'))

#### end ####

#### calc p values ####

combo = rbind(cbind(out_all, group = 'Non-gametologues'), cbind(out_gam, group = 'Gametologues'))

## compare gams (male X + Y) vs non-gams 
## not normalized

pvnongam = data.frame()
for(i in 1:length(levels(as.factor(combo$tissue)))){
  now = subset(combo, tissue == levels(as.factor(combo$tissue))[i])
  tt = t.test(abs(now$diff) ~ now$group)
  pvnongam[i,1] = levels(as.factor(combo$tissue))[i]
  pvnongam[i,2] = tt$p.value
  pvnongam[i,3] = tt$estimate[1]
  pvnongam[i,4] = tt$estimate[2]
}
colnames(pvnongam) = c('tissue','p','mean gam','mean non-gam')
pvnongam$padj = p.adjust(pvnongam$p, method = 'BH')

## compare gams (male X + Y) vs non-gams 
## normalized 

pvnongam_z = data.frame()
for(i in 1:length(levels(as.factor(combo$tissue)))){
  now = subset(combo, tissue == levels(as.factor(combo$tissue))[i])
  now$z = (now$diff - mean(now$diff)) / sd(now$diff)
  pvnongam_z[i,1] = levels(as.factor(combo$tissue))[i]
  pvnongam_z[i,2] = mean(subset(now, group == 'Gametologues')$z)
  pvnongam_z[i,3] = mean(subset(now, group == 'Non-gametologues')$z)
  tt = t.test(abs(now$z) ~ now$group)
  pvnongam_z[i,4] = tt$estimate[1]
  pvnongam_z[i,5] = tt$estimate[2]
  pvnongam_z[i,6] = tt$p.value
}
colnames(pvnongam_z) = c('tissue','mean gam','mean non','t mean gam','t mean non','p')
pvnongam_z$padj = p.adjust(pvnongam_z$p, method = 'BH')

## normalize CED across genes within X+Y and X-only
## compare normalized CED between X+Y vs X-only (per tissue)

comboX = rbind(cbind(out_all, group = 'XY'), cbind(out_gam, group = 'XY'), cbind(out_all_X, group = 'XX'), cbind(out_gam_X, group = 'XX'))
gams = unique(comboX[grepl("&", comboX$gene),'gene'])

pXvY = data.frame()
gam_test = data.frame()
for(i in 1:length(levels(as.factor(comboX$tissue)))){
  now = subset(comboX, tissue == levels(as.factor(comboX$tissue))[i])
  XY = subset(now, group == 'XY')
  XY$z = (XY$diff - mean(XY$diff)) / sd(XY$diff)
  XY = subset(XY, gene %in% gams)
  XX = subset(now, group == 'XX')
  XX$z = (XX$diff - mean(XX$diff)) / sd(XX$diff)
  XX = subset(XX, gene %in% gams)
  co = rbind(XY, XX)
  keep = names(which(table(co$gene) == 2))
  co = subset(co, gene %in% keep)
  gam_test = rbind(gam_test, co)
  tt = t.test(co$z ~ co$group, paired = T)
  pXvY[i,1] = levels(as.factor(comboX$tissue))[i]
  pXvY[i,2] = tt$p.value
  pXvY[i,3] = tt$estimate[1]
}
colnames(pXvY) = c('tissue','p','mean difference')
pXvY$padj = p.adjust(pXvY$p, method = 'BH')

## compare normalized CED between X+Y vs X-only (per gene)

pXvYG = data.frame()
gams = unique(combo[grepl("&", combo$gene),'gene'])

for(i in 1:length(levels(as.factor(gams)))){
  now = subset(gam_test, gene == levels(as.factor(gams))[i])
  keep = names(which(table(now$tissue) == 2))
  now = subset(now, tissue %in% keep)
  tt = t.test(now$z ~ now$group, paired = T)
  pXvYG[i,1] = levels(as.factor(gams))[i]
  pXvYG[i,2] = tt$p.value
  pXvYG[i,3] = tt$estimate
}
colnames(pXvYG) = c('gene','p','mean diff')
pXvYG$padj = p.adjust(pXvYG$p, method = 'BH')

## compare XY to 0 (per tissue)

pXYv0 = data.frame()
gam_test = data.frame()
for(i in 1:length(levels(as.factor(combo$tissue)))){
  now = subset(combo, tissue == levels(as.factor(combo$tissue))[i])
  meannow = mean(now$diff)
  sdnow = sd(now$diff)
  now$z = (now$diff - meannow) / sdnow
  gam_now = now$gene[which(grepl("&", now$gene))]
  now = subset(now, gene %in% gam_now)
  gam_test = rbind(gam_test, now)
  tt = t.test(now$z)
  pXYv0[i,1] = levels(as.factor(combo$tissue))[i]
  pXYv0[i,2] = tt$p.value
  pXYv0[i,3] = tt$estimate
}
colnames(pXYv0) = c('tissue','p','mean gam')
pXYv0$padj = p.adjust(pXYv0$p, method = 'BH')

## compare XY to 0 (per gene)

pXYv0G = data.frame()
gams = combo$gene[which(grepl("&", combo$gene))]
for(i in 1:length(levels(as.factor(gams)))){
  now = subset(gam_test, gene == levels(as.factor(gams))[i])
  tt = t.test(now$z)
  pXYv0G[i,1] = levels(as.factor(gams))[i]
  pXYv0G[i,2] = tt$p.value
  pXYv0G[i,3] = tt$estimate
}
colnames(pXYv0G) = c('gene','p','mean diff')
pXYv0G$padj = p.adjust(pXYv0G$p, method = 'BH')

#### end ####

#### Tables S2 S3 | Figures 2A 2B 2D ####

## p-value vectors per tissue
## XY vs non-gametologues (untransformed) = pvnongam
## XY vs non-gametologues (z) = pvnongam_z
## XY vs X only = pXvY
## XY vs 0 = pXYv0
## XY vs null (per tissue) = mean_p

colnames(mean_p)[2] = colnames(mean_p_X)[2] = 'tissue'
all_tissue_p = merge(pXvY, pvnongam_z, by = 'tissue') 
all_tissue_p = merge(all_tissue_p, pXYv0, by = 'tissue')
all_tissue_p = merge(all_tissue_p, mean_p, by = 'tissue')
all_tissue_p = merge(all_tissue_p, mean_p_X, by = 'tissue')
all_tissue_p = all_tissue_p[,c(1,4,10,13,15,17)]
colnames(all_tissue_p) = c('tissue','XvXY','XYvNONGAM','XYv0','XYvNULL','XonlyvNULL')
View(all_tissue_p)

out_all_now = out_all
out_all_X_now = out_all_X

# run

pl = data.frame()

for (i in 1:length(nospace.list)){
  
  print(nospace.list[i])
  allnow = subset(out_all_now, tissue == nospace.list[i])
  gamnow = subset(out_gam, tissue == nospace.list[i])
  combo = rbind(allnow, gamnow)  
  meannow = mean(combo$diff)
  sdnow = sd(combo$diff)
  combo$z = (combo$diff - meannow) / sdnow
  gamout = subset(combo, gene %in% gamnow$gene)
  gamout$group = 'male X+Y'
  pl = rbind(pl, gamout)
  
  allnow = subset(out_all_X_now, tissue == nospace.list[i])
  gamnow = subset(out_gam_X, tissue == nospace.list[i])
  combo = rbind(allnow, gamnow)  
  meannow = mean(combo$diff)
  sdnow = sd(combo$diff)
  combo$z = (combo$diff - meannow) / sdnow
  gamout = subset(combo, gene %in% gamnow$gene)
  gamout$group = 'male X only'
  pl = rbind(pl, gamout)
  
}

# Table S2
write.csv(pl, file = 'sex-dep-ced-per-tissue-and-gam2.csv')
write.csv(pl, file = 'sex-dep-ced-per-tissue-and-gam2_auto.csv')

pl2 = pl %>% group_by(tissue, group) %>% summarise(mean = mean(z), sd = sd(z), n = n())
pl2 = data.frame(pl2)
pl2$tissue = as.factor(pl2$tissue)
pl2$group = as.factor(pl2$group)
pl2$min = pl2$mean - pl2$sd/pl2$n
pl2$max = pl2$mean + pl2$sd/pl2$n
pl2 = merge(pl2, all_tissue_p, by = 'tissue')

# Table S3

write.csv(pl2, file = 'aced_mXY_mX.csv')
write.csv(pl2, file = 'aced_mXY_mX_auto.csv')

# Figure 2A

levels(pl2$tissue) = short.listnow
ggplot(pl2, aes(y=reorder(tissue, mean), x=mean, color = group)) + 
  geom_point()+
  geom_errorbar(aes(xmin=min, xmax=max), width=.2,
                position=position_dodge(0.05))  +
  theme_article() +
  scale_color_manual(values = c(tissue.colors[c(30,1)])) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') +
  xlab("Normalized Sex-Dependent \n Co-expression Divergence") +
  theme(axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.title = element_blank())

# Figure 2B

mxy1 = readRDS('Colon-USP9XY-male-XY.rds')
mx1 = readRDS('Colon-USP9XY-male-X.rds')
fxy1 = readRDS('Colon-USP9XY-female-X.rds')
xy1 = cbind(data.frame(mcoexp = mxy1, fcoexp = fxy1), tissue = 'Colon', group = 'XY')
rxy = round(cor.test(xy1$mcoexp, xy1$fcoexp)$estimate, 2)
fxy1 = fxy1[which(names(fxy1) %in% names(mx1))]
mx1 = mx1[which(names(mx1) %in% names(fxy1))]
x1 = cbind(data.frame(mcoexp = mx1, fcoexp = fxy1), tissue = 'Colon', group = 'X')
rx = round(cor.test(x1$mcoexp, x1$fcoexp)$estimate, 2)

pxy1 = ggplot(xy1, aes(x = mcoexp, y = fcoexp)) +
  geom_hex() +
  scale_fill_gradient2(low = "grey", high = "black") +
  annotate(geom="text", x=-0.6, y=1, label=paste("r = ",rxy,sep = "")) +
  theme_article() +
  xlim(c(-1,1)) + ylim(-1,1) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'black') +
  coord_equal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        legend.position = 'none')

px1 = ggplot(x1, aes(x = mcoexp, y = fcoexp)) +
  geom_hex() +
  scale_fill_gradient2(low = "grey", high = "black") +
  annotate(geom="text", x=-0.6, y=1, label=paste("r = ",rx,sep = "")) +
  theme_article() +
  xlim(c(-1,1)) + ylim(-1,1) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'black') +
  coord_equal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        legend.position = 'none')

mxy2 = readRDS('Adipose-USP9XY-male-XY.rds')
mx2 = readRDS('Adipose-USP9XY-male-X.rds')
fxy2 = readRDS('Adipose-USP9XY-female-X.rds')
xy2 = cbind(data.frame(mcoexp = mxy2, fcoexp = fxy2), tissue = 'Adipose', group = 'XY')
rxy = round(cor.test(xy2$mcoexp, xy2$fcoexp)$estimate, 2)
fxy2 = fxy2[which(names(fxy2) %in% names(mx2))]
mx2 = mx2[which(names(mx2) %in% names(fxy2))]
x2 = cbind(data.frame(mcoexp = mx2, fcoexp = fxy2), tissue = 'Adipose', group = 'X')
rx = round(cor.test(x2$mcoexp, x2$fcoexp)$estimate, 2)

pxy2 = ggplot(xy2, aes(x = mcoexp, y = fcoexp)) +
  geom_hex() +
  scale_fill_gradient2(low = "grey", high = "black") +
  annotate(geom="text", x=-0.6, y=1, label=paste("r = ",rxy,sep = "")) +
  theme_article() +
  xlim(c(-1,1)) + ylim(-1,1) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'black') +
  coord_equal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        legend.position = 'none')

px2 = ggplot(x2, aes(x = mcoexp, y = fcoexp)) +
  geom_hex() +
  scale_fill_gradient2(low = "grey", high = "black") +
  annotate(geom="text", x=-0.6, y=1, label=paste("r = ",rx,sep = "")) +
  theme_article() +
  xlim(c(-1,1)) + ylim(-1,1) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'black') +
  coord_equal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        legend.position = 'none')

mxy3 = readRDS('Ileum-USP9XY-male-XY.rds')
mx3 = readRDS('Ileum-USP9XY-male-X.rds')
fxy3 = readRDS('Ileum-USP9XY-female-X.rds')
xy3 = cbind(data.frame(mcoexp = mxy3, fcoexp = fxy3), tissue = 'Ileum', group = 'XY')
rxy = round(cor.test(xy3$mcoexp, xy3$fcoexp)$estimate, 2)
fxy3 = fxy3[which(names(fxy3) %in% names(mx3))]
mx3 = mx3[which(names(mx3) %in% names(fxy3))]
x3 = cbind(data.frame(mcoexp = mx3, fcoexp = fxy3), tissue = 'Ileum', group = 'X')
rx = round(cor.test(x3$mcoexp, x3$fcoexp)$estimate, 2)

pxy3 = ggplot(xy3, aes(x = mcoexp, y = fcoexp)) +
  geom_hex() +
  scale_fill_gradient2(low = "grey", high = "black") +
  annotate(geom="text", x=-0.6, y=1, label=paste("r = ",rxy,sep = "")) +
  theme_article() +
  xlim(c(-1,1)) + ylim(-1,1) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'black') +
  coord_equal() +
  theme(axis.text = element_text(size = 12),
        legend.position = 'none',
        axis.title = element_blank())

px3 = ggplot(x3, aes(x = mcoexp, y = fcoexp)) +
  geom_hex() +
  scale_fill_gradient2(low = "grey", high = "black") +
  annotate(geom="text", x=-0.6, y=1, label=paste("r = ",rx,sep = "")) +
  theme_article() +
  xlim(c(-1,1)) + ylim(-1,1) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'black') +
  coord_equal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_blank(),
        legend.position = 'none')

ggarrange(px2, pxy2, px1, pxy1, px3, pxy3, ncol = 2, nrow = 3)

pl = rbind(xy1, xy2, xy3, x1, x2, x3)
pl = pl[complete.cases(pl),]

ggplot(pl, aes(x = mcoexp, y = fcoexp)) +
  geom_hex() +
  scale_fill_gradient2(low = "grey", high = "black") +
  #annotate(geom="text", x=-0.6, y=1, label=paste("r = ",r,sep = "")) +
  facet_grid(tissue ~ group) +
  theme_article() +
  xlim(c(-1,1)) + ylim(-1,1) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'black') +
  coord_equal() +
  theme(axis.text = element_text(size = 12),
        legend.position = 'none',
        axis.title = element_text(size = 12))

# Figure 2D

arg2 = 'SmallIntestineTerminalIleum'
gam = "EIF1AX&EIF1AY"

coexp = readRDS(file = paste(arg2,'_coexp_out_norm.rds',sep=""))
coexp = coexp[,grep(gam, colnames(coexp))]
colnames(coexp) = c('Male','Female')
t = cor.test(coexp$Male, coexp$Female, method = 'pearson')
r = round(t$estimate,2)

p1 = ggplot(coexp, aes(x = Male, y = Female)) +
  geom_hex() +
  scale_fill_gradient2(low = "grey", high = "black") +
  xlab('Males (X+Y)') +
  ylab('Females (X)') +
  annotate(geom="text", x=-0.6, y=1, label=paste("r = ",r,sep = "")) +
  theme_article() +
  coord_equal() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'black') +
  xlim(c(-1,1)) + ylim(-1,1) +
  theme(axis.text = element_text(size = 12),
        legend.position = 'none',
        axis.title = element_blank())

arg2 = 'SmallIntestineTerminalIleum'
gam = "USP9X&USP9Y"

coexp = readRDS(file = paste(arg2,'_coexp_out_norm.rds',sep=""))
coexp = coexp[,grep(gam, colnames(coexp))]
colnames(coexp) = c('Male','Female')
t = cor.test(coexp$Male, coexp$Female, method = 'pearson')
r = round(t$estimate,2)

p2 = ggplot(coexp, aes(x = Male, y = Female)) +
  geom_hex() +
  scale_fill_gradient2(low = "grey", high = "black") +
  xlab('Males (X+Y)') +
  ylab('Females (X)') +
  annotate(geom="text", x=-0.6, y=1, label=paste("r = ",r,sep = "")) +
  theme_article() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'black') +
  xlim(c(-1,1)) + ylim(-1,1) +
  coord_equal() +
  theme(axis.text = element_text(size = 12),
        legend.position = 'none',
        axis.title = element_blank())

arg = "Brain - Nucleus accumbens (basal ganglia)"
arg2 = "BrainNucleusaccumbensbasalganglia"
gam = "EIF1AX&EIF1AY"

coexp = readRDS(file = paste(arg2,'_coexp_out_norm.rds',sep=""))
coexp = coexp[,grep(gam, colnames(coexp))]
colnames(coexp) = c('Male','Female')
t = cor.test(coexp$Male, coexp$Female, method = 'pearson')
r = round(t$estimate,2)

p3 = ggplot(coexp, aes(x = Male, y = Female)) +
  geom_hex() +
  scale_fill_gradient2(low = "grey", high = "black") +
  xlab('Males (X+Y)') +
  ylab('Females (X)') +
  annotate(geom="text", x=-0.6, y=1, label=paste("r = ",r,sep = "")) +
  theme_article() +
  xlim(c(-1,1)) + ylim(-1,1) +
  coord_equal() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'black') +
  theme(axis.text = element_text(size = 12),
        legend.position = 'none',
        axis.title = element_blank())

arg = "Brain - Nucleus accumbens (basal ganglia)"
arg2 = "BrainNucleusaccumbensbasalganglia"
gam = "USP9X&USP9Y"

coexp = readRDS(file = paste(arg2,'_coexp_out_norm.rds',sep=""))
coexp = coexp[,grep(gam, colnames(coexp))]
colnames(coexp) = c('Male','Female')
t = cor.test(coexp$Male, coexp$Female, method = 'pearson')
r = round(t$estimate,2)

p4 = ggplot(coexp, aes(x = Male, y = Female)) +
  geom_hex() +
  scale_fill_gradient2(low = "grey", high = "black") +
  xlab('Males (X+Y)') +
  ylab('Females (X)') +
  annotate(geom="text", x=-0.6, y=1, label=paste("r = ",r,sep = "")) +
  theme_article() +
  xlim(c(-1,1)) + ylim(-1,1) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'black') +
  coord_equal() +
  theme(axis.text = element_text(size = 12),
        legend.position = 'none',
        axis.title = element_blank())

arg = "Heart - Left Ventricle"
arg2 = "HeartLeftVentricle"
gam = "EIF1AX&EIF1AY"

coexp = readRDS(file = paste(arg2,'_coexp_out_norm.rds',sep=""))
coexp = coexp[,grep(gam, colnames(coexp))]
colnames(coexp) = c('Male','Female')
t = cor.test(coexp$Male, coexp$Female, method = 'pearson')
r = round(t$estimate,2)

p5 = ggplot(coexp, aes(x = Male, y = Female)) +
  geom_hex() +
  scale_fill_gradient2(low = "grey", high = "black") +
  xlab('Males (X+Y)') +
  ylab('Females (X)') +
  annotate(geom="text", x=-0.6, y=1, label=paste("r = ",r,sep = "")) +
  theme_article() +
  xlim(c(-1,1)) + ylim(-1,1) +
  coord_equal() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'black') +
  theme(axis.text = element_text(size = 12),
        legend.position = 'none',
        axis.title = element_blank())

arg = "Heart - Left Ventricle"
arg2 = "HeartLeftVentricle"
gam = "USP9X&USP9Y"

coexp = readRDS(file = paste(arg2,'_coexp_out_norm.rds',sep=""))
coexp = coexp[,grep(gam, colnames(coexp))]
colnames(coexp) = c('Male','Female')
t = cor.test(coexp$Male, coexp$Female, method = 'pearson')
r = round(t$estimate,2)

p6 = ggplot(coexp, aes(x = Male, y = Female)) +
  geom_hex() +
  scale_fill_gradient2(low = "grey", high = "black") +
  xlab('Males (X+Y)') +
  ylab('Females (X)') +
  annotate(geom="text", x=-0.6, y=1, label=paste("r = ",r,sep = "")) +
  theme_article() +
  xlim(c(-1,1)) + ylim(-1,1) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'black') +
  coord_equal() +
  theme(axis.text = element_text(size = 12),
        legend.position = 'none',
        axis.title = element_blank())

ggarrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 3)

#### end ####

#### compare all vs. autosomal only ####

## absolute & signed sex dep CFD (per tissue & pair)
## Figure S4

out_gam = readRDS('out_gam_norm.rds')
out_gam_auto = readRDS('out_gam_norm_auto.rds')

out_gam = readRDS('out_gam_signed_norm.rds')
out_gam_auto = readRDS('out_gam_signed_norm_auto.rds')

out_gam = readRDS('out_gam_norm_X.rds')
out_gam_auto = readRDS('out_gam_norm_X_auto.rds')

pl = cbind(out_gam, out_gam_auto$diff)
colnames(pl) = c('all_genes','tissue','pair','auto_only')
pl$tissue = as.factor(pl$tissue)
levels(pl$tissue) = short.listnow
hist(pl$all_genes - pl$auto_only)

ggplot(pl, aes(x = all_genes, y = auto_only, color = tissue)) +
  geom_point() +
  scale_color_manual(values = tissue.colors) +
  theme_article() +
  xlab('All genes') + ylab('Autosomal genes only') +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = 'lm', color = 'blue') +
  theme(legend.position = 'bottom',
        legend.title = element_blank())

ggplot(pl, aes(x = all_genes - auto_only)) +
  geom_histogram(fill = tissue.colors[30], color = tissue.colors[40]) +
  theme_article() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlab('all genes - autosomal only') + ylab('') +
  theme(legend.position = 'bottom',
        legend.title = element_blank())


out_gam = readRDS('out_gam_norm.rds')
out_gam_X = readRDS('out_gam_norm_X.rds')
comp = merge(out_gam, out_gam_X, by = c('tissue','gene'))
comp$diff = comp$diff.x - comp$diff.y
hist(comp$diff)
table(comp$diff > 0)

out_gam_auto = readRDS('out_gam_norm_auto.rds')
out_gam_auto_X = readRDS('out_gam_norm_X_auto.rds')
comp = merge(out_gam_auto, out_gam_auto_X, by = c('tissue','gene'))
comp$diff = comp$diff.x - comp$diff.y
hist(comp$diff)
table(comp$diff > 0)

#### end ####
