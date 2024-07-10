library(expss)
library(pheatmap)
library(clValid)
library(dplyr)
library(maditr)
library(dcanr)
library(corrplot)
library(ggplot2)
library(egg)
library(igraph)
library(ggnetwork)
library(tidytext)
library(tidyr)

`%!in%` = Negate(`%in%`)

#### prepare data ####

avg_diff_out_now = readRDS(file = 'out_xy.rds')

avg_diff_out_now$key = paste(avg_diff_out_now$region, avg_diff_out_now$gene.y, sep = "_")
check = data.frame(avg_diff_out_now %>% group_by(key) %>% summarise(avg = mean(diff, na.rm = TRUE)))
rm = check[is.na(check$avg),]$key
avg_diff_out_now = subset(avg_diff_out_now, key %!in% rm)

#avg_diff_out_now$diff = abs(avg_diff_out_now$x_coexp) - abs(avg_diff_out_now$y_coexp)

avg_diff_out_cast = dcast(avg_diff_out_now, formula = gene.x ~ key, value.var = 'diff')
avg_diff_out_cast = data.frame(avg_diff_out_cast)
rownames(avg_diff_out_cast) = avg_diff_out_cast$gene.x
avg_diff_out_cast = avg_diff_out_cast[,-1]
avg_diff_out_cast = t(avg_diff_out_cast)
avg_diff_out_cast = as.matrix(avg_diff_out_cast)
dim(avg_diff_out_cast)

#### end ####

#### variance partitioning | Figure 4F | Table S12 ####

library(variancePartition)
library(tidyverse)

v = t(avg_diff_out_cast)
v = v[complete.cases(v),]
dim(v) # 8116

m = data.frame(key = colnames(v))
m = separate(data = m, col = key, into = c('tissue','pair'), sep = "_", remove = FALSE)

form = ~ (1|tissue) + (1|pair)
varPart = fitExtractVarPartModel(v, form, m)
vp = sortCols(varPart)
saveRDS(vp, file = 'vp.rds')

vp = readRDS('vp.rds')

# Table S12

write.csv(vp, file = 'vp.csv')
plotVarPart(vp)

mean(vp$pair)
mean(vp$tissue)
mean(vp$Residuals)

obj = data.frame(vp)
col=c(tissue.colors[c(1,30)], "grey85")
obj$gene = rownames(obj)
data.plot = melt(obj, id="gene")
levels(data.plot$variable) = c('Pair','Tissue','Residuals')

ggplot(data=data.plot, aes(x=variable, y=value)) + 
  geom_violin( scale="width", aes(alpha = 0.8, fill = factor(variable))) +
  theme_article() +
  ylab('Proportion of variance') +
  scale_fill_manual(values=col) +
  theme(legend.position = 'none',
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 16))

#### end ####

#### collapse similar tissues | Figure S11D-E | Table S12 ####

shorts = unique(str_sub(nospace.list,1,4))
p = unique(gsub(".*\\_","",colnames(v)))
c2 = v
for(i in 1:length(shorts)){
  print(shorts[i])
  for(j in 1:length(p)){
    n = c2[,which(str_sub(colnames(c2),1,4) == shorts[i])]
    s = n[,which(gsub(".*\\_","",colnames(n)) == p[j])]
    s = data.frame(s)
    if(dim(s)[2] == 1) {colnames(s) = colnames(n)[which(gsub(".*\\_","",colnames(n)) == p[j])]} 
    if(is.null(colnames(s)) == TRUE) {print("skip")} else {
      t = rowMeans(s, na.rm = T)
      t = data.frame(t)
      c2 = c2[,which(colnames(c2) %!in% colnames(s))]
      c2 = cbind(t, c2)
      colnames(c2)[1] = paste(shorts[i], p[j],sep="_")
    }
  }
}
colnames(c2)
rem = which(c2 %>% summarise_all(funs(sum(is.nan(.)))) == dim(c2)[1])
c3 = c2[,-c(rem)]
colnames(c3)
dim(c3) # 8116

m = data.frame(key = colnames(c3))
m = separate(data = m, col = key, into = c('tissue','pair'), sep = "_", remove = FALSE)
form = ~ (1|tissue) + (1|pair)
varPart = fitExtractVarPartModel(c3, form, m)
vp = sortCols(varPart)
saveRDS(vp, file = 'vp_collapsed_tissues.rds')

vp = readRDS('vp_collapsed_tissues.rds')

# Table S12

write.csv(vp, file = 'vp_collapsed_tissues.csv')
plotVarPart(vp)

mean(vp$pair)
mean(vp$tissue)
mean(vp$Residuals)

obj = data.frame(vp)
col=c(tissue.colors[c(1,30)], "grey85")
obj$gene = rownames(obj)
data.plot = melt(obj, id="gene")
levels(data.plot$variable) = c('Pair','Tissue','Residuals')

ggplot(data=data.plot, aes(x=variable, y=value)) + 
  geom_violin( scale="width", aes(alpha = 0.8, fill = factor(variable))) +
  theme_article() +
  ylab('Proportion of variance') +
  scale_fill_manual(values=col) +
  theme(legend.position = 'none',
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 16))

#### end ####

#### Figure 4G | Table S13 ####

v = t(avg_diff_out_cast)
dim(v) # 20290

m = data.frame(key = colnames(v))
m = separate(data = m, col = key, into = c('tissue','pair'), sep = "_", remove = FALSE)

varout = data.frame()
tukout_pair = data.frame()
tukout_tissue = data.frame()

for (i in 1:length(rownames(v))){
  print(i)
  now = data.frame(diff = t(v)[,i],
                   tissue = m[,c('tissue')],
                   pair = m[,c('pair')])
  now = now[complete.cases(now),]
  if(length(unique(now$tissue)) <= 1) {
    tnow = data.frame(gene = colnames(v)[i], f_pair = NA, p_pair = NA, f_tissue = NA, p_tissue = NA)
    tukp = data.frame(diff = NA, lwr = NA, upr = NA, p.adj = NA, gene = colnames(v)[i])
    tukt = data.frame(diff = NA, lwr = NA, upr = NA, p.adj = NA, gene = colnames(v)[i])
  } else {
  t = aov(diff ~ tissue + pair, data = now)
  tuk = TukeyHSD(t)
  t = summary(t)
  tnow = data.frame(gene = rownames(v)[i],
                    f_pair = t[[1]][2,4],
                    p_pair = t[[1]][2,5],
                    f_tissue = t[[1]][1,4],
                    p_tissue = t[[1]][1,5])
  tukp = data.frame(tuk$pair, gene = rownames(v)[i])
  tukp = subset(tukp, p.adj < 0.05)
  tukt = data.frame(tuk$tissue, gene = rownames(v)[i])
  tukt = subset(tukt, p.adj < 0.05)
  
  }
  
  varout = rbind(varout, tnow)
  tukout_pair = rbind(tukout_pair, tukp)
  tukout_tissue = rbind(tukout_tissue, tukt)
  
}

saveRDS(varout, file = 'varout.rds')
saveRDS(varout, file = 'varout_collapsed.rds')

varout = readRDS('varout.rds')
varout = readRDS('varout_collapsed.rds')

table(is.na(varout$p_pair)) # 2210 genes only expressed in 1 tissue
varout = varout[complete.cases(varout),]
dim(varout) # 18080 genes remaining (expressed in at least 2 tissues)
varout$p_pair_adj = p.adjust(varout$p_pair, method = 'bonferroni')
varout$p_tissue_adj = p.adjust(varout$p_tissue, method = 'bonferroni')

# Table S13

write.csv(varout, file = 'varout.csv')
write.csv(varout, file = 'varout_collapsed.csv')

table(varout$p_pair_adj < 0.05) 
12639/18080 * 100 
table(varout$p_tissue_adj < 0.05) 
1181/18080 * 100 
table(varout$p_tissue_adj < 0.05 & varout$p_pair_adj < 0.05) 
1141/1181 * 100 
1141/18080 * 100 
12639 - 1141
(12639 - 1141)/18080 * 100 
1181 - 1141
(1181 - 1141)/18080 * 100 
table(varout$p_tissue_adj < 0.05 & varout$p_pair_adj > 0.05) 
table(varout$p_tissue_adj > 0.05 & varout$p_pair_adj > 0.05) 

# plot

groups = c('NA','pair only','tissue only','pair & tissue','neither')
counts = c(2210, 11498, 40, 1141, 5401)
sum(counts)
counts / sum(counts) * 100
pl = data.frame(groups = groups, counts = counts)
pl$groups = factor(pl$groups, levels = c("pair only", 
                                         "tissue only", 
                                         "pair & tissue", 
                                         "neither", 
                                         "NA"))
levels(pl$groups) 

ggplot(pl, aes(x="", y = counts, fill = groups)) +
  geom_bar(stat="identity", width=1, color="white", alpha = 0.8) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = tissue.colors[c(1,30,15,28,10)]) +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 18))

#### end ####

#### dimensionality reduction | Figures 4E S11A-C ####

library(Rtsne)
library(umap)
library(MASS)

out_scdced = readRDS('out_scdced.rds')
out_scdced$region = as.factor(out_scdced$region)
levels(out_scdced$region) = short.listnow
colnames(out_scdced)[6] = 'tissue'
out_scdced = out_scdced[,-1]
colnames(out_scdced)[6] = 'pair'

c = t(avg_diff_out_cast)
c = c[complete.cases(c),]
dim(c)

# pca

diff.pca = prcomp(t(c))
summary(diff.pca)
pl = data.frame(diff.pca$x[,c(1:10)])
pl$tissue = t(data.frame(str_split(rownames(pl),"_")))[,1]
pl$pair = t(data.frame(str_split(rownames(pl),"_")))[,2]
pl$tissue = as.factor(pl$tissue)
levels(pl$tissue) = short.listnow
pl$pair = gsub("\\..."," & ",pl$pair)
pl = merge(pl, out_scdced[,c(1,5,6)], by = c('tissue','pair'))
pl$scfd = ifelse(pl$diff > 0, 'X', 'Y')
head(pl)

ggplot(pl, aes(x=PC1, 
               y=PC2, 
               color = tissue, 
               shape = pair)) +
  geom_point(size=3) +
  scale_color_manual(values=tissue.colors[c(1,3,4,7,20,21,23,26,28:37,39:43)]) +
  scale_shape_manual(values=c(8:25)) +
  theme_classic() +
  theme(legend.text = element_text(size=8), 
        legend.box = 'vertical', 
        legend.title = element_blank(), 
        legend.key.size = unit(0, 'lines'), 
        legend.position = 'bottom')

#tsne

a = Rtsne(t(c))

plot.tsne = as.data.frame(a$Y)
plot.tsne$tissue = t(data.frame(str_split(colnames(c),"_")))[,1]
plot.tsne$pair = t(data.frame(str_split(colnames(c),"_")))[,2]
plot.tsne$tissue = as.factor(plot.tsne$tissue)
levels(plot.tsne$tissue) = short.listnow
plot.tsne$pair = gsub("\\..."," & ",plot.tsne$pair)
plot.tsne = merge(plot.tsne, out_scdced[,c(1,5,6)], by = c('tissue','pair'))
plot.tsne$scfd = ifelse(pl$diff > 0, 'X', 'Y')

ggplot(plot.tsne, aes(x=V1, 
                      y=V2, 
                      color = tissue,
                      shape = pair)) +
  geom_point(size=3) +
  scale_color_manual(values=tissue.colors) +
  scale_shape_manual(values=c(8:25)) +
  theme_classic() +
  theme(legend.text = element_text(size=8), 
        legend.box = 'vertical', 
        legend.title = element_blank(), 
        legend.key.size = unit(0, 'lines'), 
        legend.position = 'bottom')

# umap

set.seed(88)

a = umap(t(c), metric = 'manhattan')
plot.umap = as.data.frame(a$layout)
plot.umap$tissue = t(data.frame(str_split(colnames(c),"_")))[,1]
plot.umap$pair = t(data.frame(str_split(colnames(c),"_")))[,2]
plot.umap$tissue = as.factor(plot.umap$tissue)
levels(plot.umap$tissue) = short.listnow
plot.umap$pair = gsub("\\..."," & ",plot.umap$pair)

ggplot(plot.umap, aes(x=V1, 
                      y=V2, 
                      color = tissue,
                      shape = pair)) +
  geom_point(size=3) +
  scale_color_manual(values=tissue.colors) +
  scale_shape_manual(values=c(8:25)) +
  theme_article() +
  xlab('UMAP 1') + ylab('UMAP 2') +
  theme(legend.text = element_text(size=8), 
        legend.box = 'vertical', 
        legend.title = element_blank(), 
        legend.key.size = unit(0, 'lines'), 
        legend.position = 'bottom')

ggplot(plot.umap, aes(x=V1, 
                      y=V2, 
                      color = tissue,
                      shape = pair)) +
  geom_point(size=3) +
  scale_color_manual(values=tissue.colors) +
  scale_shape_manual(values=c(8:25)) +
  theme_article() +
  xlab('UMAP 1') + ylab('UMAP 2') +
  theme(legend.text = element_text(size=8), 
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.box = 'vertical', 
        legend.title = element_blank(), 
        legend.key.size = unit(0, 'lines'), 
        legend.position = 'none')

#### end ####

#### collapse similar tissues | Figure S11D-E ####

# collapse similar tissues
# average correlations within similar tissues
shorts = unique(str_sub(nospace.list,1,4))
p = unique(gsub(".*\\_","",colnames(c)))
c2 = c
for(i in 1:length(shorts)){
  print(shorts[i])
  for(j in 1:length(p)){
    n = c2[,which(str_sub(colnames(c2),1,4) == shorts[i])]
    s = n[,which(gsub(".*\\_","",colnames(n)) == p[j])]
    s = data.frame(s)
    if(dim(s)[2] == 1) {colnames(s) = colnames(n)[which(gsub(".*\\_","",colnames(n)) == p[j])]} 
    if(is.null(colnames(s)) == TRUE) {print("skip")} else {
    t = rowMeans(s, na.rm = T)
    t = data.frame(t)
    c2 = c2[,which(colnames(c2) %!in% colnames(s))]
    c2 = cbind(t, c2)
    colnames(c2)[1] = paste(shorts[i], p[j],sep="_")
    }
    }
}
colnames(c2)
rem = which(c2 %>% summarise_all(funs(sum(is.nan(.)))) == dim(c2)[1])
c3 = c2[,-c(rem)]
colnames(c3)

# pca

diff.pca = prcomp(t(c3))

summary(diff.pca)
pl = data.frame(diff.pca$x[,c(1:10)])
pl$tissue = t(data.frame(str_split(rownames(pl),"_")))[,1]
pl$pair = t(data.frame(str_split(rownames(pl),"_")))[,2]
pl$tissue = as.factor(pl$tissue)
levels(pl$tissue) = c('Adipose','Adrenal','Atery','Brain','Breast','Colon','Esophagus',
                      'Heart','Kidney','Liver','Lung','Salivary','Muscle','Nerve','Pancreas',
                      'Pituitary','Prostate','Skin','Ileum','Spleen','Stomach','Testes','Thyroid')
pl$pair = gsub("\\..."," & ",pl$pair)
head(pl)

ggplot(pl, aes(x=PC1, 
               y=PC2, 
               color = tissue, 
               shape = pair)) +
  geom_point(size=3) +
  scale_color_manual(values=tissue.colors[c(1,3,4,7,20,21,23,26,28:37,39:43)]) +
  scale_shape_manual(values=c(8:25)) +
  theme_classic() +
  theme(legend.text = element_text(size=8), 
        legend.box = 'vertical', 
        legend.title = element_blank(), 
        legend.key.size = unit(0, 'lines'), 
        legend.position = 'bottom')

a = Rtsne(t(c3))

plot.tsne = as.data.frame(a$Y)
plot.tsne$tissue = t(data.frame(str_split(colnames(c3),"_")))[,1]
plot.tsne$pair = t(data.frame(str_split(colnames(c3),"_")))[,2]
plot.tsne$tissue = as.factor(plot.tsne$tissue)
levels(plot.tsne$tissue) = c('Adipose','Adrenal','Atery','Brain','Breast','Colon','Esophagus',
                      'Heart','Kidney','Liver','Lung','Salivary','Muscle','Nerve','Pancreas',
                      'Pituitary','Prostate','Skin','Ileum','Spleen','Stomach','Testes','Thyroid')
plot.tsne$pair = gsub("\\..."," & ",plot.tsne$pair)

ggplot(plot.tsne, aes(x=V1, 
                      y=V2, 
                      color = tissue,
                      shape = pair)) +
  geom_point(size=3) +
  scale_color_manual(values=tissue.colors[c(1,3,4,7,20,21,23,26,28:37,39:43)]) +
  scale_shape_manual(values=c(8:25)) +
  theme_classic() +
  theme(legend.text = element_text(size=8), 
        legend.box = 'vertical', 
        legend.title = element_blank(), 
        legend.key.size = unit(0, 'lines'), 
        legend.position = 'bottom')

set.seed(88)

a = umap(t(c3), metric = 'manhattan')

plot.umap = as.data.frame(a$layout)
plot.umap$tissue = t(data.frame(str_split(colnames(c3),"_")))[,1]
plot.umap$pair = t(data.frame(str_split(colnames(c3),"_")))[,2]
plot.umap$tissue = as.factor(plot.umap$tissue)
levels(plot.umap$tissue) = c('Adipose','Adrenal','Atery','Brain','Breast','Colon','Esophagus',
                             'Heart','Kidney','Liver','Lung','Salivary','Muscle','Nerve','Pancreas',
                             'Pituitary','Prostate','Skin','Ileum','Spleen','Stomach','Testes','Thyroid')
plot.umap$pair = gsub("\\..."," & ",plot.umap$pair)

ggplot(plot.umap, aes(x=V1, 
                      y=V2, 
                      color = tissue,
                      shape = pair)) +
  geom_point(size=3) +
  scale_color_manual(values=tissue.colors[c(1,3,4,7,20,21,23,26,28:37,39:43)]) +
  scale_shape_manual(values=c(8:25)) +
  theme_article() +
  xlab('UMAP 1') + ylab('UMAP 2') +
  theme(legend.text = element_text(size=8), 
        legend.box = 'vertical', 
        legend.title = element_blank(), 
        legend.key.size = unit(0, 'lines'), 
        legend.position = 'bottom')

#### end ####

#### correlations across pair/tissue ANOVA groups ####

pair_genes = subset(varout, p_pair_adj < 0.05 & p_tissue_adj > 0.05)$gene
tissue_genes = subset(varout, p_tissue_adj < 0.05 & p_pair_adj > 0.05)$gene
pair_tissue_genes = subset(varout, p_tissue_adj < 0.05 & p_pair_adj < 0.05)$gene

M = cor.pairs(avg_diff_out_cast[,which(colnames(avg_diff_out_cast) %in% pair_genes)], cor.method = 'spearman')
dd =  as.dist(1 - M)
cl = hclust(dd, method = 'average')
cor(cophenetic(cl), dd) 
plot(cl)
saveRDS(M, file = 'M_pair_genes.rds')
saveRDS(dd, file = 'dd_pair_genes.rds')
saveRDS(cl, file = 'cl_pair_genes.rds')

M = cor.pairs(avg_diff_out_cast[,which(colnames(avg_diff_out_cast) %in% tissue_genes)], cor.method = 'spearman')
dd =  as.dist(1 - M)
cl = hclust(dd, method = 'average')
cor(cophenetic(cl), dd) 
plot(cl, cex = 0.5)
saveRDS(M, file = 'M_tissue_genes.rds')
saveRDS(dd, file = 'dd_tissue_genes.rds')
saveRDS(cl, file = 'cl_tissue_genes.rds')

M = cor.pairs(avg_diff_out_cast[,which(colnames(avg_diff_out_cast) %in% pair_tissue_genes)], cor.method = 'spearman')
dd =  as.dist(1 - M)
cl = hclust(dd, method = 'average')
cor(cophenetic(cl), dd) 
plot(cl, cex = 0.3)
saveRDS(M, file = 'M_pair_tissue_genes.rds')
saveRDS(dd, file = 'dd_pair_tissue_genes.rds')
saveRDS(cl, file = 'cl_pair_tissue_genes.rds')

#### end ####

#### find optimal N clusters ####

dd = readRDS('dd_pair_genes.rds')
set.seed(88)
clusterNo=NbClust(diss = dd, distance=NULL, min.nc=2, max.nc=15, method="average", index = "dunn")
clusterNo$Best.nc 
saveRDS(clusterNo, file = 'clusterNo_dunn_pair_genes.rds')

dd = readRDS('dd_tissue_genes.rds')
set.seed(88)
clusterNo=NbClust(diss = dd, distance=NULL, min.nc=2, max.nc=15, method="average", index = "dunn")
clusterNo$Best.nc 
saveRDS(clusterNo, file = 'clusterNo_dunn_tissue_genes.rds')

dd = readRDS('dd_pair_tissue_genes.rds')
set.seed(88)
clusterNo=NbClust(diss = dd, distance=NULL, min.nc=2, max.nc=15, method="average", index = "dunn")
clusterNo$Best.nc 
saveRDS(clusterNo, file = 'clusterNo_dunn_pair_tissue_genes.rds')

#### end ####

#### prep plots optimal N clusters | Table S13 ####

## pair predicted genes

clusterNo = readRDS('clusterNo_dunn_pair_genes.rds')
table(clusterNo$Best.partition)
clus = levels(as.factor(clusterNo$Best.partition))

write.csv(data.frame(clusterNo$Best.partition), file = 'pair_clusters.csv')

clus_out = data.frame()
for(i in 1:length(clus)){
  genes_now = names(clusterNo$Best.partition[which(clusterNo$Best.partition == clus[i])])
  diff_now = subset(avg_diff_out_now, gene.x %in% genes_now)
  means_now = diff_now %>% group_by(key) %>% summarise(mean = mean(diff, na.rm = T))
  means_now = data.frame(means_now)
  means_now$clus = clus[i]
  clus_out = rbind(clus_out, means_now)
}
saveRDS(clus_out, file = 'clus_out_pair_per_pair_tissue.rds')

clus_out = data.frame()
for(i in 1:length(clus)){
  genes_now = names(clusterNo$Best.partition[which(clusterNo$Best.partition == clus[i])])
  diff_now = subset(avg_diff_out_now, gene.x %in% genes_now)
  means_now = diff_now %>% group_by(gene.y) %>% summarise(mean = mean(diff, na.rm = T))
  means_now = data.frame(means_now)
  means_now$clus = clus[i]
  clus_out = rbind(clus_out, means_now)
}

sim_out = data.frame()
for(i in 1:length(clus)){
  print(i)
  genes_now = names(clusterNo$Best.partition[which(clusterNo$Best.partition == clus[i])])
  Ngenes = length(genes_now)
  for(j in 1:100){
    s = sample(names(clusterNo$Best.partition), Ngenes)
    d = subset(avg_diff_out_now, gene.x %in% s)
    mm = d %>% group_by(gene.y) %>% summarise(mean = mean(diff, na.rm = T))
    mm = data.frame(mm)
    mm$clus = clus[i]
    mm$it = j
    sim_out = rbind(sim_out, mm)
  }
}

for(i in 1:length(clus_out$gene.y)){
  n = clus_out[i,]
  s = subset(sim_out, gene.y == n$gene.y[1] & clus == n$clus[1])
  p = sum(abs(s$mean) > abs(n$mean[1])) / 100
  clus_out$p[i] = p
}

clus_out$padj = p.adjust(clus_out$p, method = 'BH')
table(clus_out$padj < 0.05)
clus_out$sig = ifelse(clus_out$padj < 0.05, "*", "")
saveRDS(clus_out, file = 'clus_out_pair.rds')

## tissue predicted genes

clusterNo = readRDS('clusterNo_dunn_tissue_genes.rds')
table(clusterNo$Best.partition)
clus = levels(as.factor(clusterNo$Best.partition))

write.csv(data.frame(clusterNo$Best.partition), file = 'tissue_clusters.csv')

clus_out = data.frame()
for(i in 1:length(clus)){
  genes_now = names(clusterNo$Best.partition[which(clusterNo$Best.partition == clus[i])])
  diff_now = subset(avg_diff_out_now, gene.x %in% genes_now)
  means_now = diff_now %>% group_by(key) %>% summarise(mean = mean(diff, na.rm = T))
  means_now = data.frame(means_now)
  means_now$clus = clus[i]
  clus_out = rbind(clus_out, means_now)
}
saveRDS(clus_out, file = 'clus_out_tissue_per_pair_tissue.rds')

clus_out = data.frame()
for(i in 1:length(clus)){
  genes_now = names(clusterNo$Best.partition[which(clusterNo$Best.partition == clus[i])])
  diff_now = subset(avg_diff_out_now, gene.x %in% genes_now)
  means_now = diff_now %>% group_by(region) %>% summarise(mean = mean(diff, na.rm = T))
  means_now = data.frame(means_now)
  means_now$clus = clus[i]
  clus_out = rbind(clus_out, means_now)
}

sim_out = data.frame()
for(i in 1:length(clus)){
  genes_now = names(clusterNo$Best.partition[which(clusterNo$Best.partition == clus[i])])
  Ngenes = length(genes_now)
  for(j in 1:100){
    s = sample(names(clusterNo$Best.partition), Ngenes)
    d = subset(avg_diff_out_now, gene.x %in% s)
    mm = d %>% group_by(region) %>% summarise(mean = mean(diff, na.rm = T))
    mm = data.frame(mm)
    mm$clus = clus[i]
    mm$it = j
    sim_out = rbind(sim_out, mm)
  }
}

for(i in 1:length(clus_out$region)){
  n = clus_out[i,]
  s = subset(sim_out, region == n$region[1] & clus == n$clus[1])
  p = sum(abs(s$mean) > abs(n$mean[1])) / 100
  clus_out$p[i] = p
}

clus_out$padj = p.adjust(clus_out$p, method = 'BH')
table(clus_out$padj < 0.05)
clus_out$sig = ifelse(clus_out$padj < 0.05, "*", "")
saveRDS(clus_out, file = 'clus_out_tissue.rds')

## pair + tissue predicted genes

clusterNo = readRDS('clusterNo_dunn_pair_tissue_genes.rds')
table(clusterNo$Best.partition)
clus = levels(as.factor(clusterNo$Best.partition))

write.csv(data.frame(clusterNo$Best.partition), file = 'pair_tissue_clusters.csv')

clus_out = data.frame()
for(i in 1:length(clus)){
  genes_now = names(clusterNo$Best.partition[which(clusterNo$Best.partition == clus[i])])
  diff_now = subset(avg_diff_out_now, gene.x %in% genes_now)
  means_now = diff_now %>% group_by(key) %>% summarise(mean = mean(diff, na.rm = T))
  means_now = data.frame(means_now)
  means_now$clus = clus[i]
  clus_out = rbind(clus_out, means_now)
}
saveRDS(clus_out, file = 'clus_out_pair_tissue_per_pair_tissue.rds')

clus_out = data.frame()
for(i in 1:length(clus)){
  genes_now = names(clusterNo$Best.partition[which(clusterNo$Best.partition == clus[i])])
  diff_now = subset(avg_diff_out_now, gene.x %in% genes_now)
  #means_now = diff_now %>% group_by(gene.y) %>% summarise(mean = mean(diff, na.rm = T))
  means_now = diff_now %>% group_by(region) %>% summarise(mean = mean(diff, na.rm = T))
  means_now = data.frame(means_now)
  means_now$clus = clus[i]
  clus_out = rbind(clus_out, means_now)
}

sim_out = data.frame()
for(i in 1:length(clus)){
  genes_now = names(clusterNo$Best.partition[which(clusterNo$Best.partition == clus[i])])
  Ngenes = length(genes_now)
  for(j in 1:100){
    s = sample(names(clusterNo$Best.partition), Ngenes)
    d = subset(avg_diff_out_now, gene.x %in% s)
    #mm = d %>% group_by(gene.y) %>% summarise(mean = mean(diff, na.rm = T))
    mm = d %>% group_by(region) %>% summarise(mean = mean(diff, na.rm = T))
    mm = data.frame(mm)
    mm$clus = clus[i]
    mm$it = j
    sim_out = rbind(sim_out, mm)
  }
}

#for(i in 1:length(clus_out$gene.y)){
for(i in 1:length(clus_out$region)){
  n = clus_out[i,]
  #s = subset(sim_out, gene.y == n$gene.y[1] & clus == n$clus[1])
  s = subset(sim_out, region == n$region[1] & clus == n$clus[1])
  p = sum(abs(s$mean) > abs(n$mean[1])) / 100
  clus_out$p[i] = p
}

clus_out$padj = p.adjust(clus_out$p, method = 'BH')
table(clus_out$padj < 0.05)
clus_out$sig = ifelse(clus_out$padj < 0.05, "*", "")

saveRDS(clus_out, file = 'clus_out_pair_tissue_per_pair.rds')
saveRDS(clus_out, file = 'clus_out_pair_tissue_per_tissue.rds')

#### end ####

#### Figures 4H S12 ####

clusterNo = readRDS('clusterNo_dunn_pair_genes.rds')
clus_out = readRDS('clus_out_pair.rds')

clusterNo = readRDS('clusterNo_dunn_tissue_genes.rds')
clus_out = readRDS('clus_out_tissue.rds')

clusterNo = readRDS('clusterNo_dunn_pair_tissue_genes.rds')
clus_out = readRDS('clus_out_pair_tissue_per_pair.rds')
clus_out = readRDS('clus_out_pair_tissue_per_tissue.rds')

# all
colnames(clus_out)[1] = 'group'
clus_out$group = as.factor(clus_out$group)
levels(clus_out$group)

# for tissue
levels(clus_out$group) = short.listnow

# for pair
levels(clus_out$group) = c('DDX3X/Y', 'EIF1AX/Y', 'KDM5C/D', 'UTX/Y', 'NLGN4X/Y',
                            'PCDH11X/Y','PRKX/Y', 'RPS4X/Y1', 'SOX3/SRY','TBL1X/Y',
                            'TGIF2LX/Y','TMSB4X/Y','TXLNG/Y','USP9X/Y','ZFX/Y')
clus_out = subset(clus_out, clus %in% c(1:8))                       

ca = dcast(clus_out, formula = group ~ clus, value.var = 'mean')
ca = data.frame(ca)
rownames(ca) = ca$group
ca = ca[,-1]
colnames(ca) = seq(1:length(colnames(ca)))
cat = t(ca)

dd2 = dist(ca)
cl2 = hclust(dd2, method = 'average')
plot(cl2)
ord = cl2$order
clus_out$group = factor(clus_out$group, levels = levels(clus_out$group)[ord])

dd3 =  dist(cat)
cl3 = hclust(dd3, method = 'average')
plot(cl3)
ord2 = cl3$order
clus_out$clus = as.factor(clus_out$clus)
clus_out$clus = factor(clus_out$clus, levels = levels(clus_out$clus)[ord2])

levels(clus_out$clus)
table(clusterNo$Best.partition)

for(i in 1:length(levels(clus_out$clus))){
  levels(clus_out$clus)[i] = paste(levels(clus_out$clus)[i], 
                                   " (N=", 
                                   table(clusterNo$Best.partition)[levels(clus_out$clus)[i]],
                                   ")", sep = "")}
levels(clus_out$clus)

# Figures S12 C/F

clus_out %>%
  mutate(clus = as.factor(clus),
         group = reorder_within(group, mean, clus)) %>%
  ggplot(aes(x = group, y = mean, fill = mean)) +
    geom_bar(stat = 'identity') +
    scale_fill_gradient2(low = tissue.colors[14], 
                       high = tissue.colors[7], 
                       mid = 'white', 
                       na.value = '#FDF0D5',
                       limits = c(-0.5,0.5),
                       name = bquote(mean~sCEFD[MX_MY])) +
    facet_wrap(~clus, scales = 'free_x', ncol = 1, strip.position = 'left') +
    scale_x_reordered() +
    theme_article() +
    theme(strip.placement = 'outside',
          axis.text = element_text(size = 6),
          #axis.text.x = element_text(angle = 0), # pair
          axis.text.x = element_text(angle = 90, hjust = 1), # tissue
          axis.title = element_blank(),
          legend.position = 'none')

# Figures 4H, S12 A/D/I

ggplot(clus_out, aes(y = clus, x = group, fill = mean)) +
  geom_tile() +
  scale_fill_gradient2(low = tissue.colors[14], 
                       high = tissue.colors[7], 
                       mid = 'white', 
                       na.value = '#FDF0D5',
                       limits = c(-0.5,0.5),
                       name = bquote(mean~sCEFD[MX_MY])) +
  theme_article() +
  #ylab('Cluster (tissue only)') + # tissue
  ylab('Cluster (pair only)') + #pair
  #ylab('Cluster (pair and tissue)') + #pair
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  facet_wrap(~clus, ncol = 1, scales = 'free_y') +
  geom_text(label = clus_out$sig, size = 8, vjust = 0.8) + #pair
  #geom_text(label = clus_out$sig, size = 6, vjust = 0.8) + # tissue
  theme(axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 10), # pair
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1), #pair
        #axis.text.x = element_text(size = 8, angle = 45, hjust = 1), # tissues
        #axis.text.y = element_text(size = 14), # tissue
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = 'top',
        plot.margin = margin(10, 10, 10, 20))

# Figures S12 B/E/J

#clus_out2 = readRDS('clus_out_pair_per_pair_tissue.rds')
#clus_out2 = readRDS('clus_out_tissue_per_pair_tissue.rds')
clus_out2 = readRDS('clus_out_pair_tissue_per_pair_tissue.rds')

table(clus_out2$clus)

colnames(clus_out2)[1] = 'group'

clus_out2 = separate(clus_out2, group, into = c('region','pair'), sep = "_")
clus_out2$region = factor(clus_out2$region)
levels(clus_out2$region) = short.listnow

#for pair
clus_out2 = subset(clus_out2, clus %in% c(1:8))
clus_out2$pair = as.factor(clus_out2$pair)
clus_out2$pair = factor(clus_out2$pair, levels = levels(clus_out2$pair)[ord])

#for tissue
clus_out2$region = as.factor(clus_out2$region)
clus_out2$region = factor(clus_out2$region, levels = levels(clus_out2$region)[ord])

clus_out2$clus = as.factor(clus_out2$clus)
clus_out2$clus = factor(clus_out2$clus, levels = levels(clus_out2$clus)[ord2])

ggplot(clus_out2, aes(y = region, x = pair, fill = mean)) +
  geom_tile() +
  scale_fill_gradient2(low = tissue.colors[14], 
                       high = tissue.colors[7], 
                       mid = 'white', 
                       na.value = '#FDF0D5',
                       #limits = c(-0.75,0.75),
                       limits = c(-1,1),
                       name = bquote(mean~sCEFD[MX_MY])) +
  theme_article() +
  #ylab('Cluster (tissue only)') +
  ylab('Cluster (pair only)') + #pair
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  facet_wrap(~clus, ncol = 1, scales = 'free_y') +
  #geom_text(label = clus_out$sig, size = 8, vjust = 0.8) + #pair
  #geom_text(label = clus_out$sig, size = 6, vjust = 0.8) +
  theme(axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 4), # pair
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1), #pair
        #axis.text.x = element_text(size = 8, angle = 45, hjust = 1), # tissues
        #axis.text.y = element_text(size = 14), # tissue
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = 'top',
        plot.margin = margin(10, 10, 10, 20))

ggplot(clus_out2, aes(y = pair, x = region, fill = mean)) +
  geom_tile() +
  scale_fill_gradient2(low = tissue.colors[14], 
                       high = tissue.colors[7], 
                       mid = 'white', 
                       na.value = '#FDF0D5',
                       limits = c(-1,1),
                       name = bquote(mean~sCEFD[MX_MY])) +
  theme_article() +
  ylab('Cluster (tissue only)') +
  #ylab('Cluster (pair only)') + #pair
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  facet_wrap(~clus, ncol = 1, scales = 'free_y') +
  #geom_text(label = clus_out$sig, size = 8, vjust = 0.8) + #pair
  #geom_text(label = clus_out$sig, size = 6, vjust = 0.8) +
  theme(axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 2),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = 'top',
        plot.margin = margin(10, 10, 10, 20))

#### end ####

#### Figures 4H, S12 G/H/K (enrichments) #### 

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(GOSemSim)

#clusterNo = readRDS('clusterNo_dunn_pair_genes.rds')
clusterNo = readRDS('clusterNo_dunn_tissue_genes.rds')
#clusterNo = readRDS('clusterNo_dunn_pair_tissue_genes.rds')

clus = unique(clusterNo$Best.partition)
clusgo = data.frame()
for(i in 1:length(clus)){
  print(clus[i])
  ego = enrichGO(gene         = names(clusterNo$Best.partition[which(clusterNo$Best.partition == clus[i])]),
                universe      = colnames(avg_diff_out_cast),
                keyType       = 'ENSEMBL',
                OrgDb         = org.Hs.eg.db,
                #ont           = "BP",
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.2,
                qvalueCutoff  = 0.2,
                minGSSize = 5,
                maxGSSize = 500, 
                readable = TRUE)
  ego2 = clusterProfiler::simplify(ego, measure = 'Wang', cutoff=0.5, by="pvalue", select_fun=min)
  go_res = ego@result
  if(dim(go_res)[1] == 0) {print("no enrichments")} else{
  go_res$cluster = clus[i]
  clusgo = rbind(clusgo, go_res)
}}
View(clusgo)

#saveRDS(clusgo, file = 'clusgo_pair_only_bp.rds')
#saveRDS(clusgo, file = 'clusgo_pair_only_cc.rds')
#saveRDS(clusgo, file = 'clusgo_tissue_only_bp.rds')
#saveRDS(clusgo, file = 'clusgo_tissue_only_cc.rds')
#saveRDS(clusgo, file = 'clusgo_pair_tissue_bp.rds')
saveRDS(clusgo, file = 'clusgo_pair_tissue_cc.rds')

## plot

# for pair
clusgo1 = readRDS('clusgo_pair_only_bp.rds')
clusgo2 = readRDS('clusgo_pair_only_cc.rds')

## Tables S14 S15
write.csv(clusgo1, file = 'bp_pair.csv')
write.csv(clusgo2, file = 'cc_pair.csv')

clusgo1 = subset(clusgo1, ID %!in% c("GO:0010563" , "GO:0006505", "GO:0044782"))

clusgo1 = subset(clusgo1, cluster %in% c(1:8))                        
clpl1 = clusgo1 %>% group_by(cluster) %>% slice(1:2) 
clpl1$cluster = as.factor(clpl1$cluster)
clpl1$cluster = factor(clpl1$cluster, levels = levels(clpl1$cluster)[ord2])
clusgo2 = subset(clusgo2, cluster %in% c(1:8))                        
clpl2 = clusgo2 %>% group_by(cluster) %>% slice(1:2) 
clpl2$cluster = as.factor(clpl2$cluster)
clpl2$cluster = factor(clpl2$cluster, levels = levels(clpl2$cluster)[ord2])

# for tissue
clusgo1 = readRDS('clusgo_tissue_only_bp.rds')
clusgo2 = readRDS('clusgo_tissue_only_cc.rds')

## Tables S16 S17
write.csv(clusgo1, file = 'bp_tissue.csv')
write.csv(clusgo2, file = 'cc_tissue.csv')

clpl1 = clusgo1 %>% group_by(cluster) %>% slice(1:4) 
clpl2 = clusgo2 %>% group_by(cluster) %>% slice(1:4) 

# for pair and tissue
clusgo1 = readRDS('clusgo_pair_tissue_bp.rds')
clusgo2 = readRDS('clusgo_pair_tissue_cc.rds')

## Tables S18 S19
write.csv(clusgo1, file = 'bp_pair_tissue.csv')
write.csv(clusgo2, file = 'cc_pair_tissue.csv')

clpl1 = clusgo1 %>% group_by(cluster) %>% slice(1:2)
clpl2 = clusgo2 %>% group_by(cluster) %>% slice(1:2) 

# for all
clpl1$pvalue = ifelse(-log10(clpl1$pvalue) > 20, 1e-20, clpl1$pvalue)

clpl1$Description
clpl1$Description = gsub("involved in","-",clpl1$Description)
clpl1$Description = gsub("negative regulation of ","",clpl1$Description)
clpl1$Description = gsub("positive regulation of ","",clpl1$Description)
clpl1$Description = gsub("signaling pathway","",clpl1$Description)
clpl1$Description = gsub("based on somatic recombination of immune receptors built from immunoglobulin superfamily domains","",clpl1$Description)
clpl1$Description = gsub("regulation of ","",clpl1$Description)
clpl1$Description

ggplot(clpl1, aes(y = reorder(Description, -log10(pvalue)), x = -log10(pvalue))) +
  geom_bar(stat = 'identity', fill = tissue.colors[31], alpha = 0.6) +
  facet_grid(cluster ~ ., scales = 'free') +
  theme_article() +
  xlab('-log10(p-value)') +
  scale_y_discrete(position = 'right') +
  scale_x_continuous(expand = c(0,0), limits = c(0,20)) +
  geom_vline(xintercept = -log10(0.05), linetype = 'dashed', color = 'black') +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 10),
        strip.text.y = element_blank())

clpl2$pvalue = ifelse(-log10(clpl1$pvalue) > 20, 1e-20, clpl1$pvalue)

clpl2$Description

ggplot(clpl2, aes(y = reorder(Description, -log10(pvalue)), x = -log10(pvalue))) +
  geom_bar(stat = 'identity', fill = tissue.colors[3], alpha = 0.6) +
  facet_grid(cluster ~ ., scales = 'free') +
  theme_article() +
  scale_x_continuous(expand = c(0,0), limits = c(0,20)) +
  xlab('-log10(p-value)') +
  scale_y_discrete(position = 'right') +
  geom_vline(xintercept = -log10(0.05), linetype = 'dashed', color = 'grey') +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 10),
        strip.text.y = element_blank())

#### end ####

#### find hub genes ####

hsap = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 104)
gam_bm = getBM(attributes = c("external_gene_name","ensembl_gene_id", "chromosome_name"),  mart = hsap)

clusgo = readRDS('clusgo_pair_only_bp.rds')
clusterNo = readRDS('clusterNo_dunn_pair_genes.rds')

genesnow = names(clusterNo$Best.partition[which(clusterNo$Best.partition == 4)])
genesnow = names(clusterNo$Best.partition[which(clusterNo$Best.partition == 8)])
genesnow = names(clusterNo$Best.partition[which(clusterNo$Best.partition == 7)])
genesnow = names(clusterNo$Best.partition[which(clusterNo$Best.partition == 2)])
genesnow = names(clusterNo$Best.partition[which(clusterNo$Best.partition == 3)])

m = avg_diff_out_cast[,genesnow]
c = cor.pairs(m, cor.method = 'spearman')
d = rowMeans(c, na.rm = T)
d = d[order(d, decreasing = TRUE)]
d = data.frame(d)

d$ensembl_gene_id = rownames(d)
d = merge(d, gam_bm, by = 'ensembl_gene_id')

clussub = clusgo[clusgo$Description %like% "translation", ]
clussub = clusgo[clusgo$Description %like% "ubiquitination", ]
clussub = clusgo[clusgo$Description %like% "splicing", ]
clussub = clusgo[clusgo$Description %like% "synap", ]
clussub = clusgo[clusgo$Description %like% "immun", ]

clussub = clussub %>% separate_rows(geneID, sep="/") 
colnames(clussub)[8] = "external_gene_name"
dnow = merge(d, clussub, by = 'external_gene_name')
dnow = dnow[order(dnow$d, decreasing = TRUE),]
View(dnow)

#### end ####
