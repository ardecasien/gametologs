library(dplyr)
library(rstatix)
library(rcompanion)
library(igraph)
library(ggraph)
library(stringr)

#### load data ####

## load metadata

meta = readRDS('gtex_combined_meta.rds')

# load asd data

asd = read.csv('satterstrom_cat.csv')
colnames(asd)[c(1,2)] = c('name','gene')
dim(asd)
table(asd$Functional_Category)

# load overall diff coupling data

xy = readRDS('out_xy_weight.rds')
colnames(xy)[c(1,2)] = c('gene','tissue')
clip = readRDS('xy_clip.rds')
clip$tissue = as.factor(clip$tissue)
levels(clip$tissue) = nospace.list

h = merge(xy, clip, by = c('gene','tissue'))
table(h$tissue, h$padj < 0.05)
# BrainCortex N=0
# BrainSpinalcordcervicalc1 N=0
# KidneyCortex N=0
# Testis N=0

a = merge(asd, h, by = 'gene')

# brain tissues only

tissuesnow = a$tissue[which(str_sub(unique(a$tissue),1,5) == 'Brain')]

listnow = tissue.list[which(str_sub(tissue.list,1,5) == 'Brain')]
listnow = nospace.list[c(7:19)]
sl = short.listnow[which(str_sub(tissue.list,1,5) == 'Brain')]

#### end ####

#### summarise ####

a = subset(a, tissue %in% tissuesnow)
dim(a)
write.csv(a, file = 'satterstrom-asym-clip.csv')
length(unique(a$gene)) # 100 genes
unique(asd$gene)[which(unique(asd$gene) %!in% a$gene)] # PAX5 TCF20 not in our data
xc = subset(a, padj < 0.05 & diff_xy_new > 0)
length(unique(xc$gene)) # 34 genes
table(unique(xc[,c('gene','Functional_Category2')])$Functional_Category2) # 26 GER
xcg = xc %>% group_by(gene, name, Functional_Category2) %>% summarise(n = n())
View(xcg) # CREBBP HDLBP VEZF1 CTNNB1 WAC
xcg %>% group_by(Functional_Category2) %>% summarise(mean = mean(n), max = max(n), min = min(n))
yc = subset(a, padj < 0.05 & diff_xy_new < 0)
length(unique(yc$gene)) # 36 genes
table(unique(yc[,c('gene','Functional_Category2')])$Functional_Category2) # 16 NC
ycg = yc %>% group_by(gene, name, Functional_Category2) %>% summarise(n = n())
View(ycg) # SCN2A STXBP1 SCN1A GABRB2 CACNA1E
ycg %>% group_by(Functional_Category2) %>% summarise(mean = mean(n), max = max(n), min = min(n))
unique(yc$gene)[which(unique(yc$gene) %in% xc$gene)] # N = 4 genes in both (3/4 = other)
nn = subset(a, name %!in% xc$name & name %!in% yc$name)
table(unique(nn[,c('gene','Functional_Category2')])$Functional_Category2) 

# combine with per gam asymmetric coupling

xygam = readRDS('out_xy.rds')
colnames(xygam)[9] = 'tissue'
colnames(xygam)[10] = 'gene'
pl = subset(xygam, gene %in% asd$gene)
pl = subset(pl, tissue %in% listnow)
pl = merge(pl, a, by = c('gene','tissue'))
colnames(pl)[11] = 'gams'
pl = merge(pl, gams, by = 'gams')

xyclip = readRDS(file = 'xy_clip_per_gam.rds')
xyclip$tissue = factor(xyclip$tissue)
levels(xyclip$tissue) = nospace.list
clip = subset(xyclip, tissue %in% listnow)
clip = subset(clip, gene %in% asd$gene)
colnames(clip)[7] = 'gams'
pl = merge(pl, clip, by = c('gene','gams','tissue'))
pl = subset(pl, padj.x < 0.05)
pl$overall = ifelse(pl$diff_xy_new > 0, "X", "Y")
pl$specific = ifelse(pl$diff > 0, "X", "Y")
pl$specific = ifelse(pl$padj.y > 0.05, "ns", pl$specific)
length(unique(pl$name))

check = pl %>% group_by(Functional_Category2, overall, specific, gams) %>% summarise(n = n())
check = subset(check, specific != 'ns')
View(check)

#### end ####

#### chi squared ####

aa = a
aa$sig = ifelse(aa$padj < 0.05, 'sig', 'ns')
aa %>% group_by(name) %>% summarise(minp = min(padj))
sig = subset(a, padj < 0.05)

m = matrix(c(26,6,2,12,8,16,18,6,6), nrow = 3, ncol = 3)
colnames(m) = c('X','Y','neither')
rownames(m) = c('GER','Other','NC')
sum(m)
m
ct = chisq.test(m)
ct
pairwiseNominalIndependence(m,fisher = TRUE,gtest  = FALSE,chisq  = TRUE,method = "fdr")
pairwiseNominalIndependence(t(m),fisher = TRUE,gtest  = FALSE,chisq  = TRUE,method = "fdr")

#### end ####

#### Figure 5D ####

library(wordcloud)
library(wordcloud2)

xpl = xcg
xpl = subset(xpl, name %!in% ycg$name)
table(xpl$Functional_Category2)
xpl$col = ifelse(xpl$Functional_Category2 == 'Gene expression regulation', '#1B9E77', '#666666')
xpl$col = ifelse(xpl$Functional_Category2 == 'Neuronal communication', '#D0990B', xpl$col)
table(xpl$col)
wordcloud(words = xpl$name, freq = xpl$n, min.freq = 1,           
          max.words=200, random.order=FALSE, rot.per=0.35,  scale=c(2,0.5),    
          random.color = F, ordered.colors=TRUE,          
          colors=xpl$col)

ypl = ycg
ypl = subset(ypl, name %!in% xcg$name)
table(ypl$Functional_Category2)
ypl$col = ifelse(ypl$Functional_Category2 == 'Gene expression regulation', '#1B9E77', '#666666')
ypl$col = ifelse(ypl$Functional_Category2 == 'Neuronal communication', '#D0990B', ypl$col)
table(ypl$col)
wordcloud(words = ypl$name, freq = ypl$n, min.freq = 1, scale=c(2,0.5),          
          max.words=200, random.order=FALSE, rot.per=0.35, 
          random.color = F, ordered.colors=TRUE,          
          colors=ypl$col)


# plot venn

library(ggvenn)
library(ggVennDiagram)

s = subset(a, padj < 0.05)
s$coup = ifelse(s$diff_xy_new > 0, "X", "Y")
pl = unique(s[,c('name','Functional_Category2','coup')])
table(pl$name)
pl = pl[order(pl$Functional_Category2),]
x = list(x = subset(pl, coup == 'X')$name, y = subset(pl, coup == 'Y')$name)

ggvenn(x, show_elements = T, label_sep = "\n", 
       fill_color = brewer.pal(name="Dark2",n=32), text_size = 2) +
      theme(text = element_text(face = 'italics'))

#### end ####

#### null | Table S24 ####

b = merge(asd, h, by = 'gene', all = TRUE)
b$coup = ifelse(b$padj < 0.05 & b$diff_xy_new > 0, "X", "non")
b$coup = ifelse(b$padj < 0.05 & b$diff_xy_new < 0, "Y", b$coup)
table(b$coup)

c = list()
g = list()
ts = data.frame()
null = data.frame()

for(i in 1:length(tissuesnow)){
  
  print(tissuesnow[i])
  ts[i,1] = tissuesnow[i]
  t = subset(b, tissue == tissuesnow[i])
  t$Functional_Category2[is.na(t$Functional_Category2)] = 'nonASD'
  
  table(t$Functional_Category2, t$coup)
  m = matrix(table(t$Functional_Category2, t$coup), ncol = 3, nrow = 4)
  colnames(m) = c('non','X','Y')
  rownames(m) = c('GER','NC','nonASD','Other')
  m
  c[[i]] = m
  prop.table(m, margin = 1)
  prop.table(m, margin = 2)
  #fisher.test(m)
  ct = chisq.test(m)
  ts[i,2] = ct$p.value
  #pairwise_fisher_test(m)
  pp = pairwiseNominalIndependence(m,fisher = TRUE,gtest  = FALSE,chisq  = TRUE,method = "fdr")
  ts[i,3] = pp[1,3]
  ts[i,4] = pp[2,3]
  ts[i,5] = pp[3,3]
  ts[i,6] = pp[4,3]
  ts[i,7] = pp[5,3]
  ts[i,8] = pp[6,3]
  pairwiseNominalIndependence(t(m),fisher = TRUE,gtest  = FALSE,chisq  = TRUE,method = "fdr")

  ms = t %>% group_by(Functional_Category2) %>% summarise(mean = mean(diff_xy_new))
  ms
  g[[i]] = ms
  an = aov(diff_xy_new ~ Functional_Category2, data = t)
  sa = summary(an)
  ts[i,9] = sa[[1]][1,5]
  tu = TukeyHSD(an)
  ts[i,10] = tu$Functional_Category2[1,4]
  ts[i,11] = tu$Functional_Category2[2,4]
  ts[i,12] = tu$Functional_Category2[3,4]
  ts[i,13] = tu$Functional_Category2[4,4]
  ts[i,14] = tu$Functional_Category2[5,4]
  ts[i,15] = tu$Functional_Category2[6,4]
#}
  exp = readRDS(paste(listnow[i], '_adjusted_exp_MALES.rds',sep=""))
  me = rowMeans(exp, na.rm = T)
  me = me[order(me)]
  me = data.frame(me)
  me$decile = ntile(me$me, 10)  
  
  measd = me[which(rownames(me) %in% subset(t, Functional_Category2 != 'nonASD')$gene),]
  meger = me[which(rownames(me) %in% subset(t, Functional_Category2 == 'Gene expression regulation')$gene),]
  menc = me[which(rownames(me) %in% subset(t, Functional_Category2 == 'Neuronal communication')$gene),]
  meo = me[which(rownames(me) %in% subset(t, Functional_Category2 == 'Other')$gene),]
  menonasd = me[which(rownames(me) %in% subset(t, Functional_Category2 == 'nonASD')$gene),]
  
  ito = data.frame()
  for(k in 1:1000){
        ito[k,1] = k
        
        # ger
        itger = data.frame()
        for(j in 1:length(rownames(meger))){
          itger[j,1] = sample(rownames(subset(menonasd, decile == meger$decile[j])),1)}
        colnames(itger) = 'gene'
        itgern = merge(t, itger, by = 'gene')
        ito[k,2] = mean(itgern$diff_xy_new)
        ito[k,3] = table(itgern$coup)[2]
        ito[k,4] = table(itgern$coup)[3]
        ito[k,5] = table(itgern$coup)[1]
        
        # nc
        itnc = data.frame()
        for(j in 1:length(rownames(menc))){
          itnc[j,1] = sample(rownames(subset(menonasd, decile == menc$decile[j])),1)}
        colnames(itnc) = 'gene'
        itncn = merge(t, itnc, by = 'gene')
        ito[k,6] = mean(itncn$diff_xy_new)
        ito[k,7] = table(itncn$coup)[2]
        ito[k,8] = table(itncn$coup)[3]
        ito[k,9] = table(itncn$coup)[1]
        
        # other
        itot = data.frame()
        for(j in 1:length(rownames(meo))){
          itot[j,1] = sample(rownames(subset(menonasd, decile == meo$decile[j])),1)}
        colnames(itot) = 'gene'
        itotn = merge(t, itot, by = 'gene')
        ito[k,10] = mean(itotn$diff_xy_new)
        ito[k,11] = table(itotn$coup)[2]
        ito[k,12] = table(itotn$coup)[3]
        ito[k,13] = table(itotn$coup)[1]
  }
  
  colnames(ito) = c('iteration',
                        'mean.GER','X.GER','Y.GER','non.GER',
                        'mean.NC','X.NC','Y.NC','non.NC',
                        'mean.other','X.other','Y.other','non.other')
  
  # do GER and NC genes include a high # of X or Y coupled genes
  
  null[i,1] = tissuesnow[i]
  null[i,2] = sum(ito$X.GER[complete.cases(ito$X.GER)] > m[1,2])/1000
  null[i,3] = sum(ito$Y.GER[complete.cases(ito$Y.GER)] > m[1,3])/1000
  null[i,4] = sum(ito$X.NC[complete.cases(ito$X.NC)] > m[2,2])/1000
  null[i,5] = sum(ito$Y.NC[complete.cases(ito$Y.NC)] > m[2,3])/1000
  null[i,6] = sum(ito$X.other[complete.cases(ito$X.other)] > m[4,2])/1000
  null[i,7] = sum(ito$Y.other[complete.cases(ito$Y.other)] > m[4,3])/1000
  
  # do GER and NC genes show extreme asym coupling values?
  
  null[i,8] = sum(abs(ito$mean.GER) > abs(ms$mean[1]))/1000
  null[i,9] = sum(abs(ito$mean.NC) > abs(ms$mean[2]))/1000
  null[i,10] = sum(abs(ito$mean.other) > abs(ms$mean[3]))/1000
  
  null[i,11] = sum((ito$mean.GER) > (ms$mean[1]))/1000
  null[i,12] = sum((ito$mean.NC) > (ms$mean[2]))/1000
  null[i,13] = sum((ito$mean.other) > (ms$mean[3]))/1000
  
  null[i,14] = sum((ito$mean.GER) < (ms$mean[1]))/1000
  null[i,15] = sum((ito$mean.NC) < (ms$mean[2]))/1000
  null[i,16] = sum((ito$mean.other) < (ms$mean[3]))/1000
}

colnames(null) = c('tissue','X.GER','Y.GER','X.NC','Y.NC','X.other','Y.other',
                   'mean.GER','mean.NC','mean.other',
                   'mean.GER.greater','mean.NC.greater','mean.other.greater',
                   'mean.GER.lesser','mean.NC.lesser','mean.other.lesser')
null$X.GER.adj = p.adjust(null$X.GER)
null$Y.GER.adj = p.adjust(null$Y.GER)
null$X.NC.adj = p.adjust(null$X.NC)
null$Y.NC.adj = p.adjust(null$Y.NC)
null$X.other.adj = p.adjust(null$X.other)
null$Y.other.adj = p.adjust(null$Y.other)
null$mean.GER.adj = p.adjust(null$mean.GER)
null$mean.NC.adj = p.adjust(null$mean.NC)
null$mean.other.adj = p.adjust(null$mean.other)
null$mean.GER.greater.adj = p.adjust(null$mean.GER.greater)
null$mean.NC.greater.adj = p.adjust(null$mean.NC.greater)
null$mean.other.greater.adj = p.adjust(null$mean.other.greater)
null$mean.GER.lesser.adj = p.adjust(null$mean.GER.lesser)
null$mean.NC.lesser.adj = p.adjust(null$mean.NC.lesser)
null$mean.other.lesser.adj = p.adjust(null$mean.other.lesser)
View(subset(null, tissue %!in% c('BrainCortex','BrainHippocampus','BrainSpinalcordcervicalc1')))
View(null)
saveRDS(null, 'null.rds')

names(c) = names(g) = tissuesnow

c
c$BrainCortex[,2] = c(0,0,0,0)
c$BrainCortex[,3] = c(0,0,0,0)
c$BrainSpinalcordcervicalc1[,2] = c(0,0,0,0)
c$BrainSpinalcordcervicalc1[,3] = c(0,0,0,0)
lapply(c, FUN = prop.table, margin = 1)
c2 = reshape2::melt(c)
c2

g = reshape2::melt(g)
g
g %>% group_by(Functional_Category2) %>% summarise(mean = mean(value))

colnames(ts) = c('tissue',
                 'chi.p',
                 'post.fisher.ger.nc',
                 'post.fisher.ger.non',
                 'post.fisher.ger.other',
                 'post.fisher.nc.non',
                 'post.fisher.nc.other',
                 'post.fisher.nonASD.other',
                 'anova.p',
                 'tukey.nc.ger',
                 'tukey.non.ger',
                 'tukey.other.ger',
                 'tukey.non.nc',
                 'tukey.other.nc',
                 'tukey.other.non')
ts
ts$chi.p.adj = p.adjust(ts$chi.p)
ts$anova.p.adj = p.adjust(ts$anova.p)
# tukey and fisher already p adjusted

saveRDS(ts, 'ts.rds')
saveRDS(c, 'c.rds')
saveRDS(g, 'g.rds')

write.csv(ts, 'ts.csv')
write.csv(null, 'null.csv')

#### end ####

#### ANOVA (excl non ASD) ####

tss = data.frame()

for(i in 1:length(tissuesnow)){
  
  print(tissuesnow[i])
  tss[i,1] = tissuesnow[i]
  t = subset(b, tissue == tissuesnow[i])

  an = aov(diff_xy_new ~ Functional_Category2, data = t)
  sa = summary(an)
  tss[i,2] = sa[[1]][1,5]
  tu = TukeyHSD(an)
  tss[i,3] = tu$Functional_Category2[1,4]
  tss[i,4] = tu$Functional_Category2[2,4]
  tss[i,5] = tu$Functional_Category2[3,4]
}

colnames(tss) = c('tissue',
                 'anova.p',
                 'tukey.nc.ger',
                 'tukey.other.ger',
                 'tukey.other.nc')
tss$anova.p.adj = p.adjust(tss$anova.p)
table(tss$anova.p.adj < 0.05)
table(tss$tukey.nc.ger < 0.05)
table(tss$tukey.other.ger < 0.05)
table(tss$tukey.other.nc < 0.05)
View(tss)

write.csv(tss, 'tss.csv')

#### end ####

#### Figure 5e ####

pl = b
pl$Functional_Category2[is.na(pl$Functional_Category2)] = 'nonASD'
pl = subset(pl, tissue %in% tissuesnow)
pl = subset(pl, tissue %!in% c('BrainHippocampus','BrainCerebellarHemisphere','BrainCerebellum'))
table(pl$tissue)

pl$Functional_Category2 = factor(pl$Functional_Category2,
                                 levels = c('Gene expression regulation',
                                            'nonASD','Other','Neuronal communication'))
levels(pl$Functional_Category2) = c('GER','non-ASD','other','NC')
pl = subset(pl, Functional_Category2 != 'non-ASD')
pl$tissue = as.factor(pl$tissue)
levels(pl$tissue) = sl[c(1:3,6:7,9:13)]
table(pl$tissue, pl$Functional_Category2)

ggplot(pl, aes(x = Functional_Category2, y = diff_xy_new, fill = Functional_Category2))  +
  geom_boxplot() +
  facet_wrap(~tissue, nrow = 2) +
  ylab('Overall asymmetric coupling') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_fill_manual(values = c('#1B9E77','#666666','#D0990B')) +
  theme_article() +
  theme(axis.title.x = element_blank(),
        legend.position = 'none')
  
#### end ####

#### Figure 5F  ####

gam_list = read.csv('gametologs_in_genome.csv')
gam_list$pair = as.factor(gam_list$pair)
gams = data.frame()
for(i in 1:length(levels(gam_list$pair))){
  gams[i,1] = levels(gam_list$pair)[i]
  xgene = subset(gam_list, pair == levels(gam_list$pair)[i] & Gametolog == 'X')$common_name
  ygene = subset(gam_list, pair == levels(gam_list$pair)[i] & Gametolog == 'Y')$common_name
  gams[i,2] = paste(xgene, ygene, sep = " & ")
  gams[i,3] = subset(gam_list, pair == levels(gam_list$pair)[i] & Gametolog == 'X')$ensembl_id
  gams[i,4] = subset(gam_list, pair == levels(gam_list$pair)[i] & Gametolog == 'Y')$ensembl_id
}
colnames(gams) = c('pair','gams','X','Y')

tnows = 'BrainAnteriorcingulatecortexBA24'
tnow = "Brain - Anterior cingulate cortex (BA24)"

xygam = readRDS('out_xy.rds')
colnames(xygam)[9] = 'tissue'
colnames(xygam)[10] = 'gene'
pl = subset(xygam, gene %in% asd$gene)
pl = merge(pl, asd, by = 'gene')
colnames(pl)[11] = 'gams'
pl = merge(pl, gams, by = 'gams')
pl = subset(pl, tissue == tnows)

xyclip = readRDS(file = 'xy_clip_per_gam.rds')
clip = subset(xyclip, tissue == tnow)
clip = subset(clip, gene %in% asd$gene)
colnames(clip)[7] = 'gams'
pl = merge(pl, clip, by = c('gene','gams'))
pl = subset(pl, padj < 0.05)

exp = readRDS(paste(tnows,'_adjusted_exp_MALES.rds',sep=""))
mexp = rowMeans(exp, na.rm = T)
mexp = data.frame(mexp)
mexp$name = rownames(mexp)

table(pl$x_coexp > pl$y_coexp, pl$Functional_Category2)

# only genes sig overall coupling in tissue
k = subset(h, tissue == tnows & padj < 0.05)
pl = subset(pl, gene %in% k$gene)
k = subset(k, gene %in% pl$gene)

table(pl$x_coexp > pl$y_coexp, pl$Functional_Category2)

# origin to GER NC Other Xgam Ygam
d1 = data.frame(from="origin", to=c('GER','X','NC','Other','Y'))
ap = unique(asd[,c('Functional_Category2','gene')])
colnames(ap) = c('from','to')
ap$from = ifelse(ap$from == 'Gene expression regulation', 'GER', ap$from)
ap$from = ifelse(ap$from == 'Neuronal communication', 'NC', ap$from)
gp = unique(gam_list[,c('Gametolog','ensembl_id')])
colnames(gp) = c('from','to')
d2 = rbind(ap, gp)
hierarchy <- rbind(d1, d2)

all_leaves = hierarchy$to[6:dim(hierarchy)[1]]

o = data.frame()
for(i in 1:length(ap$to)){
  genenow = ap$to[i]
  xynow = subset(pl, gene == genenow)
  xynow$con = ifelse(xynow$x_coexp > xynow$y_coexp, xynow$X, xynow$Y)
  know = subset(k, gene == genenow)
  xynow$overall = ifelse(know$diff_xy_new > 0, 'X','Y')
  o = rbind(o, xynow)
}
connect = data.frame(from = o$gene, to = o$con, diff = o$diff, overall = o$overall)
connect$gam = ifelse(connect$to %in% gams$X,1,3)
head(connect)
#connect = connect[1:10,]

vertices  =  data.frame(
  name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to)))) 
vertices$group  =  hierarchy$from[ match( vertices$name, hierarchy$to ) ]
vertices[3,1] = 'X'
vertices[4,1] = 'Other'
vertices[5,1] = 'NC'
hierarchy$from = factor(hierarchy$from, levels = c('origin','GER','X','Other','NC','Y'))
vertices$group = factor(vertices$group, levels = c('NA','origin','GER','X','Other','NC','Y'))
vertices = merge(vertices, mexp, all.x = T)
colnames(vertices)[3] = 'value'
vertices = vertices[order(vertices$group),]
head(vertices)
vertices[2,1] = 'X'
vertices[3,1] = 'Other'
vertices[4,1] = 'NC'
mygraph = graph_from_data_frame(hierarchy, vertices=vertices)
from  =  match( connect$from, vertices$name)
to  =  match( connect$to, vertices$name)

genenames = data.frame(rbind(data.frame(id = asd$gene, name = asd$name),
                             data.frame(id = gam_list$ensembl_id, name = gam_list$common_name)))
genenames = unique(genenames)
head(genenames)

for(i in 1:length(hierarchy$to)){
  tn = hierarchy$to[i]
  if(tn %in% genenames$id){hierarchy$to[i] = genenames[which(genenames$id == tn),2]} else {hierarchy$to[i] = hierarchy$to[i]}}
for(i in 1:length(vertices$name)){
  tn = vertices$name[i]
  if(tn %in% genenames$id){vertices$name[i] = genenames[which(genenames$id == tn),2]} else {vertices$name[i] = vertices$name[i]}}
for(i in 1:length(connect$from)){
  tn = connect$from[i]
  if(tn %in% genenames$id){connect$from[i] = genenames[which(genenames$id == tn),2]} else {connect$from[i] = connect$from[i]}}
for(i in 1:length(connect$to)){
  tn = connect$to[i]
  if(tn %in% genenames$id){connect$to[i] = genenames[which(genenames$id == tn),2]} else {connect$to[i] = connect$to[i]}}

vertices$id <- NA
myleaves <- which(is.na( match(vertices$name, hierarchy$from) ))
nleaves <- length(myleaves)
vertices$id[ myleaves ] <- seq(1:nleaves)
vertices$angle <- 90 - 360 * vertices$id / nleaves
vertices$hjust <- ifelse( vertices$angle < -90, 1, 0)
vertices$angle <- ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)
head(vertices)

mygraph = graph_from_data_frame(hierarchy, vertices=vertices)
from  =  match( connect$from, vertices$name)
to  =  match( connect$to, vertices$name)

cols = ifelse(connect$gam == 1, tissue.colors[7], tissue.colors[14])
al = ifelse(connect$overall == 'X' & connect$diff > 0, 0.8, 0.4)
al = ifelse(connect$overall == 'Y' & connect$diff < 0, 0.8, al)
table(al)

ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group, size=value, alpha = group)) +
  scale_color_manual(values = c('#1B9E77','#D0990B','#666666',tissue.colors[c(7,14)])) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,1,1)) +
  geom_conn_bundle(data = get_con(from = from, to = to, diff = connect$diff, gam = connect$gam), 
                   alpha = rep(al, each=100), color = rep(cols, each=100), 
                   width=0.9, tension = 0.6) +
  geom_node_text(aes(x = x*1.1, y=y*1.1, filter = leaf, label=name, angle = angle, hjust=hjust), size=1.5, alpha=1) +
  theme_void() +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))

#### end ####

#### Figure S14 ####

xygam = readRDS('out_xy.rds')
colnames(xygam)[9] = 'tissue'
colnames(xygam)[10] = 'gene'
pl = subset(xygam, gene %in% asd$gene)
pl = merge(pl, asd, by = 'gene')
colnames(pl)[11] = 'gams'
pl = merge(pl, gams, by = 'gams')
pl = subset(pl, tissue == tnows)

xyclip = readRDS(file = 'xy_clip_per_gam.rds')
clip = subset(xyclip, tissue == tnow)
clip = subset(clip, gene %in% asd$gene)
colnames(clip)[7] = 'gams'
pl = merge(pl, clip, by = c('gene','gams'))

k = subset(h, tissue == tnows & padj < 0.05)
pl = subset(pl, gene %in% k$gene)
k = subset(k, gene %in% pl$gene)

ap = unique(asd[,c('Functional_Category2','gene')])
colnames(ap) = c('from','to')
ap$from = ifelse(ap$from == 'Gene expression regulation', 'GER', ap$from)
ap$from = ifelse(ap$from == 'Neuronal communication', 'NC', ap$from)

o = data.frame()
for(i in 1:length(ap$to)){
  genenow = ap$to[i]
  xynow = subset(pl, gene == genenow)
  xynow$con = ifelse(xynow$x_coexp > xynow$y_coexp, xynow$X, xynow$Y)
  know = subset(k, gene == genenow)
  xynow$overall = ifelse(know$diff_xy_new > 0, 'X','Y')
  o = rbind(o, xynow)
}

r = o
r$diff = ifelse(r$overall == "X", 1, -1)
r$gams = 'overall'
r$code = 'overall'
o = rbind(cbind(o, code = 'each'), r)
  
ggplot(o, aes(x = gams, y = name, fill = diff)) +
  geom_tile() +
  geom_tile(data = o[which(o$padj < 0.05),], aes(x = gams, y = name), fill = "transparent", color = "black") +
  facet_grid(Functional_Category2~code, scales = 'free', space = 'free') +
  theme_article() +
  scale_fill_gradient2(high=tissue.colors[7],mid="white",
                       low=tissue.colors[14],
                       na.value="lightgrey", 
                       midpoint=0, 
                       name = bquote(sCFD[MXMY])) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x.top = element_blank(),
        axis.title = element_blank())

#### end ####
