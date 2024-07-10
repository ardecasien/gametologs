library(mashr)
library(reshape2)

#### load data ####

# load sex-biased estimates (per gene per tissue)

m = readRDS(file = 'mashr_sex.rds')

mash.beta = get_pm(m)
mash.beta = reshape2::melt(mash.beta)
mash.beta$Var2 = as.factor(mash.beta$Var2)
levels(mash.beta$Var2) = nospace.list[-c(20,36,42)]

mash.lfsr = get_lfsr(m)
mash.lfsr = reshape2::melt(mash.lfsr)
mash.lfsr$Var2 = as.factor(mash.lfsr$Var2)
levels(mash.lfsr$Var2) = nospace.list[-c(20,36,42)]

mash.res = cbind(mash.beta, mash.lfsr)
mash.res = mash.res[,c(1,2,3,6)]
colnames(mash.res) = c('gene','tissue','beta','lfsr')
mash.res$key = paste(mash.res$tissue, mash.res$gene, sep = '_')
dim(mash.res)

# Table S20

write.csv(mash.res, file = 'mash.res.csv')

# load CLIP sig diff summed X-Y coupling (per gene per tissue)

xy_clip = readRDS(file = 'xy_clip.rds')
xy_clip$tissue = as.factor(xy_clip$tissue)
levels(xy_clip$tissue) = nospace.list
xy_clip$key = paste(xy_clip$tissue, xy_clip$gene, sep = "_")
dim(xy_clip)
xy_clip$padj_all = p.adjust(xy_clip$p, method = 'BH')

# load weighted avg diff X-Y coupling (per gene per tissue)

out_xy_weight = readRDS(file = 'out_xy_weight.rds')
out_xy_weight$key = paste(out_xy_weight$region, out_xy_weight$gene.x, sep = "_")
dim(out_xy_weight)

out_xy_weight = merge(out_xy_weight, xy_clip, by = 'key')
dim(out_xy_weight)

# Table S11

write.csv(out_xy_weight, file = 'out_xy_weight.csv')

#### end ####

#### combine and plot ####

combo = merge(mash.res, out_xy_weight, by = 'key')
combo = combo[,-2]
dim(combo)

# subset and plot
# malebiased = positive beta
# xcoupled = positive diff_xy_new

library(biomaRt)

hsap = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", versio = 104)
gam_bm = getBM(attributes = c("external_gene_name","ensembl_gene_id", "chromosome_name","gene_biotype"), 
               filters = "ensembl_gene_id", values = unique(combo$gene.x), mart = hsap)
Ygenes = subset(gam_bm, chromosome_name == 'Y')
Xgenes = subset(gam_bm, chromosome_name == 'X')

pl = subset(combo, lfsr < 0.05 # sex-biased
                    & padj < 0.05 # asymmetrically coupled
                    & gene.x %!in% Ygenes$ensembl_gene_id # no Y chr
                    & gene.x %!in% Xgenes$ensembl_gene_id # no X chr
                    & gene.x %in% gam_bm$ensembl_gene_id)

cor.test(pl$beta, -pl$diff_xy_new, method = 'spearman')

pl$code = ifelse(pl$diff_xy_new > 0 & pl$beta < 0, 'concordant', 'discordant')
pl$code = ifelse(pl$diff_xy_new < 0 & pl$beta > 0, 'concordant', pl$code)
table(pl$code) / sum(table(pl$code))
cc = table(pl$diff_xy_new < 0, pl$beta > 0)
mat = as.matrix(cc)
fisher.test(mat, alternative = 'greater')

dim(pl)
length(unique(pl$gene.x))

# Figure 5a

ggplot(pl, aes(x = diff_xy_new, y = -beta)) +
  geom_point(aes(x = -diff_xy_new, y = beta, color = code),  alpha = 0.4) +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.6,0.6)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.6,0.6)) +
  scale_color_manual(values=c('#407076','#C9C5BA')) +
  ylab("Sex-biased expression") +
  xlab("Differential X-Y coupling") +
  geom_smooth(method = 'lm', se = TRUE, color = 'black', linewidth = 2) +
  theme_article() +
  annotate(geom="text", x=-0.45, y=0.45, label="67% concordant \n OR = 4.23 \n p < 0.01", size = 5) +
  annotate(geom="text", x=0.5, y=-0.5, label="rho = 0.32 \n p < 0.01", size = 5) +
  theme(axis.text = element_text(size=18),
        legend.position = 'none',
        axis.title = element_text(size=18))

# Figure S13

ggplot(pl, aes(x = diff_xy_new, y = -beta, color = tissue.x)) +
  geom_point(alpha = 0.4) +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.6,0.6)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.6,0.6)) +
  scale_color_manual(values=tissue.colors[-c(20,36,42)]) +
  ylab("Sex-biased expression") +
  xlab("Differential X-Y coupling \n (weighted mean)") +
  geom_smooth(method = 'lm', se = FALSE, alpha = 0.4) +
  theme_article() +
  theme(axis.text = element_text(size=18),
        legend.position = 'none',
        axis.title = element_text(size=18))

pl$tissue.x = droplevels(pl$tissue.x)
r = data.frame()
for (i in 1:length(levels(pl$tissue.x))){
  dat = subset(pl, tissue.x == levels(pl$tissue.x)[i])
  cc = cor.test(-dat$beta, dat$diff_xy_new, method = 'spearman')
  r[i,1] = levels(pl$tissue.x)[i]
  r[i,2] = dim(dat)[1]
  r[i,3] = cc$estimate
  r[i,4] = cc$p.value
}

r$padj = p.adjust(r$V4, method = 'BH')
r
table(r$V3 > 0, r$padj < 0.05)

# Table S21

write.csv(r, file = 'sex-bias-cor-per-tissue.csv')

#### end ####

library(biomaRt)

#### load data ####

## load weighted diff XY coupling

out_xy_weight = readRDS('out_xy_weight.rds')
out_xy_weight$region = factor(out_xy_weight$region)
levels(out_xy_weight$region) = tissue.list
out_xy_weight$key = paste(out_xy_weight$region, out_xy_weight$gene.x, sep = '_')
dim(out_xy_weight)

clip = readRDS('xy_clip.rds')
clip$padj_all = p.adjust(clip$p, method = 'fdr')
clip$key = paste(clip$tissue, clip$gene, sep = '_')
dim(clip)

out_xy_weight = merge(out_xy_weight, clip, by = 'key')

## load hartman
## SexbiasCallNr	= How often is the gene called sex-biased out of the 100 different permutations
hart = read.csv('hartman2021.csv')
hart$key = paste(hart$Tissue, hart$Ensembl, sep = '_')
table(hart$Tissue, hart$Bias)
dim(hart)

#### end ####

#### combine and plot | Figure 5b ####

combo = merge(out_xy_weight, hart, by = 'key')
dim(combo)

hsap = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 104)
gam_bm = getBM(attributes = c("external_gene_name","ensembl_gene_id", "chromosome_name"), 
               filters = "ensembl_gene_id", values = unique(combo$gene.x), mart = hsap)

colnames(combo)[2] = 'ensembl_gene_id'
combo = merge(combo, gam_bm, by = 'ensembl_gene_id', all.x = T)

## code and test across all
combo$coup = ifelse(combo$diff_xy_new > 0, 'X', 'Y')

## coupling cutoff
#combo = subset(combo, padj >= 0.05)
combo$Bias = ifelse(combo$padj > 0.05, "none", combo$Bias)

## coexp cutoff
combo = subset(combo, SexbiasCallNr >= 50)

# remove sex chr genes
combo = subset(combo, chromosome_name != "X")
combo = subset(combo, chromosome_name != "Y")
dim(combo)
length(unique(combo$gene))

# fisher tests

XandF = subset(combo, Bias == "female" & coup == 'X')
justF = subset(combo, Bias == "female")
justF = subset(justF, key %!in% XandF$key)
justX = subset(combo, coup == 'X')
justX = subset(justX, key %!in% XandF$key)
neither = subset(combo, key %!in% c(XandF$key, justF$key, justX$key))
fX = matrix(c(dim(XandF)[1], dim(justF)[1], dim(justX)[1], dim(neither)[1]), nrow = 2, ncol = 2)
fX
fisher.test(fX, alternative = 'greater')

XandM = subset(combo, Bias == "male" & coup == 'X')
justM = subset(combo, Bias == "male")
justM = subset(justM, key %!in% XandM$key)
justX = subset(combo, coup == 'X')
justX = subset(justX, key %!in% XandM$key)
neither = subset(combo, key %!in% c(XandM$key, justM$key, justX$key))
mX = matrix(c(dim(XandM)[1], dim(justM)[1], dim(justX)[1], dim(neither)[1]), nrow = 2, ncol = 2)
mX
fisher.test(mX, alternative = 'greater')

YandM = subset(combo, Bias == "male" & coup == 'Y')
justM = subset(combo, Bias == "male")
justM = subset(justM, key %!in% YandM$key)
justY = subset(combo, coup == 'Y')
justY = subset(justY, key %!in% YandM$key)
neither = subset(combo, key %!in% c(YandM$key, justM$key, justY$key))
mY = matrix(c(dim(YandM)[1], dim(justM)[1], dim(justY)[1], dim(neither)[1]), nrow = 2, ncol = 2)
mY
fisher.test(mY, alternative = 'greater')

YandF = subset(combo, Bias == "female" & coup == 'Y')
justF = subset(combo, Bias == "female")
justF = subset(justF, key %!in% YandF$key)
justY = subset(combo, coup == 'Y')
justY = subset(justY, key %!in% YandF$key)
neither = subset(combo, key %!in% c(YandF$key, justF$key, justY$key))
fY = matrix(c(dim(YandF)[1], dim(justF)[1], dim(justY)[1], dim(neither)[1]), nrow = 2, ncol = 2)
fY
fisher.test(fY, alternative = 'greater')

## code and test for each tissue

res = data.frame()

tissue.now = tissue.list[which(tissue.list %in% hart$Tissue)]

for(i in 1:length(tissue.now)){
  
  print(tissue.now[i])
  dat = subset(combo, region == tissue.now[i])

  XandF = subset(dat, Bias == "female" & coup == 'X')
  justF = subset(dat, Bias == "female")
  justF = subset(justF, key %!in% XandF$key)
  justX = subset(dat, coup == 'X')
  justX = subset(justX, key %!in% XandF$key)
  neither = subset(dat, key %!in% c(XandF$key, justF$key, justX$key))
  mX = matrix(c(dim(XandF)[1], dim(justF)[1], dim(justX)[1], dim(neither)[1]), nrow = 2, ncol = 2)
  mX
  f = fisher.test(mX, alternative = 'greater')
  res[i,1] = f$estimate
  res[i,2] = f$p.value
  
  YandM = subset(dat, Bias == "male" & coup == 'Y')
  justM = subset(dat, Bias == "male")
  justM = subset(justM, key %!in% YandM$key)
  justY = subset(dat, coup == 'Y')
  justY = subset(justY, key %!in% YandM$key)
  neither = subset(dat, key %!in% c(YandM$key, justM$key, justY$key))
  mY = matrix(c(dim(YandM)[1], dim(justM)[1], dim(justY)[1], dim(neither)[1]), nrow = 2, ncol = 2)
  mY
  f = fisher.test(mY, alternative = 'greater')
  res[i,3] = f$estimate
  res[i,4] = f$p.value
  
  XandF = subset(dat, Bias == "male" & coup == 'X')
  justF = subset(dat, Bias == "male")
  justF = subset(justF, key %!in% XandF$key)
  justX = subset(dat, coup == 'X')
  justX = subset(justX, key %!in% XandF$key)
  neither = subset(dat, key %!in% c(XandF$key, justF$key, justX$key))
  mX = matrix(c(dim(XandF)[1], dim(justF)[1], dim(justX)[1], dim(neither)[1]), nrow = 2, ncol = 2)
  mX
  f = fisher.test(mX, alternative = 'greater')
  res[i,5] = f$estimate
  res[i,6] = f$p.value
  
  YandM = subset(dat, Bias == "female" & coup == 'Y')
  justM = subset(dat, Bias == "female")
  justM = subset(justM, key %!in% YandM$key)
  justY = subset(dat, coup == 'Y')
  justY = subset(justY, key %!in% YandM$key)
  neither = subset(dat, key %!in% c(YandM$key, justM$key, justY$key))
  mY = matrix(c(dim(YandM)[1], dim(justM)[1], dim(justY)[1], dim(neither)[1]), nrow = 2, ncol = 2)
  mY
  f = fisher.test(mY, alternative = 'greater')
  res[i,7] = f$estimate
  res[i,8] = f$p.value
  
  res[i,9] = tissue.now[i]

}

colnames(res) = c('F X OR','F X P','M Y OR','M Y P','M X OR','M X P','F Y OR','F Y P')
res$FXpadj = p.adjust(res$`F X P`, method = 'BH')
res$MXpadj = p.adjust(res$`M X P`, method = 'BH')
res$FYpadj = p.adjust(res$`F Y P`, method = 'BH')
res$MYpadj = p.adjust(res$`M Y P`, method = 'BH')

View(res)

write.csv(res, file = 'hartman_comparison.csv')

## plot

library(ggplot2)
library(egg)

m = matrix(c(1.321, 0.639, 0.757, 1.565), ncol = 2, nrow = 2)
m
rownames(m) = c('female-biased','male-biased')
colnames(m) = c('X-biased','Y-biased')
m
pl = reshape2::melt(m)
pl
ggplot(pl, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() + 
  scale_x_discrete(limits=rev) +
  scale_y_discrete(limits=rev) +
  scale_fill_gradient2(midpoint = 1, high = "#E6AB02", low = 'grey', limits = c(0.6, 1.6)) +
  coord_equal() +
  theme_article() +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5))

#### end ####


