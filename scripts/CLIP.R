library(dplyr)
library(rcompanion)
library(biomaRt)

#### upload meta data ####

keep_samples = readRDS('gtex_combined_meta.rds')
keep_samples = subset(keep_samples, SMTSD %!in% remove.tissues)
table(keep_samples$SMTSD, keep_samples$SEX)
table(keep_samples$SEX)
table(keep_samples$SMTSD)
length(levels(as.factor(keep_samples$SMTSD)))

keep_samples %>% group_by(SEX) %>% summarise(mean=mean(AGE))
s = data.frame(keep_samples %>% group_by(SMTSD, SEX) %>% summarise(n=n()))
mean(subset(s, SEX == 1)$n, na.rm = T)
mean(subset(s, SEX == 2)$n, na.rm = T)

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
pairs = levels(as.factor(gam_list$pair))

#### end ####

#### CLIP per gam (X vs Y in males) ####

m = subset(pheno, SEX == 1)

xy_clip_per_gam = data.frame()
xy_clip_per_gam_subsamp = data.frame()

for (i in 1:length(tissue.list)){
  
  print(paste("Now analyzing:",tissue.list[i]))
  
  # get current data set
  m_now = subset(keep_samples, SMTSD == tissue.list[i] & SAMPID2 %in% m$SAMPID2)
  
  # limit to samples with all data available
  m_now = m_now[complete.cases(m_now[,c('AGE','SMTSISCH','SMRIN','SMNTRNRT')]),]
  
  # subsample tissues to lowest N=66 males
  subsamp = sample(m_now$SAMPID2, 66)
  m_now = subset(m_now, SAMPID2 %in% subsamp)
  
  if(length(m_now$SAMPID) == 0) {
    print('not enough samples') } else {
      
      # load adjusted expression
      exp_now = readRDS(file = paste(nospace.list[i],'_adjusted_exp_MALES.rds',sep=""))
      exp_now = exp_now[,which(colnames(exp_now) %in% m_now$SAMPID)]
      
      # scale adjusted expression
      scaled_exp = apply(exp_now, 1, scale)
      scaled_exp = t(scaled_exp)
      colnames(scaled_exp) = colnames(exp_now)
      
      # estimate products per gam
      
      for(k in 1:length(pairs)) {
      
        print(pairs[k])
        
        xgam = scaled_exp[which(rownames(scaled_exp) %in% subset(gam_list, pair == pairs[k] & Gametolog == 'X')$ensembl_id),]
        xgam = c(xgam)
      
        ygam = scaled_exp[which(rownames(scaled_exp) %in% subset(gam_list, pair == pairs[k] & Gametolog == 'Y')$ensembl_id),]
        ygam = c(ygam)
      
        nongam = scaled_exp[which(rownames(scaled_exp) %!in% gam_list$ensembl_id),]
      
        xprod = tryCatch(sweep(nongam, MARGIN=2, xgam, `*`), error=function(err) NA)
        yprod = tryCatch(sweep(nongam, MARGIN=2, ygam, `*`), error=function(err) NA)
      
        out_p = data.frame()
        
        if(is.na(xprod)[1,1] == TRUE | is.na(yprod)[1,1] == TRUE) {
          out_p[1,1] = NA
          out_p[1,2] = NA
          out_p$gene[1] = NA
          out_p$padj[1] = NA
          out_p$tissue[1] = tissue.list[i]
          out_p$pair[1] = pairs[k]
          #xy_clip_per_gam = rbind(xy_clip_per_gam, out_p)
          xy_clip_per_gam_subsamp = rbind(xy_clip_per_gam_subsamp, out_p)
          
        } else {
        for(j in 1:length(rownames(xprod))){
        
          tt = t.test(xprod[j,], yprod[j,], paired = TRUE)
          out_p[j,1] = tt$statistic
          out_p[j,2] = tt$p.value
        
      }
      
        out_p$gene = rownames(xprod)
        out_p$padj = p.adjust(out_p$V2, method = 'BH')
        out_p$tissue = tissue.list[i]
        out_p$pair = pairs[k]
        #xy_clip_per_gam = rbind(xy_clip_per_gam, out_p)
        xy_clip_per_gam_subsamp = rbind(xy_clip_per_gam_subsamp, out_p)
        
        }
      }
    }
  }

# positive t = X > Y product = X > Y correlations

xy_clip_per_gam = merge(xy_clip_per_gam, gams, by = 'pair')
xy_clip_per_gam$padj_all = p.adjust(xy_clip_per_gam$V2, method = 'BH')
colnames(xy_clip_per_gam) = c('pair','t','p','gene','padj','tissue','gene_pair','padj_all')
xy_clip_per_gam2 = data.frame()
for(i in 1:length(tissue.list)){
  dnow = subset(xy_clip_per_gam, tissue == tissue.list[i])
  dnow$padj_tissuelevel = p.adjust(dnow$p, method = 'BH')
  xy_clip_per_gam2 = rbind(xy_clip_per_gam2, dnow)
}

xy_clip_per_gam_subsamp = merge(xy_clip_per_gam_subsamp, gams, by = 'pair')
xy_clip_per_gam_subsamp$padj_all = p.adjust(xy_clip_per_gam_subsamp$V2, method = 'BH')
colnames(xy_clip_per_gam_subsamp) = c('pair','t','p','gene','padj','tissue','gene_pair','padj_all')
xy_clip_per_gam_subsamp2 = data.frame()
for(i in 1:length(tissue.list)){
  dnow = subset(xy_clip_per_gam_subsamp, tissue == tissue.list[i])
  dnow$padj_tissuelevel = p.adjust(dnow$p, method = 'BH')
  xy_clip_per_gam_subsamp2 = rbind(xy_clip_per_gam_subsamp2, dnow)
}

table(xy_clip_per_gam$tissue, xy_clip_per_gam$pair, xy_clip_per_gam$t < 0, xy_clip_per_gam$padj < 0.05)
table(xy_clip_per_gam$tissue, xy_clip_per_gam$pair, xy_clip_per_gam$t < 0, xy_clip_per_gam$padj_all < 0.05)
table(xy_clip_per_gam2$tissue, xy_clip_per_gam2$pair, xy_clip_per_gam2$t < 0, xy_clip_per_gam2$padj_tissuelevel < 0.05)

saveRDS(xy_clip_per_gam, file = 'xy_clip_per_gam.rds')
saveRDS(xy_clip_per_gam_subsamp, file = 'xy_clip_per_gam_subsamp.rds')

#### end ####

#### compare total to subsamp | Figures S9e-g ####

xy_clip_per_gam = readRDS('xy_clip_per_gam.rds')
xy_clip_per_gam_subsamp = readRDS('xy_clip_per_gam_subsamp.rds')

xy_clip_per_gam$sig = ifelse(xy_clip_per_gam$padj < 0.05, 'sig', 'ns')
g1 = xy_clip_per_gam %>% group_by(tissue, gene_pair, sig) %>% summarise(n())
g1$key = paste(g1$tissue, g1$gene_pair, g1$sig, sep = "_")
g1$key2 = paste(g1$tissue, g1$gene_pair, g1$sig, sep = "_")

xy_clip_per_gam_subsamp$sig = ifelse(xy_clip_per_gam_subsamp$padj < 0.05, 'sig', 'ns')
g2 = xy_clip_per_gam_subsamp %>% group_by(tissue, gene_pair, sig) %>% summarise(n())
g2$key = paste(g2$tissue, g2$gene_pair, g2$sig, sep = "_")

comp = merge(g1, g2, by = 'key', all = T)
comp[is.na(comp)] = 0
comp = subset(comp, sig.x == 'sig')
cor.test(comp$`n().x`, comp$`n().y`, method = 'spearman')

ggplot(comp, aes(x = `n().x`, y = `n().y`)) +
  geom_point() +
  xlab('Total (N significant)') +
  ylab('Subsample (N significant)') +
  geom_smooth(method = 'lm') +
  theme_article()

comp2 = comp %>% group_by(tissue.x) %>% summarise(x = sum(`n().x`), y = sum(`n().y`))
cor.test(comp2$x, comp2$y, method = 'spearman')
ggplot(comp2, aes(x = x, y = y)) +
  geom_point() +
  xlab('Total (N significant per tissue)') +
  ylab('Subsample (N significant per tissue') +
  geom_smooth(method = 'lm') +
  theme_article()

comp2 = comp %>% group_by(gene_pair.x) %>% summarise(x = sum(`n().x`), y = sum(`n().y`))
cor.test(comp2$x, comp2$y, method = 'spearman')
ggplot(comp2, aes(x = x, y = y)) +
  geom_point() +
  xlab('Total (N significant per pair)') +
  ylab('Subsample (N significant per pair') +
  geom_smooth(method = 'lm') +
  theme_article()

#### end ####

#### CLIP total X vs Y in males ####

m = subset(pheno, SEX == 1)

xy_clip = data.frame()

for (i in 1:length(tissue.list)){
  
  print(paste("Now analyzing:",tissue.list[i]))
  
  # get current data set
  m_now = subset(keep_samples, SMTSD == tissue.list[i] & SAMPID2 %in% m$SAMPID2)
  
  # limit to samples with all data available
  m_now = m_now[complete.cases(m_now[,c('AGE','SMTSISCH','SMRIN','SMNTRNRT')]),]
  
  if(length(m_now$SAMPID) == 0) {
    print('not enough samples') } else {
      
      # load adjusted expression
      
      exp_now = readRDS(file = paste(nospace.list[i],'_adjusted_exp_MALES.rds',sep=""))
      exp_now = exp_now[,which(colnames(exp_now) %in% m_now$SAMPID)]
      
      # scale adjusted expression
      
      scaled_exp = apply(exp_now, 1, scale)
      scaled_exp = t(scaled_exp)
      colnames(scaled_exp) = colnames(exp_now)
      
      # get summed x and y expression
      
      xgam = scaled_exp[which(rownames(scaled_exp) %in% subset(gam_list, Gametolog == 'X')$ensembl_id),]
      xgam = colSums(xgam)
      xgam = c(xgam)
      
      ygam = scaled_exp[which(rownames(scaled_exp) %in% subset(gam_list, Gametolog == 'Y')$ensembl_id),]
      ygam = colSums(ygam)
      ygam = c(ygam)
      
      nongam = scaled_exp[which(rownames(scaled_exp) %!in% gam_list$ensembl_id),]
      
      # products 
      
      xprod = sweep(nongam, MARGIN=2, xgam, `*`)
      yprod = sweep(nongam, MARGIN=2, ygam, `*`)
      
      # test difference
      
      out_p = data.frame()
      
      for(j in 1:length(rownames(xprod))){
        
        tt = t.test(xprod[j,], yprod[j,])
        out_p[j,1] = tt$statistic
        out_p[j,2] = tt$p.value
        out_p[j,3] = mean(xprod[j,])
        out_p[j,4] = mean(yprod[j,])
        out_p[j,5] = var(xprod[j,])
        out_p[j,6] = var(yprod[j,])
        
      }
      
      out_p$gene = rownames(xprod)
      out_p$padj = p.adjust(out_p$V2, method = 'BH')
      out_p$tissue = tissue.list[i]
      xy_clip = rbind(xy_clip, out_p)
      
    }}

# positive t = X > Y product = X > Y correlations

colnames(xy_clip)[1:6] = c('t','p','mean_x','mean_y','var_x','var_y')

table(xy_clip$tissue, xy_clip$t < 0, xy_clip$padj < 0.05)
table(xy_clip$tissue, xy_clip$var_x < xy_clip_mf$var_y)

saveRDS(xy_clip, file = 'xy_clip.rds')

xy_clip = readRDS('xy_clip.rds')

s = xy_clip %>% group_by(gene) %>% summarise(padj = min(padj))
table(s$padj < 0.05)

xy_clip$bias = ifelse(xy_clip$t > 0, "X", "Y")
s = xy_clip %>% group_by(bias, gene) %>% summarise(padj = min(padj))
table(subset(s, bias == "X")$padj < 0.05)
table(subset(s, bias == "Y")$padj < 0.05)
table(table(s$gene))

#### end ####

#### CLIP total XX vs XY in females vs males ####

m = subset(pheno, SEX == 1)
f = subset(pheno, SEX == 2)

xy_clip_mf = data.frame()

remove.tissues = c('Whole Blood','Breast - Mammary Tissue','Bladder','Cells - Cultured fibroblasts','Cells - EBV-transformed lymphocytes','Cells - Leukemia cell line (CML)','Cervix - Ectocervix','Cervix - Endocervix','Fallopian Tube','Kidney - Medulla','Ovary','Prostate','Testis','Uterus','Vagina')
tissue.list2 = tissue.list[which(tissue.list %!in% remove.tissues)]

for (i in 1:length(tissue.list2)){
  
  print(paste("Now analyzing:",tissue.list2[i]))
  
  # get current data set
  m_now = subset(keep_samples, SMTSD == tissue.list2[i])
  
  # limit to samples with all data available
  m_now = m_now[complete.cases(m_now[,c('AGE','SMTSISCH','SMRIN','SMNTRNRT')]),]
  
  if(length(m_now$SAMPID) == 0) {
    print('not enough samples') } else {
      
      # load male adjusted expression
      # load female adjusted expression
      
      me = readRDS(file = paste(tissue.list2[i],'_adjusted_exp_MALES.rds',sep=""))
      fe = readRDS(file = paste(tissue.list2[i],'_adjusted_exp_FEMALES.rds',sep=""))
      
      # same N in males and females
      
      dd = min(dim(me)[2], dim(fe)[2])
      fk = sample(colnames(fe), dd)
      mk = sample(colnames(me), dd)
      me = me[,mk]
      fe = fe[,fk]
      
      # scale male adjusted expression
      # get summed gam expression
      
      sme = apply(me, 1, scale)
      sme = t(sme)
      colnames(sme) = colnames(me)
      gam = sme[which(rownames(sme) %in% gam_list$ensembl_id),]
      gam = colSums(gam)
      gam = c(gam)
      ngm = sme[which(rownames(sme) %!in% gam_list$ensembl_id),]
      mprod = sweep(ngm, MARGIN=2, gam, `*`)
      
      # scale female adjusted expression
      # get summed gam expression
      
      sfe = apply(fe, 1, scale)
      sfe = t(sfe)
      colnames(sfe) = colnames(fe)
      gam = sfe[which(rownames(sfe) %in% gam_list$ensembl_id),]
      gam = colSums(gam)
      gam = c(gam)
      ngf = sfe[which(rownames(sfe) %!in% gam_list$ensembl_id),]
      fprod = sweep(ngf, MARGIN=2, gam, `*`)
      
      # same genes
      
      keep = rownames(fprod[which(rownames(fprod) %in% rownames(mprod)),])
      mprod = mprod[keep,]
      fprod = fprod[keep,]
      
      # test difference
      
      out_p = data.frame()
      
      for(j in 1:length(rownames(mprod))){
        
        tt = t.test(fprod[j,], mprod[j,])
        #boxplot(fprod[j,], mprod[j,])
        out_p[j,1] = tt$statistic
        out_p[j,2] = tt$p.value
        out_p[j,3] = mean(fprod[j,])
        out_p[j,4] = mean(mprod[j,])
        out_p[j,5] = var(fprod[j,])
        out_p[j,6] = var(mprod[j,])
      }
      
      out_p$padj = p.adjust(out_p$V2, method = 'BH')
      out_p$tissue = tissue.list2[i]
      xy_clip_mf = rbind(xy_clip_mf, out_p)
      
    }}

# positive t = F > M product = F > M correlations

table(xy_clip_mf$tissue, xy_clip_mf$V1 < 0, xy_clip_mf$padj < 0.05)
table(xy_clip_mf$tissue, xy_clip_mf$V5 < xy_clip_mf$V6)

saveRDS(xy_clip_mf, file = 'xy_clip_mf.rds')

#### end ####

#### compare CLIP t to differential X-Y coupling (per gene x pair x tissue) ####

xy_clip_per_gam = readRDS('xy_clip_per_gam.rds')
xy_clip_per_gam$tissue = factor(xy_clip_per_gam$tissue)
levels(xy_clip_per_gam$tissue) = nospace.list
xy_clip_per_gam$key = paste(xy_clip_per_gam$tissue, xy_clip_per_gam$gene_pair, xy_clip_per_gam$gene, sep = "_")

out_xy = readRDS(file = 'out_xy.rds') 
out_xy$key = paste(out_xy$region, out_xy$gene.y, out_xy$gene.x, sep = "_")

comp = merge(xy_clip_per_gam, out_xy, by = 'key')
cor.test(comp$t, comp$diff, method = 'pearson') 
cor.test(comp$t, comp$diff, method = 'spearman') 

# check
sig = subset(comp, padj_all < 0.05)
table(sig$t > 0, sig$diff > 0)

ggplot(comp, aes(x = diff, y = t)) +
  geom_point() + 
  xlab('Differential X-Y coupling') +
  ylab('CLIP t-estimate') +
  geom_smooth(method = 'lm') +
  theme_article()

ggplot(subset(comp, abs(diff) < 2.5), aes(x = diff, y = t)) +
  geom_point() + 
  xlab('Differential X-Y coupling') +
  ylab('CLIP t-estimate') +
  geom_smooth(method = 'lm') +
  theme_article()

#### end ####

#### compare CLIP to permutation approach (per gene x pair x tissue) #### 

xy_clip_per_gam = readRDS('xy_clip_per_gam.rds')
xy_clip_per_gam$tissue = factor(xy_clip_per_gam$tissue)
levels(xy_clip_per_gam$tissue) = nospace.list
xy_clip_per_gam$key = paste(xy_clip_per_gam$tissue, xy_clip_per_gam$`gene pair`, xy_clip_per_gam$gene, sep = "_")

out_xy = readRDS(file = 'out_xy.rds') 
out_xy$key = paste(out_xy$region, out_xy$gene.y, out_xy$gene.x, sep = "_")

comp = merge(xy_clip_per_gam, out_xy, by = 'key')
cor.test(comp$p, comp$padj_approach1, method = 'spearman')

ggplot(comp, aes(x = p, y = pval_approach1)) + 
  geom_point() +
  xlab('P-value (CLIP)') +
  ylab('P-value (permutation)') +
  geom_smooth(method = 'lm') +
  theme_article()

comp$sig_clip = ifelse(comp$p < 0.05, 'sig', 'ns')
comp$sig_perm = ifelse(comp$pval_approach1 < 0.05, 'sig', 'ns')
sig_clip = comp %>% group_by(tissue, gene.y, sig_clip) %>% summarise(n())
sig_clip$key = paste(sig_clip$tissue, sig_clip$gene.y, sig_clip$sig_clip)
sig_perm = comp %>% group_by(tissue, gene.y, sig_perm) %>% summarise(n())
sig_perm$key = paste(sig_perm$tissue, sig_perm$gene.y, sig_perm$sig_perm)

pl = merge(sig_clip, sig_perm, by = 'key', all = T)
pl[is.na(pl)] = 0 
pl = subset(pl, sig_clip == 'sig')

cor.test(pl$`n().x`, pl$`n().y`, method = 'spearman')
ggplot(pl, aes(x = `n().x`, y = `n().y`)) +
  geom_point() +
  xlab('CLIP (N significant)') +
  ylab('Permutation (N significant)') +
  geom_smooth(method = 'lm') +
  theme_article()

#### end ####

#### compare CLIP to differential X-Y coupling (weighted avg per gene x tissue) | Table S11 ####

xy_clip = readRDS('xy_clip.rds')
xy_clip$tissue = factor(xy_clip$tissue)
levels(xy_clip$tissue) = nospace.list
xy_clip$key = paste(xy_clip$tissue, xy_clip$gene, sep = "_")

out_xy = readRDS(file = 'out_xy_weight.rds') 
out_xy$key = paste(out_xy$region, out_xy$gene.x, sep = "_")

comp = merge(xy_clip, out_xy, by = 'key')
cor.test(comp$t, comp$diff, method = 'pearson') 
cor.test(comp$t, comp$diff, method = 'spearman') 

# Table S11

write.csv(comp, file = 'exp-weighted-avg-diff-coupling.csv')

ggplot(comp, aes(x = diff_xy_new, y = t)) +
  geom_point() + 
  xlab('Differential X-Y coupling (weighted average)') +
  ylab('CLIP t-estimate (summed)') +
  geom_smooth(method = 'lm') +
  theme_article()

#### end ####

#### load data for Figures 4C 4D S9 | Tables S8 S10 ####

xy_clip_per_gam = readRDS('xy_clip_per_gam.rds')
xy_clip_per_gam_subsamp = readRDS('xy_clip_per_gam_subsamp.rds')

library(biomaRt)
hsap = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 104)
gam_bm = getBM(attributes = c("external_gene_name","ensembl_gene_id", "chromosome_name","gene_biotype"), filters = "ensembl_gene_id", values = unique(xy_clip_per_gam$gene), mart = hsap)

Ygenes = subset(gam_bm, chromosome_name == 'Y')
Xgenes = subset(gam_bm, chromosome_name == 'X')
yclip = subset(xy_clip_per_gam_subsamp, gene %in% Ygenes$ensembl_gene_id)
xclip = subset(xy_clip_per_gam_subsamp, gene %in% Xgenes$ensembl_gene_id)
autoclip = subset(xy_clip_per_gam_subsamp, gene %!in% c(Xgenes$ensembl_gene_id,Ygenes$ensembl_gene_id))
length(unique(xclip$gene))
length(unique(yclip$gene))
length(unique(autoclip$gene))
mean(xclip$t)
mean(yclip$t)
mean(autoclip$t, na.rm = T)
table(xclip$t > 0, xclip$padj_all < 0.05)
table(yclip$t > 0, yclip$padj_all < 0.05)
table(autoclip$t > 0, autoclip$padj_all < 0.05)
xy_clip_per_gam_subsamp$chrom = ifelse(xy_clip_per_gam_subsamp$gene %in% Xgenes$ensembl_gene_id, 'X','auto')
xy_clip_per_gam_subsamp$chrom = ifelse(xy_clip_per_gam_subsamp$gene %in% Ygenes$ensembl_gene_id, 'Y',xy_clip_per_gam_subsamp$chrom)
table(xy_clip_per_gam_subsamp$chrom)
xy_clip_per_gam_subsamp$bias = ifelse(xy_clip_per_gam_subsamp$padj_all < 0.05 & xy_clip_per_gam_subsamp$t > 0, 'X','not')
xy_clip_per_gam_subsamp$bias = ifelse(xy_clip_per_gam_subsamp$padj_all < 0.05 & xy_clip_per_gam_subsamp$t < 0, 'Y',xy_clip_per_gam_subsamp$bias)
table(xy_clip_per_gam_subsamp$bias)
table(xy_clip_per_gam_subsamp$chrom, xy_clip_per_gam_subsamp$bias)
round(prop.table(table(xy_clip_per_gam_subsamp$chrom, xy_clip_per_gam_subsamp$bias),1),3)
chisq.test(table(xy_clip_per_gam_subsamp$chrom, xy_clip_per_gam_subsamp$bias))
m = table(xy_clip_per_gam_subsamp$chrom, xy_clip_per_gam_subsamp$bias)
m
pairwiseNominalIndependence(m,fisher = TRUE, gtest  = FALSE, chisq  = TRUE, method = "fdr", simulate.p.value = TRUE)
pairwiseNominalIndependence(t(m),fisher = TRUE, gtest  = FALSE, chisq  = TRUE, method = "fdr", simulate.p.value = TRUE)
perg = xy_clip_per_gam_subsamp %>% group_by(chrom, gene_pair,bias) %>% summarise(mean(t))
pergt = xy_clip_per_gam_subsamp %>% group_by(chrom, tissue,gene_pair,bias) %>% summarise(mean(t))

#### end ####

#### Tables S8 S10 ####

# Table S8

xy_clip_per_gam = xy_clip_per_gam[complete.cases(xy_clip_per_gam),]
length(unique(paste(xy_clip_per_gam$tissue, xy_clip_per_gam$`gene pair`)))
out = dcast(xy_clip_per_gam, gene ~ gene_pair + tissue, value.var = 'padj')
write.csv(out, file = 'xy_clip_per_gam.csv') 

# Table S10

xy_clip_per_gam_subsamp = xy_clip_per_gam_subsamp[complete.cases(xy_clip_per_gam_subsamp),]
length(unique(paste(xy_clip_per_gam_subsamp$tissue, xy_clip_per_gam_subsamp$gene_pair)))
out = dcast(xy_clip_per_gam_subsamp, gene ~ gene_pair + tissue, value.var = 'padj')
write.csv(out, file = 'xy_clip_per_gam_subsamp.csv') # Supplementary Table 10

#### end ####

#### Figure 4C ####

pl = xy_clip_per_gam_subsamp[complete.cases(xy_clip_per_gam_subsamp),]
pl = subset(pl, padj_all < 0.05)

pl2 = pl %>% group_by(gene_pair, tissue, t > 0) %>% summarise(tot = n())
pl2$tot = ifelse(pl2$`t > 0` == TRUE, pl2$tot, pl2$tot * -1)
pl2$`t > 0` = as.factor(pl2$`t > 0`)
levels(pl2$`t > 0`) = c('Y-coupled','X-coupled')
pl2$`t > 0` = relevel(pl2$`t > 0`, "X-coupled")
pl2$ymin = ifelse(pl2$`t > 0` == 'X-coupled',0,-5000)
pl2$ymax = ifelse(pl2$`t > 0` == 'X-coupled',5000,0)

length(unique(pl$gene))
length(unique(xy_clip_per_gam_subsamp$gene))

# select subgroup to plot
pl2$gene_pair = as.factor(pl2$gene_pair)
levels(pl2$gene_pair)
levels(pl2$gene_pair) = c('DDX3X/Y','EIF1AX/Y','KDM5C/D','UTX/Y',
                          'NLGN4X/Y','PCDH11X/Y','PRKX/Y','RPS4X/Y1',
                          'SOX3/SRY','TBL1X/Y','TGIF2LX/Y',
                          'TMSB4X/Y','TXLNG/Y','USP9X/Y','ZFX/Y')
pl_sub = subset(pl2, gene_pair %in% levels(pl2$gene_pair)[c(1:5,7)])
pl_sub = subset(pl2, gene_pair %in% levels(pl2$gene_pair)[c(8,12:15)])
View(pl2 %>% group_by(gene_pair) %>% summarise(mean(abs(tot))))
View(pl2 %>% group_by(tissue) %>% summarise(mean(abs(tot))))

ggplot(pl_sub, aes(x = reorder_within(tissue, -abs(tot), gene_pair), y = tot, group = tissue, fill = tissue)) +
  geom_bar(stat = 'identity') +
  scale_x_reordered() +
  facet_grid(`t > 0`~ gene_pair, scales = 'free', space = 'free_x') +
  scale_fill_manual(values = tissue.colors) +
  guides(fill='none') +
  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +
  theme_article() +
  ylab('Count') +
  theme(legend.position = 'none',
        strip.text.x = element_text(size = 10, angle = 45),
        strip.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank())

#### end ####

#### Figure 4D ####

pl2 = pl %>% group_by(gene, tissue,t > 0) %>% summarise(tot = n())
pl3 = pl2 %>% group_by(gene, `t > 0`) %>% summarise(tot = n())
#pl3$tot = as.character(pl3$tot)
pl4 = pl3 %>% group_by(tot, `t > 0`) %>% summarise(count = n())
pl4$count = ifelse(pl4$`t > 0` == TRUE, pl4$count, -1*pl4$count)
pl4$`t > 0` = as.factor(pl4$`t > 0`)
levels(pl4$`t > 0`) = c('Y-coupled','X-coupled')
pl4$`t > 0` = relevel(pl4$`t > 0`, "X-coupled")
pl4$ymin = ifelse(pl4$`t > 0` == 'X-coupled',0,-2500)
pl4$ymax = ifelse(pl4$`t > 0` == 'X-coupled',2500,0)

ggplot(pl4, aes(x = tot, y = count, fill = `t > 0`)) +
  geom_bar(stat = 'identity',  alpha = 0.8) +
  facet_grid(`t > 0`~ ., scales = 'free', space = 'free_x') +
  scale_fill_manual(values = c(tissue.colors[7], tissue.colors[14])) +
  guides(fill='none') +
  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +
  theme_article() +
  ylab('Count') +
  xlab('# Tissues') +
  theme(legend.position = 'none',
        strip.text.x = element_text(size = 10, angle = 45),
        strip.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 12))

#### end ####

#### Figure S9 ####

check = pl %>% group_by(gene, tissue) %>% summarise(n = n())
check2 = check %>% group_by(gene) %>% summarise(n = n())
spec = subset(check2, n == 1)$gene

pl2 = pl %>% group_by(gene, tissue,t > 0) %>% summarise(tot = n())
pl3 = pl2 %>% group_by(gene, `t > 0`) %>% summarise(tot = n())
pl3 = subset(pl3, gene %!in% spec) # remove tissue specific genes
pl4 = pl3 %>% group_by(tot, `t > 0`) %>% summarise(count = n())
pl4$count = ifelse(pl4$`t > 0` == TRUE, pl4$count, -1*pl4$count)
pl4$`t > 0` = as.factor(pl4$`t > 0`)
levels(pl4$`t > 0`) = c('Y-coupled','X-coupled')
pl4$`t > 0` = relevel(pl4$`t > 0`, "X-coupled")
pl4$ymin = ifelse(pl4$`t > 0` == 'X-coupled',0,-2500)
pl4$ymax = ifelse(pl4$`t > 0` == 'X-coupled',2500,0)

ggplot(pl4, aes(x = tot, y = count, fill = `t > 0`)) +
  geom_bar(stat = 'identity',  alpha = 0.8) +
  facet_grid(`t > 0`~ ., scales = 'free', space = 'free_x') +
  scale_fill_manual(values = c(tissue.colors[7], tissue.colors[14])) +
  guides(fill='none') +
  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +
  theme_article() +
  ylab('Count') +
  xlab('# Tissues') +
  theme(legend.position = 'none',
        strip.text.x = element_text(size = 10, angle = 45),
        strip.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 12))

#### end ####

#### Figure 4D ####

pl2 = pl %>% group_by(gene, gene_pair,t > 0) %>% summarise(tot = n())
pl3 = pl2 %>% group_by(gene, `t > 0`) %>% summarise(tot = n())
pl4 = pl3 %>% group_by(tot, `t > 0`) %>% summarise(count = n())
pl4$count = ifelse(pl4$`t > 0` == TRUE, pl4$count, -1*pl4$count)
pl4$`t > 0` = as.factor(pl4$`t > 0`)
levels(pl4$`t > 0`) = c('Y-coupled','X-coupled')
pl4$`t > 0` = relevel(pl4$`t > 0`, "X-coupled")
pl4$ymin = ifelse(pl4$`t > 0` == 'X-coupled',0,-2500)
pl4$ymax = ifelse(pl4$`t > 0` == 'X-coupled',2500,0)

ggplot(pl4, aes(x = tot, y = count, fill = `t > 0`)) +
  geom_bar(stat = 'identity',  alpha = 0.8) +
  facet_grid(`t > 0`~ ., scales = 'free', space = 'free_x') +
  scale_fill_manual(values = c(tissue.colors[7], tissue.colors[14])) +
  guides(fill='none') +
  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +
  theme_article() +
  ylab('Count') +
  xlab('# Pairs') +
  theme(legend.position = 'none',
        strip.text.x = element_text(size = 10, angle = 45),
        strip.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 12))

#### end ####

#### Figure S9 ####

## by max # tissues within pairs

#check = pl %>% group_by(gene, tissue) %>% summarise(n = n())
#check2 = check %>% group_by(gene) %>% summarise(n = n())
#spec = subset(check2, n == 1)$gene

pl2 = pl %>% group_by(gene, tissue,gene_pair,t > 0) %>% summarise(tot = n())
pl3 = pl2 %>% group_by(gene, gene_pair,`t > 0`) %>% summarise(tot = n())
#pl3 = subset(pl3, gene %!in% spec) # remove tissue specific genes
pl3 = subset(pl3, gene_pair %!in% c('RPS4X & RPS4Y2',
                                    'SOX3 & SRY',
                                    'TGIF2LX & TGIF2LY')) # remove tissue specific gametologues
pl4 = pl3 %>% group_by(gene,`t > 0`) %>% summarise(max = max(tot))
pl5 = pl4 %>% group_by(max, `t > 0`) %>% summarise(count = n())
pl5$count = ifelse(pl5$`t > 0` == TRUE, pl5$count, -1*pl5$count)
pl5$`t > 0` = as.factor(pl5$`t > 0`)
levels(pl5$`t > 0`) = c('Y-coupled','X-coupled')
pl5$`t > 0` = relevel(pl5$`t > 0`, "X-coupled")
pl5$ymin = ifelse(pl5$`t > 0` == 'X-coupled',0,-6000)
pl5$ymax = ifelse(pl5$`t > 0` == 'X-coupled',6000,0)

# how many genes biased in more than 1 tissue
length(unique(pl$gene))
bb = subset(pl4, max > 1)
length(unique(bb$gene))

# which pairs show highest consistent cross-tissue bias
pl3 %>% group_by(gene_pair) %>% summarise(mean(tot))

ggplot(pl5, aes(x = max, y = count, fill = `t > 0`)) +
  geom_bar(stat = 'identity',  alpha = 0.8) +
  facet_grid(`t > 0`~ ., scales = 'free', space = 'free_x') +
  scale_fill_manual(values = c(tissue.colors[7], tissue.colors[14])) +
  guides(fill='none') +
  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +
  theme_article() +
  ylab('Count') +
  xlab('# Tissues') +
  theme(legend.position = 'none',
        strip.text.x = element_text(size = 10, angle = 45),
        strip.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 12))

## by max # pairs within tissues

pl2 = pl %>% group_by(gene, tissue,gene_pair,t > 0) %>% summarise(tot = n())
pl3 = pl2 %>% group_by(gene, tissue,`t > 0`) %>% summarise(tot = n())
pl4 = pl3 %>% group_by(gene,`t > 0`) %>% summarise(max = max(tot))
pl5 = pl4 %>% group_by(max, `t > 0`) %>% summarise(count = n())
pl5$count = ifelse(pl5$`t > 0` == TRUE, pl5$count, -1*pl5$count)
pl5$`t > 0` = as.factor(pl5$`t > 0`)
levels(pl5$`t > 0`) = c('Y-coupled','X-coupled')
pl5$`t > 0` = relevel(pl5$`t > 0`, "X-coupled")
pl5$ymin = ifelse(pl5$`t > 0` == 'X-coupled',0,-6000)
pl5$ymax = ifelse(pl5$`t > 0` == 'X-coupled',6000,0)

ggplot(pl5, aes(x = max, y = count, fill = `t > 0`)) +
  geom_bar(stat = 'identity',  alpha = 0.8) +
  facet_grid(`t > 0`~ ., scales = 'free', space = 'free_x') +
  scale_fill_manual(values = c(tissue.colors[7], tissue.colors[14])) +
  guides(fill='none') +
  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +
  theme_article() +
  ylab('Count') +
  xlab('# Pairs') +
  theme(legend.position = 'none',
        strip.text.x = element_text(size = 10, angle = 45),
        strip.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 12))

#### end ####

#### Figure S10 ####

xy_clip_per_gam = readRDS('xy_clip_per_gam.rds')
xy_clip_per_gam$tissue = factor(xy_clip_per_gam$tissue)
levels(xy_clip_per_gam$tissue) = nospace.list
xy_clip_per_gam$key = paste(xy_clip_per_gam$tissue, xy_clip_per_gam$gene_pair, xy_clip_per_gam$gene, sep = "_")

out_xy = readRDS(file = 'out_xy.rds') 
out_xy$key = paste(out_xy$region, out_xy$gene.y, out_xy$gene.x, sep = "_")

comp = merge(xy_clip_per_gam, out_xy, by = 'key')

hsap = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 104)
gam_bm = getBM(attributes = c("external_gene_name","ensembl_gene_id", "chromosome_name","gene_biotype"), filters = "ensembl_gene_id", values = unique(xy_clip_per_gam$gene), mart = hsap)
colnames(comp)[5] = "ensembl_gene_id"
comp = merge(comp, gam_bm, by = "ensembl_gene_id", all.x = T)

comp$chrom = comp$chromosome_name
comp$chrom = ifelse(comp$chrom %!in% c('X','Y'), 'auto', comp$chrom)
comp$bias = ifelse(comp$diff > 0, 'X', 'Y')
comp$bias = ifelse(comp$padj_all > 0.05, 'none', comp$bias)
comp$bias = factor(comp$bias)
table(comp$bias)

k = comp %>% group_by(chrom, ensembl_gene_id, bias) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))
k = subset(k, bias != 'none')
k2 = k %>% group_by(chrom, ensembl_gene_id) %>% summarise(freq_bias = sum(freq))
table(k2$chrom)

k2 %>% group_by(chrom) %>% summarise(mean = mean(freq_bias))

ggplot(k2, aes(x = chrom, y = freq_bias)) +
  geom_boxplot() +
  theme_article()

mod = aov(freq_bias ~ chrom, data = k2)
summary(mod)
TukeyHSD(mod)

length(unique(comp$gene.x))
length(unique(subset(comp, padj_all < 0.05)$gene.x))
length(unique(subset(comp, padj_all < 0.05 & diff > 0)$gene.x))
length(unique(subset(comp, padj_all < 0.05 & diff < 0)$gene.x))

chisq.test(comp$bias, comp$chrom)

pl = subset(comp, padj_all < 0.05)
pl = subset(pl, abs(diff) < 2.5)
levels(pl$tissue) = short.listnow
View(subset(pl, chrom == 'Y'))

ggplot(pl, aes(x = tissue, y = diff, fill = chrom)) + 
  geom_boxplot(outlier.size = 1, linewidth = 0.4) +
  facet_grid(bias~., space = 'free', scales = 'free') +
  theme_article() +
  scale_fill_brewer(palette = 'Dark2') +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
        axis.title.x = element_blank())

tt = pl %>% group_by(chrom, bias, tissue) %>% summarise(mean = mean(diff))
ggplot(pl, aes(x = gene_pair, y = diff, fill = chrom)) + 
  geom_boxplot(outlier.size = 1) +
  facet_grid(bias~., space = 'free', scales = 'free') +
  theme_article() +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
        axis.title.x = element_blank())

hh = pl %>% group_by(chrom, bias, gene_pair) %>% summarise(mean = mean(diff))
pl$gene_pair = factor(pl$gene_pair)

out = data.frame()
for(i in 1:length(short.listnow)){
  for(j in 1:length(levels(pl$gene_pair))){
    d = subset(pl, tissue == short.listnow[i] & gene_pair == levels(pl$gene_pair)[j])
    if(dim(d)[1] == 0) {print('not analyzed')} else {
    mod = aov(abs(diff) ~ chrom, data = d)
    s = summary(mod)
    tu = TukeyHSD(mod)
    o = data.frame(s[[1]])
    o$tissue = short.listnow[i]
    o$pair = levels(pl$gene_pair)[j]
    out = rbind(out, o)
  }}
}
View(out)

out = data.frame()
out2 = data.frame()
for(i in 1:length(levels(pl$gene_pair))){
    d = subset(pl, gene_pair == levels(pl$gene_pair)[i])
    if(dim(d)[1] == 0) {print('not analyzed')} else {
      mod = aov(abs(diff) ~ chrom, data = d)
      s = summary(mod)
      tu = TukeyHSD(mod)
      o2 = data.frame(tu$chrom)
      o2$gene_pair = levels(pl$gene_pair)[i]
      o2$var = rownames(o2)
      rownames(o2) = NULL
      out2 = rbind(out2, o2)
      o = data.frame(s[[1]])
      o$gene_pair = levels(pl$gene_pair)[i]
      o$var = rownames(o)
      rownames(o) = NULL
      out = rbind(out, o)
    }
}
out = subset(out, var == out$var[1])
out$padj = p.adjust(out$Pr..F.)
View(subset(out, padj < 0.05))
View(subset(out2, p.adj < 0.05))


out = data.frame()
out2 = data.frame()
for(i in 1:length(short.listnow)){
  d = subset(pl, tissue == short.listnow[i])
  if(dim(d)[1] == 0) {print('not analyzed')} else {
    mod = aov(abs(diff) ~ chrom, data = d)
    s = summary(mod)
    tu = TukeyHSD(mod)
    o2 = data.frame(tu$chrom)
    o2$tissue = short.listnow[i]
    o2$var = rownames(o2)
    rownames(o2) = NULL
    out2 = rbind(out2, o2)
    o = data.frame(s[[1]])
    o$tissue = short.listnow[i]
    o$var = rownames(o)
    rownames(o) = NULL
    out = rbind(out, o)
  }
}
out = subset(out, var == out$var[1])
out$padj = p.adjust(out$Pr..F.)
View(subset(out, padj < 0.05))
View(subset(out2, p.adj < 0.05))

table(pl$bias)
table(pl$chrom)
chisq.test(pl$chrom, pl$bias)

#### end ####
