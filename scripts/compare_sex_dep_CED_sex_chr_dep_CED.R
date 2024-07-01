#### Figure 2E S5 | Table S5 ####

library(ggplot2)
library(egg)

`%!in%` = Negate(`%in%`)

out_gam = readRDS('out_gam_norm.rds')
out_gam2 = out_gam
out_gam2$X2_code = paste(out_gam2$gene, out_gam2$tissue, sep = "_")

diffc = readRDS('out_scdced.rds')
diffc = subset(diffc, region %!in% c('Testis','Prostate','BreastMammaryTissue'))
diffc$region = as.factor(diffc$region)
levels(diffc$region) = droplevels(diffc$region)
levels(diffc$region) = nospace.list
diffc = data.frame(diffc)
diffc$gene = str_replace_all(diffc$gene, pattern=" ", repl="")
diffc$X2_code = paste(diffc$gene, diffc$region, sep = "_")

out_gam2$X2_code[which(out_gam2$X2_code %!in% diffc$X2_code)] # sex-dep not in sex-chr-dep
diffc$X2_code[which(diffc$X2_code %!in% out_gam2$X2_code)] # sex-chr-dep not in sex-dep

combo = merge(out_gam2, diffc, by = 'X2_code')
combo$tissue = as.factor(combo$tissue)
levels(combo$tissue) = short.listnow
combo$abs_diff[which(combo$abs_diff == Inf)] = NA

ggplot(combo, aes(x = abs_diff, y = diff.x)) +
  geom_point(size = NA) +
  geom_smooth(method = 'lm', color = 'black') +
  ylab("Between-sex functional divergence (aCFD)") +
  xlab("Within-pair functional divergence (aCFD)") +
  geom_point(aes(x = abs_diff, y = diff.x, color = tissue), size = 3, alpha = 0.5) +
  geom_smooth(aes(x = abs_diff, y = diff.x, color = tissue), method = 'lm', 
              se = FALSE, alpha = 1, size = 1) +
  facet_wrap(~tissue, scales = 'free') +
  scale_color_manual(values = tissue.colors) +
  geom_point(size = NA) +
  geom_smooth(method = 'lm', color = 'black',  size = 1.5) +
  theme_article() +
  theme(axis.text = element_text(size = 8),
        legend.position = 'none',
        legend.title = element_blank(),
        strip.text.x = element_text(size = 8),
        legend.box = 'vertical',
        legend.key.size = unit(0, 'lines'),
        axis.title = element_text(size = 8))

ggplot(combo, aes(x = abs_diff, y = diff.x)) +
  geom_point(size = NA) +
  geom_smooth(method = 'lm', color = 'black') +
  ylab(bquote(aCEFD[MXY_FXX])) +
  xlab(bquote(aCEFD[MXvMY])) +
  geom_point(aes(x = abs_diff, y = diff.x, color = tissue), size = 3, alpha = 0.5) +
  geom_smooth(aes(x = abs_diff, y = diff.x, color = tissue), method = 'lm', 
              se = FALSE, alpha = 1, size = 1) +
  scale_color_manual(values = tissue.colors) +
  geom_point(size = NA) +
  geom_smooth(method = 'lm', color = 'black',  size = 1.5) +
  theme_article() +
  theme(axis.text = element_text(size = 14),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.box = 'vertical',
        legend.key.size = unit(0, 'lines'),
        axis.title = element_text(size = 14))

cor.test(combo$abs_diff, combo$diff.x)
mod = lmerTest::lmer(data = combo, formula = diff.x ~ abs_diff + (1|tissue))
summary(mod)
mod = aov(data = combo, formula = diff.x ~ abs_diff *tissue)
summary(mod)
mod = lm(data = combo, formula = diff.x ~ abs_diff *tissue)
summary(mod)

mod = data.frame()
for(i in 1:length(nospace.list)){
  nn = subset(combo, region == nospace.list[i])
  cc = cor.test(nn$abs_diff, nn$diff.x, method = 'pearson')
  mod[i,1] = nospace.list[i]
  mod[i,2] = cc$estimate
  mod[i,3] = cc$p.value
}
View(mod)
mod$padj = p.adjust(mod$V3)
table(mod$V3 < 0.05, mod$V2 > 0)
table(mod$padj < 0.05, mod$V2 > 0)

# Table S5

write.csv(mod, file = 'sdced-vs-schdepced-per-tissue-correlations.csv')

#### end ####

#### Figure S5 | Table S5 ####

out_gam = readRDS('out_gam_signed_norm.rds')
out_gam2 = out_gam
out_gam2$X2_code = paste(out_gam2$gene, out_gam2$tissue, sep = "_")

diffc = readRDS('out_scdced.rds')
diffc = subset(diffc, region %!in% c('Testis','Prostate','BreastMammaryTissue'))
diffc$region = as.factor(diffc$region)
levels(diffc$region) = droplevels(diffc$region)
levels(diffc$region) = nospace.list
diffc = data.frame(diffc)
diffc$gene = str_replace_all(diffc$gene, pattern=" ", repl="")
diffc$X2_code = paste(diffc$gene, diffc$region, sep = "_")

out_gam2$X2_code[which(out_gam2$X2_code %!in% diffc$X2_code)] # sex-dep not in sex-chr-dep
diffc$X2_code[which(diffc$X2_code %!in% out_gam2$X2_code)] # sex-chr-dep not in sex-dep

combo = merge(out_gam2, diffc, by = 'X2_code')
combo$tissue = as.factor(combo$tissue)
levels(combo$tissue) = short.listnow

ggplot(combo, aes(x = diff.x, y = -diff.y)) +
  geom_point(size = NA) +
  geom_smooth(method = 'lm', color = 'black') +
  ylab(bquote(sCEFD[MXY_FXX])) +
  xlab(bquote(sCFED[MXvMY])) +
  geom_point(aes(x = diff.x, y = -diff.y, color = tissue), size = 3, alpha = 0.5) +
  geom_smooth(aes(x = diff.x, y = -diff.y, color = tissue), method = 'lm', 
              se = FALSE, alpha = 1, size = 1) +
  scale_color_manual(values = tissue.colors) +
  geom_point(size = NA) +
  geom_smooth(method = 'lm', color = 'black',  size = 1.5) +
  theme_article() +
  theme(axis.text = element_text(size = 14),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.box = 'vertical',
        legend.key.size = unit(0, 'lines'),
        axis.title = element_text(size = 14))

cor.test(combo$diff.x, -combo$diff.y)
mod = lmerTest::lmer(data = combo, formula = -diff.y ~ diff.x + (1|tissue))
summary(mod)

mod = data.frame()
for(i in 1:length(nospace.list)){
  nn = subset(combo, region == nospace.list[i])
  cc = cor.test(nn$diff.x, -nn$diff.y, method = 'pearson')
  mod[i,1] = nospace.list[i]
  mod[i,2] = cc$estimate
  mod[i,3] = cc$p.value
}
View(mod)
mod$padj = p.adjust(mod$V3)
table(mod$V3 < 0.05, mod$V2 > 0)
table(mod$padj < 0.05, mod$V2 > 0)

# Table S5

write.csv(mod, file = 'sdced-vs-schdepced-per-tissue-correlations.csv')

#### end ####

