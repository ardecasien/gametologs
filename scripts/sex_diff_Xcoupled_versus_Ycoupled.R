library(ggplot2)
library(egg)

#### load data ####

## load CLIP sig diff summed X-Y coupling (per gene per tissue)

xy_clip = readRDS(file = 'xy_clip.rds')
xy_clip$tissue = as.factor(xy_clip$tissue)
xy_clip$padj_all = p.adjust(xy_clip$p, method = 'BH')

#### end ####

#### calculate sex diffs ####

out = data.frame()

for (i in 1:length(tissue.list)){
  
  print(tissue.list[i])
  
  Ycoup_genes = subset(xy_clip, tissue == tissue.list[i] & padj_all < 0.05 & t < 0)$gene
  Xcoup_genes = subset(xy_clip, tissue == tissue.list[i] & padj_all < 0.05 & t > 0)$gene
  
  fcoexp = readRDS(paste(nospace.list[i],'_corF_spqn_z.rds',sep=""))
  diag(fcoexp) = NA
  
  mcoexp = readRDS(paste(nospace.list[i],'_corM_spqn_z.rds',sep=""))
  diag(mcoexp) = NA
  
  keep = rownames(fcoexp)[which(rownames(fcoexp) %in% rownames(mcoexp))]
  fcoexp = fcoexp[which(rownames(fcoexp) %in% keep),which(colnames(fcoexp) %in% keep)]
  mcoexp = mcoexp[which(rownames(mcoexp) %in% keep),which(colnames(mcoexp) %in% keep)]
  
  fcoexp_XvYcoupled = fcoexp[which(rownames(fcoexp) %in% Xcoup_genes), which(colnames(fcoexp) %in% Ycoup_genes)]
  mcoexp_XvYcoupled = mcoexp[which(rownames(mcoexp) %in% Xcoup_genes), which(colnames(mcoexp) %in% Ycoup_genes)]
  
  fcoexp_XvYcoupled_avg_abs = mean(abs(fcoexp_XvYcoupled))
  mcoexp_XvYcoupled_avg_abs = mean(abs(mcoexp_XvYcoupled))
  print(fcoexp_XvYcoupled_avg_abs)
  print(mcoexp_XvYcoupled_avg_abs)
  
  fcoexp_scaled = scale(fcoexp)
  mcoexp_scaled = scale(mcoexp)
  
  fcoexp_scaled_XvYcoupled = fcoexp_scaled[which(rownames(fcoexp_scaled) %in% Xcoup_genes), which(colnames(fcoexp_scaled) %in% Ycoup_genes)]
  mcoexp_scaled_XvYcoupled = mcoexp_scaled[which(rownames(mcoexp_scaled) %in% Xcoup_genes), which(colnames(mcoexp_scaled) %in% Ycoup_genes)]
  
  fcoexp_XvYcoupled_avg = mean(fcoexp_scaled_XvYcoupled)
  mcoexp_XvYcoupled_avg = mean(mcoexp_scaled_XvYcoupled)
  print(fcoexp_XvYcoupled_avg)
  print(mcoexp_XvYcoupled_avg)
  
  out[i,1] = short.listnow[i]
  out[i,2] = length(Ycoup_genes)
  out[i,3] = length(Xcoup_genes)
  out[i,4] = fcoexp_XvYcoupled_avg_abs
  out[i,5] = mcoexp_XvYcoupled_avg_abs
  out[i,6] = fcoexp_XvYcoupled_avg
  out[i,7] = mcoexp_XvYcoupled_avg
  
  rm(fcoexp)
  rm(fcoexp_scaled)
  rm(mcoexp)
  rm(mcoexp_scaled)
  
  }

colnames(out) = c('tissue','NYcoupled','NXcoupled','f_XvY_abs','m_XvY_abs','f_XvY','m_XvY')
out$m_v_f = out$m_XvY - out$f_XvY
out$m_v_f_abs = out$m_XvY_abs - out$f_XvY_abs
View(out)  

saveRDS(out, file = 'XvY_coupled_sex_specific_coexpression.rds')

#### end ####

#### Figure 5c ####

out = readRDS('XvY_coupled_sex_specific_coexpression.rds')

out_pl = subset(out, NYcoupled > 10 & NXcoupled > 10)
length(unique(out_pl$tissue))
mean(out_pl$m_v_f)
median(out_pl$m_v_f)
min(out_pl$m_v_f)
max(out_pl$m_v_f)

ggplot(out_pl, aes(x = m_v_f)) +
  geom_histogram(data=subset(out_pl,m_v_f > 0),fill = "#C9C5BA", alpha = 0.6, color = 'black') +
  geom_histogram(data=subset(out_pl,m_v_f < 0),fill = "#407076", alpha = 0.6, color = 'black') +
  ylab("") +
  xlab('Sex differences in X vs Y coupled gene co-exp') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4.5)) +
  theme_article() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18))

out_pl$code = ifelse(out_pl$m_v_f > 0, "#C9C5BA", "#407076")
mm = data.frame(tissue = out_pl$tissue, value = out_pl$m_XvY, sex = 'Males')
ff = data.frame(tissue = out_pl$tissue, value = out_pl$f_XvY, sex = 'Females')
pl2 = rbind(mm, ff)

ggplot(pl2, aes(y=reorder(tissue, value), x=value, color = sex, group = tissue)) + 
  geom_point()+
  geom_line(color = 'grey') +
  theme_article() +
  scale_color_manual(values = c(tissue.colors[c(30,1)])) +
  xlab("Mean Scaled Co-expression \n for X-Coupled & Y-Coupled Genes") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(axis.text.y = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.15))

t.test(value ~ sex, data = pl2, paired = T)

#### end ####
