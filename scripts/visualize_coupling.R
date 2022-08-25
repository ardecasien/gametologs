# load data

library(parallel)
library(stringr)

n.cores = detectCores()-4
`%!in%` = Negate(`%in%`)

avg_diffz_out = readRDS('avg_diffz_out.rds')

write.csv(avg_diffz_out, file = 'avg_diffz_out.csv')

# plot X-Y variation across tissues and gametologues

for(i in 1:length(tissue.list)){
  
  avg_diff_plot = subset(avg_diffz_out, Region == tissue.list[i])
  
  avg_diff_plot = dcast(avg_diff_plot, formula = gene ~ variable, value.var = 'value')
  avg_diff_plot = avg_diff_plot[order(-avg_diff_plot$avg),]
  avg_diff_plot$order = c(1:length(avg_diff_plot$gene))
  avg_diff_plot=avg_diff_plot[,-1]
  avg_diff_plot = melt(avg_diff_plot, id.vars=c('order'))
  avg_diff_plot = avg_diff_plot[complete.cases(avg_diff_plot),]
  
  print(ggplot(avg_diff_plot, aes(y=value, x=order, color = variable)) + 
          theme_classic() +
          ggtitle(tissue.list[i]) +
          geom_point(alpha=0.2) +
          geom_hline(yintercept = 0, linetype = 'dashed') +
          theme(legend.position = 'bottom', legend.box="vertical", legend.margin=margin()))
  
  avg_diff_plot2 = subset(avg_diff_plot, variable != 'avg')
  print(ggplot(avg_diff_plot2, aes(x=reorder(variable,-value), y=value)) +
          geom_boxplot() +
          ggtitle(tissue.list[i]) +
          theme_classic() +
          geom_hline(yintercept = 0, linetype='dashed') +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
}

# plot median, abs median, sd X-Y within gams/tissues

library(dcanr)
library(corrplot)
library(stringr)
library(dplyr)
library(RColorBrewer)

`%!in%` = Negate(`%in%`)

avg_diff_out_now = subset(avg_diffz_out, variable != 'avg')

avg_diff_out_now$key = paste(avg_diff_out_now$Region, avg_diff_out_now$variable, sep = "_")

# plot median median (abs) sd XY difference
plotdata = avg_diff_out_now %>% group_by(key) %>% summarise(med = median(value), abs.med = median(abs(value)), sd = sd(value))
plotdata = data.frame(plotdata)
plotdata$region = as.factor(sub("_.*", "", plotdata$key))
plotdata$pair = as.factor(sub(".*_", "", plotdata$key))

levels(plotdata$region) = short.listnow
head(plotdata)
plotdata$pair = factor(plotdata$pair, levels = c("PCDH11X & PCDH11Y","TGIF2LX & TGIF2LY","PRKX & PRKY","NLGN4X & NLGN4Y","EIF1AX & EIF1AY","TMSB4X & TMSB4Y","AMELX & AMELY","ZFX & ZFY","TBL1X & TBL1Y","TXLNG & TXLNGY","USP9X & USP9Y","DDX3X & DDX3Y","KDM6A & UTY","KDM5C & KDM5D","RPS4X & RPS4Y1","RPS4X & RPS4Y2","SOX3 & SRY"))

# heatmap median X-Y

ggplot(plotdata, aes(x = region, y = pair, fill=med)) + 
  geom_tile() +
  coord_equal() +
  theme_classic() +
  scale_fill_gradient2(high="darkgreen",mid="white",low="darkred",na.value="lightgrey", midpoint=0, name = 'Median') +
  ylab("") + xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
        axis.text.y = element_text(size = 8))

# heatmap abs median X-Y

ggplot(plotdata, aes(x = region, y = pair, fill=abs.med)) + 
  geom_tile() +
  coord_equal() +
  theme_classic() +
  scale_fill_gradient2(high="darkgreen",low="white",na.value="lightgrey", name = 'Median (abs)') +
  ylab("") + xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
        axis.text.y = element_text(size = 8))

# heatmap sd X-Y

ggplot(plotdata, aes(x = region, y = pair, fill=sd)) + 
  geom_tile() +
  coord_equal() +
  theme_classic() +
  scale_fill_gradient2(high="darkgreen",low="white",na.value="lightgrey", name = 'SD') +
  ylab("") + xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
        axis.text.y = element_text(size = 8))

# correlate median and sd X-Y

library(RColorBrewer)
palette_Dark2 <- colorRampPalette(brewer.pal(17, "Dark2"))

ggplot(plotdata, aes(x=abs.med, y=sd)) +
  geom_point(aes(color=region, shape=pair)) +
  geom_smooth(method='lm') +
  #geom_abline(slope=1,intercept=0) +
  theme_classic() +
  xlab('Median (abs) X-Y coexp') + ylab('SD X-Y coexp') +
  scale_color_manual(values=palette_Dark2(44)) +
  scale_shape_manual(values=c(8:25)) +
  theme(legend.text = element_text(size=12), legend.box = 'vertical', legend.title = element_blank(), legend.key.size = unit(0, 'lines'), legend.position = 'bottom')

ggplot(plotdata, aes(x=med, y=sd, color=region, shape=pair)) +
  geom_point(size=3) +
  theme_classic() +
  xlab('Median X-Y coexp') + ylab('SD X-Y coexp') +
  scale_color_manual(values=palette_Dark2(44)) +
  scale_shape_manual(values=c(8:25)) +
  theme(legend.text = element_text(size=12), legend.box = 'vertical', legend.title = element_blank(), legend.key.size = unit(0, 'lines'), legend.position = 'bottom')

# compare abs median X-Y to coupled coexpression

coup_coexp_plot = readRDS('coup_coexp_plot_43tissues_norm.rds')
coexp_plot2 = melt(coup_coexp_plot, id = 'region')
coexp_plot2$region = as.factor(coexp_plot2$region)
levels(coexp_plot2$region) = tissue.list
coexp_plot2$key = paste(coexp_plot2$region, coexp_plot2$variable, sep = "_")

comp = merge(coexp_plot2, plotdata, by = 'key')
comp = comp[complete.cases(comp),]
ggplot(comp, aes(x=abs.med, y=value)) +
  geom_point(size=3, aes(color=region.y, shape=pair)) +
  geom_smooth(method='lm') +
  theme_classic() +
  xlab('Median (abs) X-Y coexp') + ylab('Coupled coexp') +
  scale_color_manual(values=palette_Dark2(44)) +
  scale_shape_manual(values=c(8:25)) +
  theme(legend.text = element_text(size=12), 
        legend.box = 'vertical', 
        legend.title = element_blank(), 
        legend.key.size = unit(0, 'lines'), 
        legend.position = 'bottom')
cor.test(comp$abs.med, comp$value, method = 'spearman')

# correlate XY difference across tissues/gametologues 

check = data.frame(avg_diff_out_now %>% group_by(key) %>% summarise(avg = mean(value, na.rm = TRUE)))
rm = check[is.na(check$avg),]$key
avg_diff_out_now = subset(avg_diff_out_now, key %!in% rm)

avg_diff_out_cast = dcast(avg_diff_out_now, formula = gene ~ key, value.var = 'value')
rownames(avg_diff_out_cast) = avg_diff_out_cast$gene
avg_diff_out_cast = avg_diff_out_cast[,-1]
avg_diff_out_cast = avg_diff_out_cast
M = cor(avg_diff_out_cast, use = 'pairwise.complete.obs')
#diag(M) = NA

library(expss)
library(pheatmap)

groups.key = levels(as.factor(avg_diff_out_now$key))
groups = strsplit(groups.key, "_")
groups = data.frame(do.call(rbind, groups))
colnames(groups) = c('region','variable')

# select grouping variable
groupnow = 'region'
#groupnow = 'variable'

labs = data.frame(group = levels(as.factor(groups[,groupnow])))
labs$group2 = labs$group

# group similar tissues
# labs$group2 = sub("\\ -.*", "", labs$group)

nb.cols = length(unique(labs$group2))
mycolors <- data.frame(col = colorRampPalette(brewer.pal(8, "Set2"))(nb.cols), group2 = unique(labs$group2))
labs = merge(mycolors, labs, by = 'group2')
cols = NULL
for(i in 1:length(groups[,groupnow])){
  cols[i] = vlookup(groups[i,groupnow], dict=labs, lookup_column=3, result_column=2)}
cols
ord = corrMatOrder(M, order="hclust", hclust.method = 'ward.D2')
cols = cols[ord]

corrplot(M, order = 'hclust', hclust.method = 'ward.D2', 
         addrect = 7, tl.cex = 0.1, tl.col=cols, addgrid.col = NA)

names = groups.key
names = names[ord]
tissue = data.frame(group = groups$region)
tissue = tissue[ord,]
tissue = data.frame(tissue)
rownames(tissue) = names

# group similar tissues
# tissue$tissue = sub("\\ -.*", "", tissue$tissue)

pair = data.frame(group = groups$variable)
pair = pair[ord,]
pair = data.frame(pair)
rownames(pair) = names

macolor = colorRampPalette(c("navyblue", "white", "red"))(100)
ph = pheatmap(M, 
         color = rev(macolor), 
         clustering_method = "ward.D2", 
         clustering_distance_cols = as.dist(1 - M),
         clustering_distance_rows = as.dist(1 - M),
         fontsize_row = 0.1, 
         fontsize_col = 0.1, 
         annotation_legend = F,
         #legend=F,
         annotation_row = tissue, 
         annotation_col = pair)

## compare distances

m = as.matrix(dist(M))

tissues = unique(groups$region)
out = data.frame()
for(i in 1:length(tissues)){
  mnow = m[which(sub("\\_.*", "", rownames(m)) == tissues[i]),which(sub("\\_.*", "", colnames(m)) == tissues[i])]
  diag(mnow) = NA
  avg_corr = mean(mnow, na.rm = T)
  out[i,1] = tissues[i]
  out[i,2] = avg_corr
}
View(out)
write.csv(out, file = 'dist_tissues.csv')

pairs = unique(groups$variable)
out = data.frame()
for(i in 1:length(pairs)){
  mnow = m[which(sub('.+_(.+)', '\\1', rownames(m)) == pairs[i]),which(sub('.+_(.+)', '\\1', colnames(m)) == pairs[i])]
  diag(mnow) = NA
  avg_corr = mean(mnow, na.rm = T)
  out[i,1] = pairs[i]
  out[i,2] = avg_corr
}
View(out)
write.csv(out, file = 'dist_pairs.csv')

mpl = data.frame(m)
mpl$group = rownames(mpl)
mpl = melt(mpl)
mpl$group = gsub(" - ","...",mpl$group)
mpl$group = gsub(" & ","...",mpl$group)
mpl$group = gsub("\\(",".",mpl$group)
mpl$group = gsub("\\ ",".",mpl$group)
mpl$group = gsub(")",".",mpl$group)
mpl$group = gsub("-",".",mpl$group)
mpl$tissue = ifelse(sub("\\_.*", "", mpl$group) == sub("\\_.*", "", mpl$variable),sub("\\_.*", "", mpl$variable),"mismatch")
mpl$tissue = ifelse(mpl$group == mpl$variable, 'mismatch', mpl$tissue)
mpl$pair = ifelse(sub('.+_(.+)', '\\1', mpl$group) == sub('.+_(.+)', '\\1',  mpl$variable),sub('.+_(.+)', '\\1',  mpl$variable),"mismatch")
mpl$pair = ifelse(mpl$group == mpl$variable, 'mismatch', mpl$pair)
table(mpl$tissue)

mplnow = subset(mpl, tissue != 'mismatch')
mplnow$tissue = as.factor(mplnow$tissue)
levels(mplnow$tissue) = short.listnow
ggplot(mplnow, aes(x = reorder(tissue, -value, na.rm = T, FUN=median), y = value)) +
  geom_boxplot() +
  xlab("") + ylab("Distance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

mplnow = subset(mpl, pair != 'mismatch')
mplnow$pair = as.factor(mplnow$pair)
ggplot(mplnow, aes(x = reorder(pair, -value, na.rm = T, FUN=median), y = value)) +
  geom_boxplot() +
  xlab("") + ylab("Distance") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1,size =8))

