library(stringr)
library(ggplot2)
library(egg)
library(reshape2)

`%!in%` = Negate(`%in%`)

###################
## upload meta data
###################

keep_samples = readRDS('gtex_combined_meta.rds')

gam_list = read.csv('gametologs_in_genome.csv')
gam_list$pair = as.factor(gam_list$pair)

#############################################
## load coexp and coup coexp for gametologues
#############################################

coexp_plot = data.frame()
coup_coexp_plot = data.frame()

for (j in 1:length(tissue.list)){
  
  print(tissue.list[j])
  
  coexp_now = readRDS(paste(tissue.list[j],'_male_coexp_norm.rds',sep=""))
  coup_coexp_now = readRDS(paste(tissue.list[j],'_male_coup_coexp_norm.rds',sep=""))
  
  for (i in 1:length(levels(gam_list$pair))){
    
    genes_now = subset(gam_list, pair == levels(gam_list$pair)[i])
    coexp_plot[j,i] = tryCatch(coexp_now[subset(genes_now, Gametolog=='X')$ensembl_id,subset(genes_now, Gametolog=='Y')$ensembl_id], error=function(err) NA)
    coup_coexp_plot[j,i] = tryCatch(coup_coexp_now[subset(genes_now, Gametolog=='X')$ensembl_id,subset(genes_now, Gametolog=='Y')$ensembl_id], error=function(err) NA)
  
    }
  
  rm(coexp_now)
  rm(coup_coexp_now)

  }

pairs = data.frame(pairs=c(1:17))

for (i in 1:length(levels(gam_list$pair))){
  
  genes_now = subset(gam_list, pair == levels(gam_list$pair)[i])
  pairs$pairs[i] = paste(subset(genes_now, Gametolog=='X')$common_name,"&",subset(genes_now, Gametolog=='Y')$common_name)
  
  }

pairs

#format
rownames(coexp_plot) = rownames(coup_coexp_plot) = short.listnow
colnames(coexp_plot) = colnames(coup_coexp_plot) = pairs$pairs
coexp_plot$region = rownames(coexp_plot)
coup_coexp_plot$region = rownames(coup_coexp_plot)

# save
saveRDS(coexp_plot, file = 'coexp_plot_43tissues_norm.rds')
saveRDS(coup_coexp_plot, file = 'coup_coexp_plot_43tissues_norm.rds')

########################
## plot region variation
########################

library(ggplot2)

coexp_plot = readRDS('coexp_plot_43tissues_norm.rds')
coup_coexp_plot = readRDS('coup_coexp_plot_43tissues_norm.rds')

write.csv(coexp_plot, file = 'coexp_plot.csv')
write.csv(coup_coexp_plot, file = 'coup_coexp_plot.csv')

# chose data set
coexp_plot2 = melt(coexp_plot, id = 'region')
coexp_plot2 = melt(coup_coexp_plot, id = 'region')

# remove testes specific?
testes = c('TGIF2LX & TGIF2LY','SOX3 & SRY', 'RPS4X & RPS4Y2')
coexp_plot2 = subset(coexp_plot2, variable %!in% testes)

table(coexp_plot2$value < 0)
table(coexp_plot2$region, coexp_plot2$value < 0)
table(coexp_plot2$variable, coexp_plot2$value < 0)
table(coexp_plot2$region, coexp_plot2$variable, coexp_plot2$value > 0)

n = coexp_plot2 %>% group_by(region) %>% summarise(m = mean(value, na.rm = T))
n = coexp_plot2 %>% group_by(region) %>% summarise(m = median(value, na.rm = T))
View(n)

# plot this variable
coexp_plot2$region = factor(coexp_plot2$region, levels = short.listnow)
coexp_plot2$variable = factor(coexp_plot2$variable, levels = c("PCDH11X & PCDH11Y","TGIF2LX & TGIF2LY","PRKX & PRKY","NLGN4X & NLGN4Y","EIF1AX & EIF1AY","TMSB4X & TMSB4Y","AMELX & AMELY","ZFX & ZFY","TBL1X & TBL1Y","TXLNG & TXLNGY","USP9X & USP9Y","DDX3X & DDX3Y","KDM6A & UTY","KDM5C & KDM5D","RPS4X & RPS4Y1","RPS4X & RPS4Y2","SOX3 & SRY"))
coexp_plot2$value = as.numeric(coexp_plot2$value)

ggplot(coexp_plot2, aes(x = region, y = variable, fill = value)) +
  coord_fixed() +
  geom_tile() +
  scale_fill_gradient2(high="darkgreen",mid="white",low="darkred",na.value="lightgrey", midpoint=0, name = 'Spearman correlation') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title = element_blank())

ggplot(coexp_plot2, aes(x = reorder(region,-value,na.rm=TRUE,FUN=median), y = value)) +
  geom_boxplot(outlier.colour="#999999") + 
  scale_color_manual(values = c('grey')) +
  #theme_classic() + ylab('Co-expression') +
  theme_classic() + ylab('Coupled Co-expression') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank()) + theme(legend.position = 'none') 

ggplot(coexp_plot2, aes(x = reorder(variable,-value,na.rm=TRUE,FUN=median), y = value)) +
  geom_boxplot(outlier.colour="#999999") + 
  scale_color_manual(values = c('grey')) +
  #theme_classic() + ylab('Co-expression') +
  theme_classic() + ylab('Coupled Co-expression') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank()) + theme(legend.position = 'none') 

mod=aov(value ~ region, data = coexp_plot2)
summary(mod)
tuk=TukeyHSD(mod)
View(tuk$region)

mod=aov(value ~ variable, data = coexp_plot2)
summary(mod)
tuk=TukeyHSD(mod)
tuk = data.frame(tuk$variable)
tuk$pair = rownames(tuk)
View(tuk)

########################
# plot coexp vs coupled coexp
########################

coexp_plot2 = melt(coexp_plot, id = 'region')
coup_coexp_plot2 = melt(coup_coexp_plot, id = 'region')
coexp_plot2$region = factor(coexp_plot2$region, levels = short.listnow)
coup_coexp_plot2$region = factor(coup_coexp_plot2$region, levels = short.listnow)

plot_data = merge(coexp_plot2,coup_coexp_plot2,by=c('region','variable'))
library(RColorBrewer)
palette_Dark2 <- colorRampPalette(brewer.pal(17, "Dark2"))

plot_data$region = as.factor(plot_data$region)
plot_data = plot_data[complete.cases(plot_data),]
plot_data$value.x = as.numeric(plot_data$value.x)
plot_data$value.y = as.numeric(plot_data$value.y)

plot_data[which(plot_data$value.x < 0 & plot_data$value.y > 0),]

# correlation

cor.test(x=plot_data$value.x, y=plot_data$value.y, method='spearman')

#linear plot

ggplot(plot_data, aes(x=value.x, y=value.y, color=region, shape = variable)) + 
  geom_point(size=3) +
  theme_classic() + 
  xlab('Co-expression') + ylab('Coupled Co-expression') +
  scale_color_manual(values=palette_Dark2(44)) +
  scale_shape_manual(values=c(8:25)) +
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE,aes(group=1),colour="black") +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(legend.text = element_text(size=12), legend.box = 'vertical', legend.title = element_blank(), legend.key.size = unit(0, 'lines'), legend.position = 'bottom')

mod.lm = lm(value.y ~ value.x, data = plot_data)
summary(mod.lm)
resids = data.frame(cbind(plot_data, resid = residuals(mod.lm)))
saveRDS(resids, file='coexp_coupcoexp_residuals_linear.rds')

# asymptotic 

library(drc)
library(aomisc)

mod.as = drm(value.y ~ value.x, data = plot_data, fct = DRC.asymReg(), na.action=na.omit)
summary(mod.as)
newdata <- expand.grid(coexp=seq(-0.5, 1, length=1000))
pm <- predict(mod.as, newdata=newdata, interval="confidence")
newdata$p <- pm
R2nls(mod.as)$PseudoR2

AIC(mod.lm, mod.as)
BIC(mod.lm, mod.as)

library(egg)

ggplot(plot_data, aes(x = value.x, y = value.y, color=region, shape = variable)) + 
  geom_point(size=3) +
  #theme_classic() + 
  theme_article() +
  xlab('Co-expression') + ylab('Coupled Co-expression') +
  scale_color_manual(values=palette_Dark2(44)) +
  scale_shape_manual(values=c(8:25)) +
  geom_line(data=newdata, aes(x=coexp, y=p), inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(legend.text = element_text(size=12), legend.box = 'vertical', legend.title = element_blank(), legend.key.size = unit(0, 'lines'), legend.position = 'bottom')

resids = data.frame(cbind(plot_data, resid = residuals(mod.as)))
saveRDS(resids, file='coexp_coupcoexp_residuals_asymp.rds')

# loess

model <- loess(value.y ~ value.x, data = plot_data)
summary(model)
resids = data.frame(cbind(plot_data, resid = residuals(model)))
saveRDS(resids, file='coexp_coupcoexp_residuals_loess_notestes.rds')

xrange <- range(plot_data$value.x)
xseq <- seq(from=xrange[1], to=xrange[2], length=80)
pred <- predict(model, newdata = data.frame(value.x = xseq), se=TRUE)
y = pred$fit
ci <- pred$se.fit * qt(0.95 / 2 + .5, pred$df)
ymin = y - ci
ymax = y + ci
loess.DF <- data.frame(x = xseq, y, ymin, ymax, se = pred$se.fit)

ggplot(plot_data, aes(x = value.x, y = value.y)) + 
  geom_point(size=0) +
  geom_smooth(method = 'loess') +
  geom_point(aes(color=region, shape = variable), size=3) +
  theme_article() +
  xlab('Co-expression') + ylab('Coupled Co-expression') +
  scale_color_manual(values=palette_Dark2(44)) +
  scale_shape_manual(values=c(8:25)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(legend.text = element_text(size=12), 
        legend.box = 'vertical', 
        legend.title = element_blank(), 
        legend.key.size = unit(0, 'lines'), 
        legend.position = 'bottom')

######################
# examine residuals
######################

resids = readRDS('coexp_coupcoexp_residuals_linear.rds')
resids = readRDS('coexp_coupcoexp_residuals_asymp.rds')
resids = readRDS('coexp_coupcoexp_residuals_loess.rds')

ggplot(resids, aes(x = reorder(region,-resid,na.rm=TRUE,FUN=median), y = resid, fill = region)) +
  geom_boxplot(outlier.colour="#999999") + 
  #geom_violin() + 
  scale_fill_manual(values=palette_Dark2(44)) +
  #theme_classic() + 
  theme_article() +
  ylab('Residual') +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 11.5)) +
  theme(axis.title.x = element_blank()) + theme(legend.position = 'none')
ggplot(resids, aes(x = reorder(region,-abs(resid),na.rm=TRUE,FUN=median), y = abs(resid), fill = region)) +
  geom_boxplot(outlier.colour="#999999") + 
  #geom_violin() + 
  scale_fill_manual(values=palette_Dark2(44)) +
  theme_classic() + ylab('Residual') +
  theme(axis.text.x = element_text(size = 11.5, angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank()) + theme(legend.position = 'none')

mod = aov(resid ~ region, data = resids)
summary(mod)
tuk = TukeyHSD(mod)
tuk = data.frame(tuk$region)
tuk$comp = rownames(tuk)
View(tuk)

resids2 = resids %>% group_by(region) %>% summarise(avg_resid = mean(resid, na.rm=TRUE))
resids2 = resids %>% group_by(region) %>% summarise(avg_resid = mean(abs(resid), na.rm=TRUE))
View(resids2)

ggplot(resids, aes(x = reorder(variable,-resid,na.rm=TRUE,FUN=median), y = resid)) +
  geom_boxplot(outlier.colour="#999999") + 
  theme_classic() + ylab('Residual') +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank()) + theme(legend.position = 'none')
ggplot(resids, aes(x = reorder(variable,-abs(resid),na.rm=TRUE,FUN=median), y = abs(resid))) +
  geom_boxplot(outlier.colour="#999999") + 
  theme_classic() + ylab('Residual') +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank()) + theme(legend.position = 'none')

library(dplyr)
plot_data2 = plot_data %>% group_by(region) %>% summarise(avg_coexp = mean(value.x, na.rm=TRUE), avg_coup_coexp = mean(value.y, na.rm=TRUE))
View(plot_data2)

ggplot(plot_data2, aes(x=avg_coexp, y=avg_coup_coexp, color=region)) + 
  geom_point(size=3) +
  theme_classic() + xlim(c(0.2,1)) + ylim(c(0.2,1)) +
  xlab('Co-expression') + ylab('Coupled Co-expression') +
  geom_abline(slope=1, intercept = 0) + coord_equal() +
  scale_color_manual(values=palette_Dark2(44)) +
  theme(axis.title = element_text(size=12), axis.text = element_text(size=12)) +
  theme(legend.text = element_text(size=12), legend.box = 'vertical', legend.title = element_blank(), legend.key.size = unit(0, 'lines'), legend.position = 'right')
mod = lm(avg_coup_coexp ~ avg_coexp, data = plot_data2)
res_avg = data.frame(residuals(mod))
res_avg$region = plot_data2$region
View(res_avg)

plot_data2 = plot_data %>% group_by(variable) %>% summarise(avg_coexp = mean(value.x, na.rm=TRUE), avg_coup_coexp = mean(value.y, na.rm=TRUE))
plot_data2 = data.frame(plot_data2)
plot_data2 = plot_data2[complete.cases(plot_data2),]
ggplot(plot_data2, aes(x=avg_coexp, y=avg_coup_coexp, shape=variable)) + 
  geom_point(size=3) +
  theme_classic() + xlim(c(-0.1,1)) + ylim(c(-0.1,1)) +
  geom_abline(slope=1, intercept = 0) + coord_equal() +
  xlab('Co-expression') + ylab('Coupled Co-expression') +
  scale_shape_manual(values=c(8:25)) +
  theme(axis.title = element_text(size=12), axis.text = element_text(size=12)) +
  theme(legend.text = element_text(size=12), legend.box = 'vertical', legend.title = element_blank(), legend.key.size = unit(0, 'lines'), legend.position = 'right')




