library(expss)
library(dplyr)
library(ggplot2)
library(upstartr)

#### load data ####

## gametologs

gam_list = read.csv('gametologs_in_genome.csv')
gam_list$pair = as.factor(gam_list$pair)
gams = data.frame()
for(i in 1:length(levels(gam_list$pair))){
  gams[i,1] = levels(gam_list$pair)[i]
  xgene = subset(gam_list, pair == levels(gam_list$pair)[i] & Gametolog == 'X')$common_name
  ygene = subset(gam_list, pair == levels(gam_list$pair)[i] & Gametolog == 'Y')$common_name
  xid = subset(gam_list, pair == levels(gam_list$pair)[i] & Gametolog == 'X')$ensembl_id
  yid = subset(gam_list, pair == levels(gam_list$pair)[i] & Gametolog == 'Y')$ensembl_id
  gams[i,2] = paste(xgene, ygene, sep = " & ")
  gams[i,3] = xid
  gams[i,4] = yid
}
colnames(gams) = c('pair','gene','x_id','y_id')

## load differential coupling data

out_xy = readRDS(file = 'out_xy.rds') # per gene per tissue/pair # signed diff coupling & p-value
colnames(out_xy)[11] = 'gene'
out_xy = merge(out_xy, gams, by = 'gene', all.x = T)
out_xy = out_xy[complete.cases(out_xy),]

## load adjusted gametologue expression data

gam_exp = data.frame()

for (i in 1:length(tissue.list)){
  
  exp = readRDS(paste(nospace.list[i],"_adjusted_exp_MALES.rds",sep=""))
  exp = exp[which(rownames(exp) %in% gam_list$ensembl_id),]
  meanexp = data.frame(rowMeans(exp))
  meanexp$ensembl_id = rownames(meanexp)
  rownames(meanexp) = NULL
  meanexp$region = nospace.list[i]
  gam_exp = rbind(gam_exp, meanexp)
  
}

gam_exp$meanexp = gam_exp$rowMeans.exp. + abs(min(gam_exp$rowMeans.exp.)) + 0.1
table(gam_exp$region)

#### end ####

#### expression weighted avg asymmetric coupling ####

gam_exp_merge = gam_exp[,c(2:4)]
colnames(gam_exp_merge)[1] = 'x_id'
out_xy = merge(out_xy, gam_exp_merge, by = c('x_id','region'), all.x = T)
colnames(out_xy)[15] = 'x_exp'
colnames(gam_exp_merge)[1] = 'y_id'
out_xy = merge(out_xy, gam_exp_merge, by = c('y_id','region'), all.x = T)
colnames(out_xy)[16] = 'y_exp'

out_xy$x_coexp_weight = out_xy$x_coexp * out_xy$x_exp
out_xy$y_coexp_weight = out_xy$y_coexp * out_xy$y_exp

out_xy_weight = out_xy %>% group_by(gene.x, region) %>% 
                summarise(x_coexp_weight = sum(x_coexp_weight),
                          x_coexp_denom = sum(x_exp),
                          y_coexp_weight = sum(y_coexp_weight),
                          y_coexp_denom = sum(y_exp))
out_xy_weight$x_coexp_new = out_xy_weight$x_coexp_weight / out_xy_weight$x_coexp_denom
out_xy_weight$y_coexp_new = out_xy_weight$y_coexp_weight / out_xy_weight$y_coexp_denom
out_xy_weight$diff_xy_new = out_xy_weight$x_coexp_new - out_xy_weight$y_coexp_new

saveRDS(out_xy_weight, file = 'out_xy_weight.rds')

# check similarity to non weighted averages

out_xy_nonweight = out_xy %>% group_by(gene.x, region) %>% 
                summarise(mean_diff = mean(diff))
comp = merge(out_xy_weight, out_xy_nonweight, by = c('gene.x','region'))
cor.test(comp$mean_diff, comp$diff_xy_new, method = 'spearman')
ggplot(comp, aes(x = mean_diff, y = diff_xy_new)) + 
  geom_point() + 
  geom_smooth(method = 'lm') + 
  theme_classic()
# some oddly high values likely from original subtraction of fisher Z transformed x and y coexp

#### end ####

#### compare with mean expression ####

out_xy_weight = readRDS('out_xy_weight.rds')
out = data.frame()

for (i in 1:length(tissue.list)){
  
  print(nospace.list[i])
  xy = subset(out_xy_weight, region == nospace.list[i])
  exp = readRDS(paste(tissue.list[i],"_adjusted_exp_MALES.rds",sep=""))
  meanexp = data.frame(rowMeans(exp))
  meanexp$gene.x = rownames(meanexp)
  rownames(meanexp) = NULL
  meanexp$region = nospace.list[i]
  meanexp = merge(meanexp, xy[,c('gene.x','diff_xy_new')], by = 'gene.x')
  c1 = cor.test(meanexp$rowMeans.exp., meanexp$diff_xy_new, method = 'spearman')
  plot(meanexp$rowMeans.exp., meanexp$diff_xy_new)
  c2 = cor.test(meanexp$rowMeans.exp., abs(meanexp$diff_xy_new), method = 'spearman')
  plot(meanexp$rowMeans.exp., abs(meanexp$diff_xy_new))
  out[i,1] = nospace.list[i]
  out[i,2] = c1$estimate
  out[i,3] = c1$p.value
  out[i,4] = c2$estimate
  out[i,5] = c2$p.value
  }

colnames(out) = c('tissue','rho','p','rho_abs','p_abs')
View(out)

#### end ####
