library(parallel)
library(stringr)

n.cores = detectCores()-4
`%!in%` = Negate(`%in%`)

# gametologues

gam_list = read.csv('gametologs_in_genome.csv')
colnames(gam_list)[3] = 'ensembl_gene_id'
gam_list$pair = as.factor(gam_list$pair)

pairs = data.frame(pairs=c(1:17))
for (i in 1:length(levels(gam_list$pair))){
  genes_now = subset(gam_list, pair == levels(gam_list$pair)[i])
  pairs$pairs[i] = paste(subset(genes_now, Gametolog=='X')$common_name,"&",subset(genes_now, Gametolog=='Y')$common_name)}
pairs

##################
## estimate X-Y coexpression per gam per tissue 
## save output
##################

library(clusterProfiler)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(DescTools)
library(reshape2)
library(ggplot2)
library(matrixStats)

avg_diff_out = data.frame()
avg_diffz_out = data.frame()

for (j in 1:length(tissue.list)){
  
  print(paste('working on ',tissue.list[j],sep=""))
  
  coexp_now = readRDS(paste(tissue.list[j],'_male_coexp_norm.rds',sep=""))
  genes = length(rownames(coexp_now)[rownames(coexp_now) %!in% gam_list$ensembl_gene_id])
  avg_diff = matrix(ncol=length(levels(gam_list$pair)), nrow=length(rownames(coexp_now)[rownames(coexp_now) %!in% gam_list$ensembl_gene_id]))
  avg_diffz = matrix(ncol=length(levels(gam_list$pair)), nrow=length(rownames(coexp_now)[rownames(coexp_now) %!in% gam_list$ensembl_gene_id]))
  
  print('calculating x minus y')
  
  for (k in 1:length(levels(gam_list$pair))){
    
    genes_now = subset(gam_list, pair == levels(gam_list$pair)[k])
    outnow = tryCatch(data.frame(x_coexp = coexp_now[,subset(genes_now, Gametolog=='X')$ensembl_gene_id], y_coexp = coexp_now[,subset(genes_now, Gametolog=='Y')$ensembl_gene_id]), error=function(err) NA)
    outnow = tryCatch(outnow[which(rownames(outnow) %!in% gam_list$ensembl_gene_id),], error=function(err) NA)
    outnow = tryCatch(outnow[complete.cases(outnow),], error=function(err) NA)
    diff = tryCatch(outnow$x_coexp - outnow$y_coexp, error=function(err) NA)
    
    # for normalized co-expression - some values equal 1
    if(is.na(outnow) == TRUE) {print("nope")} else {outnow$x_coexp[which(outnow$x_coexp == 1)] = 0.999999999999}
    if(is.na(outnow) == TRUE) {print("nope")} else {outnow$y_coexp[which(outnow$y_coexp == 1)] = 0.999999999999}
    
    diffz = tryCatch(FisherZ(outnow$x_coexp) - FisherZ(outnow$y_coexp), error=function(err) NA)
    avg_diff[,k] = diff
    avg_diffz[,k] = diffz
  }
  
  print('averaging across gametologs')
  
  # mean coexp difference across gametologues
  
  avg_diff = data.frame(avg_diff)
  colnames(avg_diff) = pairs[,1]
  avg_diff$avg = rowMedians(as.matrix(avg_diff), na.rm = TRUE)
  rownames(avg_diff) = rownames(coexp_now)[rownames(coexp_now) %!in% gam_list$ensembl_gene_id]
  avg_diff = avg_diff[order(-avg_diff$avg),]
  
  #avg_diff_plot = avg_diff
  #avg_diff_plot$order = c(1:length(rownames(avg_diff)))
  #avg_diff_plot = melt(avg_diff_plot, id.vars=c('order'))
  #ggplot(avg_diff_plot, aes(y=value, x=order, color = variable)) + 
  #  theme_classic() +
  #  ggtitle(tissue.list[j]) +
  #  geom_point(alpha=0.2) +
  #  theme(legend.position = 'bottom', legend.box="vertical", legend.margin=margin())
  #rm(avg_diff_plot)
  
  avg_diffz = data.frame(avg_diffz)
  colnames(avg_diffz) = pairs[,1]
  avg_diffz$avg = rowMeans(avg_diffz, na.rm = TRUE)
  rownames(avg_diffz) = rownames(coexp_now)[rownames(coexp_now) %!in% gam_list$ensembl_gene_id]
  avg_diffz = avg_diffz[order(-avg_diffz$avg),]
  
  # avg_diffz_plot = avg_diffz
  # avg_diffz_plot$order = c(1:length(rownames(avg_diffz)))
  # avg_diffz_plot = melt(avg_diffz_plot, id.vars=c('order'))
  # ggplot(avg_diffz_plot, aes(y=value, x=order, color = variable)) + 
  #   theme_classic() +
  #   ggtitle(tissue.list[j]) +
  #   geom_point(alpha=0.2) +
  #   theme(legend.position = 'bottom', legend.box="vertical", legend.margin=margin())
  # rm(avg_diffz_plot)
  
  # save difference matrices
  avg_diff$gene = rownames(avg_diff)
  avg_diff_save = melt(avg_diff, id.vars = 'gene')
  avg_diff_save$Region = tissue.list[j]
  rownames(avg_diff_save) = NULL
  avg_diff_out = rbind(avg_diff_out, avg_diff_save)
  rm(avg_diff_save)
  
  avg_diffz$gene = rownames(avg_diffz)
  avg_diffz_save = melt(avg_diffz, id.vars = 'gene')
  avg_diffz_save$Region = tissue.list[j]
  rownames(avg_diffz_save) = NULL
  avg_diffz_out = rbind(avg_diffz_out, avg_diffz_save)
  rm(avg_diffz_save)

}

saveRDS(avg_diff_out, 'avg_diff_out.rds')
saveRDS(avg_diffz_out, 'avg_diffz_out.rds')





