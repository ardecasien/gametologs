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

# load disease data

do.data = read.table('human_disease_associations.tsv',
                     sep='\t',
                     quote='',
                     col.names=c('protein_id','protein_name','do_id','do_name','z_score','confidence'),
                     stringsAsFactors=FALSE)

do.def = unique(subset(do.data,select=c('do_id','do_name')))

library(biomaRt)

hsap = useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl')
hsap.info = getBM(attributes=c('ensembl_gene_id','ensembl_peptide_id','external_gene_name'),mart = hsap)
do.ensembl = subset(do.data,grepl('^ENSP[0-9]{11}',protein_id))
do.ensembl = merge(do.ensembl,hsap.info,by.x='protein_id','ensembl_peptide_id',all.x=FALSE,all.y=FALSE)
do.ensembl$ensembl_peptide_id = do.ensembl$protein_id
do.proteinname = subset(do.data,!grepl('^ENSP[0-9]{11}',protein_id))
do.proteinname = merge(do.proteinname,hsap.info,by.x='protein_name',by.y='external_gene_name',all.x=FALSE,all.y=FALSE)
do.proteinname$external_gene_name = do.proteinname$protein_name
do.all = rbind(do.ensembl[intersect(names(do.ensembl),names(do.proteinname))],do.proteinname[intersect(names(do.ensembl),names(do.proteinname))])

saveRDS(do.all, file = 'do.all.rds')

##################
## estimate X-Y coexpression per gam per tissue 
## save output
## DO + GO for ranked X-Y per pair/tissue combination and averaged within tissue
## save output
##################

library(clusterProfiler)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(DescTools)
library(reshape2)
library(ggplot2)
library(matrixStats)

diseases = data.frame()
biol = data.frame()
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
    
    avg_diffz = data.frame(avg_diffz)
    colnames(avg_diffz) = pairs[,1]
    avg_diffz$avg = rowMeans(avg_diffz, na.rm = TRUE)
    rownames(avg_diffz) = rownames(coexp_now)[rownames(coexp_now) %!in% gam_list$ensembl_gene_id]
    avg_diffz = avg_diffz[order(-avg_diffz$avg),]
    
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
    
    print('GO and DO analysis')
    
    # select gene list 
    
    #gene_list_now = avg_diff[,-19]
    gene_list_now = avg_diffz[,-19]
    
    #RUN ALL
    for(i in 1:length(colnames(gene_list_now))){
    
        
      gene_list = gene_list_now[,i]
      names(gene_list) = rownames(gene_list_now)
      gene_list = gene_list[order(gene_list, decreasing = T)]
      
      if(is.na(gene_list[1]) == TRUE) {print(paste('skipping',colnames(gene_list_now)[i],sep=" "))} else {
    
      # gene ontology
    
      gse <- gseGO(geneList=gene_list, 
                 ont ="BP", 
                 keyType = "ENSEMBL", 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Hs.eg.db, 
                 pAdjustMethod = "BH")
    
      gse_results = gse@result
      rownames(gse_results) = NULL
      gse_results$tissue = tissue.list[j]
      gse_results$pair = colnames(gene_list_now)[i]
      biol = rbind(biol, gse_results)
    
      # disease enrichment
    
      gene_list = data.frame(effect = gene_list, ensembl_gene_id = names(gene_list))
      all.region.do = merge(do.all, gene_list, by = 'ensembl_gene_id')
      all.region.do.pass = subset(all.region.do,do_id %in% names(which(table(subset(all.region.do,confidence >= 0)$do_id) >= 10)))
      all.region.do.split = split(all.region.do.pass,all.region.do.pass$do_id)
    
      all.region.do.test = do.call(rbind,mclapply(names(all.region.do.split),function(i) {
      x = all.region.do.split[[i]]
      inc.kst.test = ks.test(x$effect,subset(all.region.do.pass,do_id != i)$effect,alternative='less')
      dec.kst.test = ks.test(x$effect,subset(all.region.do.pass,do_id != i)$effect,alternative='greater')
      
      data.frame(
        do_id = unique(x$do_id),
        inc.kst.score = inc.kst.test$statistic,
        dec.kst.score = dec.kst.test$statistic,
        inc.kst.pval = inc.kst.test$p.value,
        dec.kst.pval = dec.kst.test$p.value
      )
    },mc.cores=n.cores))
    
    all.region.do.test$inc.p.adj = p.adjust(all.region.do.test$inc.kst.pval)
    all.region.do.test$dec.p.adj = p.adjust(all.region.do.test$dec.kst.pval)
    
    all.region.do.results = merge(all.region.do.test,do.def,by='do_id')
    all.region.do.results$tissue = tissue.list[j]
    all.region.do.results$pair = colnames(gene_list_now)[i]
    diseases = rbind(diseases, all.region.do.results)
    }
  }
    
    # done
    
    rm(coexp_now)
}

saveRDS(avg_diff_out, 'avg_diff_out.rds')
saveRDS(avg_diffz_out, 'avg_diffz_out.rds')
saveRDS(biol, "go_enrichment_all_tissues_z.rds")
saveRDS(diseases, "disease_enrichment_all_tissues_z.rds")
