#################
# load data
#################

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

# disease data

library(biomaRt)

hsap = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl')
hsap.info = getBM(attributes=c('ensembl_gene_id','ensembl_peptide_id','external_gene_name'),mart = hsap)

do.data = read.table('human_disease_associations.tsv',
                     sep='\t',
                     quote='',
                     col.names=c('protein_id','protein_name','do_id','do_name','z_score','confidence'),
                     stringsAsFactors=FALSE)

do.def = unique(subset(do.data,select=c('do_id','do_name')))
do.ensembl = subset(do.data,grepl('^ENSP[0-9]{11}',protein_id))
do.ensembl = merge(do.ensembl,hsap.info,by.x='protein_id','ensembl_peptide_id',all.x=FALSE,all.y=FALSE)
do.ensembl$ensembl_peptide_id = do.ensembl$protein_id
do.proteinname = subset(do.data,!grepl('^ENSP[0-9]{11}',protein_id))
do.proteinname = merge(do.proteinname,hsap.info,by.x='protein_name',by.y='external_gene_name',all.x=FALSE,all.y=FALSE)
do.proteinname$external_gene_name = do.proteinname$protein_name
do.all = rbind(do.ensembl[intersect(names(do.ensembl),names(do.proteinname))],do.proteinname[intersect(names(do.ensembl),names(do.proteinname))])

# layer and disease gene lists

kk = read.csv('all_gene_lists.csv')
colnames(kk)[1] = 'external_gene_name'
hsap.info2 = unique(hsap.info[,c('ensembl_gene_id','external_gene_name')])
kk2 = merge(kk, hsap.info2, by = 'external_gene_name')
head(kk2)

##################
## disease + GO for ranked X-Y per gam/tissue combination and averaged within tissue
## save output
##################

library(clusterProfiler)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(DescTools)
library(reshape2)
library(ggplot2)
library(matrixStats)

biol = data.frame()
diseases = data.frame()
diseases_kk = data.frame()

avg_diffz_out = readRDS('avg_diffz_out.rds')

gene_list_now = dcast(avg_diffz_out, formula = gene + Region ~ variable, value.var = 'value')

for (j in 1:length(tissue.list)){
  
  print(paste('working on ',tissue.list[j],sep=""))
  
  gene_list_now2 = subset(gene_list_now, Region == tissue.list[j])

  for(i in 1:(length(colnames(gene_list_now2))-2)){
    
    print(colnames(gene_list_now2)[i+2])
    gene_list = gene_list_now2[,i+2]
    names(gene_list) = gene_list_now2$gene
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
      gse_results$pair = colnames(gene_list_now2)[i+2]
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
      
      all.region.do.test$inc.p.adj = p.adjust(all.region.do.test$inc.kst.pval, method = 'BH')
      all.region.do.test$dec.p.adj = p.adjust(all.region.do.test$dec.kst.pval, method = 'BH')
      
      all.region.do.results = merge(all.region.do.test,do.def,by='do_id')
      all.region.do.results$tissue = tissue.list[j]
      all.region.do.results$pair = colnames(gene_list_now2)[i+2]
      diseases = rbind(diseases, all.region.do.results)
      
      ## layer and disease enrichments
      
      kk_out = data.frame()
      for(m in 1:(length(colnames(kk2))-2)){
        gene_list_now = subset(gene_list, ensembl_gene_id %in% kk2[which(kk2[,(m+1)] == 'True'),]$ensembl_gene_id)
        inc.kst.test = ks.test(gene_list_now$effect,gene_list$effect,alternative='less')
        dec.kst.test = ks.test(gene_list_now$effect,gene_list$effect,alternative='greater')
        kk_out[m,1] = colnames(kk2)[m+1] 
        kk_out[m,2] = inc.kst.test$statistic
        kk_out[m,3] = inc.kst.test$p.value
        kk_out[m,4] = dec.kst.test$statistic
        kk_out[m,5] = dec.kst.test$p.value
      }
      
      kk_out$tissue = tissue.list[j]
      kk_out$pair = colnames(gene_list_now2)[i+2]
      colnames(kk_out) = c('list','higherX_D','higherX_P','higherY_D','higherY_P','tissue','pair')
      kk_out$higherX_Padj = p.adjust(kk_out$higherX_P, method = 'BH')
      kk_out$higherY_Padj = p.adjust(kk_out$higherY_P, method = 'BH')
      
      diseases_kk = rbind(diseases_kk, kk_out)
    }
  }
}

saveRDS(biol, "go_enrichment_all_tissues_z.rds")
saveRDS(diseases, "disease_enrichment_all_tissues_z.rds")
saveRDS(diseases_kk, "disease_layer_enrichment_all_tissues_z.rds")
