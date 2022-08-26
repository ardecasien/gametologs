library(parallel)
library(stringr)
library(DescTools)

#install.packages('maditr')
library(maditr)

`%!in%` = Negate(`%in%`)

print('loading data')  

attrib = read.delim('GTExAnalysisv8AnnotationsSampleAttributesDS.txt')
head(attrib)
keep_samples = subset(attrib, SMAFRZE == "RNASEQ")
table(keep_samples$SMTSD)
tissue.list = levels(as.factor(attrib$SMTSD))
remove.tissues = c('Whole Blood','Bladder','Cells - Cultured fibroblasts','Cells - EBV-transformed lymphocytes','Cells - Leukemia cell line (CML)','Cervix - Ectocervix','Cervix - Endocervix','Fallopian Tube','Kidney - Medulla','Ovary','Uterus','Vagina')
listnow = tissue.list[which(tissue.list %!in% remove.tissues)]
short.listnow = c('Adipose(Sub)','Adipose(Vis)','Adrenal','Artery(Aor)','Artery(Cor)','Artery(Tib)','Brain(Amy)','Brain(BA24)','Brain(Caud)','Brain(Cblm1)','Brain(Cblm2)','Brain(Cort)','Brain(BA9)','Brain(Hip)','Brain(Hyp)','Brain(NAc)','Brain(Put)','Spinal(C1)','Brain(SBN)','Mammary','Colon(Sig)','Colon(Trans)','Esoph(Gas)','Esoph(Muc)','Esoph(Musc)','Heart(Atr)','Heart(Ven)','Kidney(Cor)','Liver','Lung','Salivary','Muscle','Nerve(Tib)','Pancreas','Pituitary','Prostate','Skin(NoSun)','Skin(Sun)','Ileum','Spleen','Stomach','Testes','Thyroid')

gam_list = read.csv('gametologsingenome.csv')
colnames(gam_list)[3] = 'ensembl_gene_id'
gam_list$pair = as.factor(gam_list$pair)

pairs = data.frame(pairs=c(1:17))
for (i in 1:length(levels(gam_list$pair))){
  genes_now = subset(gam_list, pair == levels(gam_list$pair)[i])
  pairs$pairs[i] = paste(subset(genes_now, Gametolog=='X')$common_name,"&",subset(genes_now, Gametolog=='Y')$common_name)}

arg = commandArgs(trailingOnly=TRUE)

coexp_now = readRDS(paste(arg,'_male_coexp_norm.rds',sep=""))  
genes = length(rownames(coexp_now)[rownames(coexp_now) %!in% gam_list$ensembl_gene_id])
avg_diffz = matrix(ncol=length(levels(gam_list$pair)), nrow=length(rownames(coexp_now)[rownames(coexp_now) %!in% gam_list$ensembl_gene_id]))

print('analyzing each gametologue pair')

diffz_out_all = data.frame()
diffz2_out_all = data.frame()

for (k in 1:length(levels(gam_list$pair))){
    
    genes_now = subset(gam_list, pair == levels(gam_list$pair)[k])
    print(paste(genes_now[1,1], genes_now[2,1], sep = " "))
    outnow = tryCatch(data.frame(x_coexp = coexp_now[,subset(genes_now, Gametolog=='X')$ensembl_gene_id], y_coexp = coexp_now[,subset(genes_now, Gametolog=='Y')$ensembl_gene_id]), error=function(err) NA)
    outnow = tryCatch(outnow[which(rownames(outnow) %!in% gam_list$ensembl_gene_id),], error=function(err) NA)
    outnow = tryCatch(outnow[complete.cases(outnow),], error=function(err) NA)
    
    print('calc actual difference')
    diffz = tryCatch(FisherZ(outnow$x_coexp) - FisherZ(outnow$y_coexp), error=function(err) NA)
    
    print('randomize 100x')
    outnow2 = data.frame()
    for(i in 1:100){
      outnow2_it = tryCatch(data.frame(x_coexp = sample(outnow$x_coexp, size = length(outnow$x_coexp)), y_coexp = sample(outnow$y_coexp, size = length(outnow$y_coexp))) , error=function(err) NA)
      outnow2 = rbind(outnow2, outnow2_it)}
    
    print('calc random difference')
     diffz2 = tryCatch(FisherZ(outnow2$x_coexp) - FisherZ(outnow2$y_coexp), error=function(err) NA)
    if(colnames(outnow2)[1] == 'NA.') {diffz2 = NA} else {diff2 = diffz2}
    diffz2_out = data.frame(diff = diffz2, pair = levels(gam_list$pair)[k])
    diffz2_out$it = 1:length(diffz2_out$diff)
    
    print('calc p values')
    p_diffz = data.frame()
    for(i in 1:length(diffz)) {
      p_diffz[i,1] = sum(abs(diffz2) > abs(diffz[i])) / length(diffz2) }
    diffz_out = data.frame(diff = diffz, p = p_diff[,1])
    diffz_out$padj = p.adjust(diffz_out$p, method = 'BH')
    diffz_out$region = arg
    diffz_out$pair = levels(gam_list$pair)[k]
    for(i in 1:length(diffz_out$diff)){
      if(is.null(rownames(outnow)) == TRUE) {diffz_out$gene[i] = 'NA'} else {diffz_out$gene[i] = rownames(outnow)[i]}
    }    
    
    print('rbind outputs')
    diffz_out_all = rbind(diffz_out_all, diffz_out)
    diffz2_out_all = rbind(diffz2_out_all, diffz2_out)
    
}
    
print('saving outputs')

print('XY difference per gene per gam pair and p value (adjusted within each gam pair) for this region')

saveRDS(diffz_out_all, file = paste(arg,'diffz_out_all.rds'))

print('random iterations per gametolog for this region')

saveRDS(diffz2_out_all, file = paste(arg,'diffz2_out_all.rds'))

print('estimating average across gametologs')

avgz_out = data.frame()

# average across gametolog pairs for this region
    
diffz_out_all_diff = diffz_out_all[,c(1,5,6)]
diffz_out_all_avg = dcast(diffz_out_all_diff, gene ~ pair, value.var = 'diff')
diffz_out_all_avg2 = rowMeans(diffz_out_all_avg[,c(2:18)], na.rm = TRUE)
names(diffz_out_all_avg2) = diffz_out_all_avg$gene
  
diffz2_out_all_diff = diffz2_out_all
diffz2_out_all_avg = dcast(diffz2_out_all_diff, it ~ pair, value.var = 'diff')
diffz2_out_all_avg2 = rowMeans(diffz2_out_all_avg[,c(2:18)], na.rm = TRUE)
names(diffz2_out_all_avg2) = diffz2_out_all_avg$it
  
p_diffz_avg = c()
for(i in 1:length(diffz_out_all_avg2)) {
    p_diffz_avg[i] = sum(abs(diffz2_out_all_avg2) > abs(diffz_out_all_avg2[i])) / length(diffz2_out_all_avg2)}

avgz_out = data.frame(gene = names(diffz_out_all_avg2), diff = diffz_out_all_avg2, p = p_diffz_avg, tissue = arg)
avgz_out$padj = p.adjust(avgz_out$p, method = 'BH')
  
print('saving outputs')

print('XY difference per gene averaged across gam pair and p value (adjusted within each gam pair) for this region')

saveRDS(avgz_out, file = paste(arg,'XminusY_sig_avgz_all.rds')) 
