library(parallel)
library(stringr)
library(DescTools)
library(dplyr)
library(maditr)
library(spqn)
library(mosaic)

`%!in%` = Negate(`%in%`)

print('loading data')  

attrib = read.delim('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
head(attrib)
keep_samples = subset(attrib, SMAFRZE == "RNASEQ")
table(keep_samples$SMTSD)
tissue.list = levels(as.factor(attrib$SMTSD))
remove.tissues = c('Whole Blood','Bladder','Cells - Cultured fibroblasts','Cells - EBV-transformed lymphocytes','Cells - Leukemia cell line (CML)','Cervix - Ectocervix','Cervix - Endocervix','Fallopian Tube','Kidney - Medulla','Ovary','Uterus','Vagina')
listnow = tissue.list[which(tissue.list %!in% remove.tissues)]
short.listnow = c('ADs','ADv','ADR','ARa','ARc','ARt','AMY','BA24','CAU','CBH','CBL','COR','BA9','HIP','HYP','NAC','PUT','SPC','SBN','MAM','COs','COt','Ega','Em1','Em2','HEa','HEv','KID','LIV','LNG','SAL','MUS','TIB','PAN','PIT','PRO','SKn','SKe','ILM','SPL','STO','TES','THY')
short.listnow = c('Adipose(Sub)','Adipose(Vis)','Adrenal','Artery(Aor)','Artery(Cor)','Artery(Tib)','Brain(Amy)','Brain(BA24)','Brain(Caud)','Brain(Cblm1)','Brain(Cblm2)','Brain(Cort)','Brain(BA9)','Brain(Hip)','Brain(Hyp)','Brain(NAc)','Brain(Put)','Spinal(C1)','Brain(SBN)','Mammary','Colon(Sig)','Colon(Trans)','Esoph(Gas)','Esoph(Muc)','Esoph(Musc)','Heart(Atr)','Heart(Ven)','Kidney(Cor)','Liver','Lung','Salivary','Muscle','Nerve(Tib)','Pancreas','Pituitary','Prostate','Skin(NoSun)','Skin(Sun)','Ileum','Spleen','Stomach','Testes','Thyroid')

gam_list = read.csv('gametologs_in_genome.csv')
colnames(gam_list)[3] = 'ensembl_gene_id'
gam_list$pair = as.factor(gam_list$pair)

pairs = data.frame(pairs=c(1:17))
for (i in 1:length(levels(gam_list$pair))){
  genes_now = subset(gam_list, pair == levels(gam_list$pair)[i])
  pairs$pairs[i] = paste(subset(genes_now, Gametolog=='X')$common_name,"&",subset(genes_now, Gametolog=='Y')$common_name)}

arg = commandArgs(trailingOnly=TRUE)

coexp_now = readRDS(paste('norm_coexp/',arg[1],'malecoexpnorm.rds',sep="")) 
m_exp =  readRDS(paste('adj_exp/',arg[1],'adjustedexpMALES.rds',sep="")) 

genes = length(rownames(coexp_now)[rownames(coexp_now) %!in% gam_list$ensembl_gene_id])
avg_diffz = matrix(ncol=length(levels(gam_list$pair)), nrow=length(rownames(coexp_now)[rownames(coexp_now) %!in% gam_list$ensembl_gene_id]))

print('analyzing each gametologue pair')

diff_coupling_res = data.frame()
xy_ced_res = data.frame()

for (k in 1:length(levels(gam_list$pair))){
    
    genes_now = subset(gam_list, pair == levels(gam_list$pair)[k])
    print(paste(genes_now[1,1], genes_now[2,1], sep = " "))
    outnow = tryCatch(data.frame(x_coexp = coexp_now[,subset(genes_now, Gametolog=='X')$ensembl_gene_id], y_coexp = coexp_now[,subset(genes_now, Gametolog=='Y')$ensembl_gene_id]), error=function(err) NA)
    outnow = tryCatch(outnow[which(rownames(outnow) %!in% gam_list$ensembl_gene_id),], error=function(err) NA)
    outnow = tryCatch(outnow[complete.cases(outnow),], error=function(err) NA)
    
    print('calc actual difference')
    
    diffz = tryCatch(FisherZ(outnow$x_coexp) - FisherZ(outnow$y_coexp), error=function(err) NA)
    
    # combine results
    
    outnow$diff = tryCatch(diffz, error=function(err) NA)
    
    print('calc significance of differential coupling (approach 1)')
    
    # scramble x and y coupling 1000x
    outnow2 = data.frame()
    for(i in 1:1000){
      outnow2_it = tryCatch(data.frame(x_coexp = sample(outnow$x_coexp, size = length(outnow$x_coexp)), y_coexp = sample(outnow$y_coexp, size = length(outnow$y_coexp))) , error=function(err) NA)
      outnow2 = rbind(outnow2, outnow2_it)}
    
    # calculate difference among random iterations
    diffz2 = tryCatch(FisherZ(outnow2$x_coexp) - FisherZ(outnow2$y_coexp), error=function(err) NA)
    if(dim(outnow2)[1] == 0) {diffz2 = NA} else {diffz2 = diffz2}
    diffz2_out = data.frame(diff = diffz2, pair = levels(gam_list$pair)[k])
    diffz2_out$it = 1:length(diffz2_out$diff)
    
    # calculate p value per gene (compare to distribution)
    p_diffz = data.frame()
    for(i in 1:length(diffz)) {
      p_diffz[i,1] = sum(abs(diffz2) > abs(diffz[i])) / length(diffz2) }
    
    # combine results
    outnow$pval_approach1 = p_diffz$V1
    outnow$padj_approach1 = p.adjust(outnow$pval_approach1, method = 'BH')
         
    print('calc significance of differential coupling (approach 2)')
    
    # get mean coexp per gene
    
    mean_coexp = rowMeans(coexp_now, na.rm = T)
	mean_coexp = mean_coexp[order(mean_coexp)]
	mean_coexp = data.frame(mean_coexp)
	mean_coexp$decile = ntile(mean_coexp$mean_coexp, 10)  

    # get decile groups (genes in same mean coexp deciles as x and y gams)
    
    xdec = mean_coexp[which(rownames(mean_coexp) == subset(genes_now, Gametolog == 'X')$ensembl_gene_id),]$decile
  	xgroup = subset(mean_coexp, decile == xdec)
  	xgroup = rownames(xgroup[which(rownames(xgroup) %!in% genes_now$ensembl_gene_id),])
  
    ydec = mean_coexp[which(rownames(mean_coexp) == subset(genes_now, Gametolog == 'Y')$ensembl_gene_id),]$decile
  	ygroup = subset(mean_coexp, decile == ydec)
  	ygroup = rownames(ygroup[which(rownames(ygroup) %!in% genes_now$ensembl_gene_id),])
 
	# get random differences based on decile genes
	
	if(is.null(dim(outnow))) {
    
    print('no values for this gametologue pair') 
    outnow = data.frame(x_coexp = NA, y_coexp = NA, diff = NA, 
    					pval_approach1 = NA, padj_approach1 = NA, 
    					pval_approach2 = NA, padj_approach2 = NA, 
    					region = arg, pair = levels(gam_list$pair)[k])
    } else {
    
    outnow2 = data.frame(NA_col = rep(NA, dim(outnow)[1]))
        
      for(i in 1:1000){
        xsim = sample(xgroup, 1)
        ysim = sample(ygroup, 1)
        outnow2_it = tryCatch(data.frame(x_coexp = coexp_now[,xsim], y_coexp = coexp_now[,ysim]), error=function(err) NA)
        outnow2_it = tryCatch(outnow2_it[which(rownames(outnow2_it) %!in% gam_list$ensembl_gene_id),], error=function(err) NA)
        outnow2_it = tryCatch(outnow2_it[complete.cases(outnow2_it),], error=function(err) NA)
        diffz_it = tryCatch(FisherZ(outnow2_it$x_coexp) - FisherZ(outnow2_it$y_coexp), error=function(err) NA)
        outnow2[,i] = diffz_it
      }
      
    # calc p value per gene 
    
    p_diffz = data.frame()
  	for(i in 1:length(diffz)) {
    	dist_now = abs(as.numeric(as.vector(outnow2[i,])))
    	p_diffz[i,1] = sum(abs(dist_now) > abs(diffz[i])) / 1000 }
    	
    # combine results
    
    outnow$pval_approach2 = p_diffz$V1
    outnow$padj_approach2 = p.adjust(outnow$pval_approach2, method = 'BH')
    outnow$region = arg
    outnow$pair = levels(gam_list$pair)[k]
    
    }
    
    print('combining per gene diff coupling results with other gametologues')
  
  	diff_coupling_res = rbind(diff_coupling_res, outnow)

    print('calc significance of mean sex-chr-dep CED')
        
    outnow3 = mutate_all(outnow2, function(x) as.numeric(as.character(x)))
    outnow3 = as.matrix(outnow3)
	outnow3[which(outnow3 == Inf)] = NA
	outnow3[which(outnow3 == -Inf)] = NA

    mean_observed = mean(diffz, na.rm = T)
    abs_mean_observed = mean(abs(diffz), na.rm = T)

    mean_iterations = colMeans(outnow3, na.rm = TRUE)
    abs_mean_iterations = colMeans(abs(outnow3), na.rm = T)
    
    p_mean = sum(abs(mean_iterations) > abs(mean_observed))/ 1000
    p_abs_mean = sum(abs(abs_mean_iterations) > abs(abs_mean_observed))/ 1000
    
    ced_out = data.frame(diff = mean_observed, p_diff = p_mean, abs_diff = abs_mean_observed, p_abs_diff = p_abs_mean)
    ced_out$region = arg
    ced_out$pair = levels(gam_list$pair)[k]
    
    xy_ced_res = rbind(xy_ced_res, ced_out)

    }
    
	print('calc significance of differential coupling (resampling approach)')
    
    resamp = data.frame()
    for (m in 1:100){
		print(m)
		msamp = sample(colnames(m_exp), replace = TRUE)
		mnow = m_exp[,msamp]
    	mco = cor(t(mnow), use= 'everything', method = 'spearman')
    	mavg_exp = rowMeans(m_exp)
    	mnormco = normalize_correlation(mco, ave_exp=mavg_exp, ngrp=20, size_grp=1000, ref_grp=18)
   		rownames(mnormco) = colnames(mnormco) = rownames(mco)
   		
		for (k in 1:length(levels(gam_list$pair))){
			
			genes_now = subset(gam_list, pair == levels(gam_list$pair)[k])
    		print(paste(genes_now[1,1], genes_now[2,1], sep = " "))
			mout = tryCatch(data.frame(x_coexp = mnormco[,subset(genes_now, Gametolog=='X')$ensembl_gene_id], y_coexp = mnormco[,subset(genes_now, Gametolog=='Y')$ensembl_gene_id]), error=function(err) NA)
    		mout = tryCatch(mout[which(rownames(mout) %!in% gam_list$ensembl_gene_id),], error=function(err) NA)
    		mout = tryCatch(mout[complete.cases(mout),], error=function(err) NA)    
    		mdiff = tryCatch(FisherZ(mout$x_coexp) - FisherZ(mout$y_coexp), error=function(err) NA)
    		meanit = mean(mdiff, na.rm = T)
    		meanabsit = mean(abs(mdiff), na.rm = T)
    		outnow = data.frame(tissue = arg[1], pair = levels(gam_list$pair)[k], meanit = meanit, meanabsit = meanabsit, it = m)
    		resamp = rbind(resamp, outnow)
    }}
    
saveRDS(diff_coupling_res, file = paste(arg[1],'_xy_coupling_res.rds',sep=""))
saveRDS(xy_ced_res, file = paste(arg[1],'_sexchr_dep_ced_res.rds',sep=""))
saveRDS(resamp, file = paste(arg[1],'_xy_resamp.rds',sep=""))
    
print('done')
