library(limma)
library(matrixStats)
library(tidyverse)
library(data.table)
library(readr)
library(stringr)
library(plyr)
library(dplyr)
library(Hmisc)
library(DescTools)
library(spqn)

`%!in%` = Negate(`%in%`)

# set options

attrib = read.delim('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
head(attrib)
keep_samples = subset(attrib, SMAFRZE == "RNASEQ")
tissue.list = levels(as.factor(attrib$SMTSD))

pheno = read.csv('gtex_meta_edit.csv')
colnames(pheno)[2] = 'SAMPID2'
head(pheno)
table(pheno$SEX)

library(RColorBrewer)
palette_Dark2 <- colorRampPalette(brewer.pal(17, "Dark2"))

gam_list = read.csv('gametologs_in_genome.csv')
gam_list$pair = as.factor(gam_list$pair)

remove.tissues = c('Whole Blood','Breast - Mammary Tissue','Bladder','Cells - Cultured fibroblasts','Cells - EBV-transformed lymphocytes','Cells - Leukemia cell line (CML)','Cervix - Ectocervix','Cervix - Endocervix','Fallopian Tube','Kidney - Medulla','Ovary','Prostate','Testis','Uterus','Vagina')
tissue.list = tissue.list[which(tissue.list %!in% remove.tissues)]
short.listnow = c('Adipose(Sub)','Adipose(Vis)','Adrenal','Artery(Aor)','Artery(Cor)','Artery(Tib)','Brain(Amy)','Brain(BA24)','Brain(Caud)','Brain(Cblm1)','Brain(Cblm2)','Brain(Cort)','Brain(BA9)','Brain(Hip)','Brain(Hyp)','Brain(NAc)','Brain(Put)','Spinal(C1)','Brain(SBN)','Colon(Sig)','Colon(Trans)','Esoph(Gas)','Esoph(Muc)','Esoph(Musc)','Heart(Atr)','Heart(Ven)','Kidney(Cor)','Liver','Lung','Salivary','Muscle','Nerve(Tib)','Pancreas','Pituitary','Skin(NoSun)','Skin(Sun)','Ileum','Spleen','Stomach','Thyroid')

## upload genes

gene2chrom = readRDS('gene2chrom.rds')
auto = subset(gene2chrom, chromosome_name %in% c(1:22))

## upload meta data

keep_samples = readRDS('gtex_combined_meta.rds')
keep_samples$short_tissue = keep_samples$SMTSD
keep_samples$short_tissue = str_remove_all(keep_samples$short_tissue, pattern = fixed(" "))
keep_samples$short_tissue = str_remove_all(keep_samples$short_tissue, pattern = fixed("_"))
keep_samples$short_tissue = str_remove_all(keep_samples$short_tissue, pattern = fixed("-"))
keep_samples$short_tissue = str_remove_all(keep_samples$short_tissue, pattern = fixed("("))
keep_samples$short_tissue = str_remove_all(keep_samples$short_tissue, pattern = fixed(")"))

gam_list = read.csv('gametologs_in_genome.csv')
pairs = levels(as.factor(gam_list$pair))
pairs = pairs[-14]

# get X co-expression IN MALES 
# get X co-expression IN FEMALES

m = subset(pheno, SEX == 1)
f = subset(pheno, SEX == 2)

out_all = data.frame()
out_gam = data.frame()
mean_p = data.frame()
pergam_p = data.frame()
coexp_out = data.frame()

arg = commandArgs(trailingOnly=TRUE)

print('loading male data')

  m_now = subset(keep_samples, short_tissue == arg[1] & SAMPID2 %in% m$SAMPID2)
  m_now = m_now[complete.cases(m_now[,c('AGE','SMTSISCH','SMRIN','SMNTRNRT')]),]

  if(length(m_now$SAMPID) == 0) {
    print('not enough samples') } else {

      # load adjusted expression
      	exp_M = readRDS(file = paste('/data/DNU/alex/gtex/remapped_lengthScaledTPM/adj_exp/',arg[1],'_adjusted_exp_MALES.rds',sep=""))
      exp_M = exp_M[,which(colnames(exp_M) %in% m_now$SAMPID)]

# print('subset to autosomal genes')

# exp_M = exp_M[which(rownames(exp_M) %in% c(auto$ensembl_gene_id, gam_list$ensembl_id)),]

print("getting male X")

      exp_now2 = exp_M[which(rownames(exp_M) %!in% gam_list$ensembl_id),]
      gam_exp = exp_M[which(rownames(exp_M) %in% gam_list$ensembl_id),]
      gam_exp_combo = data.frame(matrix(NA, ncol=dim(exp_now2)[2], nrow = 16))
      colnames(gam_exp_combo) = colnames(gam_exp)

      for(j in 1:length(pairs)){

        print(pairs[j])
        xnow = subset(gam_list, pair == pairs[j] & Gametolog == 'X')
        ynow = subset(gam_list, pair == pairs[j] & Gametolog == 'Y')

        xcount = sum(rownames(gam_exp) %in% xnow$ensembl_id)
        ycount = sum(rownames(gam_exp) %in% ynow$ensembl_id)

        gam_tot = data.frame(gam_exp[which(rownames(gam_exp) %in% c(xnow$ensembl_id)),])

        if(sum(xcount, ycount) == 0) {
          
          print('Both X/Y members not available') } else if (xcount == 0 && ycount == 1) {
            print('Only Y member available') } else {
              
              print('X member available')
              gam_tot2 = gam_tot[,1]
              gam_mat = tryCatch(c(gam_tot2), error=function(err) NA) 
              gam_exp_combo[j,] = gam_mat
              rownames(gam_exp_combo)[j] = paste(xnow$common_name, ynow$common_name, sep = "&") 
              } 
      }}

      gam_exp_combo = gam_exp_combo[rowSums(is.na(gam_exp_combo)) != ncol(gam_exp_combo), ]
      exp_now2 = rbind(exp_now2, gam_exp_combo)
      corM = cor(t(exp_now2), use= 'everything', method = 'spearman')

print('normalizing male correlation')

      avg_exp = rowMeans(exp_now2)
      corM_spqn = normalize_correlation(corM, ave_exp=avg_exp,
                                         ngrp=20, size_grp=1000, ref_grp=18)
      rownames(corM_spqn) = colnames(corM_spqn) = rownames(corM)
      corM_spqn[which(corM_spqn == 1)] = 0.999999999999
      rm(corM)

      corM_spqn_z = apply(corM_spqn, 1, FisherZ)
      rm(corM_spqn)
      diag(corM_spqn_z) = NA

print('loading female data')

      f_now = subset(keep_samples, short_tissue == arg[1] & SAMPID2 %in% f$SAMPID2)
      f_now = f_now[complete.cases(f_now[,c('AGE','SMTSISCH','SMRIN','SMNTRNRT')]),]

      if(length(f_now$SAMPID) == 0) {
        print('not enough samples') } else {

          # load adjusted expression
      	exp_F = readRDS(file = paste('/data/DNU/alex/gtex/remapped_lengthScaledTPM/adj_exp/',arg[1],'_adjusted_exp_FEMALES.rds',sep=""))
          exp_F = exp_F[,which(colnames(exp_F) %in% f_now$SAMPID)]

# print('subset to autosomal genes')

# exp_F = exp_F[which(rownames(exp_F) %in% c(auto$ensembl_gene_id, gam_list$ensembl_id)),]

print("getting female X")

          exp_now2 = exp_F[which(rownames(exp_F) %!in% gam_list$ensembl_id),]
          gam_exp = exp_F[which(rownames(exp_F) %in% gam_list$ensembl_id),]
          gam_exp_combo = data.frame(matrix(NA, ncol=dim(exp_now2)[2], nrow = 16))
          colnames(gam_exp_combo) = colnames(gam_exp)

          for(j in 1:length(pairs)){

            print(pairs[j])
            xnow = subset(gam_list, pair == pairs[j] & Gametolog == 'X')
            ynow = subset(gam_list, pair == pairs[j] & Gametolog == 'Y')

            gam_tot = data.frame(gam_exp[which(rownames(gam_exp) %in% c(xnow$ensembl_id)),])

            if(dim(gam_tot)[1] == 0) {

              print('X member not available') } else {

                gam_mat = tryCatch(c(gam_tot[,1]), error=function(err) NA)
                gam_exp_combo[j,] = gam_mat
                rownames(gam_exp_combo)[j] = paste(xnow$common_name, ynow$common_name, sep = "&")
              }
          }}

      gam_exp_combo = gam_exp_combo[rowSums(is.na(gam_exp_combo)) != ncol(gam_exp_combo), ]
      exp_now2 = rbind(exp_now2, gam_exp_combo)
      corF = cor(t(exp_now2), use= 'everything', method = 'spearman')

print('normalizing female correlation')

      avg_exp = rowMeans(exp_now2)
      corF_spqn = normalize_correlation(corF, ave_exp=avg_exp,
                                        ngrp=20, size_grp=1000, ref_grp=18)
      rownames(corF_spqn) = colnames(corF_spqn) = rownames(corF)
      corF_spqn[which(corF_spqn == 1)] = 0.999999999999
      rm(corF)

      corF_spqn_z = apply(corF_spqn, 1, FisherZ)
      rm(corF_spqn)
      diag(corF_spqn_z) = NA

      keep = colnames(corM_spqn_z)[which(colnames(corM_spqn_z) %in% colnames(corF_spqn_z))]
      corM_spqn_z = corM_spqn_z[which(rownames(corM_spqn_z) %in% keep), which(colnames(corM_spqn_z) %in% keep)]
      corF_spqn_z = corF_spqn_z[which(rownames(corF_spqn_z) %in% keep), which(colnames(corF_spqn_z) %in% keep)]

print('saving sex-specific gametologue co-expression fingerprints')

saveRDS(corM_spqn_z, file = paste("remapped_lengthScaledTPM/",arg,'_corM_spqn_z_X.rds',sep=""))
saveRDS(corF_spqn_z, file = paste("remapped_lengthScaledTPM/",arg,'_corF_spqn_z_X.rds',sep=""))

      gams = colnames(corM_spqn_z)[which(grepl("&", colnames(corM_spqn_z)) == TRUE)]
      corM_spqn_z_gam = corM_spqn_z[,which(colnames(corM_spqn_z) %in% gams)]
      colnames(corM_spqn_z_gam) = paste(colnames(corM_spqn_z_gam),".MALE.",arg[1], sep = "")
      corF_spqn_z_gam = corF_spqn_z[,which(colnames(corF_spqn_z) %in% gams)]
      colnames(corF_spqn_z_gam) = paste(colnames(corF_spqn_z_gam),".FEMALE.",arg[1], sep = "")
      coexp_out = merge(coexp_out, corM_spqn_z_gam, by = 'row.names',all = T)
      rownames(coexp_out) = coexp_out$Row.names
      coexp_out = coexp_out[,-1]
      coexp_out = merge(coexp_out, corF_spqn_z_gam, by = 'row.names',all = T)
      rownames(coexp_out) = coexp_out$Row.names
      coexp_out = coexp_out[,-1]

print('calculating absolute sex difference')

      cor_diff = corM_spqn_z - corF_spqn_z
      cor_diff = abs(cor_diff)
      cor_diff_mean = rowMeans(cor_diff, na.rm = T)
      gams = names(cor_diff_mean)[which(grepl("&", names(cor_diff_mean)) == TRUE)]
      gams_cor_diff_mean = cor_diff_mean[gams]
      cor_diff_mean = cor_diff_mean[which(names(cor_diff_mean) %!in% gams)]

      cor_diff_mean = data.frame(diff = cor_diff_mean, tissue = arg[1])
      cor_diff_mean$gene = rownames(cor_diff_mean)
      rownames(cor_diff_mean) = NULL
      out_all = rbind(out_all, cor_diff_mean)
      gams_cor_diff_mean = data.frame(diff = gams_cor_diff_mean, tissue = arg[1])
      gams_cor_diff_mean$gene = rownames(gams_cor_diff_mean)
      rownames(gams_cor_diff_mean) = NULL
      out_gam = rbind(out_gam, gams_cor_diff_mean)

print('estimating p values')

      corM_mean = rowMeans(corM_spqn_z, na.rm = T)
      corF_mean = rowMeans(corF_spqn_z, na.rm = T)
      cor_mean = cbind(corM_mean, corF_mean)
      cor_mean = data.frame(rowMeans(cor_mean))
      cor_mean$decile = ntile(cor_mean$rowMeans.cor_mean., 10)
      cor_mean_nogams = cor_mean[which(rownames(cor_mean) %!in% gams),]

      itmean = data.frame()
      for(b in 1:1000){
        itgenes = data.frame()
        for(q in 1:length(gams)){
        decnow = cor_mean[which(rownames(cor_mean) == gams[q]),]$decile
        genesnow = rownames(subset(cor_mean_nogams, decile == decnow))
        matchgene = sample(genesnow, 1)
        itgenes[q,1] = matchgene
        itgenes[q,2] = subset(cor_diff_mean, gene == matchgene)$diff
        }
      itmean[b,1] = mean(abs(itgenes$V2))
      }

      mean_p[1,1] = sum(abs(itmean$V1) > mean(abs(gams_cor_diff_mean$diff))) / length(itmean$V1)
      mean_p[1,2] = arg[1]

      itmean = data.frame()
      for(q in 1:length(gams)){
          itgenes = data.frame()
          for(b in 1:1000){
          decnow = cor_mean[which(rownames(cor_mean) == gams[q]),]$decile
          genesnow = rownames(subset(cor_mean_nogams, decile == decnow))
          matchgene = sample(genesnow, 1)
          itgenes[b,1] = matchgene
          itgenes[b,2] = subset(cor_diff_mean, gene == matchgene)$diff
          }

        itmean[q,1] = gams_cor_diff_mean$gene[q]
        itmean[q,2] = arg[1]
        itmean[q,3] = sum(abs(itgenes$V2) > abs(gams_cor_diff_mean[q,]$diff)) / length(itgenes$V1)
      }

      pergam_p = itmean
      
print('resampling males and females')
      
      out_gam_resamp = data.frame()
    	
    	for (m in 1:100){
		
		print(m)
		
		msamp = sample(colnames(exp_M), length(colnames(exp_M)), replace = TRUE)
		mnow = exp_M[,msamp]
    	mco = cor(t(mnow), use= 'everything', method = 'spearman')
    	mavg_exp = rowMeans(mnow)
    	mnormco = normalize_correlation(mco, ave_exp=mavg_exp, ngrp=20, size_grp=1000, ref_grp=18)
   		rownames(mnormco) = colnames(mnormco) = rownames(mco)
   		
   		fsamp = sample(colnames(exp_F), length(colnames(exp_F)), replace = TRUE)
		fnow = exp_F[,fsamp]
    	fco = cor(t(fnow), use= 'everything', method = 'spearman')
    	favg_exp = rowMeans(fnow)
    	fnormco = normalize_correlation(fco, ave_exp=favg_exp, ngrp=20, size_grp=1000, ref_grp=18)
   		rownames(fnormco) = colnames(fnormco) = rownames(fco)
   		
   		keep = colnames(mnormco)[which(colnames(mnormco) %in% colnames(fnormco))]
   		mnormco = mnormco[which(rownames(mnormco) %in% keep), which(colnames(mnormco) %in% keep)]
        fnormco = fnormco[which(rownames(fnormco) %in% keep), which(colnames(fnormco) %in% keep)]
   		cor_diff_now = mnormco - fnormco
      	cor_diff_now = abs(cor_diff_now)
   		cor_diff_mean_now = rowMeans(cor_diff_now, na.rm = T)
   		
   		xgams = subset(gam_list, Gametolog == 'X')$ensembl_id
      	gams_cor_diff_mean_now = cor_diff_mean_now[which(names(cor_diff_mean_now) %in% xgams)]
      	gams_cor_diff_mean_now = data.frame(diff = gams_cor_diff_mean_now, tissue = arg[1])
      	gams_cor_diff_mean_now$ensembl_id = rownames(gams_cor_diff_mean_now)
      	rownames(gams_cor_diff_mean_now) = NULL
      	gams_cor_diff_mean_now$it = m
      	gams_cor_diff_mean_now = merge(gams_cor_diff_mean_now, gam_list[,c('ensembl_id','common_name')], by = 'ensembl_id')
      	out_gam_resamp = rbind(out_gam_resamp, gams_cor_diff_mean_now)

		}

saveRDS(out_all, file = paste("remapped_lengthScaledTPM/",arg,'_out_all_norm_X.rds',sep=""))
saveRDS(out_gam, file = paste("remapped_lengthScaledTPM/",arg,'_out_gam_norm_X.rds',sep=""))
saveRDS(mean_p, file = paste("remapped_lengthScaledTPM/",arg,'_mean_p_norm_X.rds',sep=""))
saveRDS(pergam_p, file = paste("remapped_lengthScaledTPM/",arg,'_pergam_p_norm_X.rds',sep=""))
saveRDS(coexp_out, file = paste("remapped_lengthScaledTPM/",arg,'_coexp_out_norm_X.rds',sep=""))
saveRDS(out_gam_resamp, file = paste("remapped_lengthScaledTPM/",arg,'_out_gam_resamp_norm_X.rds',sep=""))

print('done')

