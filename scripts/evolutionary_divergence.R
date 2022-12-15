setwd("~/Desktop/Papers/GTEX/v8_remapped")

`%!in%` = Negate(`%in%`)

#######################
# load gametologue info
#######################

gam_list = read.csv('gametologs_in_genome.csv')
colnames(gam_list)[3] = 'ensembl_gene_id'

# strata from Pandey et al 2013

library(biomaRt)
hsap = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
gam_bm = getBM(attributes = c("external_gene_name","ensembl_gene_id", "ensembl_transcript_id_version", "chromosome_name","transcript_biotype","transcript_gencode_basic"), filters = "ensembl_gene_id", values = gam_list$ensembl_gene_id, mart = hsap)
gam_list = merge(gam_list, gam_bm, by = 'ensembl_gene_id')
head(gam_list)
table(gam_list$transcript_gencode_basic, gam_list$transcript_biotype)

library(ape)
library(Biostrings)

gam_list_unique = unique(gam_list[,c(1,12,4,14)])
gam_list_unique_new = data.frame(pair = unique(gam_list_unique$pair))
for(i in 1:length(levels(as.factor(gam_list_unique_new$pair)))){
  gam_list_unique_new$x_gene[i] = subset(gam_list_unique, pair == levels(as.factor(gam_list_unique_new$pair))[i] & chromosome_name == 'X')$external_gene_name
  gam_list_unique_new$x_gene_id[i] = subset(gam_list_unique, pair == levels(as.factor(gam_list_unique_new$pair))[i] & chromosome_name == 'X')$ensembl_gene_id
  gam_list_unique_new$y_gene[i] = subset(gam_list_unique, pair == levels(as.factor(gam_list_unique_new$pair))[i] & chromosome_name == 'Y')$external_gene_name
  gam_list_unique_new$y_gene_id[i] = subset(gam_list_unique, pair == levels(as.factor(gam_list_unique_new$pair))[i] & chromosome_name == 'Y')$ensembl_gene_id
}
gam_list_unique_new

#####################
# promoter similarity
#####################

# get promoter regions and align

for (i in 1:length(levels(as.factor(gam_list_unique_new$pair)))){
  
  pairnow = gam_list_unique_new[i,]
  xseq = getSequence(id = pairnow$x_gene_id, type="ensembl_gene_id",seqType="gene_flank",upstream=2000,mart=hsap) 
  yseq = getSequence(id = pairnow$y_gene_id, type="ensembl_gene_id",seqType="gene_flank",upstream=2000,mart=hsap) 
  gam_list_unique_new$seq_X[i] = xseq$gene_flank
  gam_list_unique_new$seq_Y[i] = yseq$gene_flank
  
  align = pairwiseAlignment(xseq$gene_flank, yseq$gene_flank, type = "global")
  seq <- c(alignedPattern(align), alignedSubject(align))
  
  gam_list_unique_new$align1[i] = as.character(seq[1])
  gam_list_unique_new$align2[i] = as.character(seq[2])
  gam_list_unique_new$align_len[i] = nchar(align)
  indel = nindel(align)
  gam_list_unique_new$nindel[i] = indel@insertion[2] + indel@deletion[2]
  gam_list_unique_new$nmatch[i] = nmatch(align)
  gam_list_unique_new$nmismatch[i] = nmismatch(align)
  gam_list_unique_new$score[i] = align@score
}

# estimate promoter similarity measures

library(reshape2)

gam_list_unique_new$prop_match_overall = gam_list_unique_new$nmatch / gam_list_unique_new$align_len
gam_list_unique_new$prop_match_nongapped = gam_list_unique_new$nmatch / (gam_list_unique_new$nmatch + gam_list_unique_new$nmismatch)
gam_list_unique_new$pair_genes = paste(gam_list_unique_new$x_gene, gam_list_unique_new$y_gene, sep = " & ")
gam_list_unique_new_plot = melt(gam_list_unique_new[,c("pair_genes","prop_match_overall","prop_match_nongapped")], id.vars = 1)

# plot promoter similarity measures

library(ggplot2)

ggplot(gam_list_unique_new_plot, aes(x = reorder(pair_genes,-value), y = value, fill = variable)) +
  geom_bar(stat = "identity",position = "dodge") + theme_classic() +
  ylab('Promoter Similarity') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9"), name=NULL, breaks=c("prop_match_overall", "prop_match_nongapped"), labels=c("overall", "non-gapped"))

############################
# DNA and protein similarity
############################

library(ggplot2)

# get DNA and protein similarity from Skaletsy et al. 2003

skalet = read.csv('Skaletsky_comparison.csv', check.names = F)
skalet = subset(skalet, Gametologue == "Y")
head(skalet)
ggplot(skalet, aes(x=`Protein divergence perc`, y = `DNA divergence perc`)) + 
  geom_point(size=3) +
  theme_classic() + 
  geom_smooth(method='lm') +
  xlab('Protein Divergence') + ylab('DNA divergence') +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18))
cor.test(skalet$`Protein divergence perc`, skalet$`DNA divergence perc`, method='pearson')

skalet$DNA_sim = (100 - skalet$`DNA divergence perc`)/100
skalet$prot_sim = (100 - skalet$`Protein divergence perc`)/100
skalet_plot = melt(skalet[,c("pair_genes","DNA_sim","prot_sim")], id.vars = 1)
skalet_plot = skalet_plot[which(skalet_plot$pair_genes %in% gam_list_unique_new_plot$pair_genes),]

# plot similarity measures

combo_plot = rbind(gam_list_unique_new_plot,skalet_plot)

ggplot(combo_plot, aes(x = reorder(pair_genes,-value), y = value, fill = variable)) +
  geom_bar(stat = "identity",position = "dodge") + theme_classic() +
  ylab('Similarity') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9","#999999","darkolivegreen4"), 
                    name=NULL, 
                    breaks=c("prop_match_overall", "prop_match_nongapped","DNA_sim","prot_sim"), 
                    labels=c("promoter (overall)", "promoter (non-gapped)", "DNA sequence", "protein sequence"))

## compare similarity measures

library(corrplot)
corp = dcast(data = combo_plot,formula = pair_genes~variable,fun.aggregate = sum,value.var = "value")
colnames(corp) = c('pair_genes','prom-all','prom-ng','DNA-sim','pro-sim')
m = cor(corp[,c(2:5)])
corrplot(m, addCoef.col = 'white', 
         order = 'AOE', 
         cl.pos ='b', 
         tl.pos = 'd',
         #diag = FALSE, 
         tl.col = "black")

####################
# isoform similarity
####################

gam_list2 = read.csv('gametologs_in_genome.csv')
colnames(gam_list2)[3] = 'ensembl_gene_id'

# get isoform sequences

nucl = getSequence(id = gam_list$ensembl_transcript_id_version, type = 'ensembl_transcript_id_version', seqType = 'cdna', mart = hsap)
seq = getSequence(id = gam_list$ensembl_transcript_id_version, type = 'ensembl_transcript_id_version', seqType = 'coding', mart = hsap)
pep = getSequence(id = gam_list$ensembl_transcript_id_version, type = 'ensembl_transcript_id_version', seqType = 'peptide', mart = hsap)

gam_list = merge(gam_list, nucl, by = 'ensembl_transcript_id_version')
gam_list = merge(gam_list, seq, by = 'ensembl_transcript_id_version')
gam_list = merge(gam_list, pep, by = 'ensembl_transcript_id_version')

# remove bad sequences

# the ones without an asterisk don't have a stop codon because they are 3' incomplete
# the coding sequences that don't start with ATG are 5' incomplete
# only include complete sequences and coding sequences for isoform similarity

protein_coding = subset(gam_list, transcript_biotype == 'protein_coding')

library(stringr)

# incomplete protein coding transcripts
gam_list_bad_5 = subset(protein_coding, str_sub(protein_coding$coding,1,3) != "ATG")$ensembl_transcript_id_version
gam_list_bad_3 = subset(protein_coding, str_sub(protein_coding$peptide,-1) != "*")$ensembl_transcript_id_version
gam_list_bad = unique(c(gam_list_bad_5, gam_list_bad_3)) #153 (all) #24 (protein coding)

# good transcripts
gam_list2 = protein_coding[which(protein_coding$ensembl_transcript_id_version %!in% gam_list_bad),] # 156
table(gam_list2$external_gene_name, gam_list2$pair) 
check = data.frame(pair = c(1:17), count = c(1:17))
for(i in 1:length(levels(as.factor(gam_list2$pair)))){
  check$pair[i] = levels(as.factor(gam_list2$pair))[i]
  check$count[i] = length(unique(subset(gam_list2, pair == levels(as.factor(gam_list2$pair))[i])$external_gene_name))}
good_pairs = subset(check, count == 2)$pair 
gam_list2 = subset(gam_list2, pair %in% good_pairs)

gam_list2 = gam_list2[!duplicated(gam_list2[,c('coding','peptide')]),]
table(gam_list2$external_gene_name)

# plot the # of isoforms across gametologues

library(ggplot2)
library(RColorBrewer)
library(forcats)

sum = data.frame(table(gam_list2$external_gene_name))
colnames(sum) = c('external_gene_name','count')
sum = merge(sum, gam_list_unique, by = 'external_gene_name')
sum = sum[order(sum$pair),]
sum$pair2 = as.factor(sum$pair)
nb.cols <- length(levels(sum$pair2))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
ggplot(sum, aes(x=reorder(external_gene_name,pair), y = count, fill = pair2)) + theme_classic() + 
  geom_bar(stat='identity') + scale_fill_manual(values=mycolors) + coord_flip() +
  theme(legend.position = 'none', axis.title.y = element_blank())

# plot the # of isoforms vs the mean sequence length across gametologues

library(dplyr)

for(i in 1:length(gam_list2$pep)){
  gam_list2$length[i] = nchar(gam_list2$pep[i])}
gam_comp = data.frame(gam_list2 %>% group_by(external_gene_name) %>% summarise(n = n(), len = mean(length)))
colnames(gam_comp) = c('external_gene_name','count','mean AA length')
gam_comp = merge(gam_comp, gam_list_unique, by = 'external_gene_name')
gam_comp = gam_comp[order(gam_comp$pair),]
gam_comp$pair2 = as.factor(gam_comp$pair)
nb.cols <- length(levels(gam_comp$pair2))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
cor.test(x = gam_comp$`mean AA length`, y = gam_comp$count)
ggplot(gam_comp, aes(x = `mean AA length`, y = count, color = pair2)) + 
  geom_point(size=4, alpha=0.8) + theme_classic() + scale_color_manual(values=mycolors) +
  theme(legend.position = 'none') 

# get isoforms and align all pairs

scores_mat = matrix(0, nrow = length(gam_list2$ensembl_transcript_id_version), ncol = length(gam_list2$ensembl_transcript_id_version))
rownames(scores_mat) = colnames(scores_mat) = gam_list2$ensembl_transcript_id_version

for(i in 1:length(colnames(scores_mat))){
  for(j in 1:length(rownames(scores_mat))){
    s1 = AAString(c(subset(gam_list2, ensembl_transcript_id_version == colnames(scores_mat)[i])$peptide))
    s2 = AAString(c(subset(gam_list2, ensembl_transcript_id_version == rownames(scores_mat)[j])$peptide))
    align = try(pairwiseAlignment(s1, s2, substitutionMatrix = 'BLOSUM50', scoreOnly = FALSE, type = "global"))
    #align = try(pairwiseAlignment(s1, s2, substitutionMatrix = 'PAM30', scoreOnly = FALSE, type = "global"))
    sc = try(align@score)
    len = nchar(align)
    if("try-error" %in% class(align)) sc=len=NA
    scores_mat[i,j] = sc/len
  }
}

diag(scores_mat) = NA
scores_df = melt(replace(scores_mat, lower.tri(scores_mat, TRUE), NA), na.rm = TRUE)
tnames = gam_list2[c(13,1)]
colnames(tnames) = c('external_gene_name','Var1')
scores_df = merge(scores_df, tnames, by = 'Var1')
colnames(tnames) = c('external_gene_name','Var2')
scores_df = merge(scores_df, tnames, by = 'Var2')
scores_df = scores_df[,c(2,1,4,5,3)]
colnames(scores_df) = c('t1.id', 't2.id', 't1.name', 't2.name', 'score')
head(scores_df)

saveRDS(scores_df, file = 'scores_df_BLOSUM50.rds')
#saveRDS(scores_df, file = 'scores_df_PAM30.rds')

scores_df = readRDS('scores_df_BLOSUM50.rds')
#scores_df = readRDS('scores_df_PAM30.rds')

# plot isoform similarity scores across all pairs of isoforms

matches = c(paste(gam_list_unique_new$x_gene, gam_list_unique_new$y_gene),paste(gam_list_unique_new$y_gene, gam_list_unique_new$x_gene))
matches
scores_df$key = ifelse(paste(scores_df$t1.name,scores_df$t2.name) %in% matches, 'pair','not pair') 
scores_df$key = ifelse(scores_df$t1.name == scores_df$t2.name, 'self', scores_df$key)
table(scores_df$key)

ggplot(scores_df,aes(x=score)) + theme_classic() + 
  xlab("Score") + ylab("Count") +
  geom_histogram(data=subset(scores_df,key == 'self'),fill = "red", alpha = 0.2, bins = 60) +
  geom_histogram(data=subset(scores_df,key == 'pair'),fill = "blue", alpha = 0.2, bins = 60) +
  geom_histogram(data=subset(scores_df,key == 'not pair'),fill = "green", alpha = 0.2, bins = 60) +
  theme(axis.title = element_text(size=18), axis.text = element_text(size=18))

scores_df %>% group_by(key) %>% summarise(avg = mean(score))

# get isoform similarity scores within gametologues

library(dplyr)
head(scores_df)
pair_sim = subset(scores_df, key == 'pair')

xgenes = subset(gam_list_unique, chromosome_name == "X")$external_gene_name
ygenes = subset(gam_list_unique, chromosome_name == "Y")$external_gene_name
for (i in 1:length(pair_sim$key)){
  if(pair_sim$t2.name[i] %in% xgenes) {
    pair_sim$pair_genes[i] = paste(pair_sim$t2.name[i], pair_sim$t1.name[i], sep = " & ")
  } else {
    pair_sim$pair_genes[i] = paste(pair_sim$t1.name[i], pair_sim$t2.name[i], sep = " & ") 
  }}

# calculate non-weighted mean isoform similarity
pair_sim_mean = pair_sim %>% group_by(pair_genes) %>% summarise(mean_score = mean(score))
pair_sim_mean

# calculate expression weighted mean isoform similarity
trans_exp = read.csv('avg_gam_transcript_exp.csv')
colnames(trans_exp) = c('t1.id','t1.id_mean_exp')
pair_sim2 = merge(pair_sim, trans_exp, by = 't1.id')
colnames(trans_exp) = c('t2.id','t2.id_mean_exp')
pair_sim2 = merge(pair_sim2, trans_exp, by = 't2.id')
pair_sim2$mean_exp_both = rowSums(pair_sim2[,c(8,9)])
head(pair_sim2)
pair_sim_mean2 = pair_sim2 %>% group_by(pair_genes) %>% summarise(weighted_mean_score = weighted.mean(score,mean_exp_both))
pair_sim_mean2

comp = merge(pair_sim_mean2, pair_sim_mean, by = 'pair_genes')
ggplot(comp, aes(x = mean_score, y = weighted_mean_score)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_classic()

#scale mean similarity scores
pair_sim_mean$scale = (pair_sim_mean$mean_score - min(pair_sim$score)) / (max(pair_sim$score) - min(pair_sim$score))
pair_sim_mean

pair_sim_mean2$scale = (pair_sim_mean2$weighted_mean_score - min(pair_sim2$score)) / (max(pair_sim2$score) - min(pair_sim2$score))
pair_sim_mean2

## plot similarity measures

pair_sim_mean_plot = data.frame(pair_genes = pair_sim_mean2$pair_genes, variable = "mean_iso_sim", value = pair_sim_mean2$scale)

colnames(combo_plot) = colnames(pair_sim_mean_plot)
combo_plot2 = rbind(combo_plot, pair_sim_mean_plot)
combo_plot2

simout = dcast(combo_plot2, pair_genes ~ variable, value.var = 'value')

write.csv(simout, file = 'similarity_measures.csv')

library(egg)
combo_plot2$pair = factor(combo_plot2$pair, levels = c("PCDH11X & PCDH11Y","TGIF2LX & TGIF2LY","PRKX & PRKY","NLGN4X & NLGN4Y","EIF1AX & EIF1AY","TMSB4X & TMSB4Y","AMELX & AMELY","ZFX & ZFY","TBL1X & TBL1Y","TXLNG & TXLNGY","USP9X & USP9Y","DDX3X & DDX3Y","KDM6A & UTY","KDM5C & KDM5D","RPS4X & RPS4Y1","RPS4X & RPS4Y2","SOX3 & SRY"))
a = c("black","black","red","black","black","black","black","black","black","red","black","black","black","black","black","black")
ggplot(combo_plot2, aes(x = pair, y = value, fill = variable)) +
  geom_bar(stat = "identity",position = "dodge") + theme_article() +
  ylab('Similarity') +
  theme(axis.title.y = element_text(size = 12), axis.title.x = element_blank(), axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1, color = a)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9","#999999","darkolivegreen4","tomato3"), 
                    name=NULL, 
                    breaks=c("prop_match_overall", "prop_match_nongapped","DNA_sim","prot_sim","mean_iso_sim"), 
                    labels=c("Promoter (overall)", "Promoter (non-gapped)", "DNA Sequence", "Protein","Isoform (mean)"))

