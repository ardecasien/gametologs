`%!in%` = Negate(`%in%`)

library(biomaRt)
library(ape)
library(Biostrings)
library(reshape2)
library(ggplot2)
library(egg)
library(corrplot)

#### load gametolog info ####

gam_list = read.csv('gametologs_in_genome.csv')
colnames(gam_list)[3] = 'ensembl_gene_id'

hsap = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 104)
gam_bm = getBM(attributes = c("external_gene_name","ensembl_gene_id", "ensembl_transcript_id_version", "chromosome_name","transcript_biotype","transcript_gencode_basic"), filters = "ensembl_gene_id", values = gam_list$ensembl_gene_id, mart = hsap)
gam_list = merge(gam_list, gam_bm, by = 'ensembl_gene_id')
head(gam_list)
table(gam_list$transcript_gencode_basic, gam_list$transcript_biotype)

gam_list_unique = unique(gam_list[,c(1,12,4,14)])
gam_list_unique_new = data.frame(pair = unique(gam_list_unique$pair))
for(i in 1:length(levels(as.factor(gam_list_unique_new$pair)))){
  gam_list_unique_new$x_gene[i] = subset(gam_list_unique, pair == levels(as.factor(gam_list_unique_new$pair))[i] & chromosome_name == 'X')$external_gene_name
  gam_list_unique_new$x_gene_id[i] = subset(gam_list_unique, pair == levels(as.factor(gam_list_unique_new$pair))[i] & chromosome_name == 'X')$ensembl_gene_id
  gam_list_unique_new$y_gene[i] = subset(gam_list_unique, pair == levels(as.factor(gam_list_unique_new$pair))[i] & chromosome_name == 'Y')$external_gene_name
  gam_list_unique_new$y_gene_id[i] = subset(gam_list_unique, pair == levels(as.factor(gam_list_unique_new$pair))[i] & chromosome_name == 'Y')$ensembl_gene_id
}
gam_list_unique_new

#### end ####

##### estimate promoter similarity ####

for (i in 1:length(levels(as.factor(gam_list_unique_new$pair)))){
  
  pairnow = gam_list_unique_new[i,]
  xseq = getSequence(id = pairnow$x_gene_id, type="ensembl_gene_id",seqType="gene_flank",upstream=2000,mart=hsap) 
  yseq = getSequence(id = pairnow$y_gene_id, type="ensembl_gene_id",seqType="gene_flank",upstream=2000,mart=hsap) 
  gam_list_unique_new$seq_X[i] = xseq$gene_flank
  gam_list_unique_new$seq_Y[i] = yseq$gene_flank
  
  align = pairwiseAlignment(xseq$gene_flank, yseq$gene_flank, type = "global")
  seq = c(alignedPattern(align), alignedSubject(align))
  
  gam_list_unique_new$align1[i] = as.character(seq[1])
  gam_list_unique_new$align2[i] = as.character(seq[2])
  gam_list_unique_new$align_len[i] = nchar(align)
  indel = nindel(align)
  gam_list_unique_new$nindel[i] = indel@insertion[2] + indel@deletion[2]
  gam_list_unique_new$nmatch[i] = nmatch(align)
  gam_list_unique_new$nmismatch[i] = nmismatch(align)
  gam_list_unique_new$score[i] = align@score
}

gam_list_unique_new$prop_match_overall = gam_list_unique_new$nmatch / gam_list_unique_new$align_len
gam_list_unique_new$prop_match_nongapped = gam_list_unique_new$nmatch / (gam_list_unique_new$nmatch + gam_list_unique_new$nmismatch)
gam_list_unique_new$pair_genes = paste(gam_list_unique_new$x_gene, gam_list_unique_new$y_gene, sep = " & ")
saveRDS(gam_list_unique_new, 'gam_list_unique_new_promoter.rds')

gam_list_unique_new_plot = melt(gam_list_unique_new[,c("pair_genes","prop_match_overall","prop_match_nongapped")], id.vars = 1)

ggplot(gam_list_unique_new_plot, aes(x = reorder(pair_genes,value), y = 1-value, fill = variable)) +
  geom_bar(stat = "identity",position = "dodge") + theme_classic() +
  ylab('Promoter Dissimilarity') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9"), name=NULL, breaks=c("prop_match_overall", "prop_match_nongapped"), labels=c("overall", "non-gapped"))

#### end ####

#### DNA and protein similarity ####

skalet = read.csv('Skaletsky_comparison.csv', check.names = F) # from Skaletsy et al. 2003
skalet = subset(skalet, Gametologue == "Y")

skalet$DNA_sim = (100 - skalet$`DNA divergence perc`)/100
skalet$prot_sim = (100 - skalet$`Protein divergence perc`)/100
skalet_plot = melt(skalet[,c("pair_genes","DNA_sim","prot_sim")], id.vars = 1)
skalet_plot = skalet_plot[which(skalet_plot$pair_genes %in% gam_list_unique_new_plot$pair_genes),]
combo_plot = rbind(gam_list_unique_new_plot,skalet_plot)
saveRDS(combo_plot, file = 'combo_plot_promoter_DNA_protein.rds')

simout = dcast(combo_plot, pair_genes ~ variable, value.var = 'value')

# Table S5

write.csv(simout, file = 'similarity_measures.csv')

ggplot(combo_plot, aes(x = reorder(pair_genes,value), y = 1-value, fill = variable)) +
  geom_bar(stat = "identity",position = "dodge") + theme_classic() +
  ylab('Dissimilarity') +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9","#999999","darkolivegreen4"), 
                    name=NULL, 
                    breaks=c("prop_match_overall", "prop_match_nongapped","DNA_sim","prot_sim"), 
                    labels=c("promoter (overall)", "promoter (non-gapped)", "DNA sequence", "protein sequence"))

#### end ####

#### Figure S6A | Table S6 ####

combo_plot = readRDS(file = 'combo_plot_promoter_DNA_protein.rds')
combo_plot$value = 1 - combo_plot$value
corp = dcast(data = combo_plot,formula = pair_genes~variable,fun.aggregate = sum,value.var = "value")
colnames(corp) = c('pair_genes','prom-all','prom-ng','DNA-sim','pro-sim')
corp
m = cor(corp[,c(2:5)])
m
corrplot.mixed(m, 
               order = 'AOE', 
               tl.col = 'black',
               tl.pos = 'd')
cor.mtest(corp[,c(2:5)])

#### end ####

#### Figure 3B ####

groups = read.csv('gam_groups.csv')
colnames(combo_plot)[1] = 'pair'
combo_plot2 = merge(combo_plot, groups, by = 'pair')
combo_plot2$pair = factor(combo_plot2$pair, levels = c("SOX3 & SRY", "RPS4X & RPS4Y2", "RPS4X & RPS4Y1", "KDM5C & KDM5D", "KDM6A & UTY", "DDX3X & DDX3Y", "USP9X & USP9Y", "TXLNG & TXLNGY", "TBL1X & TBL1Y", "ZFX & ZFY", "AMELX & AMELY", "TMSB4X & TMSB4Y", "EIF1AX & EIF1AY", "NLGN4X & NLGN4Y", "PRKX & PRKY", "TGIF2LX & TGIF2LY", "PCDH11X & PCDH11Y"))
a = c("black","black","black","black","black","black","black","grey","black","black","black","black","black","black","grey","black","black")

ggplot(combo_plot2, aes(x = pair, y = 1-value, fill = variable)) +
  geom_bar(stat = "identity",position = "dodge") + theme_article() +
  ylab('Dissimilarity') +
  ylim(c(0,1)) +
  facet_grid(~group, scale = 'free_x', space = 'free') +
  scale_fill_manual(values=tissue.colors[c(1,10,20,30,40)], 
                    name=NULL, 
                    breaks=c("prop_match_overall", "prop_match_nongapped","DNA_sim","prot_sim"), 
                    labels=c("Promoter (overall)", "Promoter (non-gapped)", "DNA", "Protein")) +
  theme(axis.title.y = element_text(size = 12), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1, color = a),
        #axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1, color = a, face = 'italic'),
        legend.position = 'bottom',
        legend.text = element_text(size = 12)) +
  guides(fill=guide_legend(ncol=2,byrow=TRUE)) 

#### end ####

#### Figure S6B-D ####

simout = read.csv('similarity_measures.csv')
simout = simout[,-1]
simout = reshape2::melt(simout)

gam_list = read.csv('gametologs_in_genome.csv')
gam_list$pair = as.factor(gam_list$pair)
gams = data.frame()
for(i in 1:length(levels(gam_list$pair))){
  gams[i,1] = levels(gam_list$pair)[i]
  xgene = subset(gam_list, pair == levels(gam_list$pair)[i] & Gametolog == 'X')$common_name
  ygene = subset(gam_list, pair == levels(gam_list$pair)[i] & Gametolog == 'Y')$common_name
  gams[i,2] = paste(xgene, ygene, sep = " & ")
}
colnames(gams) = c('pair','pair_genes')

simout = merge(simout, gams, by = 'pair_genes')
simout = merge(simout, unique(gam_list[,c('pair','Strata2')]), by = 'pair')
simout$Strata2 = ifelse(simout$Strata2 == "", "XTR", simout$Strata2)

m = unique(simout$variable)
for(i in 1:length(m)){
  print(m[i])
  n = subset(simout, variable == m[i])
  mod = aov((1-value) ~ Strata2, data = n)
  print(summary(mod))
  print(TukeyHSD(mod))
}

# Figure S6B

simout$variable = factor(simout$variable)
levels(simout$variable) = c("Promoter (overall)", "Promoter (non-gapped)", "DNA", "Protein")
ggplot(simout, aes(x = Strata2, y = 1-value, fill = variable)) +
  geom_boxplot() +
  ylim(c(0,1)) +
  scale_fill_manual(values=tissue.colors[c(1,10,20,30)]) +
  facet_wrap(~variable) +
  theme_article() +
  ylab('Dissimilarity') + xlab("") +
  theme(legend.position = 'none')

# Figure S6C

p = unique(simout[,c('pair_genes','Strata2')])
colnames(p)[1] = 'gene'
out_scdced = readRDS(file = 'out_scdced.rds')
out_scdced = merge(out_scdced, p, by = 'gene')

ggplot(out_scdced, aes(x = Strata2, y = abs_diff)) +
  geom_boxplot() +
  ylim(c(0,1)) +
  theme_article() +
  ylab('aCFD') + xlab("") +
  theme(legend.position = 'none')

# Figure S6D

ggplot(out_scdced, aes(x = Strata2, y = abs_diff, fill = gene)) +
  geom_boxplot(position = position_dodge(preserve = "single")) +
  ylim(c(0,1)) +
  theme_article() +
  ylab('aCFD') + xlab("") +
  theme(legend.position = 'bottom',
        legend.title = element_blank())

#### end ####
