###################
## upload meta data
###################

# sample attributes (N=22951)
attrib = read.delim('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
head(attrib)
keep_samples = subset(attrib, SMAFRZE == "RNASEQ")

# individual phenotypes (N=980)
pheno = read.csv('gtex_meta_edit.csv')
colnames(pheno)[2] = 'SAMPID2'
head(pheno)
table(pheno$SEX)

# merge meta data
library(stringr)
x = sub(".*GTEX-", "", keep_samples$SAMPID) 
xx = str_sub(x,1,6)
xxx = sub("-.*","",xx)
keep_samples$SAMPID2 = paste("GTEX-",xxx,sep="")
keep_samples = merge(keep_samples, pheno, by = 'SAMPID2')
table(keep_samples$SMTSD, keep_samples$SEX)
dim(keep_samples)

# save
saveRDS(keep_samples, file = 'gtex_combined_meta.rds')
