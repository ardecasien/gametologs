library(reshape2)

attrib = read.delim('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
keep_samples = subset(attrib, SMAFRZE == "RNASEQ")
tissue.list = levels(as.factor(attrib$SMTSD))
tissue.list = gsub(" ","",tissue.list)
tissue.list = gsub("[(]","",tissue.list)
tissue.list = gsub("[)]","",tissue.list)
tissue.list = gsub("-","",tissue.list)

for (i in 1:length(tissue.list)){
	line = paste("Rscript /data/DNU/alex/gtex/coupling_sig.R", tissue.list[i], sep = " ")
	write(line, file = "XminusY.swarm", append = TRUE, sep = "\n")
	}
  
