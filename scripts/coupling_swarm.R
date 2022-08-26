library(reshape2)

attrib = read.delim('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
keep_samples = subset(attrib, SMAFRZE == "RNASEQ")
tissue.list = levels(as.factor(attrib$SMTSD))
remove.tissues = c('Whole Blood','Bladder','Cells - Cultured fibroblasts','Cells - EBV-transformed lymphocytes','Cells - Leukemia cell line (CML)','Cervix - Ectocervix','Cervix - Endocervix','Fallopian Tube','Kidney - Medulla','Ovary','Uterus','Vagina')
tissue.list = tissue.list[which(tissue.list %!in% remove.tissues)]

tissue.list = gsub(" ","",tissue.list)
tissue.list = gsub("[(]","",tissue.list)
tissue.list = gsub("[)]","",tissue.list)
tissue.list = gsub("-","",tissue.list)
tissue.list

for (i in 1:length(tissue.list)){
	line = paste("Rscript /data/DNU/alex/gtex/coupling_sig.R", tissue.list[i], sep = " ")
	write(line, file = "XminusY.swarm", append = TRUE, sep = "\n")
	}
  
