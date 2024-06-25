library(stringr)

`%!in%` = Negate(`%in%`)

# sample attributes (N=22951)
attrib = read.delim('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
head(attrib)
keep_samples = subset(attrib, SMAFRZE == "RNASEQ")
tissue.list = levels(as.factor(attrib$SMTSD))

# remove cells/blood + sex-specific/differentiated tissues (PRO, TES, MAM) (N = 40 tissues)

remove.tissues = c('Whole Blood','Breast - Mammary Tissue','Bladder','Cells - Cultured fibroblasts','Cells - EBV-transformed lymphocytes','Cells - Leukemia cell line (CML)','Cervix - Ectocervix','Cervix - Endocervix','Fallopian Tube','Kidney - Medulla','Ovary','Prostate','Testis','Uterus','Vagina')
tissue.list = tissue.list[which(tissue.list %!in% remove.tissues)]
tissue.list = str_remove_all(tissue.list, pattern = fixed(" "))
tissue.list = str_remove_all(tissue.list, pattern = fixed("_"))
tissue.list = str_remove_all(tissue.list, pattern = fixed("-"))
tissue.list = str_remove_all(tissue.list, pattern = fixed("("))
tissue.list = str_remove_all(tissue.list, pattern = fixed(")"))

for (i in 1:length(tissue.list)) {
		line = paste("Rscript /checkpoint/sex_dependent_CFD_MXY_FXX.R", tissue.list[i], sep = " ")
		write(line, file = "sex_dependent_CFD_MXY_FXX.swarm", append = TRUE, sep = "\n")
	}
	
for (i in 1:length(tissue.list)) {
		line = paste("Rscript /checkpoint/sex_dependent_CFD_MX_FXX.R", tissue.list[i], sep = " ")
		write(line, file = "sex_dependent_CFD_MX_FXX.swarm", append = TRUE, sep = "\n")
	}
