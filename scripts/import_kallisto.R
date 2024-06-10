## import .h5 files into R

library(tximport)
library(rhdf5)
library(biomaRt)
library(stringr)
library(data.table)

`%!in%` = Negate(`%in%`)

samplesF=scan('samplesF.txt', what = "")
samplesM=scan('samplesM.txt', what = "")

filesF = file.path('/data/NIMH_scratch/GTEx/process/kallisto_ssc/females',samplesF,'abundance.h5')
names(filesF) = samplesF

filesM = file.path('/data/NIMH_scratch/GTEx/process/kallisto_ssc/males',samplesM,'abundance.h5')
names(filesM) = samplesM

txi.kallistoF = tximport(filesF, type = 'kallisto', txOut = TRUE, countsFromAbundance="lengthScaledTPM")
txi.kallistoM = tximport(filesM, type = 'kallisto', txOut = TRUE, countsFromAbundance="lengthScaledTPM")

## get which Y chr genes are not in females file
Yadd = rownames(txi.kallistoM$counts)[which(rownames(txi.kallistoM$counts) %!in% rownames(txi.kallistoF$counts))]

## add Y chr to females with 0 count
YaddF = data.frame(matrix(0, nrow = length(Yadd), ncol = length(colnames(txi.kallistoF$counts))))
rownames(YaddF) = Yadd
colnames(YaddF) = colnames(txi.kallistoF$counts)
txi.kallistoF$abundance = rbind(txi.kallistoF$abundance, YaddF)
txi.kallistoF$counts = rbind(txi.kallistoF$counts, YaddF)
txi.kallistoF$length = rbind(txi.kallistoF$length, YaddF)

## reorder data frames
idx = match(rownames(txi.kallistoF$length), rownames(txi.kallistoM$length))
txi.kallistoM$abundance = txi.kallistoM$abundance[idx,]
txi.kallistoM$counts = txi.kallistoM$counts[idx,]
txi.kallistoM$length = txi.kallistoM$length[idx,]

## check
table(rownames(txi.kallistoF$abundance) == rownames(txi.kallistoM$abundance))
table(rownames(txi.kallistoF$counts) == rownames(txi.kallistoM$counts))
table(rownames(txi.kallistoF$length) == rownames(txi.kallistoM$length))

## merge males and females
txi.kallisto = list(abundance=as.matrix(cbind(txi.kallistoF$abundance,txi.kallistoM$abundance)),counts=as.matrix(cbind(txi.kallistoF$counts,txi.kallistoM$counts)),length=as.matrix(cbind(txi.kallistoF$length,txi.kallistoM$length)),countsFromAbundance=c("no"))

## save transcript-level 
saveRDS(txi.kallisto, "gtex_kallisto_transcripts.rds")

## import 
txi.kallisto = readRDS("gtex_kallisto_transcripts.rds")
idx = match(rownames(txi.kallisto$abundance), rownames(txi.kallisto$length))
txi.kallisto$length = txi.kallisto$length[idx,]

## convert to counts per gene
hsap = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl') 

## get a file that matches ensembl transcripts (with version appended) to ensembl genes
tx2gene = getBM(attributes=c('ensembl_transcript_id_version', 'ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position'), filters = 'ensembl_transcript_id_version', values = rownames(txi.kallisto$abundance), mart = hsap)

## summarize your kallisto mapped data to the gene level
txi.gene = summarizeToGene(txi.kallisto, tx2gene) # transcripts missing from tx2gene: 3651

## save gene-level 
saveRDS(txi.gene, file = "gtex_kallisto_genes.rds")
