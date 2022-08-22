###############################
## download human transcriptome
###############################

mkdir transcriptomes

cd transcriptomes

H38_URL_cdna=ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

wget $H38_URL_cdna

gunzip -c Homo_sapiens.GRCh38.cdna.all.fa.gz > h38_cdna.fa

###############################
## remove Y for females
###############################

grep ENST h38_cdna.fa | grep -v "chromosome:GRCh38:Y" | cut -f1 -d " " | sed -e "s/>//g" > transcripts_nonY

module load samtools/1.13

samtools faidx h38_cdna.fa

samtools faidx -r transcripts_nonY h38_cdna.fa > h38_cdna_noY.fa

###############################
## remove PAR for males 
###############################

# no need to edit for males since PAR should not be on Y in Ensembl
# https://useast.ensembl.org/info/genome/genebuild/human_PARS.html

# chromosome:GRCh38:Y:10001 - 2781479 is shared with X: 10001 - 2781479 (PAR1)

grep ENST h38_cdna.fa | grep "chromosome:GRCh38:Y" | awk -F':' '$4 > 10001 && $5 < 2781479' 
# only 2 transcripts for XGY2 gene (spans PAR1 boundary 2752083 - 2854641)
# ENST00000683447.1 2752296 - 2774739
# ENST00000680009.1 2752296 - 2774996
# but XGY2 is pseudogene consisting of the first 3 exons of XG in PAR1
# https://www.exeley.com/exeley/journals/immunohematology/36/1/pdf/10.21307_immunohematology-2020-035.pdf

grep ENST h38_cdna.fa | grep "XGY2"
# all are in PAR or span PAR boundary

# confirmed

###############################
## index transcriptomes
###############################

mkdir indexed_t

module load kallisto/0.46.2

kallisto index -i indexed_t/h38.idx transcriptomes/h38_cdna.fa
kallisto index -i indexed_t/h38_noY.idx transcriptomes/h38_cdna_noY.fa
