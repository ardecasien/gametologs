#!/bin/bash

module load kallisto/0.46.2

idx_path=/data/indexed_t/h38.idx

fastq_path=/data/fastq

output_path=/data/kallisto_ssc/males

sample=sed -n ${SLURM_ARRAY_TASK_ID}p samplesM.txt

kallisto quant -bias -i ${idx_path} -t 24 -o ${output_path}/${sample} ${fastq_path}/${sample}/${sample}_1.fastq.gz ${fastq_path}/${sample}/${sample}_2.fastq.gz

print('done')
