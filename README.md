# Evolution and disease implications of human X/Y gametologue genes

#### Repository for the analysis of human gametologue co-expression and coupled co-expression

This repository contains scripts used in the analysis of human gametologue co-expression and coupled co-expression. 

Note that we ran most steps on the NIH ([Biowulf](https://hpc.nih.gov/)) high-performance computing cluster. We have aimed to generalize the code here by removing system-specific references to installed software and modules. Instead, we document required software and version numbers below (excluding standard Unix programs and R). For HPC systems, the required scripts and binaries must be in the PATH. The easiest way to do this is to use an existing module or to install your own. In these cases, the modules should be loaded prior to running the appropriate code below.

As Biowulf uses the [slurm](https://slurm.schedmd.com/documentation.html) scheduler, most code below should run on slurm systems with little or no modification. For non-slurm HPC systems, slurm scripts and environmental variables will need to be adjusted, though hopefully without too much hassle.

We ran most analysis steps using [R](https://cran.r-project.org/) (v4.1.2). We recommend the following utility or visualization packages to extend base R's functionality.

# Inputs

The following files are expected:

* GTEx v8 fastq files should be placed in the ```fastq/``` folder with the naming convention ```${sample}/${sample}_1.fastq.gz``` (read 1) and ```${sample}/${sample}_2.fastq.gz``` (read 2).

* The GTEx phenotype file should placed in ```data/gtex_meta_edit.csv```

* The GTEx attributes file should placed in ```data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt```

* A .csv file including gametologue information should be placed in ```data/gametologs_in_genome.csv```
  
# Pipeline
  
### Combine metadata

* **Key libraries:** stringr

```
# Import metadata and combine
scripts/combine_meta.R
```
### Map reads with pseudoalignment software

* **Required software**: kallisto (v0.46.2), samtools (v1.13)

```
# Create sex specific transcriptomes and index
scripts/kallisto_transcriptomes.sh

# Count transcripts
sbatch --array=1-$(wc -l checkpoints/samplesM.txt | cut -d ' ' -f 1) scripts/run_kallisto_M.sh
sbatch --array=1-$(wc -l checkpoints/samplesF.txt | cut -d ' ' -f 1) scripts/run_kallisto_F.sh
```

### Import and combine expression data

* **Key libraries:** tximport, rdf5, biomaRt

```
# Import male and kallisto results into R and combine
scripts/import_kallisto.R
```

### Normalize, filter, and adjust gene expression

* **Key libraries:** edgeR, limma

```
# Apply filters to gene expression dataset
# Normalize data
# Adjust expression for age + technical effects
scripts/normalize_adjust.R
```

### Calculate co-expression (males & females), apply spatial quantile normalization, and calculate coupled co-expression (males) 

* **Key libraries:** spqn

```
# Apply filters to gene expression dataset
scripts/coexpression_coupled.R
```

### Visualize co-expression and coupled co-expression in males 

* **Key libraries:** ggplot2

```
# Visualize results
scripts/plot_coexpression.R
```

### Compare to previous results 

* **Key libraries:** ggplot2

```
# Visualize results
scripts/compare_to_previous.R
```

### Calculate differential coupling 

* **Key libraries:** clusterProfiler, reshape2

```
# Calculate difference in X vs. Y co-expression per gene
scripts/calc_coupling.R
```

### Estimate significance of differential coupling 

* **Key libraries:** parallel, stringr

```
# Visualize results
scripts/coupling_sig.sh
```

### Visualize differential coupling 

* **Key libraries:** ggplot2

```
# visualize results
scripts/visualize_GO_DO.R
```

### Perform GO/DO enrichments 

* **Key libraries:** ggplot2

```
# GO and DO analyses
scripts/GO_DO_coupling.R
```

### Visualize GO and DO results 

* **Key libraries:** ggplot2

```
# GO and DO analyses
scripts/visualize_GO_DO.R
```

### Calculate sex effects on gene expression 

* **Key libraries:** limma, mashr

```
# Visualize results
scripts/calc_sex_biased_expression.R
```


