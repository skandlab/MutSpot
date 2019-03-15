## Non-coding MUTation hotSPOT dectection in cancer genomes

This package systematically scans cancer genomes to detect small, focal regions with higher than expected mutation recurrence (mutation hotspots) through binomial logistic modelling. Such regions display signatures of positive selection and therefore are likely to be causal to cancer development.

------------------------------------------------------------------------------------

## Installation

```{r}
install.packages("devtools")
library(devtools)
install_github("changmmm/MutSpot", subdir="MutSpot_Rpackage")
```

----------------------------------------------------------------------------------

## MutSpot analysis workflow
The full MutSpot workflow include the following 9 steps:

1. Sample non-mutated sites as negative examples for logistic regression.

2. Calculate local mutation rates in 100kb bins across the whole genome.

3. Select sequence features using LASSO logistic regression.

4. Select epigenetic features using LASSO logistic regression.

5. Compute feature matrix for model fitting based on selected features for all sites in whole genome.

6. Fit logistic regression model of background mutation probabilities.

7. Predict mutation hotspots.

8. Annotate hotspots.

9. Generate figures.

By default, the MutSpot() function runs the entire workflow. However, it is possible to run specific steps of the workflow by specifiying the *run.to* parameter (see full documentation).

## Usage example
All intermediate and output files will be saved in the output directory specified by the user. Package should be run in the same directory.

Install the package using the following commands.
```r
library("MutSpot")
output.dir = "./"
setwd(output.dir)
```

Download the test data sets from https://github.com/skandlab/MutSpot/tree/master/test-data into your working directory.

Run the analysis using the following commands:

Identify SNV and indel hotspots genome-wide

*Whole genome analysis will take up to 1 day using 2 cores*

```r
MutSpot(snv.mutations = "subset_snv_mutations_sid.MAF", indel.mutations = "subset_indel_mutations_sid.MAF", genomic.features = "genomic_features_genome.txt", fit.sparse = TRUE, min.count.snv = 3, min.count.indel = 3)
```

Identify SNV hotspots in CTCF binding site only, including clinical subtype and cosmic signatures as sample specific features.

*CTCF analysis will take up to 2 hours using 1 core*

```r
MutSpot(snv.mutations = "subset_snv_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 9, genomic.features = "genomic_features_ctcf.txt",
sample.snv.features = "sample_features_table_snv.txt", drop = TRUE)
```

## Main arguments

 Parameter                                    | Description
--------------------------------------------- | -----------------------------------------
 snv.mutations                                | List of SNVs in MAF format
 indel.mutations                              | List of indels in MAF format
 genomic.features                             | File paths of all potential genomic features for LASSO selection for background mutation model. E.g. Replication timing profile, histone modification profiles
 sample.snv.features/sample.indel.features    | Tab delimited file of sample specific features. E.g. clinical subtype
 min.count                                    | Minimum number of mutated samples in each hotspot (default=2)
 fit.sparse                                   | To fit background model using sparse matrix with GLMNET (default=FALSE)
 region.of.interest                           | Restrict hotspot analysis to regions in the given bed file

## Input files
##### 1. Mutations
Mutation files contain all SNVs or all indels of all tumors in the study in the MAF format. MAF file should be tab delimited with exactly 6 columns: chromosome, start position (1-based), end position (1-based), reference allele, alternate allele, and sample ID. There is no header row in MAF files.
Example MAF file:

|      |          |          |   |   |         |
|------|----------|----------|---|---|---------|
| chr1 | 16265287 | 16265287 | G | C | patient1|
| chr1 | 17320166 | 17320166 | C | T | patient2|
| chr1 | 19497536 | 19497536 | G | C | patient3|
| ...  | ...      | ...      | ..| ..| ...     |


##### 2. Genomic features
Genomic features can be continuous or binary. Continuous features, such as replication timing profile, are input as bigwig files. Binary features, such as peak calls of histone modifications, are input as bed files. All continuous features will be discretized into n bins (n is specified by the user). The the logistic regression model will be fit from a frequency table of the counts of mutated and non-mutated sites for all combinations of the covariates. It is recommended for the user to input genomic covariates as binary features where possible to reduce the memory usage of the function.

The genomic features are input as a tab delimited file with 4 columns:

  1. feature name
  2. file path of genomic feature (either bigwig or bed format).
  3. Binary value indicating if the feature is continuous or binary (1 for continuous, 0 for binary)
  4. Number of bins to discretize continuous feature into (max=10, NA for binary features).

Example format:

 feature_name  | file_path                                               | feature_type | nbins
-------------- | ------------------------------------------------------- | ------------ | -----
 mean_rep_time | ./features/wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig | 1            | 8     
 E094-DNase    | ./features/E094-DNase.bed                               | 0            | NA    
 E094-H3K27ac  | ./features/E094-H3K27ac.bed                             | 0            | NA    
 ...           | ...                                                     | ...          | ...   

A binary feature bed file should include the following columns:

  1. Chromosome
  2. Start position
  3. End position

*Genomic regions in the bed file are considered 1 for the binary feature and regions not in the bed file are considered 0 for the binary feature*


##### 3. Sample specific features
The user can choose to include sample specific features in the background mutation model, such as the clinical subtype of the tumor. Note that sample specific features will not undergo LASSO feature selection and will be automatically included in the final model. Sample specific features are to be supplied as a tab delimited file where each row corresponds to a sample and each column corresponds to a feature.

Example format:

 SampleID | subtype | feature1 | feature2
 -------- | ------- | -------- | --------
 patient1 | EBV     | 0.32     | 0.6      
 patient2 | MSI     | 0.41     | 1.5      
 patient3 | GS      | 0.18     | -0.3     
 ...      | ...     | ...      | ...      


##### 4. Region of interest
Instead of finding mutation hotspots genome-wide, the user could restrict the hotspot analysis to certain regions of interest, such as promoters, enhancers, or UTRs, by supplying a bed file with the following columns:

  1. Chromosome
  2. Start position
  3. End position

-----------------------------------------------------------------------------------

## Adjusting threshold of LASSO feature selection

For a more/less stringent nucleotide selection, users may choose to redefine frequency threshold by re-running step 3.2. This can be done by specifying the *run.to* parameter and the new selection threshold has to be defined by *cutoff.nucleotide.new*.

```{r}
MutSpot(run.to = 3.2, snv.mutations = "subset_snv_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cutoff.nucleotide.new = 0.98, cores = 9, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt", drop = TRUE)
```


Similarly, for a more/less stringent epigenetic features selection, users may choose to redefine frequency threshold by running step 4.2. This can be done by specifying the *run.to* parameter and the new selection thresholds for SNVs and indels can be defined by *cutoff.nucleotide.new.snv* and *cutoff.features.new.indel* respectively.

```{r}
MutSpot(run.to = 4.2, snv.mutations = "subset_snv_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 9, genomic.features = "genomic_features_ctcf.txt",
sample.snv.features = "sample_features_table_snv.txt", cutoff.features.new.snv = 0.8, drop = TRUE)
```

-----------------------------------------------------------------------------------

## Output files

#### Hotspot summary file
MutSpot outputs a TSV file of all hotspot regions found, ordered by significance of recurrence. Overlapping hotspot regions are merged and annotations columns are added to indicate if the hotspots are in gene promoters, UTRs or introns. The output file has the following fields:

1. Chromosome
2. Start position
3. End position
4. p-value
5. Number of bases
6. Background p-value
7. Number of mutated samples
8. FDR
9. Promoters
10. 3'UTR
11. 5'UTR

-----------------------------------------------------------------------------------

#### Figures
At the end of the analysis, 3 figures will be generated by MutSpot:

- Bar plot of feature importance in the background mutation model
- Manhattan plot of hotspots across the genome
- Distribution of mutations in the top n hotspots (default n=3, see documentation on changing the number hotspots to plot)

-------------------------------------------------------------------------------------
