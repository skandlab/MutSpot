MutSpot
====================================

## Non-coding MUTation hotSPOT dectection in cancer genomes
The MutSpot R package systematically and unbiasedly scans cancer whole genomes to detect mutation hotspots. MutSpot first builds a background mutation model that corrects for covariates of mutation probability, such as local nucleotide context, replication timing and epigenomic features. Then MutSpot evaluates the mutation recurrence of focal DNA regions using a Poisson binomial model to account for varying mutation rates across different tumors. Mutation hotspots identified have significantly higher mutation recurrence compared to the background genomic mutation rate, suggesting positive selection in cancer and involvement in tumorigenesis.

Reference: [Guo et al., Nature Communications, 2018](https://www.nature.com/articles/s41467-018-03828-2)

------------------------------------------------------------------------------------

## Contents

[Installation](#installation)
<br/>[MutSpot analysis workflow](#workflow)
<br/>[Usage example](#usage)
<br/>[Main arguments](#arguments)
<br/>[Input files](#input)
<br/>[Adjusting threshold of LASSO feature selection](#adjusting-threshold)
<br/>[Output files](#output)

------------------------------------------------------------------------------------

<a name="installation"></a>

## Installation

MutSpot runs on R (requires at least 3.2.0). Install the package from Github using the following R commands.

```{r}
install.packages("devtools")
library(devtools)
install_github("skandlab/MutSpot", subdir="MutSpot_Rpackage")
```

Alternatively, the package may downloaded from Github and installed in R:
```{r}
# Clone/download MutSpot into the current working dirctory with the following command: *git clone https://github.com/skandlab/MutSpot.git* 
library(devtools)
install("my/current/directory/MutSpot/MutSpot_Rpackage")
```

*The alternative method takes a longer time as it downloads all test data sets at the same time*

----------------------------------------------------------------------------------

<a name="workflow"></a>

## MutSpot analysis workflow
The full MutSpot workflow includes the following 7 steps:

1. Sample non-mutated sites as negative examples for logistic regression.

2. Calculate local mutation rates in 100kb bins across the whole genome.

3. Select sequence features using LASSO logistic regression.

4. Select epigenetic features using LASSO logistic regression.

5. Compute feature matrix for model fitting based on selected features for all sites in whole genome.

6. Fit logistic regression model of background mutation probabilities.

7. Predict mutation hotspots.

By default, the *MutSpot()* function runs the entire workflow. However, it is possible to run specific steps of the workflow by specifiying the *run.to* parameter (see full documentation).

-----------------------------------------------------------------------------------

<a name="usage"></a>

## Usage example
By default, MutSpot runs in the current working directory unless specified by the user. All intermediate and output files will be saved in the *results* folder created by MutSpot in the working directory. MutSpot runs genome-wide. However, if the user provides a BED file containing the coordinates of a specific region under the *region.of.interest* parameter, MutSpot runs on the user-specified region (e.g. CTCF binding sites) to find only the hotspots in the given region.

```r
library("MutSpot")
```

Download the test data sets from https://github.com/skandlab/MutSpot/tree/master/test-data into your working directory.

Run the analysis using the following commands:


Identify SNV hotspots genome-wide.
(*Whole genome analysis will take less than 1 day using 2 cores*)
```r
MutSpot(snv.mutations = "subset_snv_mutations_sid.MAF", cores = 2, cutoff.nucleotide.new = 1, genomic.features = "genomic_features_genome.txt")
```

Identify SNV hotspots in CTCF binding sites only, including clinical subtype and cosmic signatures as sample specific features.
(*CTCF analysis will take about 30 minutes using 2 cores*)
```r
MutSpot(snv.mutations = "subset_snv_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, genomic.features = "genomic_features_ctcf.txt",
sample.snv.features = "sample_features_table_snv.txt")
```

----------------------------------------------------------------------------------

<a name="arguments"></a>

## Main arguments

 Parameter                                    | Description
--------------------------------------------- | -----------------------------------------
 snv.mutations                                | List of SNVs in MAF format
 indel.mutations                              | List of indels in MAF format
 genomic.features                             | File paths of all potential genomic features for LASSO selection for background mutation model. E.g. Replication timing profile, histone modification profiles
 sample.snv.features/sample.indel.features    | Tab delimited file of sample specific features. E.g. clinical subtype
 min.count                                    | Minimum number of mutated samples in each hotspot (default = 2)
 region.of.interest                           | Restrict hotspot analysis to regions in the given BED file
 cores                                        | Number of cores (default = 1) [The maximum number of cores for Windows users is 1.]

----------------------------------------------------------------------------------

<a name="input"></a>

## Input files

#### 1. Mutations
Mutation files contain all SNVs or indels of all tumors in the study in MAF format. MAF file should be tab delimited with the following 6 columns:

  1. Chromosome
  2. Start position (1-based)
  3. End position (1-based)
  4. Reference allele
  5. Alternate allele
  6. Sample ID

Example MAF file:

|      |          |          |   |   |          |
|------|----------|----------|---|---|----------|
| chr1 | 16265287 | 16265287 | G | C | patient1 |
| chr1 | 17320166 | 17320166 | C | T | patient2 |
| chr1 | 19497536 | 19497536 | G | C | patient3 |
| ...  | ...      | ..       |.. |.. | ...      |    

*There is no header row in a MAF file.*


#### 2. Genomic features
Genomic features can be continuous or binary. Continuous features, such as replication timing profile, are input as bigwig files. Binary features, such as peak calls of histone modifications, are input as BED files. All continuous features will be discretized into *n* bins (*n* is specified by the user). The logistic regression model will be fitted from a frequency table of the counts of mutated and non-mutated sites for all combinations of the covariates. It is recommended for the user to input genomic covariates as binary features where possible to reduce the memory usage of the function.

The genomic features are input as a tab delimited file with 4 columns:

  1. Feature name
  2. File path of genomic feature (either BigWig or BED format).
  3. Binary value indicating if the feature is continuous or binary (1 for continuous, 0 for binary)
  4. Number of bins to discretize continuous feature into (max = 10, NA for binary features).

Example format:

 feature_name  | file_path                                               | feature_type | nbins
-------------- | ------------------------------------------------------- | ------------ | -----
 mean_rep_time | ./features/wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig | 1            | 10     
 E094-DNase    | ./features/E094-DNase.bed                               | 0            | NA    
 E094-H3K27ac  | ./features/E094-H3K27ac.bed                             | 0            | NA    
 ...           | ...                                                     | ...          | ...   

A binary feature BED file should include the following columns:

  1. Chromosome
  2. Start position (0-based)
  3. End position (0-based)

*For binary features, genomic regions that are found in the feature BED file are assigned value of 1, else value of 0*

A list of genomic feature files (Transcription factors, DNA secondary structure, Replication timing) can be downloaded from https://github.com/skandlab/MutSpot/tree/master/features into the *features* folder in your working directory. The user may choose to run the analysis using these features by specifying *genomic.features = "./features/genomic_features_genome_default.txt"* in the *MutSpot()* function, else he/she may create a similar text file containing desired/other features.


#### 3. Sample specific features (optional)
The user may choose to include sample specific features in the background mutation model, such as the clinical subtype of the tumor. Note that sample specific features will not undergo LASSO feature selection and will be automatically included in the final model. Sample specific features are to be supplied as a tab delimited file where each row corresponds to a sample and each column corresponds to a feature.

Example format:

 SampleID | subtype | feature1 | feature2
 -------- | ------- | -------- | --------
 patient1 | EBV     | 0.32     | 0.6      
 patient2 | MSI     | 0.41     | 1.5      
 patient3 | GS      | 0.18     | -0.3     
 ...      | ...     | ...      | ...      


#### 4. Region of interest (optional)
Instead of finding mutation hotspots genome-wide, the user may restrict the hotspot analysis to certain regions of interest, such as promoters, enhancers, or UTRs, by supplying a BED file with the following columns:

  1. Chromosome
  2. Start position (0-based)
  3. End position (0-based)

  Example BED file:

  |      |          |          |
  |------|----------|----------|
  | chr1 | 15786447 | 16265287 |
  | chr1 | 27891466 | 28456878 |
  | chr1 | 42456878 | 45468785 |
  | ...  | ...      | ..       |   

  *There is no header row in a BED file.*

-----------------------------------------------------------------------------------

<a name="adjusting-threshold"></a>

## Adjusting threshold of LASSO feature selection

For a more/less stringent nucleotide selection, users may choose to re-define frequency threshold by re-running step 3.2. This can be done by specifying *run.to = 3.2* and the new selection threshold as the *cutoff.nucleotide.new* parameter. Users may also choose to select the top *n* features based on the mean coefficients by specifying the number of contexts to select as the *top.nucleotide* parameter.

```{r}
MutSpot(run.to = 3.2, snv.mutations = "subset_snv_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cutoff.nucleotide.new = 0.98, top.nucleotide = 3, cores = 2, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt")
```


Similarly, for a more/less stringent epigenetic features selection, users may choose to re-define frequency threshold by running step 4.2. This can be done by specifying *run.to = 4.2* and the new selection thresholds for SNVs and indels as the *cutoff.nucleotide.new.snv* and *cutoff.features.new.indel* parameters respectively. Users may also choose to select the top *n* features based on the mean coefficients by specifying the number of features to select as the *top.features* parameter.

```{r}
MutSpot(run.to = 4.2, snv.mutations = "subset_snv_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, genomic.features = "genomic_features_ctcf.txt",
sample.snv.features = "sample_features_table_snv.txt", cutoff.features.new.snv = 0.8, top.features = 11)
```

-----------------------------------------------------------------------------------

<a name="output"></a>

## Output files

#### Hotspot summary file
MutSpot outputs a TSV file of all hotspot regions found, ordered by significance of recurrence. Overlapping hotspot regions are merged and annotations columns are added to indicate if the hotspots are located in gene promoters or UTRs. The output file has the following fields:

1. Chromosome
2. Start position (1-based)
3. End position (1-based)
4. *P*-value
5. Length of hotspot (bp)
6. Mean background mutation probability
7. Number of mutated samples
8. FDR
9. Transcripts overlapping hotspot in their promoters
10. Transcripts overlapping hotspot in their 3'UTRs
11. Transcripts overlapping hotspot in their 5'UTRs


#### Figures
At the end of the analysis, 3 figures will be generated by MutSpot:

- Bar plot of feature importance of the background mutation model
- Manhattan plot of hotspots across the genome
- Distribution of mutations in the top *n* hotspots (default *n* = 3)

---------------------------------------------------------------
