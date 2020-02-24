MutSpot
====================================

## Non-coding MUTation hotSPOT dectection in cancer genomes
The MutSpot R package systematically and unbiasedly scans cancer whole genomes to detect mutation hotspots. MutSpot first builds a background mutation model that corrects for known covariates of mutation probability, such as local nucleotide context, replication timing and epigenomic features. Then MutSpot evaluates the mutation recurrence of focal DNA regions using a Poisson binomial model to account for varying mutation rates across different tumors. Mutation hotspots identified have significantly higher mutation recurrence compared to the background genomic mutation rate, suggesting positive selection in cancer and involvement in tumorigenesis.

Reference: [Guo et al., Nature Communications, 2018](https://www.nature.com/articles/s41467-018-03828-2)

------------------------------------------------------------------------------------

## Contents

[Arguments](#arguments)
<br/>[Step 1: Sample sites](#sample-sites)
<br/>[Step 2: Compute local mutation rate](#compute-local-mutation-rate)
<br/>[Step 3: Select nucleotide contexts](#select-nucleotide-contexts)
<br/>[Step 4: Select genomic features](#select-genomic-features)
<br/>[Step 5: Compute feature matrix](#compute-feature-matrix)
<br/>[Step 6: Fit prediction model](#fit-prediction-model)
<br/>[Step 7: Predict hotspots](#predict-hotspots)

------------------------------------------------------------------------------------

<a name="arguments"></a>

## Arguments

| Parameter | Type | Description | Default |
|-----------|------|---------|-------------|
| run.to | Numeric | Steps to run | 1,2,3.1,3.2,4.1,4.2,5.1,5.2,5.3,5.4,5.5,6,7 |
| working.dir | String | Working directory | current working directory |
| chromosomes | String | Chromosomes to run |  1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X |
| snv.mutations | String | URL of MAF file containing SNV mutations | NULL |
| indel.mutations | String | URL of MAF file containing indel mutations | NULL |
| mask.regions.file | String | URL of RDS file containing genomic regions to be masked during analysis e.g. non-mappable sites| mask_regions.RDS |
| all.sites.file | String | URL of RDS of all genomic sites | all_sites.RDS |
| region.of.interest | String | URL of BED file containing specific region of analysis e.g. promoter region | NULL |
| sample | Logical | To sample for non-mutated sites or not | TRUE |
| cores | Numeric | Number of cores to use during analysis | 1 |
| cutoff.nucleotide | Numeric | Frequency threshold for nucleotide context lasso selection | 0.90 |
| cutoff.nucleotide.new | Numeric | Adjusted frequency threshold for nucleotide context lasso selection | NULL |
| top.nucleotide | Numeric | Number of top nucleotide contexts to select | NULL |
| genomic.features.snv | String | URL of text file containing list of potential SNV genomic features | NULL |
| genomic.features.indel | String | URL of text file containing list of potential indel genomic features | NULL |
| genomic.features | String | URL of text file containing list of potential SNV and indel genomic features (same) | NULL |
| genomic.features.fixed.snv | String | URL of text file containing list of SNV genomic features that have to be in the model | NULL |
| genomic.features.fixed.indel | String | URL of text file containing list of indel genomic features that have to be in the model | NULL |
| genomic.features.fixed | String | URL of text file containing list of SNV and indel genomic features that have to be in the model (same) | NULL |
| sample.snv.features | String | URL of text file containing SNV sample specific features | NULL |
| sample.indel.features | String | URL of text file containing indel sample specific features | NULL |
| cutoff.features | Numeric | Frequency threshold for SNV and indel genomic feature lasso selection | 0.75 |
| cutoff.features.new.snv | Numeric | Adjusted frequency threshold for SNV genomic feature lasso selection | NULL |
| cutoff.features.new.indel | Numeric | Adjusted frequency threshold for indel genomic feature lasso selection | NULL |
| top.features | Numeric | Number of top genomic features to select | NULL |
| fit.sparse | Logical | To fit model using GLM or GLMNET | FALSE |
| drop | Logical | To drop variables that are not significant from the fitted model | FALSE |
| min.count.snv | Numeric | Minimum number of mutated samples in each SNV hotspot | 2 |
| min.count.indel | Numeric | Minimum number of mutated samples in each indel hotspot | 2 |
| hotspot.size | Numeric | Size of each hotspot | 21 |
| genome.size | Numeric | Genome size | 2533374732 |
| hotspots | Logical | To run hotspot analysis or region-based analysis | TRUE |
| collapse.regions | Logical | To collapse region of interest or not | FALSE |
| promoter.file | String | URL of BED file containing promoter regions | Ensembl75.promoters.coding.bed |
| utr3.file | String | URL of BED file containing 3'UTR regions | Ensembl75.3UTR.coding.bed |
| utr5.file | String | URL of BED file containing 5'UTR regions | Ensembl75.5UTR.coding.bed |
| other.annotations | String | URL of text file containing list of additional regions to be annotated | NULL |
| fdr.cutoff | Numeric | FDR cutoff | 0.05 |
| color.line | String | Cutoff line color in manhattan plot | red |
| color.dots | String | Color significant points in manhattan plot | maroon |
| merge.hotspots | Logical | To merge overlapping hotspots in figures or not | TRUE |
| color.muts | String | Color points in individual hotspot plot | orange |
| top.no | Numeric | Number of top hotspots to visualize individually | 3 |
| debug | Logical | To keep intermediate output files or not | FALSE |

----------------------------------------------------------------------------------

## Details of intermediate functions and file formats

| Step | run.to |
|------|--------|
| 1 | 1 |
| 2 | 2 |
| 3 | 3.1, 3.2 |
| 4 | 4.1, 4.2 |
| 5 | 5.1, 5.2, 5.3, 5.4, 5.5 |
| 6 | 6 |
| 7 | 7 |

------------------------------------------------------------------------------------

<a name="sample-sites"></a>

#### Step 1 - Sample sites

As the whole genome/region of interest may be too big, we sample for non-mutated sites across the whole genome/region of interest so as to make analysis more feasible. The resulting list of mutated sites and sampled non-mutated sites will be used as the positive and negative responses respectively in fitting the prediction model. If the number of mutations exceed 2,000,000, it downsamples the number of mutations before sampling an equal proportion of non-mutated sites. If the number of mutations is below 4000, it samples a higher proportion of non-mutated sites. The number of non-mutated sites sampled will always be at least 10,000.

```r
MutSpot(run.to = 1, snv.mutations = "subset_snv_mutations_sid.MAF", indel.mutations = "subset_indel_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt", sample.indel.features = "sample_features_table_indel.txt")
```


Output generated in Step 1:

| Filename | Description | Required fields |
|----------|-------------|-----------------|
| SNV_region.MAF | Lists positions of SNV mutations that are in the region of interest, outputs only when user-specified region provided | Chromosome, Start position, End position, Reference allele, Alternate allele, Sample name e.g. <br/> *chr1 3407955 3407955 A G HK-pfg160 <br/> chr1 5570109 5570109 A C TCGA-D7-6527 <br/> chr1 5570110 5570110 A C TCGA-F1-6875 <br/> chr1 5570127 5570127 A G HK-pfg060* |
| sampled.sites.snv.RDS | List of sampled non-mutated and mutated SNV sites in a GRanges object | Chromosome, Start position, End position, Mutation status (0 encodes non-mutation, 1 encodes mutation) e.g. <br/> *chr1 840142 840142 0 <br/> chr1 856644 856644 0 <br/> chrX 154069348 154069348 0 <br/> chrX 154069371 154069371 1* |
| indel_region.MAF | Lists positions of indel mutations that are in the region of interest, outputs only when user-specified region provided | Chromosome, Start position, End position, Reference allele, Alternate allele, Sample name e.g. <br/> *chr1 156426530 156426530 A AG HK-pfg222 <br/> chr12 95537679 95537683 CCTTT C HK-pfg038 <br/> chr14 61928665 61928665 C CATT HK-pfg160 <br/> chr15 95567861 95567865 GTTAC G HK-pfg030* |
| sampled.sites.indel.RDS | List of sampled non-mutated and mutated indel sites in a GRanges object | Chromosome, Start position, End position, Mutation status *(0 encodes non-mutation, 1 encodes mutation)* e.g. <br/> *chr1 1365965 1365965 0 <br/> chr1 12244480 12244480 0 <br/> chrX 103810639 103810639 0 <br/> chrX 109245627 109245627 1* |
| downsampled.sites.snv.RDS | Lists positions of downsampled SNV mutations in a GRanges object, outputs only when number of mutations exceeds threshold | Chromosome, Start Position, End position, Reference allele, Alternate allele, Sample name e.g. <br/> *chr1 24365487 24365487 T G TCGA-BR-1542 <br/> chr2 897987 897987 C T TCGA-D5-4567 <br/> chr2 54587 54587 A T TCGA-F1-6875 <br/> chr1 15458 15458 C G HK-pfg568* |
| downsampled.sites.indel.RDS | Lists positions of downsampled indel mutations in a GRanges object, outputs only when number of mutations exceeds threshold | Chromosome, Start Position, End position, Reference allele, Alternate allele, Sample name e.g. <br/> *chr1 265487 265488 TA G TCGA-BR-1542 <br/> chr2 98754 98754 C TTGCG TCGA-D5-4567 <br/> chr2 574578 574578 AT TCG TCGA-F1-6875 <br/> chr1 45878 45878 A TCAGT HK-pfg568* |

--------------------------------------------------------------------------------------

<a name="compute-local-mutation-rate"></a>

#### Step 2 - Compute local mutation rate

As local mutation rate is a potential feature used in the background model, we calculate 100kb binned mutation rates across the whole genome for SNV and indel separately.

```r
MutSpot(run.to = 2, snv.mutations = "subset_snv_mutations_sid.MAF", indel.mutations = "subset_indel_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt", sample.indel.features = "sample_features_table_indel.txt")
```


Output generated in Step 2:

| Filename | Description | Required fields |
|----------|-------------|-----------------|
| localmutrate_snv.bed | Lists 100kb binned mutation rates across whole genome | Chromosome, Start position, End position, Mutation rate e.g. <br/> *chr1 0 100000 2.632983e-07 <br/> chr1 100000 200000 2.632983e-07 <br/> chr1 200000 300000 2.632983e-07 <br/> chr1 300000 400000 2.632983e-07* |
| localmutrate_indel.bed | Lists 100kb binned mutation rates across whole genome | Chromosome, Start position, End position, Mutation rate e.g. <br/> *chr1 0 100000 0 <br/> chr1 100000 200000 0 <br/> chr1 200000 300000 0 <br/> chr1 300000 400000 0* |

-------------------------------------------------------------------------------------

<a name="select-nucleotide-contexts"></a>

#### Step 3 - Select nucleotide contexts

Part 1

Nucleotide/sequence contexts (e.g. 1mer, 3mer, 5mer) may be associated with mutations, hence we use LASSO to select for few nucleotide contexts that are most strongly associated with SNV mutations through a stability test. We only consider one in each reverse complement pairs. For indels, we filter away indels that lie in poly A/T/C/G (at least 8 repeated bases) regions.

An example:
**CTAGT**, where A is the mutated site

| Context name     | Context | Reverse complement | e.g.         |
| -----------------|---------|--------------------|--------------|
| 1mer             | A       | T                  | oneA         |
| 3mer             | TAG     | CTA                | threeTAG     |
| 3mer right flank | AG      | CT                 | threeRightAG |
| 3mer left flank  | TA      | TA                 | threeLeftTA  |
| 5mer             | CTAGT   | ACTAG              | fiveCTAGT    |
| 5mer right flank | AGT     | ACT                | fiveRightAGT |
| 5mer left flank  | CTA     | TAG                | fiveLeftCTA  |

```r
MutSpot(run.to = 3.1, snv.mutations = "subset_snv_mutations_sid.MAF", indel.mutations = "subset_indel_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt", sample.indel.features = "sample_features_table_indel.txt")
```


Output generated in Step 3 part 1:

| Filename | Description | Required fields |
|----------|-------------|-----------------|
| nucleotide_stabs_freq.RDS | List containing dataframe of frequency each nucleotide context is selected and dataframe of mean coefficients of each nucleotide context out of 100 LASSO runs| Context name, Frequency e.g. <br/> *oneA 1 <br/> threeAAA 0.05 <br/> threeAAC 0.01* <br/> Context name, Coefficient e.g. <br/> *oneA 0.548 <br/> threeAAA 0.0125* |
| nucleotide_selected.RDS | Character vector containing nucleotide contexts that passed the given frequency threshold | e.g. *oneA three.leftAA five.rightAGA oneG three.rightAG five.rightAGT* |
| indel_polyAT.MAF | Lists positions of indel mutations after removing indels that are in poly A/T/C/G regions | Chromosome, Start position, End position, Reference allele, Alternate allele, Sample name e.g. <br/> *chr1 156426530 156426530 A AG HK-pfg222 <br/> chr12 95537679 95537683 CCTTT C HK-pfg038 <br/> chr14 61928665 61928665 C CATT HK-pfg160 <br/> chr15 95567861 95567865 GTTAC G HK-pfg030* |


Part 2

For a more/less stringent threshold, re-define threshold (i.e. *cutoff.nucleotide.new*) and run part 2 of step 3 to adjust threshold for nucleotide context selection. We may further select the top features (i.e. *top.nucleotide*) based on their frequency and coefficients.

```r
MutSpot(run.to = 3.2, snv.mutations = "subset_snv_mutations_sid.MAF", indel.mutations = "subset_indel_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, cutoff.nucleotide.new = 0.98, top.nucleotide = 3, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt", sample.indel.features = "sample_features_table_indel.txt")
```


Output generated in Step 3 part 2:

| Filename | Description | Required fields |
|----------|-------------|-----------------|
| nucleotide_selected.RDS | Character vector containing nucleotide contexts that passed the new frequency threshold | e.g. *oneA three.leftAA five.rightAGA* |

------------------------------------------------------------------------------------------

<a name="select-genomic-features"></a>

#### Step 4

Part 1

We use LASSO to select genomic features that are most strongly associated with SNV/indel mutations through a stability test. Possible features include local mutation rate, replication timing, histone modification, transcription factors and DNA secondary structures. Features provided under *genomic.features.fixed.snv*/*genomic.features.fixed.indel*/*genomic.fetaures.fixed* do not undergo LASSO selection and will automatically be incorporated as features into the background model.

```r
MutSpot(run.to = 4.1, snv.mutations = "subset_snv_mutations_sid.MAF", indel.mutations = "subset_indel_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt", sample.indel.features = "sample_features_table_indel.txt")
```


Output generated in Step 4 part 1:

| Filename | Description | Required fields |
|----------|-------------|-----------------|
| features_stabs_snv.RDS | List containing dataframe of frequency each SNV genomic feature is selected and dataframe of mean coefficients of each feature out of 100 LASSO runs | Feature name, Frequency e.g. <br/> *mean_rep_time 0.82 <br/> local_mutrate 1 <br/> E094-DNase 0.02 <br/> E094-H3K27ac 1* <br/> Feature name, Coefficient e.g. <br/> *mean_rep_time 0.125 <br/> local_mutrate 16.54* |
| continuous_features_selected_snv_url.txt | Lists names and URLs of SNV continuous features that passed the given frequency threshold | Feature name, Filename of feature e.g. <br/> *mean_rep_time ./features/mean_rep_time.bed <br/> local_mutrate ./features/localmutrate_snv.bed* |
| discrete_features_selected_snv_url.txt | Lists names and URLs of SNV discrete features that passed the given frequency threshold | Feature name, Filename of feature e.g. <br/> *E094-H3K27ac ./features/E094-H3K27ac.bed <br/> POL2 ./features/POL2.bed* |
| features_stabs_indel.RDS | List containing dataframe of frequency each indel genomic feature is selected and dataframe of mean coefficients of each feature out of 100 LASSO runs | Feature name, Frequency e.g. <br/> *mean_rep_time 0 <br/> local_mutrate 0* <br/> Feature name, Coefficient e.g. <br/> *mean_rep_time 0 <br/> local_mutrate 0* |
| continuous_features_selected_indel_url.txt | Lists names and URLs of indel continuous features that passed the given frequency threshold | Feature name, Filename of feature |
| discrete_features_selected_indel_url.txt | Lists names and URLs of indel discrete features that passed the given frequency threshold | Feature name, Filename of feature |
| features_sds.RDS | List containing standard deviations of each genomic feature, separately for SNV and indel | Feature name, Standard deviation e.g. <br/> *mean_rep_time 0.012 <br/> local_mutrate 0.124* |



Part 2

For a more/less stringent threshold, re-define threshold (i.e. *cutoff.features.new.snv*, *cutoff.features.new.indel*) and run part 2 of step 4 to adjust threshold for genomic feature selection. We may further select the top features (i.e. *top.features*) based on their frequency and coefficients.

```r
MutSpot(run.to = 4.2, snv.mutations = "subset_snv_mutations_sid.MAF", indel.mutations = "subset_indel_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt", cutoff.features.new.snv = 0.8, cutoff.features.new.indel = 0.8, top.features = 11)
```


Output generated in Step 4 part 2:

| Filename | Description | Required fields |
|----------|-------------|-----------------|
| continuous_features_selected_snv_url.txt | Lists names and URLs of SNV continuous features that passed the given frequency threshold | Feature name, Filename of feature e.g. <br/> *mean_rep_time ./features/mean_rep_time.bed <br/> local_mutrate ./features/localmutrate_snv.bed* |
| discrete_features_selected_snv_url.txt | Lists names and URLs of SNV discrete features that passed the given frequency threshold | Feature name, Filename of feature e.g. <br/> *E094-H3K27ac ./features/E094-H3K27ac.bed <br/> POL2 ./features/POL2.bed* |
| continuous_features_selected_indel_url.txt | Lists names and URLs of indel continuous features that passed the given frequency threshold | Feature name, Filename of feature |
| discrete_features_selected_indel_url.txt | Lists names and URLs of indel discrete features that passed the given frequency threshold | Feature name, Filename of feature |

------------------------------------------------------------------------------------------------------

<a name="compute-feature-matrix"></a>

#### Step 5

Part 1

Compute feature matrix based on features that were selected in preceding steps for all sites in the whole genome/region of interest.

```r
# Chromosome 1
MutSpot(run.to = 5.1, chromosomes = "1", snv.mutations = "subset_snv_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt")


# Chromosome 13
MutSpot(run.to = 5.1, chromosomes = "13", snv.mutations = "subset_snv_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt")


# Chromosomes 1 to X
MutSpot(run.to = 5.1, snv.mutations = "subset_snv_mutations_sid.MAF", indel.mutations = "subset_indel_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt", sample.indel.features = "sample_features_table_indel.txt")
```


Output generated in Step 5 part 1:

For each chromosome, a list is saved in mutCovariate_chr*.RDS and indel_mutCovariate_chr*.RDS for SNV and indel respectively. The elements in the list are:

| Element | Required fields |
|---------|-----------------|
| Covariates matrix for mutated sites | Features columns (these columns may vary depending on the selected features ), Frequency (number of mutated sites having the features scores combination) e.g. <br/> *sid mean_rep_time local_mutrate E094-H3K27ac POL2 SMC3 ZNF143 oneMerA threeLeftAA fiveRightAGA freq <br/> 4771 61.85170 2.574273e-06 0 0 1 1 1 0 0 1 <br/> 6803 49.86691 3.829896e-06 0 0 1 1 0 0 0 1 <br/> 7050 29.20818 7.166831e-06 0 0 1 0 1 0 0 1 <br/> 7114 74.40081 3.065796e-06 1 0 1 1 1 1 0 1* |
| Covariates matrix for all sites | Features columns (these columns may vary depending on the selected features e.g. sid, local_mutrate, E094_H3K36me3, oneMerA)*, Frequency *(number of mutated sites having the features scores combination) e.g. <br/> *mean_rep_time local_mutrate E094-H3K27ac POL2 SMC3 ZNF143 oneMerA threeLeftAA fiveRightAGA x <br/> 43.48897 2.632983e-07 0 0 0 0 0 0 0 16 <br/> 67.89957 2.632983e-07 0 0 0 0 0 0 0 61 <br/> 74.40081 2.632983e-07 0 0 0 0 0 0 0 31 <br/> 49.86691 1.734790e-06 0 0 0 0 0 0 0 13* |



Part 2

```r
MutSpot(run.to = 5.2, snv.mutations = "subset_snv_mutations_sid.MAF", indel.mutations = "subset_indel_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt", sample.indel.features = "sample_features_table_indel.txt")
```


Output generated in Step 5 part 2:

| Filename | Description | Required fields |
|----------|-------------|-----------------|
| mutCovariate-compile-part1.RDS | Full SNV covariates matrix with mutation and non-mutation counts for all chromosomes in a dataframe | Features columns (these columns may vary depending on the selected features), Mutation count (number of mutated sites having the features scores combination), Total count (total number of sites having the features scores combination), Non-mutation count (number of non-mutated sites having the features scores combination) e.g. <br/> *sid mean_rep_time local_mutrate E094-H3K27ac POL2 SMC3 ZNF143 oneMerA threeLeftAA fiveRightAGA mut.count tot.count nonmut.count <br/> 1018 14.08371 2.632983e-07 0 0 0 0 0 0 0 0 33 33 <br/> 1018 14.08371 2.632983e-07 0 0 0 0 1 0 0 0 20 20 <br/> 1018 14.08371 2.632983e-07 0 0 0 0 1 0 1 0 1 1 <br/> 1018 14.08371 2.632983e-07 0 0 0 0 1 1 0 0 4 4* |



Part 3

Convert feature matrices to sparse matrices for SNV.

```r
MutSpot(run.to = 5.3, snv.mutations = "subset_snv_mutations_sid.MAF", indel.mutations = "subset_indel_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt", sample.indel.features = "sample_features_table_indel.txt")
```


Output generated in Step 5 part 3:

| Filename | Description | Required fields |
|----------|-------------|-----------------|
| mutCovariate-sparse-p1.RDS | SNV covariates sparse matrix | Features columns (these columns may vary depending on the selected features) e.g. <br/> *1018 14.08371 2.632983e-07 . . . . . . . <br/> 1018 14.08371 2.632983e-07 . . . . 1 . . <br/> 1018 14.08371 2.632983e-07 . . . . 1 . 1 <br/> 1018 14.08371 2.632983e-07 . . . . 1 1 .* |
| mutCovariate-sparse-p2.RDS | SNV response matrix | Mutation count, Non-mutation count e.g. <br/> *mut.count nonmut.count <br/> 0 33 <br/> 0 20 <br/> 0 1 <br/> 0 4* |



Part 4

Compile indel feature matrices of all chromosomes.

```r
MutSpot(run.to = 5.4, snv.mutations = "subset_snv_mutations_sid.MAF", indel.mutations = "subset_indel_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt", sample.indel.features = "sample_features_table_indel.txt")
```


Output generated in Step 5 part 4:

| Filename | Description | Required fields |
|----------|-------------|-----------------|
| mutCovariate-indel-compile-part1.RDS | Full indel covariates matrix with mutation and non-mutation counts for all chromosomes in a dataframe | Features columns (these columns may vary depending on the selected features), Mutation count (number of mutated sites having the features scores combination), Total count (total number of sites having the features scores combination), Total length (length of indels having the features scores combination), Non-mutation count (number of non-mutated sites having the features scores combination) |



Part 5

Convert feature matrices to sparse matrices for indel.

```r
MutSpot(run.to = 5.5, snv.mutations = "subset_snv_mutations_sid.MAF", indel.mutations = "subset_indel_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt", sample.indle.features = "sample_features_table_indel.txt")
```

Output generated in Step 5 part 5:

| Filename | Description | Required fields |
|----------|-------------|-----------------|
| indel-mutCovariate-sparse-p1.RDS | Indel covariates sparse matrix | Features columns (these columns may vary depending on the selected features) |
| indel-mutCovariate-sparse-p2.RDS | Indel response matrix | Mutation count, Non-mutation count |

-------------------------------------------------------------------------------------------------------

<a name="fit-prediction-model"></a>

#### Step 6

Fit logistic regression model based on covariates matrix from preceding step. This step plots the feature importance plot of the fitted background model.

```r
MutSpot(run.to = 6, snv.mutations = "subset_snv_mutations_sid.MAF", indel.mutations = "subset_indel_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt", sample.indel.features = "sample_features_table_indel.txt")
```


Output generated in Step 6:

| Filename | Description | Required fields |
|----------|-------------|-----------------|
| snv-LRmodel | Fitted SNV prediction model | e.g. *sid 5.211938e-05 <br/> mean_rep_time -1.020669e-02 <br/> local_mutrate 1.335094e+05 <br/> -4.706059e-01 <br/> POL2 -4.112153e-01 <br/> SMC3 9.520664e-01 <br/> ZNF143 8.923972e-01 <br/> oneMerA 1.532452e+00 <br/> threeLeftAA 8.528978e-01 <br/> fiveRightAGA 3.414982e-01* |
| SNV_features.pdf | Feature importance barplot of fitted model | |
| indel-LRmodel | Fitted indel prediction model | | |
| indel_features.pdf | Feature importance barplot of fitted model | |

----------------------------------------------------------------------------------------------------

<a name="predict-hotspots"></a>

#### Step 7

Predict hotspot mutations using fitted prediction models from preceding step. This step annotates hotspots and outputs the manhattan plot and mutated sample count/position plot for each significant hotspot.

```r
MutSpot(run.to = 7, snv.mutations = "subset_snv_mutations_sid.MAF", indel.mutations = "subset_indel_mutations_sid.MAF", region.of.interest = "gastric_ctcf_motif.bed", cores = 2, genomic.features = "genomic_features_ctcf.txt", sample.snv.features = "sample_features_table_snv.txt", sample.indel.features = "sample_features_table_indel.txt")
```


Output generated in Step 7:

| Filename | Description | Required fields |
|----------|-------------|-----------------|
| snv_hotspots.tsv | Predicted SNV hotspot mutation recurrence | Chromosome, Start position, End position, P-value, Length of hotspot, Background p-value, Number of mutated samples, FDR e.g. <br/> *chrom start end pval length p.bg k fdr <br/> Mutation1378 chr8 71000992 71001012 7.718045e-10 21 0.0009075422 6 0.0009869712 <br/> Mutation459 chr14 70285576 70285596 2.693087e-08 21 0.0007888790 5 0.0094519870 <br/> Mutation460 chr14 70285578 70285598 2.956555e-08 21 0.0008039330 5 0.0094519870 <br/> Mutation461 chr14 70285579 70285599 2.956555e-08 21 0.0008039330 5 0.0094519870* |
| snv_hotspots_merged.tsv | Merged predicted SNV hotspot mutation recurrence | Chromosome, Start position, End position, P-value, Length of hotspot, Background p-value, Number of mutated samples, FDR e.g. <br/> *chrom start end pval length p.bg k fdr <br/> Mutation1378 chr8 71000900 71001012 7.718045e-10 21 0.0009075422 8 0.0009869712* |
| snv_hotspots_annotated.tsv | SNV hotspots with annotated regions | Chromosome, Start position, End position, P-value, Length of hotspot, Background of p-value, Number of mutated samples, FDR, Promoter, 3'UTR, 5'UTR e.g. <br/> *chrom start end pval length p.bg k fdr promoter 3UTR 5UTR <br/> Mutation1378 chr8 71000992 71001012 7.718045e-10 21 0.0009075422 6 0.0009869712 NA NA NA* |
| snv_manhattan.pdf | Each point represents a hotspot. Colored ones are significant. Line represents FDR cutoff. Y-axis shows log p-value and x-axis shows genomic position | |
| snv_hotspot_*.pdf | Each point represents a hotspot. Colored ones are significant. Line represents FDR cutoff. Y-axis shows log p-value and x-axis shows genomic position | |
| indel_hotspots.tsv | Predicted indel hotspot mutation recurrence | Chromosome, Start position, End position, P-value, Length of hotspot, Background p-value, Number of mutated samples, FDR |
| indel_hotspots_merged.tsv | Predicted indel hotspot mutation recurrence | Chromosome, Start position, End position, P-value, Length of hotspot, Background p-value, Number of mutated samples, FDR |
| indel_hotspots_annotated.tsv | Indel hotspots with annotated regions | Chromosome, Start position, End position, P-value, Length of hotspot, Background p-value, Number of mutated samples, FDR, Promoter, 3'UTR, 5'UTR |
| indel_manhattan.pdf | Each point represents a hotspot. Colored ones are significant. Line represents FDR cutoff. Y-axis shows log p-value and x-axis shows genomic position | | |
| indel_hotspot_*.pdf | Each point represents a hotspot. Colored ones are significant. Line represents FDR cutoff. Y-axis shows log p-value and x-axis shows genomic position | | |

----------------------------------------------------------------------------------------------------
