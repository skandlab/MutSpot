###
### Workflow for hotspot analysis
###
1.1 Calculate local mutation rate in 100kb non-overlapping bins
	mask: CDS regions, nonmappable regions, ig loci
	calculate local mutation rate in 100kb bins -> write precalculated rate in bed or wg format
local_mut_rate.R


1.2 discretize replication timing into n bins and calculate the mean value for each bin
RepTime_binMeans.R	

1.3 sample sites for feature selection
sample_sites.R

2. create feature matrix with epigenetic features
mutrate_gastric_site_current.R

3. feature selection for regional and site models
lasso_feature_selection.R

4. nucleotide context feature selection
MotifLassoSelection.R

5 create frequency table for logistic regression
mutCovariate_hotspot.R

6. fit logistic regression
mutCovariate_LRfit_hotspot.R

7. calculate mutation recurrence for regions of interest from predicted background mutation rates
mutrec_logistic_hotspot.R
