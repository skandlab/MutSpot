###
### Workflow for hotspot analysis
###
1.1 Calculate local mutation rate in 100kb non-overlapping bins<br />
mask: CDS regions, nonmappable regions, ig loci<br />
calculate local mutation rate in 100kb bins -> write precalculated rate in bed or wg format<br />
local_mut_rate.R<br /><br />
1.2 discretize replication timing into n bins and calculate the mean value for each bin<br />
RepTime_binMeans.R<br /><br />	
1.3 sample sites for feature selection<br />
sample_sites.R<br /><br />
2. create feature matrix with epigenetic features<br />
mutrate_gastric_site_current.R<br /><br />
3. feature selection for regional and site models<br />
lasso_feature_selection.R<br /><br />
4. nucleotide context feature selection<br />
MotifLassoSelection.R<br /><br />
5 create frequency table for logistic regression<br />
mutCovariate_hotspot.R<br /><br />
6. fit logistic regression<br />
mutCovariate_LRfit_hotspot.R<br /><br />
7. calculate mutation recurrence for regions of interest from predicted background mutation rates<br />
mutrec_logistic_hotspot.R<br /><br />
