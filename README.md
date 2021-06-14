
# Data appendix for Bland & Rosokha (2021) "Learning Under Uncertainty with Multiple Priors: Experimental Investigation"

This data appendix includes all data and code (implimented in Matlab) needed to replicate the figures and tables in this paper. The files are described below, organized into several classes:

* Data files: data from the experiment
* Estimation files: Code used to estimate the models (posterior simulations)
* Postestimation files: Code used to transform posterior simulations into plots, tables, and numbers mentioned in the text.
* Other files: Stand-alone plots and tables, etc.

You will need to change some variables called "fpath" in order to point Matlab in the direction of folders that exist on your machine. The posterior simulation files are somewhat large (up to about 500 Mb), so it may not be advantagous to (for example) not set this to a folder that is synchronized with Github.

Set your workibng directory to this folder. All tables and figures will be output into the child folder "figures". 

## Data files

* Exp1.mat - Data from the main experiment
* Expt2.mat - Data from Moreno & Rosokha (2016)

## Estimation files

 * est_MIX_BetaPriors.m - Beta priors mixture model reported in:
    + Table 1(b)
    + Table 2
    + Table 3
    + Figure 6
    + Figure 7
    + Figure C-3
  
* est_MIX_SimplexPriors_general.m - simplex priors mixture model reported in  Reported in Table 1(a)
  
* est_MIX_BetaPriorsBoth.m - Beta priors mixture model incorporating data from Bland & Rosokha (2021) and Moreno & Rosokha (2016) reported in Tables B-1(a) and B-2

* est_MIX_BetaPriorsRestricted.m - Beta priors mixture model restricting (truncating) beta priors to the interval (0.25, 0.75). Reported in Tables B-1(b) and B-3

* est_MIX_BetaPriorsLambdaConstant.m - Beta priors mixture model restricting the choice precision parameter lambda to be constant across subjects. Reported in Figure 7 and Table C-4
  
### Dependencies of estimation files

**Likelihood functions**

  * A_loglike_BAYES_Beta.m: Log-likelihood function for Bayesian types with Beta priors (A- and C-tasks)
  * R_loglike2.m: Log-likelihood for the R-task
  * A_loglike_MPBeta_black.m: Log-likelihood for MP-type with Beta priors
  * A_loglike_BAYES_Simplex: log-likelihood for Bayesian simplex priors type
  * A_loglike_MP_Simplex3.m: log-likelihood for MP simplex type
  * A_loglike_BAYES_BetaRestricted.m: Log-likelihood of Bayesian beta prior type, truncating prior to (0.25,0.75)
  * A_loglike_MPBetaRestricted_black.m: Log-likelihood of MP beta prior type, truncating prior to (0.25,0.75)
  
**Other** 

  * discretesampleJB.m: Samples from the categorical distribution. 
 
## Postestimation files

Run the estimation files first.

* MIX_postestimation_CompareLambdaConstant.m: Compares the results of the unrestricted Beta priors mixture model (Table 1) and the model imposing the restrion that choice precision (lambda) is constant across subjects. Produces Figure 7

* MIX_postestimation.m: Main file for producing estimates summaries and some parameter plots. Make sure to adjust the file path to the place where you have saves your simulation results. 

* Mix_post_compareMix.m: Produces tables like Table 1 and Table 2.

* MIX_comparison.m: Model evaluation and selection
 
## Other files
 
 * SimplexExamplePlots.m - Produces plots for Figures 4 and 5.
 * whitespace.m: removes the white space around plots
 * plotscript.m: needed to produce the plot showing the individual-level parameters (e.g Figure 6)
