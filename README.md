# DPTM

Dynamic Panel Multiple Threshold Model with Fixed Effects

Q1: Why you need to use it?

1.  The DPTM can supply the estimation and test for the dynamic panel threshold model with multiple thresholds while there is no software or tools can doing this. 
2.  The DPTM is based on maximum likelihood (ML) estimation while most of dynamic panel threshold model used IVs or GMM, which brings a better performance and a more convenient in use.
3.  The DPTM is based on Markov Chain Monte Carlo (MCMC) method while most of threshold model are based grid search method, which brings a better performance and avoids the selection of step size of grid.
4.  The DPTM not only allows fixed effects, but also allows time trend term or time fixed effects, which is more suitable for reality and application.
5.  Various forms of threshold models are available。

Q2: How can you use it?

Please see the specific example in Help Pages of functions in the DPTM.

Q3: How can you get it?

1.  Use "devtools::install_github("HujieBai/DPTM")" for installation.

2.  The R package available for local installation is "DPTM_x.x.x.tar.gz."
