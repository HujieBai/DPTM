# DPTM

Dynamic Panel Multiple Threshold Model with Fixed Effects

Q1: Why do you need to use it?

1.  The DPTM can supply the estimation and test for the dynamic panel threshold model with multiple thresholds while there is no software or tools that can do this. 
2.  The DPTM is based on maximum likelihood (ML) estimation, while most dynamic panel threshold models use IVs or GMM, which improves performance and is more convenient to use.
3.  The DPTM is based on the Markov Chain Monte Carlo (MCMC) method, while most threshold models use a grid search method, which improves performance and avoids selecting the grid's step size.
4.  The DPTM not only allows fixed effects but also allows the time fixed effects, which is more suitable for reality and application.
5.  Various forms of threshold models are available。

Q2: How can you use it?

Please see the specific example in the Help Pages of functions in the DPTM.

Q3: How can you get it?

1.  Use "devtools::install_github("HujieBai/DPTM")" for installation.

2.  The R package available for local installation is "DPTM_x.x.x.tar.gz."

#2025.3.22 (Version 3.0.2)

1. Rewrote the entire package in OOP (R6 class).

2. Detailed instructions and examples have been added to make it more user-friendly.

3. The formula interface of the R language is exploited, which is more standardized and easy to use.

Please use?DPTM::DPTS to see more useful details.

(I'm a PhD candidate at Lingnan College of Sun Yat-sen University, thank you very much for your recognition of my work!)
