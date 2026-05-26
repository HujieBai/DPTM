
\name{DPTS}
\alias{DPTS}
\title{Dynamic panel multiple threshold model with fixed effects.}
\usage{
DPTS(formula = NULL, formula_cv = NULL, data, index=NULL, Th = 1, q, timeFE = FALSE, 
NoY = FALSE, y1 = NULL, iterations = 2000, sro = 0.1, r0x = NULL, r1x = NULL,
grid_search = FALSE, grids = 100, grid_search_type = c("jointly", "sequential"),
grid_search_iter = 1,
...)

## S6 method for class 'DPTM' 
#print(...)
}
\arguments{
\item{formula}{formula of the covariates with threshold effects; If a setting is not provided, defaults (no covariates with threshold effects) will be used. Defaults to `NULL`.}

\item{formula_cv}{formula of the covariates without threshold effects; If a setting is not provided, defaults (no covariates without threshold effects) will be used. Defaults to `NULL`.}

\item{data}{data frame of the observed data.}

\item{index}{variable names of individuals and period; If a setting is not provided, defaults (the first variables in data will be as "id", while the second will be "year") will be used.Defaults to `NULL`.}

\item{Th}{number of thresholds; Defaults to `1`.}

\item{q}{threshold variable.}

\item{timeFE}{logicals. If TRUE the time fixed effects will be allowed. Defaults to `FALSE`.}

\item{NoY}{logicals. If TRUE the lags of dependent variables will be without threshold effects. Defaults to `FALSE`.}

\item{y1}{lags of dependent variables; If a setting is not provided, defaults (the first-order lag) will be used. Defaults to `NULL`.}

\item{iterations}{MCMC iterations (50\% used for burnining). Defaults to `2000`.}

\item{sro}{regime (subsample) proportion; If a setting is not provided, defaults (10\%) will be used. Defaults to `0.1`.}

\item{r0x}{lower bound of threshold parameter space; If a setting is not provided, defaults (15\% quantile of threshold variable) will be used.}

\item{r1x}{upper bound of threshold parameter space; If a setting is not provided, defaults (85\% quantile of threshold variable) will be used.}

\item{grid_search}{logicals. If TRUE, thresholds are searched by deterministic grid search instead of MCMC. Defaults to `FALSE`.}

\item{grids}{number of grid points (threshold vectors) searched when \code{grid_search = TRUE}. Defaults to `100`.}

\item{grid_search_type}{type of grid search for multiple thresholds. \code{"jointly"} searches the full Cartesian grid (\code{grids^Th}); \code{"sequential"} adds thresholds one at a time (\code{grids * Th}). Defaults to \code{"jointly"}.}

\item{grid_search_iter}{number of coordinate-refinement cycles after the initial grid search when \code{grid_search = TRUE}. For two thresholds, one cycle re-estimates the first threshold conditional on the second and then re-estimates the second conditional on the first. Set to \code{0} to skip this refinement. Defaults to \code{1}.}

\item{...}{additional arguments to be passed to the settings of MCMC (see BayesianTools::applySettingsDefault)}
}
\value{
DPTS returns an object of class "DPTM".
The function \code{print} are used to obtain and print a print of the results.   
An object of class "DPTM" is a list containing at least the following components:
\item{coefficients}{a named vector of coefficients}
\item{NNLL}{the negative log-likelihood function value}
\item{Zvalues}{a vector of t statistics}
\item{Ses}{a vector of standard errors}
\item{covariance_matrix}{a covariance matrix}
\item{Th}{the number of thresholds}
\item{thresholds}{a named vector of thresholds}
\item{threshold_search}{threshold-search method, either \code{"MCMC"} or \code{"grid"}}
\item{grid_points}{number of searched grid points when \code{grid_search = TRUE}}
\item{grid_search_type}{grid-search type when \code{grid_search = TRUE}}
\item{grid_search_iter}{number of coordinate-refinement cycles when \code{grid_search = TRUE}}
}
\description{
Use a MCMC-MLE or grid-search-MLE based two-step procedure to estimate the dynamic panel multiple threshold model with fixed effects.
}
\details{
\code{DPTS} can fit the dynamic panel threshold model with fixed effects proposed 
by Ramírez-Rondán (2020), and also allow a multiple threshold model by setting 
\code{Th} > 1.

Given the diverse forms and versatile applications of threshold models, we advocate for aligning model selection with specific research objectives, thereby granting users autonomy in specifying the model structure.

Take the model with one threshold (Ramírez-Rondán, 2020) as example.

For a standard threshold model
\deqn{\begin{aligned}y_{i t} &=\left(\rho_1 y_{i t-1}+\beta_1 x_{i t}\right) I\left(q_{i t}\leq \gamma\right)+\left(\rho_2 y_{i t-1}+\beta_2 x_{i t}\right) I\left(q_{i t}> \gamma\right)\\&+\alpha_i+u_{i t},\end{aligned}},
can use \code{DPTS(y~x,data = data, q = q, Th = 1)}.

For a threshold model who has regressors with threshold effects (\eqn{x}) and regressors without threshold effects (\eqn{z})
\deqn{\begin{aligned}y_{i t} &=\left(\rho_1 y_{i t-1}+\beta_1 x_{i t}\right) I\left(q_{i t}\leq \gamma\right)+\left(\rho_2 y_{i t-1}+\beta_2 x_{i t}\right) I\left(q_{i t}> \gamma\right)\\&+\theta z_{i t}+\alpha_i+u_{i t},\end{aligned}},
can use \code{DPTS(y~x,y~z,data = data, q = q, Th = 1)}.


If user only cares about the regressors with threshold effects (thus hopes there is no threshold effects in the lag of dependent variable \eqn{y_1}), like
\deqn{\begin{aligned}y_{i t} &= \rho y_{i t-1}+ \beta_1 x_{i t} I\left(q_{i t}\leq \gamma\right)+\beta_2 x_{i t} I\left(q_{i t}> \gamma\right)\\&+\theta z_{i t}+\alpha_i+u_{i t},\end{aligned}}, 
can use \code{DPTS(y~x,y~z,data = data, q = q, Th = 1, NoY = TRUE)}.

And, the threshold model with the following form 
\deqn{\begin{aligned}y_{i t} &=\rho_1 y_{i t-1}I\left(q_{i t}\leq \gamma\right)+\rho_2 y_{i t-1}I\left(q_{i t}> \gamma\right)+\beta x_{i t}\\&+\theta z_{i t}+\alpha_i+u_{i t},\end{aligned}},
is also allowed by \code{DPTS(NULL,y~x+z,data = data, q = q, Th = 1)}.

In addition, a special threshold model having the following form 
\deqn{\begin{aligned}\Delta y_{i t} &=\left(\rho_1 y_{i t-1}+\beta_1 x_{i t}\right) I\left(q_{i t}\leq \gamma\right)+\left(\rho_2 y_{i t-1}+\beta_2 x_{i t}\right) I\left(q_{i t}> \gamma\right)\\&+\theta z_{i t}+\alpha_i+u_{i t},\end{aligned}},
can use \code{DPTS(dy~x,dy~z,data = data, q = q, Th = 1)} with \code{y1}\eqn{= y_{it-1}}.

The MCMC we used is based on \pkg{BayesianTools}, and the default method is "DREAMzs" (see Vrugt et al., 2009).
If user wants to use other MCMC, can use \code{...} (see BayesianTools::applySettingsDefault).
As for the length of iterations, it can be set by \code{iterations} (50\% used for burnining) and default length is 2000.
The trace plot and the Gelman and Rubin's convergence diagnostic are supplied by \code{DPTS} (\code{print}) to test the convergence of MCMC sample.

Alternatively, set \code{grid_search = TRUE} to use deterministic grid search for the unknown thresholds.
In that case, \code{grids} controls the number of threshold vectors evaluated; for example, \code{grids = 100} searches 100 candidate points in the ordered threshold-parameter space.
For \code{Th > 1} and \code{grid_search_type = "jointly"}, the full Cartesian grid is searched, so \code{Th = 2} and \code{grids = 100} evaluates \code{100 * 100} candidate threshold pairs.
For \code{grid_search_type = "sequential"}, thresholds are added one at a time, so \code{Th = 2} and \code{grids = 100} evaluates about \code{200} candidates in total.
When \code{Th > 1}, \code{grid_search_iter} controls an additional coordinate-refinement step after the initial grid search. With two thresholds and \code{grid_search_iter = 1}, \code{DPTS} re-estimates the first threshold conditional on the second, then re-estimates the second conditional on the first, and sorts the final threshold estimates.
Grid search does not produce MCMC threshold confidence intervals or Gelman-Rubin diagnostics.


Additionally, we assume the exogenous regressor \eqn{x} is weakly exogenous, and 
thus the first period after difference is given by
\deqn{\Delta y_{i1}=\delta_0 + {\boldsymbol\delta}'_1 \Delta {\bf x}_{i1}+ v_{i1},}
where \eqn{E(v_{i1}| \Delta {\bf x}_{i1} )=0}. \eqn{E(v_{i1}^2)=\sigma^2_v}, 
\eqn{E(v_{i1}\Delta u_{i2})=-\sigma^2_u} and \eqn{E(v_{i1} \Delta u_{it})=0}
for \eqn{t=3,4,...,T} and \eqn{i=1,...,N}.
For more details, see Hsiao et al. (2002).

Finally, we solve the log-likelihood function by \code{stats::nlm} who uses \code{iterlim}
to set the maximum number of iterations, and thus \code{iterlim} is allowed by \code{...} in \code{DPTS}.
}
\examples{
data(d1)

# single threshold

## standard form 
#Model1_1 <- DPTS(y~x,data = d1, index = c('id','year'), q = d1$q, Th = 1, 
#iterations = 1000)
#print(Model1_1)

## deterministic grid search
#Model1_grid <- DPTS(y~x,data = d1, index = c('id','year'), q = d1$q, Th = 1,
#grid_search = TRUE, grids = 100)
#print(Model1_grid)

\donttest{
### Examples elapsed time > 15s
## with x \& z
#Model2_1 <- DPTS(y~x,y~z,data = d1, index = c('id','year'), q = d1$q, Th = 1, 
#iterations = 1000)
#print(Model2_1)

## with x \& z (y1 no threshold effects)
#Model3_1 <- DPTS(y~x,y~z,data = d1, index = c('id','year'), q = d1$q, Th = 1,
#NoY = TRUE, iterations = 1000)
#print(Model3_1)

## only y1 with threshold effects
#Model4_1 <- DPTS(NULL,y~x+z,data = d1, index = c('id','year'), q = d1$q, Th = 1, 
#iterations = 1000)
#print(Model4_1)

# two thresholds (Th = 2)
## with x \& z
#Model2_2 <- DPTS(y~x,y~z,data = d1, index = c('id','year'), q = d1$q, Th = 2, 
#iterations = 1000)
#print(Model2_2)

# Adding time fixed effects (timeFE = TRUE)
#Model2_2T <- DPTS(y~x,y~z,data = d1, index = c('id','year'), q = d1$q, Th = 2,
#timeFE = TRUE, iterations = 1000)
#print(Model2_2T)
}
}
\references{
Ramírez-Rondán, N. R. (2020). Maximum likelihood estimation of dynamic panel threshold models. Econometric Reviews, 39(3), 260-276.

Vrugt, Jasper A., et al. (2009)."Accelerating Markov chain Monte Carlo simulation by differential evolution with self-adaptive randomized subspace sampling." International Journal of Nonlinear Sciences and Numerical Simulation 10.3: 273-290.

Hsiao, C., Pesaran, M. H., & Tahmiscioglu, A. K. (2002).
Maximum likelihood estimation of fixed effects dynamic panel data models covering short time periods. Journal of econometrics, 109(1), 107-150.
}
\author{
Hujie Bai
}
