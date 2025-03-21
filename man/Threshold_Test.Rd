
\name{Threshold_Test}
\alias{Threshold_Test}
\title{Tests for multiple thresholds.}
\usage{
Threshold_Test(formula = NULL, formula_cv = NULL, data, index=NULL, Th = 1, q, 
timeFE = FALSE, bt = 100,NoY = FALSE, y1 = NULL, iterations = 2000, sro = 0.1,
r0x = NULL, r1x = NULL, parallel=TRUE, seed = NULL,...)
}
\arguments{
\item{formula}{formula of the covariates with threshold effects; If a setting is not provided, defaults (no covariates with threshold effects) will be used. Defaults to `NULL`.}

\item{formula_cv}{formula of the covariates without threshold effects; If a setting is not provided, defaults (no covariates without threshold effects) will be used. Defaults to `NULL`.}

\item{data}{data frame of the observed data.}

\item{index}{variable names of individuals and period; If a setting is not provided, defaults (the first variables in data will be as "id", while the second will be "year") will be used.Defaults to `NULL`.}

\item{q}{threshold variable.}

\item{timeFE}{logicals. If TRUE the time fixed effects will be allowed. Defaults to `FALSE`.}

\item{bt}{the number of bootstrap; If a setting is not provided, defaults (bt = 100) will be used. Defaults to `100`.}

\item{NoY}{logicals. If TRUE the lags of dependent variables will be without threshold effects. Defaults to `FALSE`.}

\item{y1}{lags of dependent variables; If a setting is not provided, defaults (the first-order lag) will be used. Defaults to `NULL`.}

\item{iterations}{MCMC iterations (50\% used for burnining). Defaults to `2000`.}

\item{sro}{regime (subsample) proportion; If a setting is not provided, defaults (10\%) will be used. Defaults to `0.1`.}

\item{r0x}{lower bound of threshold parameter space; If a setting is not provided, defaults (15\% quantile of threshold variable) will be used.}

\item{r1x}{upper bound of threshold parameter space; If a setting is not provided, defaults (85\% quantile of threshold variable) will be used.}

\item{parallel}{logicals. If TRUE test will run in parallel for saving time. Defaults to `TRUE`.}

\item{seed}{set seeds to guarantee the replication the test (see set.seed);}

\item{...}{additional arguments to be passed to the settings of MCMC (see BayesianTools::applySettingsDefault)}

\item{Th}{number of thresholds; Defaults to `1`.}
}
\value{
A list with class "htest" containing the following components:
\item{statistic}{the value of the F-statistic.}
\item{parameter}{the degrees of freedom for the F-statistic.}
\item{p.value}{the p-value for the test.}
\item{null.value}{the specified hypothesized value of the null hypothesis.}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{method}{a character string indicating what type of test was performed.}
\item{data.name}{a character string giving the name(s) of the data.}
\item{estimate}{the critical value of the statistic (5\% significance level).}
\item{LRs}{a vector of statistics from bootstrap.}
}
\description{
Tests for models with different thresholds, using bootstrap method.
}
\details{
\code{Threshold_Test} can run the Test for multiple thresholds (\code{Th} is H1).
The statistic is
\deqn{F_s=\frac{S\left(\hat{\gamma}_{s-1}\right)-S\left(\hat{\gamma}_s\right)}{S\left(\hat{\gamma}_s\right) / N(T-1)},}
where \eqn{s} is the number of thresholds in H1, \eqn{S\left(\hat{\gamma}_{s-1}\right)=-\ln L\left(\hat{\gamma}_{s-1}\right)} and \eqn{S\left(\hat{\gamma}_s\right)=-\ln L\left(\left(\hat{\gamma}_{s-1}^{\prime}, \hat{\gamma}_s\right)^{\prime}\right)}.
And the p-value is computed by bootstrap method (see Ramírez-Rondán, 2020).

Take the two threshold model as example.
User must set \code{Th} = 1 firstly to reject the null hypothesis of no threshold effects;
Then he should set \code{Th} = 2 to reject the null hypothesis of only one threshold;
Lastly, set \code{Th} = 3 to accept the null hypothesis of two thresholds.
In other words, p-values of the first test (\code{Th} = 1) and the second test (\code{Th} = 1) should be less than significant level while the third test (\code{Th} = 3) is not.

\code{Threshold_Test} contains all augments in \code{DPTS}, but with three new augments: \code{bt}, \code{parallel} and \code{seed}.
\code{bt} is the number of bootstrap (by default is 100);
\code{parallel} can allow user to run test in parallel to save time;
\code{seed} is used to guarantee the replication of tests.

It is worthy noting that the test shrinks to the so-called threshold existence test when \code{Th} = 1.
}
\examples{
\donttest{
### Examples elapsed time > 15s

#data(d1)

# H0: no threshold effects (no threshold)
#test0 <- Threshold_Test(y~x,y~z,data = d1, index = c('id','year'), q = d1$q, Th = 1,
#bt = 50, iterations = 500)
#test0

# H0: one threshold 
#test1 <- Threshold_Test(y~x,y~z,data = d1, index = c('id','year'), q = d1$q, Th = 2,
#bt = 50, iterations = 500)
#test1

# H0: two threshold 
#test2 <- Threshold_Test(y~x,y~z,data = d1, index = c('id','year'), q = d1$q, Th = 3,
#bt = 50, iterations = 500)
#test2

}


}
\references{
Ramírez-Rondán, N. R. (2020). Maximum likelihood estimation of dynamic panel threshold models. Econometric Reviews, 39(3), 260-276.
}
\author{
Hujie Bai
}
