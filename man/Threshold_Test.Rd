% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tests.R
\name{Threshold_Test}
\alias{Threshold_Test}
\title{The test for the number of thresholds.}
\usage{
Threshold_Test(
  y,
  y1 = NULL,
  x = NULL,
  q,
  cvs = NULL,
  time_trend = FALSE,
  time_fix_effects = FALSE,
  x1 = NULL,
  tt,
  nn,
  Th = 0,
  ms = 1000,
  burnin = 1000,
  types = "DREAMzs",
  ADs = FALSE,
  r0x = NULL,
  r1x = NULL,
  NoY = FALSE,
  restart = FALSE,
  Only_b = FALSE,
  w = NULL,
  var_u = NULL,
  nCR = 3,
  autoburnin = TRUE,
  bt = 100,
  parallel = TRUE,
  sro = 0.1,
  display = TRUE
)
}
\arguments{
\item{y}{the dependent variable; vector type input.}

\item{y1}{the lag dependent variable; vector type input; By default, y1 is NULL,
and then y1 will be computed by y automatically.}

\item{x}{the independent variable; matrix type input.}

\item{q}{the threshold variable; vector type input.}

\item{cvs}{the set of control variables; matrix type input;By default, cvs is NULL.}

\item{time_trend}{the time trend; By default, it is FALSE.}

\item{time_fix_effects}{the time fixed effects; By default, it is FALSE.}

\item{x1}{the initial values of independent variable; matrix type input.
By default, x1 is NULL, and thus x1 will be computed by x automatically.}

\item{tt}{the length of time period.}

\item{nn}{the number of individuals.}

\item{Th}{the number of thresholds.}

\item{ms}{the length of MCMC chains after burn-in.}

\item{burnin}{the length of burn-in.}

\item{types}{the type of MCMC used; More details see BayesianTools::runMCMC.}

\item{ADs}{the options for MCMC; More details see BayesianTools::runMCMC.}

\item{r0x}{the lower bound of thresholds; By default, r0x is NULL,
and thus r0x will be computed by q automatically.}

\item{r1x}{the upper bound of thresholds; By default, r0x is NULL,
and thus r1x will be computed by q automatically.}

\item{NoY}{the option of threshold effects on the lag dependent variable;
By default, NoY is False, and thus there will be threshold effects on y1.}

\item{restart}{the option of iterations; By default, restart is FALSE,
if encounters iteration failure, please set restart as TRUE.}

\item{Only_b}{the option of initial equation;By default, Only_b is FALSE, and if Only_b is TRUE, initial delta y will be a constant C.
Please see Hsiao (2002) and Ramírez-Rondán (2020) for more details.}

\item{w}{the variance ratio; By default, is NULL; It must be greater than 1.}

\item{var_u}{the option of variance of error term; By default, is NULL; It must be
greater than 0; When meet relevant ERROR, please change the var_u.}

\item{nCR}{parameter determining the number of cross-over proposals of DREAM MCMC. If nCR = 1
all parameters are updated jointly.}

\item{autoburnin}{a logical flag indicating of the Gelman and Rubin's convergence diagnostic,
whether variables in x should be transformed to improve the normality of the
distribution.  If set to TRUE, a log transform or logit transform, as appropriate,
will be applied.}

\item{bt}{the number of bootstrap.}

\item{parallel}{the option of parallel; By default, parallel is FALSE, when parallel is TRUE, this test will run in parallel.}

\item{sro}{the least ratio of sample in regimes.}

\item{display}{the option of whether to print the messages of estimated results; By default, the display is TRUE.}
}
\value{
A list containing the following components:
\item{ps}{  the p-value of test}
\item{crit}{  the crit value of test}
\item{LR}{  the statistic}
\item{LRs}{  a vector of statistics in bootstrap}
}
\description{
Threshold_Test This is a test for the numer of thresholds, and it is noted that 
when under H0 the number of Thresholds is 0, this test is the so called threshold existence test.
}
\examples{
data("data", package = "DPTM")
y <- data$data_test$y
q <- data$data_test$q
x <- as.matrix(data$data_test$x)
z <- as.matrix(data$data_test$z)
tt <- data$data_test$tt
nn <- data$data_test$nn
\donttest{
### Examples elapsed time > 5s
m1 <- Threshold_Test(y=y,x=x,q=q,cvs=z,tt=tt,nn=nn,Th=0,ms = 500,burnin=500,
bt=10,parallel=FALSE)
m1$ps
}


}
\references{
Ramírez-Rondán, N. R. (2020). Maximum likelihood estimation
of dynamic panel threshold models. Econometric Reviews, 39(3), 260-276.

Hsiao, C., Pesaran, M. H., & Tahmiscioglu, A. K. (2002).
Maximum likelihood estimation of fixed effects dynamic panel data models covering short time periods. Journal of econometrics, 109(1), 107-150.
}
\author{
Hujie Bai
}
