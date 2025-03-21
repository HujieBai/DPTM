#'@title Dynamic panel model with fixed effects.
#'
#'@description
#'Use a MLE procedure to estimate the dynamic panel model with fixed effects.
#'
#'@param formula formula of the covariates with threshold effects.
#'@param data data frame of the observed data.
#'@param index variable names of individuals and period; If a setting is not provided, defaults (the first variables in data will be as "id", while the second will be "year") will be used. Defaults to `NULL`.
#'@param timeFE logicals. If TRUE the time fixed effects will be allowed. Defaults to `FALSE`.
#'@param y1 lags of dependent variables; If a setting is not provided, defaults (the first-order lag) will be used. Defaults to `NULL`.
#'@param ... additional arguments, see\code{stats::nlm}.
#'
#'@usage DPML(formula, data, index=NULL, timeFE = FALSE, y1 = NULL,...)
#'
#'## S6 method for class 'DPTM' 
#'#print(...)
#'
#'@details
#'\code{DPML} can fit the dynamic panel model with fixed effects proposed by Hsiao et al. (2002), which is based on the first difference and the maximum likelihood (MLE) method.
#'
#'For a classical dynamic panel model with fixed effects having following form:
#'\deqn{y_{it}=\rho y_{it-1}+\beta_1x_{1,it}+\beta_2x_{2,it}+\alpha_i+u_{it}},
#'can use \code{y~x1+x2}.
#'
#'For a special dynamic panel model with fixed effects having the following form:
#'\deqn{\Delta y_{it}=\rho y_{it-1}+\beta_1x_{1,it}+\beta_2x_{2,it}+\alpha_i+u_{it}},
#'can use \code{ dy~x1+x2} with \code{y1}\eqn{= y_{it-1}}.
#'
#'We assume the exogenous regressor \eqn{x} is weakly exogenous, and 
#'thus the first period after difference is given by
#'\deqn{\Delta y_{i1}=\delta_0 + {\boldsymbol\delta}'_1 \Delta {\bf x}_{i1}+ v_{i1},}
#'where \eqn{E(v_{i1}| \Delta {\bf x}_{i1} )=0}. \eqn{E(v_{i1}^2)=\sigma^2_v}, 
#'\eqn{E(v_{i1}\Delta u_{i2})=-\sigma^2_u} and \eqn{E(v_{i1} \Delta u_{it})=0}
#'for \eqn{t=3,4,...,T} and \eqn{i=1,...,N}.
#'For more details, see Hsiao et al. (2002).
#'
#'
#'In addition, we solve the log-likelihood function by \code{stats::nlm} who uses \code{iterlim}
#'to set the maximum number of iterations, and thus \code{iterlim} is allowed by \code{...} in \code{DPML}.
#'
#'
#'
#'@references Hsiao, C., Pesaran, M. H., & Tahmiscioglu, A. K. (2002).
#' Maximum likelihood estimation of fixed effects dynamic panel data models covering short time periods. Journal of econometrics, 109(1), 107-150.
#'@author Hujie Bai
#'
#'@return DPML returns an object of class "DPTM".
#'The function \code{print} are used to obtain and print a print of the results.   
#'An object of class "DPTM" is a list containing at least the following components:
#'\item{coefficients}{a named vector of coefficients}
#'\item{NNLL}{the negative log-likelihood function value}
#'\item{Zvalues}{a vector of t statistics}
#'\item{Ses}{a vector of standard errors}
#'\item{covariance_matrix}{a covariance matrix}
#'\item{Th}{the number of thresholds}
#'\item{thresholds}{a named vector of thresholds}
#'
#'@examples
#'data(d1)
#'
#'# No time fixed effects
#'model1 <- DPML(y~x+z, data = d1)
#'print(model1)
#'
#'# With time fixed effects
#'model2 <- DPML(y~x+z, data = d1, timeFE = TRUE)
#'print(model2)
#'
#'@export
DPML <- function(formula,data,index = NULL,timeFE = FALSE,y1 = NULL,...){
  
  model_dy <- DPTM$new(data=data,index=index,...)
  
  #capture input &initx
  model_dy$capture_input(formula=formula,timeFE=timeFE,y1=y1)
  
  #MLE
  model_dy$MLE()
  
  return(model_dy)
}

