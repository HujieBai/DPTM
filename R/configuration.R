#'@useDynLib DPTM
#'@import Rcpp
#'@import BayesianTools
#'@import stats
#'@import parabar
#'@import R6
#'@importFrom purrr map_dbl map_chr
#'@importFrom parallel detectCores
#'@importFrom MASS ginv
#'@importFrom coda gelman.diag
#'@importFrom utils capture.output


MAP2 <- function (bayesianOutput, ...)
{
  samples = BayesianTools::getSample(bayesianOutput, parametersOnly = F,
                                     ...)
  if ("mcmcSamplerList" %in% class(bayesianOutput))
    nPars <- bayesianOutput[[1]]$setup$numPars
  else nPars = bayesianOutput$setup$numPars
  best = which(samples[, nPars + 1] == max(samples[, nPars + 1]))
  cdb = length(best)
  if(cdb == 1){
    parametersMAP = samples[best, 1:nPars]
    valuesMAP = samples[best,(nPars + 1):(nPars + 3)]
  }else{
    samples2 = colMeans(samples[best,])
    parametersMAP = samples2[1:nPars]
    valuesMAP = samples2[(nPars + 1):(nPars + 3)]
  }
  
  return(list(parametersMAP = parametersMAP, valuesMAP = valuesMAP))
}


DPML0 <- function(formula = NULL,formula_cv = NULL,data,index = NULL,timeFE = FALSE,y1 = NULL,...){
  
  model_dy <- DPTM$new(data=data,index=index,...)
  
  #capture input &initx
  model_dy$capture_input(formula=formula,formula_cv = formula_cv,timeFE=timeFE,y1=y1)
  
  #MLE
  model_dy$MLE()
  
  return(model_dy)
}