#'@title Tests for multiple thresholds.
#'
#'@description
#'Tests for models with different thresholds, using bootstrap method.
#'
#'@param formula formula of the covariates with threshold effects; If a setting is not provided, defaults (no covariates with threshold effects) will be used. Defaults to `NULL`.
#'@param formula_cv formula of the covariates without threshold effects; If a setting is not provided, defaults (no covariates without threshold effects) will be used. Defaults to `NULL`.
#'@param data data frame of the observed data.
#'@param index variable names of individuals and period; If a setting is not provided, defaults (the first variables in data will be as "id", while the second will be "year") will be used.Defaults to `NULL`.
#'@param Th number of thresholds; Defaults to `1`.
#'@param q threshold variable.
#'@param timeFE logicals. If TRUE the time fixed effects will be allowed. Defaults to `FALSE`.
#'@param bt the number of bootstrap; If a setting is not provided, defaults (bt = 100) will be used. Defaults to `100`.
#'@param NoY logicals. If TRUE the lags of dependent variables will be without threshold effects. Defaults to `FALSE`.
#'@param y1 lags of dependent variables; If a setting is not provided, defaults (the first-order lag) will be used. Defaults to `NULL`.
#'@param iterations MCMC iterations (50\% used for burnining). Defaults to `2000`.
#'@param sro regime (subsample) proportion; If a setting is not provided, defaults (10\%) will be used. Defaults to `0.1`.
#'@param r0x lower bound of threshold parameter space; If a setting is not provided, defaults (15\% quantile of threshold variable) will be used.
#'@param r1x upper bound of threshold parameter space; If a setting is not provided, defaults (85\% quantile of threshold variable) will be used.
#'@param parallel logicals. If TRUE test will run in parallel for saving time. Defaults to `TRUE`.
#'@param seed set seeds to guarantee the replication the test (see set.seed);
#'@param ... additional arguments to be passed to the settings of MCMC (see BayesianTools::applySettingsDefault)
#'
#'@usage Threshold_Test(formula, formula_cv, data, index=NULL, Th = 0, q, 
#'timeFE = FALSE, bt = 100,NoY = FALSE, y1 = NULL, iterations = 2000, sro = 0.1,
#'r0x = NULL, r1x = NULL, parallel=TRUE, seed = NULL,...)
#'
#'@details
#'\code{Threshold_Test} can run the Test for multiple thresholds (\code{Th} is H1).
#'The statistic is
#'\deqn{F_s=\frac{S\left(\hat{\gamma}_{s-1}\right)-S\left(\hat{\gamma}_s\right)}{S\left(\hat{\gamma}_s\right) / N(T-1)},}
#'where \eqn{s} is the number of thresholds in H1, \eqn{S\left(\hat{\gamma}_{s-1}\right)=-\ln L\left(\hat{\gamma}_{s-1}\right)} and \eqn{S\left(\hat{\gamma}_s\right)=-\ln L\left(\left(\hat{\gamma}_{s-1}^{\prime}, \hat{\gamma}_s\right)^{\prime}\right)}.
#'And the p-value is computed by bootstrap method (see Ramírez-Rondán, 2020).
#'
#'Take the two threshold model as example.
#'User must set \code{Th} = 1 firstly to reject the null hypothesis of no threshold effects;
#'Then he should set \code{Th} = 2 to reject the null hypothesis of only one threshold;
#'Lastly, set \code{Th} = 3 to accept the null hypothesis of two thresholds.
#'In other words, p-values of the first test (\code{Th} = 1) and the second test (\code{Th} = 1) should be less than significant level while the third test (\code{Th} = 3) is not.
#'
#'\code{Threshold_Test} contains all augments in \code{DPTS}, but with three new augments: \code{bt}, \code{parallel} and \code{seed}.
#'\code{bt} is the number of bootstrap (by default is 100);
#'\code{parallel} can allow user to run test in parallel to save time;
#'\code{seed} is used to guarantee the replication of tests.
#'
#'It is worthy noting that the test shrinks to the so-called threshold existence test when \code{Th} = 1.
#'
#'@references Ramírez-Rondán, N. R. (2020). Maximum likelihood estimation of dynamic panel threshold models. Econometric Reviews, 39(3), 260-276.
#'@author Hujie Bai
#'
#'
#'@returns A list with class "htest" containing the following components:
#'\item{statistic}{the value of the F-statistic.}
#'\item{parameter}{the degrees of freedom for the F-statistic.}
#'\item{p.value}{the p-value for the test.}
#'\item{null.value}{the specified hypothesized value of the null hypothesis.}
#'\item{alternative}{a character string describing the alternative hypothesis.}
#'\item{method}{a character string indicating what type of test was performed.}
#'\item{data.name}{a character string giving the name(s) of the data.}
#'\item{estimate}{the critical value of the statistic (5\% significance level).}
#'\item{LRs}{a vector of statistics from bootstrap.}
#'
#'@examples
#'\donttest{
#'### Examples elapsed time > 15s
#'
#'#data(d1)
#'
#'# H0: no threshold effects (no threshold)
#'#test0 <- Threshold_Test(y~x,y~z,data = d1, index = c('id','year'), q = d1$q, Th = 1,
#'#bt = 50, iterations = 500)
#'#test0
#'
#'# H0: one threshold 
#'#test1 <- Threshold_Test(y~x,y~z,data = d1, index = c('id','year'), q = d1$q, Th = 2,
#'#bt = 50, iterations = 500)
#'#test1
#'
#'# H0: two threshold 
#'#test2 <- Threshold_Test(y~x,y~z,data = d1, index = c('id','year'), q = d1$q, Th = 3,
#'#bt = 50, iterations = 500)
#'#test2
#'
#'}
#'
#'@noRd
#'@export
Threshold_Test <- function(formula = NULL, formula_cv = NULL, data, index = NULL,Th = 1
                           ,q,timeFE = FALSE, bt = 100,NoY = FALSE,y1 = NULL,
                           iterations = 2000,sro = 0.1,r0x = NULL,r1x = NULL,
                           parallel = TRUE, seed = NULL,...){
  
  cat("\n","Test for the number of Thresholds","\n")
  cat(" It is noted that when under H0 the number of Thresholds is 0, this test is the so called threshold existence test.","\n")
  cat(" H0: There are ",Th - 1," thresholds","\n")
  
  cat(" H1: There are ",Th," thresholds","\n")
  cat("---------------------------------------------------","\n")
  
  seed <- ifelse(is.null(seed),2024,seed)
  set.seed(seed)
  bt.seed <- sample(1:bt^2,bt,replace = TRUE)
  
  if(Th < 1){
    stop("Th does not be less than 1!")
  }
    
  if(Th == 1){
    
    set.seed(seed)
    m0 <- DPML0(formula = formula,formula_cv = formula_cv,data = data,
                index = index,timeFE = timeFE,y1 = y1,...)
    
  }else{
    
    set.seed(seed)
    m0 <- DPTS(formula = formula,formula_cv = formula_cv,data = data,
               index = index,Th = Th - 1,q = q,timeFE = timeFE,
               NoY = NoY,y1 = y1,iterations = iterations,sro = sro,
               r0x = r0x,r1x = r1x,...)
    
  }
  
  tt <- m0$.__enclos_env__$private$tt
  nn <- m0$.__enclos_env__$private$nn
  
  du <- matrix(m0$duit,tt - 1)
  dy <- matrix(m0$dy0,tt - 1)
  dfit <- dy - du
  m0s <- as.vector(m0$NNLL)
  
  set.seed(seed)
  m1 <- DPTS(formula = formula,formula_cv = formula_cv,data = data,
             index = index,Th = Th,q = q,timeFE = timeFE,
             NoY = NoY,y1 = y1,iterations = iterations,sro = sro,
             r0x = r0x,r1x = r1x,...)
  
  m1s <- as.vector(m1$NNLL)
  
  LRs = (m0s - m1s)/(m1s/((tt-2)*nn))
  
  btprocedure = function(j){
    
    
    set.seed(bt.seed[j])
    cy = sample(1:nn,nn,replace = TRUE)
    dub = du[,cy]
    dyb = matrix(dfit + dub,ncol = 1)
    
    if(Th == 1){
      
      set.seed(bt.seed[j])
      m0b <- suppressMessages(try(DPML0(formula = formula,formula_cv = formula_cv,data = data,
                                        index = index,timeFE = timeFE,y1 = y1,delty0=dyb,...),silent = TRUE))
      
      
    }else{
      
      set.seed(bt.seed[j])
      m0b <- suppressMessages(try(DPTS(formula = formula,formula_cv = formula_cv,data = data,
                                       index = index,Th = Th - 1,q = q,timeFE = timeFE,
                                       NoY = NoY,y1 = y1,iterations = iterations,sro = sro,
                                       r0x = r0x,r1x = r1x,delty0=dyb,...),silent = TRUE))
      
      
    }
    
    set.seed(bt.seed[j])
    m1b <- suppressMessages(try(DPTS(formula = formula,formula_cv = formula_cv,data = data,
                                     index = index,Th = Th,q = q,timeFE = timeFE,
                                     NoY = NoY,y1 = y1,iterations = iterations,sro = sro,
                                     r0x = r0x,r1x = r1x,delty0=dyb,...),silent = TRUE))
    
    while("try-error" %in% class(m0b) | "try-error" %in% class(m1b)){
      
      set.seed(bt.seed[j]+1)
      cy = sample(1:nn,nn,replace = TRUE)
      dub = du[,cy]
      dyb = matrix(dfit + dub,ncol = 1)
      
      if(Th == 1){
        
        set.seed(bt.seed[j])
        m0b <- suppressMessages(try(DPML0(formula = formula,formula_cv = formula_cv,data = data,
                                          index = index,timeFE = timeFE,y1 = y1,delty0=dyb,...),silent = TRUE))
        
        
      }else{
        
        set.seed(bt.seed[j])
        m0b <- suppressMessages(try(DPTS(formula = formula,formula_cv = formula_cv,data = data,
                                         index = index,Th = Th - 1,q = q,timeFE = timeFE,
                                         NoY = NoY,y1 = y1,iterations = iterations,sro = sro,
                                         r0x = r0x,r1x = r1x,delty0=dyb,...),silent = TRUE))
        
        
      }
      
      set.seed(bt.seed[j])
      m1b <- suppressMessages(try(DPTS(formula = formula,formula_cv = formula_cv,data = data,
                                       index = index,Th = Th,q = q,timeFE = timeFE,
                                       NoY = NoY,y1 = y1,iterations = iterations,sro = sro,
                                       r0x = r0x,r1x = r1x,delty0=dyb,...),silent = TRUE))
      
    }
    
    m0bs <- as.vector(m0b$NNLL)
    
    m1bs <- as.vector(m1b$NNLL)
    
    
    
    
    LRsb = (m0bs - m1bs)/(m1bs/((tt-2)*nn))
    pbt = "No"
    if(LRsb>=LRs){
      pbt = "exceeded"
    }
    
    cat("\n",j,"/",bt,pbt,"\n")
    return(LRsb)
  }
  
  if(!isTRUE(parallel)){
    FS = as.numeric(na.omit(purrr::map_dbl(1:bt,btprocedure)))
  }else{
    
    Btimes <- bt
    
    parabar::set_option("progress_track", TRUE)
    cores_number <- parallel::detectCores()
    backend <- parabar::start_backend(cores = cores_number, cluster_type = "psock", backend_type = "async")
    
    parabar::configure_bar(type = "basic", style = 3)
    
    FS <- parabar::par_sapply(backend, 1:Btimes, function(j) {
      bt_p = btprocedure(j)
      return(bt_p)
    })
    
    stop_backend(backend)
  }
  
  
  
  ps = mean(FS>=LRs)
  crit = stats::quantile(FS,probs=0.95)

  my_test <- list(
    statistic = c(LR = LRs),
    parameter = c(df = (tt-2)*nn),
    p.value = ps,
    null.value = c("number of thresholds" = (Th - 1)),
    alternative = "greater",
    method = "Multiple Threshold Test",
    data.name = paste("Bootstrap",bt),
    estimate = c('Critical Value' = crit),
    LRs = FS
  )
  
  class(my_test) <- "htest"

  
  return(my_test)
}

