#'@title The test for the number of thresholds.
#'@param y the dependent variable; vector type input.
#'@param x the independent variable; matrix type input.
#'@param y1 the lag dependent variable; vector type input; By default y1 is NULL,
#'and then y1 will be computed by y automatically.
#'@param q the threshold variable; vector type input.
#'@param cvs the set of control variables; matrix type input;By default cvs is NULL.
#'@param time_trend the time trend; By default it is FALSE.
#'@param time_fix_effects the time fixed effects; By default it is FALSE.
#'@param x1 the initial values of independent variable; matrix type input.
#'By default x1 is NULL, and thus x1 will be computed by x automatically.
#'@param tt the length of time period.
#'@param nn the number of individuals.
#'@param Th the number of thresholds.
#'@param ms the length of MCMC chains after burn-in.
#'@param burnin the length of burn-in.
#'@param types the type of MCMC used; More details see BayesianTools::runMCMC.
#'@param ADs the options for MCMC; More details see BayesianTools::runMCMC.
#'@param r0x the lower bound of thresholds; By default r0x is NULL,
#'and thus r0x will be computed by q automatically.
#'@param r1x the upper bound of thresholds; By default r0x is NULL,
#'and thus r1x will be computed by q automatically.
#'@param NoY the option of threshold effects on the lag dependent variable;
#'By default NoY is False, and thus there will be threshold effects on y1.
#'@param assumption the option of assumption; By default assumption is 1, and it can be 2;
#'More details see Hsiao (2002).
#'@param restart the option of iterations; By default restart is FALSE,
#'if encounters iteration failure, please set restart as TRUE.
#'@param Only_b the option of initial equation;By default Only_b is FALSE, and if Only_b is TRUE, initial delta y will be a constant C.
#'More details please see Hsiao (2002) and Ramirez-Rondan (2020).
#'@param w the variance ratio; By default is NULL; It must be greater than 1, and
#'only works when assumption is 2.
#'@param var_u the option of variance of error term; By default is NULL; It must be
#'greater than 0; When meet relevant ERROR, please change the var_u.
#'@param nCR parameter determining the number of cross-over proposals of DREAM MCMC. If nCR = 1
#'all parameters are updated jointly.
#'@param autoburnin a logical flag indicating of the Gelman and Rubin's convergence diagnostic,
#'whether variables in x should be transformed to improve the normality of the
#'distribution.  If set to TRUE, a log transform or logit transform, as appropriate,
#'will be applied.
#'@param sro the least ratio of sample in regimes.
#'@param bt the number of bootstrap.
#'@param parallel the option of parallel; By default parallel is FALSE, when parallel is TRUE, this test will run in parallel.
#'@param seeds the random seed.
#'@references Ramirez-Rondan, N. R. (2020). Maximum likelihood estimation
#' of dynamic panel threshold models. Econometric Reviews, 39(3), 260-276.
#'@references Hsiao, C., Pesaran, M. H., & Tahmiscioglu, A. K. (2002).
#' Maximum likelihood estimation of fixed effects dynamic panel data models covering short time periods. Journal of econometrics, 109(1), 107-150.
#'@author Hujie Bai
#'@examples
#'data("data", package = "DPTM")
#'y <- data$data_test$y
#'q <- data$data_test$q
#'x <- as.matrix(data$data_test$x)
#'z <- as.matrix(data$data_test$z)
#'tt <- data$data_test$tt
#'nn <- data$data_test$nn
#'\dontrun{
#'### Examples elapsed time > 5s
#'m1 <- Threshold_Test(y=y,x=x,q=q,cvs=z,tt=tt,nn=nn,Th=0,ms = 500,burnin=500,
#'assumption = 1,bt=100,parallel=TRUE)
#'m1$ps
#'}
#'
#'
#'@description
#'DPTS This is a dynamic panel threshold model with fixed effects, which
#'allows multiple thresholds, time trend term or time fixed effects.
#'@export
Threshold_Test <- function(y,y1=NULL,x=NULL,q,cvs=NULL,time_trend =FALSE,time_fix_effects=FALSE
                           ,x1=NULL,tt,nn,Th=0,ms = 1000,burnin=1000,types = "DREAMzs",
                           ADs = FALSE,r0x=NULL,r1x=NULL,NoY = FALSE,assumption = 1,
                           restart = FALSE,Only_b = FALSE,w=NULL,var_u = NULL,
                           nCR = 3,autoburnin=TRUE,bt=100,parallel=TRUE,seeds = 2024,sro =0.1){
  cat("\n","Test for the number of Thresholds","\n")
  cat("\n","It is noted that when under H0 the number of Thresholds is 1, this test is the so called threshold existence test.","\n")
  cat("\n","---------------------------------------------------","\n")
  cat("\n","H0: There are ",Th," thresholds","\n")

  cat("\n","H1: There are ",Th+1," thresholds","\n")
  cat("\n","---------------------------------------------------","\n")

  if(Th == 0){
    set.seed(seeds)
    m0 <- DPML(y=y,y1=y1,x=cbind(x,cvs),w=w,var_u = var_u,tt,nn,assumption = assumption,
               time_trend = time_trend,time_fix_effects=time_fix_effects,restart = restart,
               x1=x1,Only_b = Only_b)

    du <- matrix(m0$duit,tt - 1)
    dy <- matrix(m0$dy0,tt - 1)
    dfit <- dy - du
    m0s <- as.vector(m0$ssemin)
  }else{
    set.seed(seeds)
    m0 <- DPTS(y=y,y1=y1,x=x,q=q,cvs=cvs,time_trend =time_trend,time_fix_effects=time_fix_effects
               ,x1=x1,tt=tt,nn=nn,Th=Th,ms = ms,burnin=burnin,types = types,
               ADs = ADs,r0x=r0x,r1x=r1x,NoY = NoY,assumption = assumption,
               restart = restart,Only_b = Only_b,w=w,var_u = var_u,
               nCR = nCR,autoburnin=autoburnin,sro = sro)

    du <- matrix(m0$model$duit,tt - 1)
    dy <- matrix(m0$model$dy0,tt - 1)
    dfit <- dy - du
    m0s <- as.vector(m0$ssemin)
  }

  set.seed(seeds)
  m1 <- DPTS(y=y,y1=y1,x=x,q=q,cvs=cvs,time_trend =time_trend,time_fix_effects=time_fix_effects
             ,x1=x1,tt=tt,nn=nn,Th=Th+1,ms = ms,burnin=burnin,types = types,
             ADs = ADs,r0x=r0x,r1x=r1x,NoY = NoY,assumption = assumption,
             restart = restart,Only_b = Only_b,w=w,var_u = var_u,
             nCR = nCR,autoburnin=autoburnin,sro = sro)

  m1s <- as.vector(m1$ssemin)

  cat("\n","Bootstrap:","\n")
  cat("\n","Parallel: ",parallel,"\n")

  LRs = (m0s - m1s)/(m1s/((tt-2)*nn))



  btprocedure = function(j){

    set.seed(j)
    cy = sample(1:nn,nn,replace = TRUE)
    dub = du[,cy]
    dyb = matrix(dfit + dub,ncol = 1)

    if(Th == 0){
      set.seed(j)
      ins <- utils::capture.output(m0b <- try(DPML(y=y,y1=y1,x=cbind(x,cvs),w=w,var_u = var_u,tt,nn,assumption = assumption,
                                            time_trend = time_trend,time_fix_effects=time_fix_effects,restart = restart,
                                            x1=x1,Only_b = Only_b,delty0=dyb)))


    }else{
      set.seed(j)
      ins <- utils::capture.output(m0b <- try(DPTS(y=y,y1=y1,x=x,q=q,cvs=cvs,time_trend =time_trend,time_fix_effects=time_fix_effects
                                            ,x1=x1,tt=tt,nn=nn,Th=Th,ms = ms,burnin=burnin,types = types,
                                            ADs = ADs,r0x=r0x,r1x=r1x,NoY = NoY,assumption = assumption,
                                            restart = restart,Only_b = Only_b,w=w,var_u = var_u,
                                            nCR = nCR,autoburnin=autoburnin,delty0=dyb,sro = sro)))


    }
    set.seed(j)
    ins <- utils::capture.output(m1b <- try(DPTS(y=y,y1=y1,x=x,q=q,cvs=cvs,time_trend =time_trend,time_fix_effects=time_fix_effects
                                          ,x1=x1,tt=tt,nn=nn,Th=Th+1,ms = ms,burnin=burnin,types = types,
                                          ADs = ADs,r0x=r0x,r1x=r1x,NoY = NoY,assumption = assumption,
                                          restart = restart,Only_b = Only_b,w=w,var_u = var_u,
                                          nCR = nCR,autoburnin=autoburnin,delty0=dyb,sro = sro)))
    jj = j
    while("try-error" %in% class(m0b) | "try-error" %in% class(m1b)){
      jj = jj + 1000
      set.seed(jj)
      cy = sample(1:nn,nn,replace = TRUE)
      dub = du[,cy]
      dyb = matrix(dfit + dub,ncol = 1)

      if(Th == 0){
        set.seed(jj)
        ins <- utils::capture.output(m0b <- try(DPML(y=y,y1=y1,x=cbind(x,cvs),w=w,var_u = var_u,tt,nn,assumption = assumption,
                                              time_trend = time_trend,time_fix_effects=time_fix_effects,restart = restart,
                                              x1=x1,Only_b = Only_b,delty0=dyb)))


      }else{
        set.seed(jj)
        ins <- utils::capture.output( m0b <- try(DPTS(y=y,y1=y1,x=x,q=q,cvs=cvs,time_trend =time_trend,time_fix_effects=time_fix_effects
                                               ,x1=x1,tt=tt,nn=nn,Th=Th,ms = ms,burnin=burnin,types = types,
                                               ADs = ADs,r0x=r0x,r1x=r1x,NoY = NoY,assumption = assumption,
                                               restart = restart,Only_b = Only_b,w=w,var_u = var_u,
                                               nCR = nCR,autoburnin=autoburnin,delty0=dyb,sro = sro)))


      }
      set.seed(jj)
      ins <- utils::capture.output(m1b <- try(DPTS(y=y,y1=y1,x=x,q=q,cvs=cvs,time_trend =time_trend,time_fix_effects=time_fix_effects
                                            ,x1=x1,tt=tt,nn=nn,Th=Th+1,ms = ms,burnin=burnin,types = types,
                                            ADs = ADs,r0x=r0x,r1x=r1x,NoY = NoY,assumption = assumption,
                                            restart = restart,Only_b = Only_b,w=w,var_u = var_u,
                                            nCR = nCR,autoburnin=autoburnin,delty0=dyb,sro = sro)))

    }

    m0bs <- as.vector(m0b$ssemin)

    m1bs <- as.vector(m1b$ssemin)




    LRsb = (m0bs - m1bs)/(m1bs/((tt-2)*nn))
    pbt = "No"
    if(LRsb>=LRs){
      pbt = "exceeded"
    }
    cat("\n")
    cat(j,"/",bt,pbt)
    cat("\n")
    return(LRsb)
  }

  if(!isTRUE(parallel)){
    FS = as.numeric(na.omit(purrr::map_dbl(1:bt,btprocedure)))
  }else{

    Btimes <- bt

    set_option("progress_track", TRUE)
    
    backend <- start_backend(cores = 4, cluster_type = "psock", backend_type = "async")
    
    configure_bar(type = "basic", style = 3)
    
    FS <- parabar::par_sapply(backend, 1:Btimes, function(j) {
      bt_p = btprocedure(j)
      return(bt_p)
    })
    
    stop_backend(backend)
  }



  ps = mean(FS>=LRs)
  crit = stats::quantile(FS,probs=0.95)

  js = list(ps = ps,crit = crit,LR=LRs,LRs = FS)
  cat("P-value = ",ps)
  return(js)


}


