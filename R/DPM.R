
MLE <- function(y,x=NULL,delty0 =NULL,x1=NULL,cvs=NULL,ny=1,w=NULL,var_u = NULL,tt,nn,assumption = 1,
                restart = FALSE,stages = 1){

  if(is.null(delty0)){
    delty0 <- lag_transform(as.matrix(y),tt,nn,1,FALSE) -  lag_transform(as.matrix(y),tt,nn,1,TRUE)
  }


  ininxc <- matrix(0,nrow = (tt-1),ncol = nn)
  ininxc[1,] <- rep(1,nn)
  ininxc <- matrix(ininxc,ncol = 1)

  XX = x

  if(!is.null(cvs)){
    XX = cbind(XX,cvs)
  }

  if(!is.null(x1)){
    deltxi1 =  lag_transform(x1,tt,nn,1,FALSE) -  lag_transform(x1,tt,nn,1,TRUE)

    deltxi1 = matrix(deltxi1,nrow = tt-1)

    deltxi1[2:(tt-1),] = 0

    deltxi1 = matrix(deltxi1,nrow =  (tt-1)*nn)

    ininxc<- cbind(ininxc,deltxi1)
  }

  deltxs =  lag_transform(XX,tt,nn,1,FALSE) -  lag_transform(XX,tt,nn,1,TRUE)
  deltxs = matrix(deltxs,nrow = tt-1)
  deltxs[1,] = 0
  deltxs = matrix(deltxs,nrow =  (tt-1)*nn)


  x1int = matrix(matrix(ininxc,nrow = tt-1)[1,],nrow = nn)
  yint = matrix(matrix(delty0,nrow = (tt-1))[1,],ncol = 1)

  if(qr(x1int)$rank != ncol(x1int)){
    stop("\n","The inital deltaX1 Matrix is singular!","\n",
         "Please change the parameter x1!","\n")
  }
  delts <- qr.solve(x1int,as.vector(yint))

  omega = diag(2,nrow = (tt-1),ncol = (tt-1))
  diag(omega[-1,-(tt-1)]) = -1
  diag(omega[-(tt-1),-1]) = -1

  deltxx <- cbind(deltxs,ininxc)


  dy0_iv <-  lag_transform(delty0,tt-1,nn,2,FALSE)

  dyy <- as.matrix(deltxs[,1:ny])
  dxx <- deltxs[,(ny+1):ncol(deltxs)]

  dx3 <-  lag_transform(dxx,tt-1,nn,2,FALSE)
  dy1 <- cbind( lag_transform(dyy,tt-1,nn,2,FALSE),dx3)
  dy2 <- cbind( lag_transform( lag_transform(dyy,tt-1,nn,1,FALSE),tt-2,nn,1,TRUE),dx3)

  leftiv <- t(dy2)%*%dy1
  if(qr(leftiv)$rank != ncol(leftiv)){
    stop("\n","The IV Matrix (t(Y_2)Y_1) is singular!","\n",
         "Please check the input y and x!","\n")
  }

  betas_iv <- c(as.vector(inverse_cpp(leftiv)%*%(t(dy2)%*%dy0_iv)),delts)

  deltu <- as.vector(delty0 - deltxx%*%betas_iv)
  if(!is.null(var_u)){
    varu <- var_u
  }else{
    varu <-  sum(deltu^2)/(2*length(deltu))
  }


  pars <- (c(betas_iv,varu))

  varv1 <- (sum(yint^2)/nn)


  if(assumption != 1){

    if(!is.null(w)){
      what <- w
    }else{
      what <- varv1/varu
    }

    if(what <= 1){
      stop("\n","When assumption == 2, initial w must > 1, thus please input a new w.","\n")
    }

    pars <- c(pars,what)
    cd = length(pars)

    result_in <- suppressWarnings(try(stats::nlm(three_two,pars,delty0=delty0,evs=deltxx,
                                    omega=omega,cd=cd,tt=tt,nn=nn,iterlim = 500)
                         ,silent = TRUE))

    if("try-error" %in% class(result_in)){
      if(isTRUE(restart)){
        pars[1:ny] = runif(ny,-1/ny,1/ny)
        pars[cd-1] = runif(1,2,5)
        pars[cd] = pars[cd-1]/varv1

        result_in = suppressWarnings(try(stats::nlm(three_two,pars,delty0=delty0,evs=deltxx,
                                                    omega=omega,cd=cd,tt=tt,nn=nn,iterlim = 500)
                                         ,silent = TRUE))
      }else{
        stop("\n","Encounter iteration failure!","\n",
             "Please set restart = True or check other inputs!","\n")
      }
    }



  }else{

    losf <- ifelse(stages == 2,three_oneb,three_one)

    cd = length(pars)

    result_in = suppressWarnings(try(stats::nlm(losf,pars,delty0=delty0,evs=deltxx,
                                                omega=omega,cd=cd,tt=tt,nn=nn,
                                                ny=ny,varv1=varv1,iterlim = 500)
                                     ,silent = TRUE))

    if("try-error" %in% class(result_in)){
      if(isTRUE(restart)){
        pars[1:ny] = runif(ny,-1/ny,1/ny)
        pars[cd] = varv1/2
        result_in = suppressWarnings(try(stats::nlm(losf,pars,delty0=delty0,evs=deltxx,
                                                    omega=omega,cd=cd,tt=tt,nn=nn,
                                                    ny=ny,varv1=varv1,iterlim = 500)
                                         ,silent = TRUE))
      }else{
        stop("\n","Encounter iteration failure!","\n",
             "Please set restart = True or check other inputs!","\n")
      }
    }



  }

  pars <- result_in$estimate
  ssemin <- result_in$minimum

  evs  <- length(betas_iv)

  bets <- pars[1:evs]

  duit <- (delty0 - deltxx%*%bets)

  cms = matrix(0,nrow = evs,ncol = evs)

  ww <- ifelse(assumption == 1,pars[cd-1],varv1/pars[cd])

  xomega = omega
  xomega[1,1] = ww
  xomega = xomega*pars[cd]
  xomega = inverse_cpp(xomega)

  for (i in 1:nn) {
    cms = cms + t(deltxx[((tt-1)*(i-1)+1):((tt-1)*i),])%*%xomega%*%(deltxx[((tt-1)*(i-1)+1):((tt-1)*i),])
  }
  cms = try(MASS::ginv(cms),silent = TRUE)

  if("try-error" %in% class(cms)){
    ses = "singular, and please check inputs"
    zvalues = "singular, and please check inputs"
    cms = "singular, and please check inputs"
  }else{
    ses = sqrt(abs(diag(cms)))
    zvalues = bets/ses
  }

  jieguo = list(ssemin = ssemin,coefs = pars[1:evs],pars = pars,duit = duit,dy0=delty0,xx=deltxx
                ,covariance_matrix = cms,Ses= ses,Zvalues = zvalues,ccd = evs - length(delts))


  return(jieguo)

}


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

#'@title The dynamic panel linear model with fixed effects
#'@param y the dependent variable; vector type input.
#'@param x the independent variable; matrix type input.
#'@param y1 the lag dependent variable; vector type input; By default y1 is NULL,
#'and then y1 will be computed by y automatically.
#'@param time_trend the time trend; By default it is FALSE.
#'@param time_fix_effects the time fixed effects; By default it is FALSE.
#'@param x1 the initial values of independent variable; matrix type input.
#'By default x1 is NULL, and thus x1 will be computed by x automatically.
#'@param tt the length of time period.
#'@param nn the number of individuals.
#'@param assumption the option of assumption; By default assumption is 1, and it can be 2;
#'More details see Hsiao (2002).
#'@param restart the option of iterations; By default restart is FALSE,
#'if encounters iteration failure, please set restart as TRUE.
#'@param Only_b the option of initial equation;By default Only_b is FALSE, and if Only_b is TRUE, initial delta y will be a constant C.;
#'More details please see Hsiao (2002) and Ramirez-Rondan (2020).
#'@param w the variance ratio; By default is NULL; It must be greater than 1, and
#'only works when assumption is 2.
#'@param var_u the option of variance of error term; By default is NULL; It must be
#'greater than 0; When meet relevant ERROR, please change the var_u.
#'@param delty0 the option of delta_y; By default delty0 is NULL; Pleas do not change delty0.
#'@param Only_b the option of initial equation;By default Only_b is FALSE, and if Only_b is TRUE, initial delta y will be a constant C.
#'More details please see Hsiao (2002) and Ramirez-Rondan (2020).
#'@references Ramirez-Rondan, N. R. (2020). Maximum likelihood estimation
#' of dynamic panel threshold models. Econometric Reviews, 39(3), 260-276.
#'@references Hsiao, C., Pesaran, M. H., & Tahmiscioglu, A. K. (2002).
#' Maximum likelihood estimation of fixed effects dynamic panel data models covering short time periods. Journal of econometrics, 109(1), 107-150.
#'@author Hujie Bai
#'@examples
#'data("data", package = "DPTM")
#'y <- data$data_test_linear$y
#'q <- data$data_test_linear$q
#'x <- as.matrix(data$data_test_linear$x)
#'z <- as.matrix(data$data_test_linear$z)
#'tt <- data$data_test_linear$tt
#'nn <- data$data_test_linear$nn
#'xx <- cbind(x,z)
#'m1 <- DPML(y=y,x=xx,tt=tt,nn=nn,assumption = 1)
#'m1$Coefs
#'@describeIn DPML This is a dynamic panel linear model with fixed effects, which
#'allows time trend term or time fixed effects.
#'@returns A List of estimate results.
#'@useDynLib DPTM
#'@import Rcpp
#'@import RcppEigen
#'@import BayesianTools
#'@import Matrix
#'@import foreach
#'@import purrr
#'@import MASS
#'@import stats
#'@import snow
#'@import doSNOW
#'@import utils
#'@import parallel
#'@import coda
#'@export
DPML <- function(y,y1=NULL,x=NULL,w=NULL,var_u = NULL,tt,nn,assumption = 1,
                 time_trend =FALSE,time_fix_effects=FALSE,restart = FALSE,
                 x1=NULL,delty0=NULL,Only_b = FALSE){
  ny=1
  time_shifts <- as.matrix(rep(1:tt,nn))
  time_effects <- kronecker(rep(1,nn),diag(tt))[,-c(1,2,3)]
  if(all(c(time_trend,time_fix_effects))){
    stop("\n","time_fix_effects or time_shifts, that both are TRUE can not be accepted! ","\n")
  }


  if(is.null(y1)){
    y1 = as.vector(matrix(rbind(rep(0,nn),matrix(y,tt,nn)[-tt,]),ncol = 1))
  }else{
    y1 = as.vector(y1)
  }

  if(isTRUE(Only_b)){
    x1 = NULL
  }else{
    if(is.null(x1)){
      x1 = x
    }
  }

  cvs0 = NULL
  if(isTRUE(time_trend)){
    cvs0 = cbind(cvs0,time_shifts)
  }

  if(isTRUE(time_fix_effects)){
    cvs0 = cbind(cvs0,time_effects)
  }

  fit_model <- MLE(y=y,x=cbind(y1,x),x1=x1,cvs=cvs0,ny=ny,w=w,var_u = var_u,tt=tt,nn=nn,
             assumption = assumption,restart = restart,delty0=delty0)



  coefs <- round(fit_model$coefs,3)
  Zvalues <- fit_model$Zvalues
  alpha_values <- round(abs(qt(c(0.05,0.025,0.005),1e7)),3)

  coefs_names <- c()
  xzx <- c()


  for (i in 1:length(coefs)) {

    if(abs(Zvalues[i])>= alpha_values[3]){
      xzx[i] <- "***"
    }else{
      if(abs(Zvalues[i])< alpha_values[1]){
        xzx[i] <- ""
      }else{
        if(abs(Zvalues[i])<= alpha_values[2]){
          xzx[i] <- "*"
        }else{
          xzx[i] <- "**"
        }

      }

    }

    if(i == 1){
      coefs_names <- c(coefs_names,paste("y", i, sep = ""))

    }else{
      coefs_names <- c(coefs_names,paste("x_", i, sep = ""))
    }

  }

  jgs <- cbind(coefs,xzx)
  rownames(jgs) <-  coefs_names
  colnames(jgs) <- c("Coefs","Significance")

  jgs <- jgs[1:fit_model$ccd,]

  cat("\n","This is a Dynamic panel modedl with fixed effects.","\n",
          "It follows Assumption ",assumption,", and Only_b =",Only_b,"!\n")
  cat("\n","---------------------------------------------------","\n")
  cat("\n","Time Fixed Effects: ",time_fix_effects," !\n")
  cat("\n","---------------------------------------------------","\n")
  cat("\n","Time Shifts: ",time_trend," !\n")
  cat("\n","---------------------------------------------------","\n")


  cat("\n","The coefs are: ","\n")
  print(jgs)



  fit_model$Coefs <- jgs

  return(fit_model)
}


#'@title The Dynamic panel threshold model with multiple thresholds
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
#'@param Only_b the option of initial equation;By default Only_b is FALSE, and if Only_b is TRUE, initial delta y will be a constant C.;
#'More details please see Hsiao (2002) and Ramirez-Rondan (2020).
#'@param w the variance ratio; By default is NULL; It must be greater than 1, and
#'only works when assumption is 2.
#'@param var_u the option of variance of error term; By default is NULL; It must be
#'greater than 0; When meet relevant ERROR, please change the var_u.
#'@param delty0 the option of delta_y; By default delty0 is NULL; Pleas do not change delty0.
#'@param nCR parameter determining the number of cross-over proposals of DREAM MCMC. If nCR = 1
#'all parameters are updated jointly.
#'@param autoburnin a logical flag indicating of the Gelman and Rubin's convergence diagnostic,
#'whether variables in x should be transformed to improve the normality of the
#'distribution.  If set to TRUE, a log transform or logit transform, as appropriate,
#'will be applied.
#'@param sro the least ratio of sample in regimes.
#'@references Ramirez-Rondan, N. R. (2020). Maximum likelihood estimation
#' of dynamic panel threshold models. Econometric Reviews, 39(3), 260-276.
#'@references Hsiao, C., Pesaran, M. H., & Tahmiscioglu, A. K. (2002).
#' Maximum likelihood estimation of fixed effects dynamic panel data models covering short time periods. Journal of econometrics, 109(1), 107-150.
#'@author Hujie Bai
#'@examples
#'data("data", package = "DPTM")
#'y <- data$data_test$y
#'q <-data$data_test$q
#'x <- as.matrix(data$data_test$x)
#'z <- as.matrix(data$data_test$z)
#'tt <- data$data_test$tt
#'nn <- data$data_test$nn
#'m1 <- DPTS(y=y,q=q,x=x,cvs = z,tt=tt,nn=nn,Th=1,assumption = 1)
#'m1$Ths
#'m1$Ths_IC
#'m1$Coefs
#'m1$MCMC_Convergence_Diagnostic
#'plot(m1$MCMC)
#'@describeIn DPTS This is a dynamic panel threshold model with fixed effects, which
#'allows multiple thresholds, time trend term or time fixed effects.
#'@returns A List of estimate results.
#'@export
DPTS <- function(y,y1=NULL,x=NULL,q,cvs=NULL,time_trend =FALSE,time_fix_effects=FALSE
                 ,x1=NULL,tt,nn,Th=1,ms = 1000,burnin=1000,types = "DREAMzs",
                 ADs = FALSE,r0x=NULL,r1x=NULL,NoY = FALSE,assumption = 1,
                 restart = FALSE,Only_b = FALSE,w=NULL,var_u = NULL,delty0=NULL,
                 nCR = 3,autoburnin=TRUE,sro =0.1){

  if(Th < 1){
    stop("\n","Th must be greater than 0 !","\n")
  }

  ny= Th+1

  time_shifts <- as.matrix(rep(1:tt,nn))
  time_effects <- kronecker(rep(1,nn),diag(tt))[,-c(1,2,3)]
  if(all(c(time_trend,time_fix_effects))){
    stop("\n","time_fix_effects or time_shifts, that both are TRUE can not be accepted ! ","\n")
  }


  r0 <- stats::quantile(q, probs = 0.15)
  if(!is.null(r0x)){
    r0=r0x
  }
  r1 <- stats::quantile(q, probs = 0.85)
  if(!is.null(r1x)){
    r1=r1x
  }

  if(is.null(x1)){
    x1 = cbind(x,cvs)
  }

  if(is.null(y1)){
    y1 = as.vector(matrix(rbind(rep(0,nn),matrix(y,tt,nn)[-tt,]),ncol = 1))
  }else{
    y1 = as.vector(y1)
  }

  xx = NULL
  nx = 0
  if(!is.null(x)){
    nx=ncol(x)
    for (i in 1:ny) {
      xx <- rbind(xx,x)
    }
    xx = matrix(xx,nrow = nrow(x))
  }

  cvs0 = cvs
  if(isTRUE(time_trend)){
    cvs0 = cbind(cvs0,time_shifts)
  }

  if(isTRUE(time_fix_effects)){
    cvs0 = cbind(cvs0,time_effects)
  }

  nyy <- ifelse(isTRUE(NoY),1,ny)

  mm0 <- MLE(y=y,x=cbind(y1,x),x1=x1,cvs=cvs0,ny=1,w=w,var_u = var_u,tt=tt,nn=nn,
                   assumption = assumption,restart = restart,delty0=delty0)
  sse0x = mm0$ssemin

  ybl <- length(y1)

  yy_1 = matrix(y1,nrow = length(y1),ncol = Th)

  msx = function(ga){

    if(all(ga == sort(ga))){

      rts <- matrix(q,nrow = nn*tt,ncol = (Th)) - matrix(ga,nrow  = nn*tt,ncol = Th,byrow = TRUE)
      rts[which(rts>0,arr.ind = TRUE)] <- 1
      rts[which(rts<=0,arr.ind = TRUE)] <- 0
      rts0 <- cbind(rep(1,nn*tt),rts)

      if(any(colSums(rts)<sro*ybl)){
        sse0 = sse0x
      }else{


        for (i in 1:(ncol(rts0)-1)) {
          rts0[,i] <- rts0[,i] - rts0[,i+1]
        }

        xxx <- NULL
        rts0x <- NULL
        if(!is.null(x)){

          for (i in 1:nx) {
            rts0x <- cbind(rts0x,rts0)
          }
          xxx = xx*rts0x
        }

        if(isFALSE(NoY)){
          xxx = cbind(y1,yy_1*rts,xxx)
        }else{
          xxx = cbind(y1,xxx)
        }

        if(qr(xxx)$rank != ncol(xxx)){
          sse0 = sse0x
        }else{

          mx <- try(MLE(y=y,x=xxx,x1=x1,cvs=cvs0,ny=nyy,w=w,var_u = var_u,tt=tt,nn=nn,
                     assumption = assumption,restart = restart,delty0=delty0),silent = TRUE)

          if("try-error" %in% class(mx)){
            sse0 = sse0x
          }else{
            sse0 = mx$ssemin
          }

        }
      }
    }else{
      sse0 = sse0x
    }

    return(sse0)
  }

  pf1 = function(ga){
    s1 = msx(ga)
    return(-s1)
  }

  priors <- BayesianTools::createUniformPrior(lower = c(rep(r0,Th)), upper = c(rep(r1,Th)))
  bayesianSetup <- BayesianTools::createBayesianSetup(pf1, prior = priors)
  if(types == "Metropolis"){
    settings = list(iterations = (ms+burnin),burnin=burnin,adapt = ADs,adaptationNotBefore = (ms+burnin), nCR = nCR)
  }else{
    settings = list(iterations = (ms+burnin),burnin=burnin, nCR = nCR)
  }


  out <- BayesianTools::runMCMC(bayesianSetup = bayesianSetup,
                                sampler = types,
                                settings = settings)

  samples = BayesianTools::getSample(out,parametersOnly = F)
  mres = MAP2(out)
  gamma0 = mres[[1]]

  ssemin = -mres[[2]][2]

  rts <- matrix(q,nrow = nn*tt,ncol = Th) - matrix(gamma0,nrow  = nn*tt,ncol = Th,byrow = TRUE)
  rts[which(rts>0,arr.ind = TRUE)] <- 1
  rts[which(rts<=0,arr.ind = TRUE)] <- 0
  rts0 <- cbind(rep(1,nn*tt),rts)

  for (i in 1:(ncol(rts0)-1)) {
    rts0[,i] <- rts0[,i] - rts0[,i+1]
  }

  xxx <- NULL
  rts0x <- NULL
  if(!is.null(x)){

    for (i in 1:nx) {
      rts0x <- cbind(rts0x,rts0)
    }
    xxx = xx*rts0x
  }

  if(isFALSE(NoY)){
    xxx = cbind(cbind(y1,yy_1)*rts0,xxx)
  }else{
    xxx = cbind(y1,xxx)
  }

  if(qr(xxx)$rank != ncol(xxx)){
    stop("\n","The independent variable matrix with threshold effects is singular!","\n",
         "Please check x and other inputs!")
  }

  mx <- try(MLE(y=y,x=xxx,x1=x1,cvs=cvs0,ny=nyy,w=w,var_u = var_u,tt=tt,nn=nn,
                assumption = assumption,restart = restart,delty0=delty0,stages = 2),silent = TRUE)

  if("try-error" %in% class(mx)){
    stop("There is an Error after given thresholds, please check any inputs!")
  }

  IC = matrix(0,3,Th*2)

  colnames(IC) <- rep(c("Lower","Upper"),Th)
  rownames(IC) <- c("90%","95%","99%")

  for (i in 1:Th) {
    IC[1,(2*i-1):(2*i)] = c(stats::quantile(samples[,i],0.05),stats::quantile(samples[,i],0.95))
    IC[2,(2*i-1):(2*i)] = c(stats::quantile(samples[,i],0.025),stats::quantile(samples[,i],0.975))
    IC[3,(2*i-1):(2*i)] = c(stats::quantile(samples[,i],0.005),stats::quantile(samples[,i],0.995))
  }

  chains <- out$chain

  for (i in 1:length(chains)) {
    chains[[i]] <- chains[[i]][,1:Th]
  }
  MCMC_Convergence_Diagnostic <- try(coda::gelman.diag(chains,autoburnin=autoburnin))
  if("try-error" %in% class(MCMC_Convergence_Diagnostic)){
    MCMC_Convergence_Diagnostic <- coda::gelman.diag(chains, autoburnin= FALSE)
  }

  betas <- mx$coefs[1:mx$ccd]


  nxx <- ifelse(is.null(x),0,ncol(x)*ncol(rts0))
  ncc <- ifelse(is.null(cvs),0,ncol(cvs))

  Coefs0 <-  betas[1:nyy]
  Coefs1 <- as.vector(t( matrix(betas[(nyy+1):(nyy+nxx)],ncol(x),byrow = TRUE)))
  Coefs2 <-  NULL
  if(!is.null(cvs)){
    Coefs2 <-  betas[(nyy+nxx+1):(nyy+nxx+ncc)]
  }




  Zvalues <- mx$Zvalues[1:(nyy+nxx+ncc)]
  alpha_values <- round(abs(qt(c(0.05,0.025,0.005),1e7)),3)

  coefs_names <- c()
  xzx <- c()


  for (i in 1:length(Zvalues)) {

    if(abs(Zvalues[i])>= alpha_values[3]){
      xzx[i] <- "***"
    }else{
      if(abs(Zvalues[i])< alpha_values[1]){
        xzx[i] <- ""
      }else{
        if(abs(Zvalues[i])<= alpha_values[2]){
          xzx[i] <- "*"
        }else{
          xzx[i] <- "**"
        }

      }

    }


  }

  if(isTRUE(NoY)){
    coefs_names <- c(coefs_names,paste("y", 1, sep = ""))
  }else{

    for (i in 1:nyy) {
      coefs_names <- c(coefs_names,paste("y_regime_", i, sep = ""))
    }

  }

  if(!is.null(x)){
    for (i in 1:ncol(x)) {

      for (j in 1:ny) {
        coefs_names <- c(coefs_names,paste("x",i,"_regime_", j, sep = ""))
      }

    }
  }

 if(!is.null(cvs)){
   for (i in 1:ncol(cvs)) {
     coefs_names <- c(coefs_names,paste("Control_", i, sep = ""))
   }
 }


  jgs <- cbind(round(c(Coefs0,Coefs1,Coefs2),3),xzx)
  rownames(jgs) <-  coefs_names
  colnames(jgs) <- c("Coefs","Significance")

  cat("\n","This is a Dynamic panel threshold modedl with fixed effects.","\n",
          "It follows Assumption ",assumption,", and Only_b =",Only_b,"!\n")
  cat("\n","---------------------------------------------------","\n")
  cat("\n","Time Fixed Effects: ",time_fix_effects," !\n")
  cat("\n","---------------------------------------------------","\n")
  cat("\n","Time Shifts: ",time_trend," !\n")
  cat("\n","---------------------------------------------------","\n")
  cat("\n","The number of threshold is ",Th," !\n")
  cat("\n","---------------------------------------------------","\n")
  cat("\n","The estiamtes of thresholds: ", "\n")
  print(round(gamma0,3))
  cat("\n","Their confidence intervals are : "," !\n")
  print(round(IC,3))
  cat("\n","---------------------------------------------------","\n")
  cat("\n","The coefs are: ","\n")
  print(jgs)
  cat("\n","---------------------------------------------------","\n")
  cat("\n"," The Gelman and Rubin Convergence Diagnostic is ","\n")
  print(MCMC_Convergence_Diagnostic)
  cat("\n","If the  Upper C.I. are not close to 1, please set a longer burn-in or ms !","\n")
  cat("\n","If there any Inf or NaN, please set a longer burn-in or ms, or set nCR as 1 !","\n")
  cat("\n","---------------------------------------------------","\n")




  mx$Coefs <- jgs

  return(list(Ths = gamma0, Ths_IC = IC, Coefs = jgs,ssemin = ssemin,
              MCMC_Convergence_Diagnostic = MCMC_Convergence_Diagnostic,model = mx,MCMC = out))
}


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
#'@param delty0 the option of delta_y; By default delty0 is NULL; Pleas do not change delty0.
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
#'m1 <- Threshold_Test(y=y,x=x,q=q,cvs=z,tt=tt,nn=nn,Th=0,ms = 500,burnin=500,
#'assumption = 1,bt=10,parallel=FALSE)
#'m1$ps
#'@describeIn DPTS This is a dynamic panel threshold model with fixed effects, which
#'allows multiple thresholds, time trend term or time fixed effects.
#'@returns A List of estimate results.
#'@export
Threshold_Test <- function(y,y1=NULL,x=NULL,q,cvs=NULL,time_trend =FALSE,time_fix_effects=FALSE
                           ,x1=NULL,tt,nn,Th=0,ms = 1000,burnin=1000,types = "DREAMzs",
                           ADs = FALSE,r0x=NULL,r1x=NULL,NoY = FALSE,assumption = 1,
                           restart = FALSE,Only_b = FALSE,w=NULL,var_u = NULL,
                           nCR = 3,autoburnin=TRUE,bt=100,parallel=FALSE,seeds = 2024,sro =0.1){
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
      ins <- capture.output(m0b <- try(DPML(y=y,y1=y1,x=cbind(x,cvs),w=w,var_u = var_u,tt,nn,assumption = assumption,
                 time_trend = time_trend,time_fix_effects=time_fix_effects,restart = restart,
                 x1=x1,Only_b = Only_b,delty0=dyb)))


    }else{
      set.seed(j)
      ins <- capture.output(m0b <- try(DPTS(y=y,y1=y1,x=x,q=q,cvs=cvs,time_trend =time_trend,time_fix_effects=time_fix_effects
                 ,x1=x1,tt=tt,nn=nn,Th=Th,ms = ms,burnin=burnin,types = types,
                 ADs = ADs,r0x=r0x,r1x=r1x,NoY = NoY,assumption = assumption,
                 restart = restart,Only_b = Only_b,w=w,var_u = var_u,
                 nCR = nCR,autoburnin=autoburnin,delty0=dyb,sro = sro)))


    }
    set.seed(j)
    ins <- capture.output(m1b <- try(DPTS(y=y,y1=y1,x=x,q=q,cvs=cvs,time_trend =time_trend,time_fix_effects=time_fix_effects
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
        ins <- capture.output(m0b <- try(DPML(y=y,y1=y1,x=cbind(x,cvs),w=w,var_u = var_u,tt,nn,assumption = assumption,
                        time_trend = time_trend,time_fix_effects=time_fix_effects,restart = restart,
                        x1=x1,Only_b = Only_b,delty0=dyb)))


      }else{
        set.seed(jj)
        ins <- capture.output( m0b <- try(DPTS(y=y,y1=y1,x=x,q=q,cvs=cvs,time_trend =time_trend,time_fix_effects=time_fix_effects
                        ,x1=x1,tt=tt,nn=nn,Th=Th,ms = ms,burnin=burnin,types = types,
                        ADs = ADs,r0x=r0x,r1x=r1x,NoY = NoY,assumption = assumption,
                        restart = restart,Only_b = Only_b,w=w,var_u = var_u,
                        nCR = nCR,autoburnin=autoburnin,delty0=dyb,sro = sro)))


      }
      set.seed(jj)
      ins <- capture.output(m1b <- try(DPTS(y=y,y1=y1,x=x,q=q,cvs=cvs,time_trend =time_trend,time_fix_effects=time_fix_effects
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
    xc= parallel::detectCores()
    xc = xc -1
    Btimes <- bt
    pb <- utils::txtProgressBar(min=1, max=Btimes, style=3)
    progress <- function(n) utils::setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    cl <- snow::makeSOCKcluster(xc)
    doSNOW::registerDoSNOW(cl)

    FS <- foreach::foreach(i=1:Btimes, .packages = c("purrr","SparseM","Rcpp","BayesianTools"), .options.snow=opts,
                           .combine=c,.errorhandling = "remove") %dopar%{
                             bt_p = btprocedure(i)
                             return(bt_p)
                           }

    close(pb)
    snow::stopCluster(cl)
  }



  ps = mean(FS>=LRs)
  crit = stats::quantile(FS,probs=0.95)

  js = list(ps = ps,crit = crit,LR=LRs,LRs = FS)
  cat("P-value = ",ps)
  return(js)


}


