#'@title The Dynamic panel threshold model with multiple thresholds
#'@param y the dependent variable; vector type input.
#'@param x the independent variable; matrix type input.
#'@param y1 the lag dependent variable; vector type input; By default, y1 is NULL,
#'and then y1 will be computed by y automatically.
#'@param q the threshold variable; vector type input.
#'@param cvs the set of control variables; matrix type input;By default, cvs is NULL.
#'@param time_trend the time trend; By default, it is FALSE.
#'@param time_fix_effects the time fixed effects; By default, it is FALSE.
#'@param x1 the initial values of independent variable; matrix type input.
#'By default, x1 is NULL, and thus x1 will be computed by x automatically.
#'@param tt the length of time period.
#'@param nn the number of individuals.
#'@param Th the number of thresholds.
#'@param ms the length of MCMC chains after burn-in.
#'@param burnin the length of burn-in.
#'@param types the type of MCMC used; More details see BayesianTools::runMCMC.
#'@param ADs the options for MCMC; More details see BayesianTools::runMCMC.
#'@param r0x the lower bound of thresholds; By default, r0x is NULL,
#'and thus r0x will be computed by q automatically.
#'@param r1x the upper bound of thresholds; By default, r0x is NULL,
#'and thus r1x will be computed by q automatically.
#'@param NoY the option of threshold effects on the lag dependent variable;
#'By default, NoY is False, and thus there will be threshold effects on y1.
#'@param restart the option of iterations; By default, restart is FALSE,
#'if encounters iteration failure, please set restart as TRUE.
#'@param Only_b the option of initial equation;By default, Only_b is FALSE, and if Only_b is TRUE, initial delta y will be a constant C.;
#'Please see Hsiao (2002) and Ramírez-Rondán (2020) for more details.
#'@param w the variance ratio; By default, is NULL; It must be greater than 1.
#'@param var_u the option of variance of error term; By default, is NULL; It must be
#'greater than 0; When meet relevant ERROR, please change the var_u.
#'@param delty0 the option of delta_y; By default, delty0 is NULL; Please do not change delty0.
#'@param nCR parameter determining the number of cross-over proposals of DREAM MCMC. If nCR = 1
#'all parameters are updated jointly.
#'@param autoburnin a logical flag indicating of the Gelman and Rubin's convergence diagnostic,
#'whether variables in x should be transformed to improve the normality of the
#'distribution.  If set to TRUE, a log transform or logit transform, as appropriate,
#'will be applied.
#'@param sro the least ratio of sample in regimes.
#'@param display the option of whether to print the messages of estimated results; By default, the display is TRUE.
#'@references Ramírez-Rondán, N. R. (2020). Maximum likelihood estimation
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
#'m1 <- DPTS(y=y,q=q,x=x,cvs = z,tt=tt,nn=nn,Th=1,ms = 100,burnin = 100)
#'m1$Ths
#'m1$Ths_IC
#'m1$Coefs
#'m1$MCMC_Convergence_Diagnostic
#'plot(m1$MCMC)
#'@description
#'DPTS This is a dynamic panel threshold model with fixed effects, which
#'allows multiple thresholds, time trend term or time fixed effects.
#'@returns A list containing the following components:
#'\item{ssemin}{  the negaive log-likelihood function value}
#'\item{Ths}{  a vector of multiple thresholds in order}
#'\item{Ths_IC}{  a matrix of confidence intervals of all thresholds}
#'\item{Coefs}{  parameter estimates containing t-values}
#'\item{MCMC_Convergence_Diagnostic}{  the Gelman and Rubin's convergence diagnostic results 
#'of MCMC sample}
#'\item{model}{  a list of results of DMPL}
#'\item{MCMC}{  an object of class mcmcSampler (if one chain is run) or mcmcSamplerList, 
#'more details see BayesianTools::runMCMC}
#'@export
DPTS <- function(y,y1=NULL,x=NULL,q,cvs=NULL,time_trend =FALSE,time_fix_effects=FALSE
                 ,x1=NULL,tt,nn,Th=1,ms = 1000,burnin=1000,types = "DREAMzs",
                 ADs = FALSE,r0x=NULL,r1x=NULL,NoY = FALSE,
                 restart = FALSE,Only_b = FALSE,w=NULL,var_u = NULL,delty0=NULL,
                 nCR = 3,autoburnin=TRUE,sro =0.1,display = TRUE){

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
             restart = restart,delty0=delty0)
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
                        restart = restart,delty0=delty0),silent = TRUE)

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
                restart = restart,delty0=delty0),silent = TRUE)

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
  alpha_values <- round(abs(stats::qt(c(0.05,0.025,0.005),1e7)),3)

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

  Zvalues <- round(Zvalues[1:length(xzx)],3)
  jgs <- cbind(round(c(Coefs0,Coefs1,Coefs2),3),xzx,Zvalues)
  rownames(jgs) <-  coefs_names
  colnames(jgs) <- c("Coefs","Significance","t-value")

  if(display == TRUE){
    cat("\n","This is a Dynamic panel threshold modedl with fixed effects.","!\n")
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
  }
  




  mx$Coefs <- jgs

  return(list(Ths = gamma0, Ths_IC = IC, Coefs = jgs,ssemin = ssemin,
              MCMC_Convergence_Diagnostic = MCMC_Convergence_Diagnostic,model = mx,MCMC = out))
}


