
MLE <- function(y,x=NULL,delty0 =NULL,x1=NULL,cvs=NULL,ny=1,w=NULL,var_u = NULL,tt,nn,
                restart = FALSE){

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
  
  
  
  
  if(sum(x1int) == 0){
    delts <- NULL
  }else{
    delts <- qr.solve(x1int,as.vector(yint))
  }

  omega = diag(2,nrow = (tt-1),ncol = (tt-1))
  diag(omega[-1,-(tt-1)]) = -1
  diag(omega[-(tt-1),-1]) = -1


  deltxx <- cbind(deltxs,ininxc)


  dy0_iv <-  lag_transform(delty0,tt-1,nn,2,FALSE)

  dyy <- as.matrix(deltxs[,1:ny])
  ny1 <- ifelse(ncol(deltxs)==1,1,ny+1)
  dxx <- deltxs[,ny1:ncol(deltxs)]

  dx3 <-  lag_transform(as.matrix(dxx),tt-1,nn,2,FALSE)
  dy1 <- cbind( lag_transform(dyy,tt-1,nn,2,FALSE),dx3)
  dy2 <- cbind( lag_transform( lag_transform(dyy,tt-1,nn,1,FALSE),tt-2,nn,1,TRUE),dx3)

  leftiv <- t(dy2)%*%dy1

  betas_iv <- c(as.vector(MASS::ginv(leftiv)%*%(t(dy2)%*%dy0_iv)),delts)

  deltu <- as.vector(delty0 - deltxx%*%betas_iv)
  if(!is.null(var_u)){
    varu <- var_u
  }else{
    varu <-  sum(deltu^2)/(2*length(deltu))
  }


  pars <- (c(betas_iv,varu))

  varv1 <- (sum(yint^2)/nn)


  if(!is.null(w)){
    what <- w
  }else{
    what <- varv1/varu
  }
  
  if(what <= 1){
    stop("\n","The initial w must > 1, thus please input a new w.","\n")
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
  
  if(result_in$code == 1 & isTRUE(restart)){
    
    pars[1:ny] = runif(ny,-1/ny,1/ny)
    pars[cd-1] = runif(1,2,5)
    pars[cd] = pars[cd-1]/varv1
    
    result_in = suppressWarnings(try(stats::nlm(three_two,pars,delty0=delty0,evs=deltxx,
                                                omega=omega,cd=cd,tt=tt,nn=nn,iterlim = 500)
                                     ,silent = TRUE))
    
  }
  
  if("try-error" %in% class(result_in)){
    stop("\n","Encounter iteration failure!","\n",
         "Please set restart = True or check other inputs!","\n")
  }
  

  pars <- result_in$estimate
  ssemin <- result_in$minimum

  evs  <- length(betas_iv)

  bets <- pars[1:evs]

  duit <- (delty0 - deltxx%*%bets)

  cms = matrix(0,nrow = evs,ncol = evs)

  ww <- pars[cd]

  xomega = omega
  xomega[1,1] = ww
  xomega = xomega*pars[cd-1]
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
#'@param y1 the lag dependent variable; vector type input; By default, y1 is NULL,
#'and then y1 will be computed by y automatically.
#'@param time_trend the time trend; By default, it is FALSE.
#'@param time_fix_effects the time fixed effects; By default, it is FALSE.
#'@param x1 the initial values of independent variable; matrix type input.
#'By default, x1 is NULL, and thus x1 will be computed by x automatically.
#'@param tt the length of time period.
#'@param nn the number of individuals.
#'@param restart the option of iterations; By default, restart is FALSE,
#'if encounters iteration failure, please set restart as TRUE.
#'@param w the variance ratio; By default, is NULL; It must be greater than 1.
#'@param var_u the option of variance of error term; By default, is NULL; It must be
#'greater than 0; When meet relevant ERROR, please change the var_u.
#'@param delty0 the option of delta_y; By default, delty0 is NULL; Please do not change delty0.
#'@param Only_b the option of initial equation;By default, Only_b is FALSE, and if Only_b is TRUE, initial delta y will be a constant C.
#'Please see Hsiao (2002) and Ramírez-Rondán (2020) for more details.
#'@param display the option of whether to print the messages of estimated results; By default, the display is TRUE.
#'@references Ramírez-Rondán, N. R. (2020). Maximum likelihood estimation
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
#'m1 <- DPML(y=y,x=xx,tt=tt,nn=nn)
#'m1$Coefs
#'@describeIn DPML This is an dynamic panel linear model with fixed effects, which
#'allows time trend term or time fixed effects.
#'@returns A list containing the following components:
#'\item{ssemin}{  the negaive log-likelihood function value}
#'\item{Coefs}{  parameter estimates containing Z-values}
#'\item{pars}{  iterated results for all parameters}
#'\item{duit}{  the first-difference form of residuals}
#'\item{dy0}{  the first-difference form of dependent variable}
#'\item{xx}{  the independent variables and their initial values}
#'\item{covariance_matrix}{  the covariance matrix}
#'\item{Ses}{  the standard errors of coefs}
#'\item{Zvalues}{  the values of the statistic}
#'\item{ccd}{  the number of independent variables}
#'\item{coefs}{  parameter estimates containing their initial valuess}
#'@useDynLib DPTM
#'@import Rcpp
#'@import BayesianTools
#'@import stats
#'@import parabar
#'@importFrom purrr map_dbl
#'@importFrom MASS ginv
#'@importFrom coda gelman.diag
#'@importFrom utils capture.output
#'@export
DPML <- function(y,y1=NULL,x=NULL,w=NULL,var_u = NULL,tt,nn,
                 time_trend =FALSE,time_fix_effects=FALSE,restart = FALSE,
                 x1=NULL,delty0=NULL,Only_b = FALSE,display = TRUE){
  ny=1
  time_shifts <- as.matrix(rep(1:tt,nn))
  time_effects <- kronecker(rep(1,nn),diag(tt))[,-c(1,2)]
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
             restart = restart,delty0=delty0)



  coefs <- round(fit_model$coefs,3)
  Zvalues <- fit_model$Zvalues
  alpha_values <- round(abs(stats::qnorm(c(0.05,0.025,0.005))),3)

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
  
  Ses <- round(fit_model$Ses[1:length(coefs)],3)
  Zvalues <- round(Zvalues[1:length(coefs)],3)
  pvalues <- round(2*(1-pnorm(abs(Zvalues))),3)
  jgs <- cbind(coefs,Ses,Zvalues,pvalues,xzx)
  rownames(jgs) <-  coefs_names
  colnames(jgs) <- c("Estimate","Std. Error","Z-value","Pr(>|z|)","Significance")

  jgs <- jgs[1:fit_model$ccd,]
  
  if(display == TRUE){
    cat("\n","This is an dynamic panel linear model with fixed effects","\n")
    cat("\n","---------------------------------------------------","\n")
    cat("\n","Time Fixed Effects: ",time_fix_effects," !\n")
    cat("\n","---------------------------------------------------","\n")
    cat("\n","Time Shifts: ",time_trend," !\n")
    cat("\n","---------------------------------------------------","\n")
    
    
    cat("\n","The coefs are: ","\n")
    print(jgs)
    
  }




  fit_model$Coefs <- jgs

  return(fit_model)
}


