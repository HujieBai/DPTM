#'@title Dynamic Panel Multiple Threshold Model with Fixed Effects (DPTM)
#'
#'@keywords internal
#'
#'@description
#'Use a MCMC-MLE based on two-step procedure to estimate the dynamic panel multiple threshold model with fixed effects.
#'
#'
#'@format [R6::R6Class] object.
DPTM <- R6::R6Class(
  classname = "DPTM",
  public = list(
    #'@field coefficients a named vector of coefficients
    coefficients = NULL, 
    #'@field NNLL the negative log-likelihood function value
    NNLL = NULL, 
    #'@field Zvalues a vector of t statistics
    Zvalues = NULL, 
    #'@field Ses a vector of standard errors
    Ses = NULL, 
    #'@field covariance_matrix a covariance matrix
    covariance_matrix = NULL, 
    #'@field duit a vector of residuals after difference
    duit = NULL, 
    #'@field dy0 a vector of dependent variable after difference
    dy0 = NULL,

    #'@field Th the number of thresholds
    Th = NULL,
    #'@field thresholds a named vector of thresholds
    thresholds = NULL,
    
    #'@description
    #'initialize Initializing method
    #'@param data data.frame used
    #'@param index variable names of individuals and period; If a setting is not provided, defaults (the first variables in data will be as "id", while the second will be "year") will be used
    #'@param Th number of thresholds; If a setting is not provided, defaults (Th = 0) will be used
    #'@param iterations MCMC iterations (50\% used for burnining)
    #'@param sro regime (subsample) proportion; If a setting is not provided, defaults (10\%) will be used
    #'@param w variances ratio initial value; If a setting is not provided, defaults (automatic calculation) will be used
    #'@param var_u variances (T>=2) initial value; If a setting is not provided, defaults (automatic calculation) will be used
    #'@param iterlim the maximum number of iterations; If a setting is not provided, defaults (iterlim = 500) will be used
    #'@param restart logicals. If MLE fails, set it as TRUE
    #'@param delty0 a vector of dependent variable after difference
    initialize = function(data,index = NULL,Th = NULL, iterations = NULL,
                          sro = NULL,w = NULL,var_u = NULL,iterlim = NULL,
                          restart = FALSE,delty0 = NULL) {
      private$xtset_b(data=data,index=index)
      
      if(!is.null(Th)){
        if(is.character(Th) | Th <= 0){
          stop("Check the input Th.")
        }
        self$Th <- as.integer(Th)
      }
      
      if(!is.null(iterations)){
        private$iterations <- iterations
      }
      
      if(!is.null(sro)){
        private$sro <- sro
      }
      
      if(!is.null(w)){
        private$w <- w
      }
      
      if(!is.null(var_u)){
        private$var_u <- var_u
      }
      
      if(!is.null(delty0)){
        private$delty0 <- delty0
      }
      
      if(!is.null(iterlim)){
        private$iterlim <- iterlim
      }
      
      if(!isTRUE(restart)){
        private$restart <- restart
      }
      
      
    },
    #'@description
    #'Identify and capturing inputs
    #'@param formula formula of the covariates with threshold effects;If a setting is not provided, defaults (no covariates with threshold effects) will be used
    #'@param formula_cv formula of the covariates without threshold effects;If a setting is not provided, defaults (no covariates without threshold effects) will be used
    #'@param timeFE logicals. If TRUE the time fixed effects will be allowed
    #'@param y1 lags of dependent variables; If a setting is not provided, defaults (the first-order lag) will be used
    #'@param q threshold variable
    #'@param r0x lower bound of threshold parameter space; If a setting is not provided, defaults (15\% quantile of threshold variable) will be used
    #'@param r1x upper bound of threshold parameter space; If a setting is not provided, defaults (85\% quantile of threshold variable) will be used
    #'@param NoY logicals. If TRUE the lags of dependent variables will be without threshold effects
    capture_input = function(formula = NULL,formula_cv = NULL,timeFE,y1 = NULL,
                             q = NULL,r0x = NULL,r1x = NULL,NoY = FALSE){
      
      if(is.null(self$Th)){
        
        if(is.null(formula)&is.null(formula_cv)){
          stop("The formula is empty!")
        }
        
        X1X <- NULL
        if(!is.null(formula)){
          mf <- model.frame(formula, private$data_used)
          X1X <- as.matrix(model.matrix(formula, mf)[,-1])
          
          if(qr(X1X)$rank != ncol(X1X)){
            stop('Regressors are singular!')
          }
        }
        
        X2X <- NULL
        if(!is.null(formula_cv)){
          mf <- model.frame(formula_cv, private$data_used)
          X2X <- as.matrix(model.matrix(formula_cv, mf)[,-1])
          
          if(qr(X2X)$rank != ncol(X2X)){
            stop('Regressors are singular!')
          }
        }
        
        X <- cbind(X1X,X2X)
        
        private$x1 <- X
        
        X_names <- colnames(X)
        if(is.null(X_names)){
          X_names <- rep("X",ncol(X))
        }
        
        if(isTRUE(timeFE)){
          timeFE <- kronecker(rep(1,private$nn),diag(private$tt))[,-c(1,2)]
          X <- cbind(X,timeFE)
          private$TFE <- TRUE
        }
        
        y <- model.response(mf)
        
        if(is.null(y1)){
          y1 <- as.vector(matrix(rbind(rep(0,private$nn),matrix(y,private$tt,private$nn)[-private$tt,]),ncol = 1))
        }else{
          y1 <- as.vector(y1)
        }
        
        private$coefs_names <- c("L1.Y",X_names)
        
        private$x <- cbind(y1,X)
        
      }else{
        private$r0 <- stats::quantile(q, probs = 0.15)
        private$r1 <- stats::quantile(q, probs = 0.85)
        if(!is.null(r1x)){
          private$r1 = r1x
        }
        if(!is.null(r0x)){
          private$r0 = r0x
        }
        
        private$q = q
        
        private$ny = self$Th+1
        
        if(is.null(formula)&is.null(formula_cv)){
          stop("Needs regressors!")
        }
        
        x1 <- NULL
        x_names <- NULL
        if(!is.null(formula)){
          mf <- model.frame(formula, private$data_used)
          X <- as.matrix(model.matrix(formula, mf)[,-1])
          if(qr(X)$rank != ncol(X)){
            stop('Regressors of formula are singular!')
          }
          
          if(is.null(colnames(X))){
            for (i in 1:private$ny) {
              x_names <- c(x_names,paste(1:ncol(X),'regime',i))
            }
          }else{
            for (i in 1:private$ny) {
              x_names <- c(x_names,paste(colnames(X),'regime',i))
            }
          }
          x_names <- as.vector(matrix(x_names,private$ny,ncol(X),byrow = TRUE))
          
          x1 <- cbind(x1,X)
          xx <- NULL
          for (i in 1:private$ny) {
            xx <- rbind(xx,X)
          }
          private$xx <- matrix(xx,nrow = nrow(X))
        }
        
        cvs_names <- NULL
        if(!is.null(formula_cv)){
          mf <- model.frame(formula_cv, private$data_used)
          X <- as.matrix(model.matrix(formula_cv, mf)[,-1])
          if(qr(X)$rank != ncol(X)){
            stop('Regressors of formula_cv are singular!')
          }
          if(is.null(colnames(X))){
            cvs_names <- paste("cv",1:ncol(X))
          }else{
            cvs_names <- colnames(X)
          }
          
          x1 <- cbind(x1,X)
          private$cvs <- X
        }
        
        y <- model.response(mf)
        
        if(is.null(y1)){
          y1 <- as.vector(matrix(rbind(rep(0,private$nn),matrix(y,private$tt,private$nn)[-private$tt,]),ncol = 1))
        }else{
          y1 <- as.vector(y1)
        }
        
        private$y1 <- y1
        
        private$x1 <- x1
        
        if(isTRUE(timeFE)){
          timeFE <- kronecker(rep(1,private$nn),diag(private$tt))[,-c(1,2)]
          private$cvs <- cbind(private$cvs,timeFE)
          private$TFE <- TRUE
        }
        
        if(isFALSE(NoY)){
          private$yy_1 = matrix(y1,nrow = length(y1),ncol = private$ny)
          private$coefs_names <- c(purrr::map_chr(1:private$ny,function(x){paste("L1.Y",'regime',x)}),x_names,cvs_names)
        }else{
          private$coefs_names <- c("L1.Y",x_names,cvs_names)
        }
        
      }
      
      private$y <- y
      
      private$initx()
      
      private$omega = diag(2,nrow = (private$tt-1),ncol = (private$tt-1))
      diag(private$omega[-1,-(private$tt-1)]) = -1
      diag(private$omega[-(private$tt-1),-1]) = -1
      
      invisible(self)
    },
    
    #'@description
    #'Maximum likelihood estimation method
    #'@param ny the number of regimes
    MLE = function(ny=1) {
      
      if(is.null(private$delty0)){
        delty0 <- lag_transform(as.matrix(private$y),private$tt,private$nn,1,FALSE) -  
          lag_transform(as.matrix(private$y),private$tt,private$nn,1,TRUE)
      }else{
        delty0 <- private$delty0
      }
      
      deltxs =  lag_transform(private$x,private$tt,private$nn,1,FALSE) -  lag_transform(private$x,private$tt,private$nn,1,TRUE)
      deltxs = matrix(deltxs,nrow = private$tt-1)
      deltxs[1,] = 0
      deltxs = matrix(deltxs,nrow =  (private$tt-1)*private$nn)
      
      x1int = matrix(matrix(private$ininxc,nrow = private$tt-1)[1,],nrow = private$nn)
      yint = matrix(matrix(delty0,nrow = (private$tt-1))[1,],ncol = 1)
      
      
      #init
      delts <- qr.solve(x1int,as.vector(yint))
      
      deltxx <- cbind(deltxs,private$ininxc)
      
      dy0_iv <-  lag_transform(delty0,private$tt-1,private$nn,2,FALSE)
      
      dyy <- as.matrix(deltxs[,1:ny])
      ny1 <- ifelse(ncol(deltxs)==1,1,ny+1)
      dxx <- deltxs[,ny1:ncol(deltxs)]
      
      dx3 <-  lag_transform(as.matrix(dxx),private$tt-1,private$nn,2,FALSE)
      dy1 <- cbind( lag_transform(dyy,private$tt-1,private$nn,2,FALSE),dx3)
      dy2 <- cbind( lag_transform( lag_transform(dyy,private$tt-1,private$nn,1,FALSE),private$tt-2,private$nn,1,TRUE),dx3)
      
      leftiv <- t(dy2)%*%dy1
      
      betas_iv <- c(as.vector(MASS::ginv(leftiv)%*%(t(dy2)%*%dy0_iv)),delts)
      
      deltu <- as.vector(delty0 - deltxx%*%betas_iv)
      
      if(!is.null(private$var_u)){
        varu <- private$var_u
      }else{
        varu <-  sum(deltu^2)/(2*length(deltu))
      }
      
      pars <- (c(betas_iv,varu))
      
      varv1 <- (sum(yint^2)/private$nn)
      
      if(!is.null(private$w)){
        what <- private$w
      }else{
        what <- varv1/varu
      }
      
      pars <- c(pars,what)
      cd = length(pars)
      
      result_in <- suppressWarnings(try(stats::nlm(three_two,pars,delty0=delty0,evs=deltxx,
                                                   omega=private$omega,cd=cd,tt=private$tt,nn=private$nn,iterlim = private$iterlim)
                                        ,silent = TRUE))
      if("try-error" %in% class(result_in)){
        if(isTRUE(private$restart)){
          seeds0 <- sample(1:1000,1)
          set.seed(seeds0)
          pars[1:ny] = runif(ny,-1/ny,1/ny)
          pars[cd-1] = runif(1,2,5)
          pars[cd] = pars[cd-1]/varv1
          
          result_in = suppressWarnings(try(stats::nlm(three_two,pars,delty0=delty0,evs=deltxx,
                                                      omega=private$omega,cd=cd,tt=private$tt,nn=private$nn,iterlim = private$iterlim)
                                           ,silent = TRUE))
        }else{
          stop("\n","Encounter iteration failure!","\n",
               "Please set restart = True or check other inputs!","\n")
        }
      }
      
      if(result_in$code == 1 & isTRUE(private$restart)){
        
        seeds0 <- sample(1:1000,1)
        set.seed(seeds0)
        pars[1:ny] = runif(ny,-1/ny,1/ny)
        pars[cd-1] = runif(1,2,5)
        pars[cd] = pars[cd-1]/varv1
        
        result_in = suppressWarnings(try(stats::nlm(three_two,pars,delty0=delty0,evs=deltxx,
                                                    omega=private$omega,cd=cd,tt=private$tt,nn=private$nn,iterlim = private$iterlim)
                                         ,silent = TRUE))
        
      }
      
      if("try-error" %in% class(result_in)){
        stop("\n","Encounter iteration failure!","\n",
             "Please set restart = True or check other inputs!","\n")
      }
      
      private$pars <- result_in$estimate
      self$NNLL <- result_in$minimum
      
      evs  <- length(betas_iv)
      private$ccd <- evs - length(delts)
      if(isTRUE(private$TFE)){
        private$ccd <- evs - length(delts) - length(private$year) + 2
      }
      private$coefs <- private$pars[1:evs]
      
      bets <- private$pars[1:evs]
      
      self$duit <- (delty0 - deltxx%*%bets)
      self$dy0 <- delty0
      cms = matrix(0,nrow = evs,ncol = evs)
      
      ww <- private$pars[cd]
      
      xomega = private$omega
      xomega[1,1] = ww
      xomega = xomega*(private$pars[cd-1])
      xomega = inverse_cpp(xomega)
      
      for (i in 1:private$nn) {
        cms = cms + t(deltxx[((private$tt-1)*(i-1)+1):((private$tt-1)*i),])%*%xomega%*%(deltxx[((private$tt-1)*(i-1)+1):((private$tt-1)*i),])
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
      
      self$Zvalues = zvalues
      self$Ses= ses
      self$covariance_matrix = cms
      # private$dy0=delty0
      # private$xx=deltxx
      
      self$coefficients = private$coefs[1:length(private$coefs_names)]
      names(self$coefficients) <- private$coefs_names
      
      # 返回self以实现链式调用
      invisible(self)
    },
    
    #'@description
    #'Compute coefficients given thresholds
    #'@param ga thresholds
    TModel_fit = function(ga){
      
      rts <- matrix(private$q,nrow = private$nn*private$tt,ncol = (self$Th)) - matrix(ga,nrow  = private$nn*private$tt,ncol = self$Th,byrow = TRUE)
      rts[which(rts>0,arr.ind = TRUE)] <- 1
      rts[which(rts<=0,arr.ind = TRUE)] <- 0
      rts0 <- cbind(rep(1,private$nn*private$tt),rts)
      
      for (i in 1:(ncol(rts0)-1)) {
        rts0[,i] <- rts0[,i] - rts0[,i+1]
      }
      
      xxx <- NULL
      rts0x <- NULL
      if(!is.null(private$xx)){
        
        for (i in 1:(ncol(private$xx)/private$ny)) {
          rts0x <- cbind(rts0x,rts0)
        }
        xxx = private$xx*rts0x
      }
      
      xxx <- cbind(xxx,private$cvs)
      
      if(is.null(private$yy_1)){
        xxx <- cbind(private$y1,xxx)
      }else{
        xxx <- cbind(private$yy_1*rts0,xxx)
      }
      
      private$x <- xxx
      
      self$MLE(ny=private$ny)
      invisible(self)
    },
    
    #'@description
    #'Use MCMC to compute thresholds
    #'@param proportion the proportion of burning in the whole iterations
    #'@param types the type of MCMC, see BayesianTools::runMCMC
    #'@param ADs the parameter of MCMC, see BayesianTools::runMCMC
    #'@param nCR the parameter of MCMC, see BayesianTools::runMCMC
    #'@param ... the settings of MCMC, see BayesianTools::applySettingsDefault
    MCMC_process = function(proportion = 0.5,types = "DREAMzs",ADs = FALSE,nCR = 3,...){
      
      priors <- BayesianTools::createUniformPrior(lower = c(rep(private$r0,self$Th)), upper = c(rep(private$r1,self$Th)))
      
      bayesianSetup <- BayesianTools::createBayesianSetup(private$pf, prior = priors)
      
      if(types == "Metropolis"){
        settings = list(iterations = private$iterations,
                        burnin=ceiling(proportion*private$iterations),
                        adapt = ADs,adaptationNotBefore = private$iterations, nCR = nCR,...)
      }else{
        settings = list(iterations = private$iterations,
                        burnin=ceiling(proportion*private$iterations), nCR = nCR,...)
      }
      
      
      out <- BayesianTools::runMCMC(bayesianSetup = bayesianSetup,
                                    sampler = types,
                                    settings = settings)
      
      private$chains <- out
      private$MCMC_sample = BayesianTools::getSample(out,parametersOnly = F)
      mres = MAP2(out)
      self$thresholds = mres[[1]]
      names(self$thresholds) <- paste("Threshold",1:length(self$thresholds))
      
      self$NNLL = -mres[[2]][2]
      invisible(self)
      
    },
    
    #'@description
    #'print and print estimated results
    #'@param ... DPTM object
    print = function(...){
      coefs <- round(private$coefs,3)
      Zvalues <- self$Zvalues
      alpha_values <- round(abs(stats::qt(c(0.05,0.025,0.005),1e7)),3)
      
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
        
      }
      
      Zvalues <- round(Zvalues[1:length(coefs)],3)
      Ses <- round(self$Ses[1:length(coefs)],3)
      jgs <- cbind(coefs,Ses,Zvalues,xzx)
      rownames(jgs) <-  private$coefs_names[1:length(coefs)]
      colnames(jgs) <- c("Coefs","Std.Error","t-value","Significance")
      
      jgs <- jgs[1:private$ccd,]
      
      if(is.null(self$Th)){
        
        cat("\n","This is a Dynamic panel model with fixed effects (Linear).","\n")
        cat("---------------------------------------------------","\n")
        cat(" The coefficients are ","\n")
        print(jgs)
        
        cat("\n","Number of individuals:",private$nn,".\n")
        cat(" Number of periods:",private$tt,".\n")
        cat(" Time Fixed Effects:",private$TFE,".\n")
        cat(" Non-Negative log-likelihood:",self$NNLL,".\n")
        
      }else{
        IC = matrix(0,3,self$Th*2)
        
        colnames(IC) <- rep(c("Lower","Upper"),self$Th)
        rownames(IC) <- c("90%","95%","99%")
        
        
        for (i in 1:self$Th) {
          IC[1,(2*i-1):(2*i)] = round(c(stats::quantile(private$MCMC_sample[,i],0.05),stats::quantile(private$MCMC_sample[,i],0.95)),3)
          IC[2,(2*i-1):(2*i)] = round(c(stats::quantile(private$MCMC_sample[,i],0.025),stats::quantile(private$MCMC_sample[,i],0.975)),3)
          IC[3,(2*i-1):(2*i)] = round(c(stats::quantile(private$MCMC_sample[,i],0.005),stats::quantile(private$MCMC_sample[,i],0.995)),3)
          
        }
        
        chain <- private$chains$chain
        for (i in 1:length(chain)) {
          chain[[i]] <- chain[[i]][,1:self$Th]
        }
        MCMC_Convergence_Diagnostic <- try(coda::gelman.diag(chain,autoburnin=TRUE))
        if("try-error" %in% class(MCMC_Convergence_Diagnostic)){
          MCMC_Convergence_Diagnostic <- coda::gelman.diag(chain, autoburnin= FALSE)
        }
        
        cat("\n","This is a Dynamic panel multiple threshold model with fixed effects.","\n")
        cat("---------------------------------------------------","\n")
        cat(" The number of thresholds is",self$Th,".\n")
        cat(" The thresholds are","\n")
        print(self$thresholds)
        cat(" The CIs of thresholds are ","\n")
        print(IC)
        cat("\n","The coefficients are ","\n")
        print(jgs)
        
        cat("\n","Number of individuals:",private$nn,".\n")
        cat(" Number of periods:",private$tt,".\n")
        cat(" Time Fixed Effects:",private$TFE,".\n")
        cat(" Non-Negative log-likelihood:",self$NNLL,".\n")
        
        cat("\n"," The Gelman-Rubin Convergence Diagnostic is ","\n")
        print(MCMC_Convergence_Diagnostic)
        plot(private$chains)
        
      }
      
    }
  ),
  private = list(
    x1 = NULL,
    x = NULL,
    y = NULL,
    ininxc = NULL,
    id = NULL,
    year = NULL,
    omega = NULL,
    data_used = NULL, 
    ccd = NULL, 
    xx = NULL, 
    pars = NULL, 
    TFE = FALSE,
    coefs_names = NULL, 
    tt = NULL, 
    nn = NULL, 
    coefs = NULL,
    w=NULL,
    var_u = NULL,
    iterlim = 500,
    delty0 =NULL,
    restart = FALSE,
    #threshold
    r0 = NULL,
    r1 = NULL,
    ny = NULL,
    cvs = NULL,
    yy_1 = NULL,
    ssemin0 = Inf,
    q = NULL,
    y1 = NULL,
    iterations = 2000,
    sro = 0.1,
    MCMC_sample = NULL,
    chains = NULL,
    
    msx = function(ga){
      if(all(ga == sort(ga))){
        
        rts <- matrix(private$q,nrow = private$nn*private$tt,ncol = (self$Th)) - matrix(ga,nrow  = private$nn*private$tt,ncol = self$Th,byrow = TRUE)
        rts[which(rts>0,arr.ind = TRUE)] <- 1
        rts[which(rts<=0,arr.ind = TRUE)] <- 0
        rts0 <- cbind(rep(1,private$nn*private$tt),rts)
        
        for (i in 1:(ncol(rts0)-1)) {
          rts0[,i] <- rts0[,i] - rts0[,i+1]
        }
        
        if(any(colSums(rts0) < private$sro*(private$tt*private$nn))){
          sse0 = private$ssemin0
        }else{
          
          xxx <- NULL
          rts0x <- NULL
          if(!is.null(private$xx)){
            
            for (i in 1:(ncol(private$xx)/private$ny)) {
              rts0x <- cbind(rts0x,rts0)
            }
            xxx = private$xx*rts0x
          }
          
          xxx <- cbind(xxx,private$cvs)
          
          if(is.null(private$yy_1)){
            xxx <- cbind(private$y1,xxx)
          }else{
            xxx <- cbind(private$yy_1*rts0,xxx)
          }
          
          if(qr(xxx)$rank != ncol(xxx)){
            sse0 <- private$ssemin0
          }else{
            
            private$x <- xxx
            
            self$MLE(ny=private$ny)
            
            sse0 <- self$NNLL
            
          }
        }
      }else{
        sse0 = private$ssemin0
      }
      
      return(sse0)
    },
    
    pf = function(ga){
      mx <- try(private$msx(ga),silent = TRUE)
      NNLV <- ifelse("try-error" %in% class(mx),-private$ssemin0,-mx)
      return(NNLV)
    },
    
    initx = function(){
      ininxc <- matrix(0,nrow = (private$tt-1),ncol = private$nn)
      ininxc[1,] <- rep(1,private$nn)
      ininxc <- matrix(ininxc,ncol = 1)
      
      deltxi1 <-  lag_transform(private$x1,private$tt,private$nn,1,FALSE) -  lag_transform(private$x1,private$tt,private$nn,1,TRUE)
      deltxi1 <- matrix(deltxi1,nrow = private$tt-1)
      deltxi1[2:(private$tt-1),] <- 0
      deltxi1 <- matrix(deltxi1,nrow =  (private$tt-1)*private$nn)
      ininxc0 <- cbind(ininxc,deltxi1)
      
      x1int = matrix(matrix(ininxc0,nrow = private$tt-1)[1,],nrow = private$nn)
      if(qr(x1int)$rank == ncol(x1int)){
        ininxc <- ininxc0
      }else{
        QRininxc0 <- qr(as.matrix(ininxc0))
        ininxc <- ininxc0[,QRininxc0$pivot[1:(QRininxc0$rank)]]
      }
      private$ininxc <- ininxc
      
      invisible(self)
    },
    #xtset_b
    xtset_b = function(data,index=NULL){
      
      data <- as.data.frame(data)
      
      if(!is.null(index)){
        ids <- which(colnames(data) == index[1])[1]
        years <- which(colnames(data) == index[2])[1]
        
        if(length(ids) == 0 | length(years) == 0|is.na(ids)|is.na(years)){
          stop("Check the index!")
        }
        
        private$id <- unique(as.vector(data[,ids]))
        private$year <- sort(unique(as.numeric(as.vector(data[,years]))))
        
        private$data_used <- data[order(data[,ids],data[,years]),]
        rownames(private$data_used) <- NULL
        
      }else{
        private$id <- unique(as.vector(data[,1]))
        private$year <- as.numeric(unique(as.vector(data[,2])))
        private$data_used <- data
      }
      private$nn <- length(private$id)
      private$tt <- length(private$year)
      
      if(length(as.vector(data[,1])) != (private$nn*private$tt)){
        stop("Please check the index.")
      }
      
      invisible(self)
    }
  )
)
