### =================================================================================================== ###
###                                                                                                     ###
###   TVGJR:  All the code in this functions file uses custom-defined list structures for TV & GARCH    ###
###           To use this file effectively, please read the TV & GARCH object definitions carefully.    ###
###                                                                                                     ###
###           The main functions in this file are:                                                      ###
###           * EstimateTV()                                                                            ###
###           * EstimateGARCH()                                                                         ###
###           * CalculateTestStat_TV()                                                                  ###
###           * CalcProbabilityDist()                                                                   ###
###           Think of these as 'public methods'.  They provide a reasonable level of error-handling.   ###
###                                                                                                     ###
###                                                                                                     ###
###           Latest Ver.4.0.0                                                                          ###
###                                                                                                     ###
###     Change Log:  Ver.2.x.x
###         
###     27-Feb-2017 - List Objects changed, to remove need for position vectors                        
###     03-Mar-2017 - $var_target added to TV & GARCH objects
###                 - EstimateGARCH() modified to use compiled LogLik Fn for better performance   
###     10-Mar-2017 - Cut-down version of this file saved as: 'functions_tvgjr_WIG20.r'                 ###
###     18-Mar-2017 - EstimateTVGARCH() function finished & tested                                      ###
###
###     Change Log:  Ver.4.x.x
###         
###     20-Feb-2018 - CalculateTest_TV() changed to include a 'useRcpp' parameter
###                 - CalcProbabilityDist() changed to fix errors with length of runSimRow
###                 - dg.dt, myG & calculate_g re-written to improve performance
###                 - TV$linparsMode removed as it is not needed.
###                 - TV$Tobs (Total number of observations in sample) added.
###                 - Tests added for CCC_STCC correlation...  Work in progress
###
###                                                                                                     ###
### =================================================================================================== ###

### ======================================== Initialisation =========================================== ###
#

require(compiler)
setCompilerOptions(suppressAll = TRUE, optimize = 3)
enableJIT(3)
# We are using cmpfun() on various slow functions to speed them up!
#
library(moments)
library(Matrix)
#library(matrixcalc)
#
### ======================================== Initialisation =========================================== ###


### ======================================================== ###
###  Global Variables:
###
err_output <- -1e10

### ======================================================== ###


Main_Object_Definitions <- function() {
  #This is not a function! - Just making it visible in the list!
  
  ### ======================================================= ###
  ###        Data-Series list object definition               ###
  ### ======================================================= ###
  
  # Series <- list()
  # Series$tv <- list()       # TV object as defined below
  # Series$garch <- list()    # GARCH object as defined below
  # Series$var_target <- logical()
  # Series$tvparlimit <- vector()
  # Series$garchparlimit <- vector()
  # Series$value
  # Series$error
  
  ### ============================================== ###
  ###        TV list object definition               ###
  ### ============================================== ###
  
  # TV <- list()
  # TV$delta0 <- vector("numeric",1L)   # (MANDATORY) delta0 - single element vector (not integer) since we need to use it in vector functions 
  # TV$pars <- vector() or NULL         # delta1, speed1, loc1.1, <loc1.2>=NA, delta2, speed2, loc2.1, <loc2.2>,... loc2.n MUST provide NA if no value
  #										                     #    Note: All code expects a 4-elements-per-transition vector, so don't forget to include NA where there is no second location!
  # TV$linpars <- vector() or NULL      # [l1,l2,l3,l4] = constant, linear, quadratic, cubic terms. MUST provide NA if no value.  All code expects a 4-element vector!
  # TV$linparsMode <- numeric() or NULL # 1: Only the $linpars are free when optimising.  2: All TV$delta0, TV$pars & T$linpars are free.  NULL => $linpars are ignored
  #  TODO: remove all references to TV$linparsMode
  # TV$var_target <- logical()          # Variance Targeting switch. BOOL = FALSE => delta0 will NOT be used in estimation!!!  Always set = TRUE for this code.
  # TV$optimcontrol <- list()           # Control list passed to the optim() function.  Must be set for ALL parameters - Estimate.. functions will handle var targetting.
  # TV$speedoption <- integer()         # speedoptions are: 1=gamma, 2=gamma/sd(st), 3=exp(eta), 4=1/lambda^2
  # TV$shape <- vector() or NULL        # vector of integers describing the transitions (G's)
  #                                     #    E.g.  c(1,2,3) => 3 transitions  Must be NULL for TV Order0, as 'is.na(TV$shape)' is used to test for Order=0.
  #                                     #    1=Single location, 2=Double location, 3=Single location, with Squared transition
  # TV$st <- vector()                   # (MANDATORY) smooth transition variable.  Typically = seq(1,Tobs)/Tobs, i.e. linear transition
  #					    	                      #    This is mandatory, because the code tests for 'length(TV$st)' to determine the number of data observations
  #									                    #    Not a great idea, but too late to refactor now!  TODO: replace length(TV$st) with TV$Tobs
  # TV$condvars <- vector()             # vector of conditional variances, 'g'
  ##
  ## -- Elements returned from the EstimateTV() function: -- ##
  ##
  # TV$value <- numeric()               # numeric value returned from the optim() function for the log-liklihood value
  # TV$hess <- matrix()                 # Hessian matrix of TV$pars, returned from optim() or optimHess()
  # TV$error <- logical()               # True/False.  Set to TRUE when optim fails, False otherwise.
  # TV$stderr <- matrix()               # Residual Std Err matrix calculated from the returned Hessian matrix - Often causes a numerical error
  
  
  ### ============================================== ###
  ###       GARCH list object definition             ###
  ### ============================================== ###
  
  # GARCH <- list()
  # GARCH$type <- numeric()          # 0: Constant GARCH, 1:GARCH, 2:GJR-GARCH, 3: GJR-GARCH with omega=0
  # GARCH$pars <- vector()           # omega, alpha, beta, delta  (Single-order Garch only implemented)
  #										               #    Note: All the code relies on all 4 parameters being provided, so don't forget to include NA where needed (Garch Types 0,1 & 3)!!
  # GARCH$var_target <- logical()    # Variance Targeting switch. BOOL = FALSE => All $pars will be used in estimation. (TRUE => omega will be calculated, not free)
  #										               #    Note: Variance targetting is not used in this code, but comes into play for TV with GARCH estimation.
  # GARCH$optimcontrol <- list()     # Control list passed to the optim() function.  Must be set for ALL parameters - Estimate.. functions will handle var targetting.
  ##
  ## -- Elements returned from the EstimateGARCH() function: -- #
  ##
  # GARCH$value <- numeric()         # numeric value returned from the optim() function for the log-liklihood value
  # GARCH$condvars <- vector()       # vector of conditional variances 'h'
  # GARCH$error <- logical()         # True/False.  Set to TRUE when optim fails, False otherwise
  # GARCH$hess <- matrix()           # Hessian matrix of GARCH$pars, returned from optim() or optimHess()
  # GARCH$stderr <- matrix()         # Residual Std Err matrix calculated from the returned Hessian matrix - Often causes a numerical error
  
  
  ### ============================================== ###
  ###         CORR list object definition            ###
  ### ============================================== ###
  
  # CORR <- list()                 # This object is just a wrapper container for all the specific correlation models
  # CORR$type <- vector(char)      # "CCC", "STCC", "STEC", ...
  # CORR$CCC <- list()
  # CORR$STCC <- list()
  # CORR$STEC <- list()
  #
  
  ### ============================================== ###
  ###         STCC list object definition            ###
  ### ============================================== ###
  
  # STCC <- list()
  # STCC$type <- as.numeric()             # 0:"CCC", 1:"STCC",...
  # STCC$P1pars <- vector()
  # STCC$P2pars <- vector()
  # STCC$pars <- as.vector()              # A vector of the remaining pars: Speed, Loc1, <Loc2>  
  # STCC$speedoption <- as.integer()      # speedoptions are: 1=gamma 2=gamma/sd(st) 3=exp(eta) 4=...
  # STCC$st <- as.vector()                # transition variable.  Default = seq(1,Tobs)/Tobs, i.e. linear
  # STCC$shape <- as.vector()             # vector of integers describing the transitions (G's)
  #                                       # e.g. c(NA) => No transitions.  c(1,2,3) => 3 transitions
  #                                       # 1=Single location, 2=Double location, 3=Single location, with Squared transition
  # STCC$condcorrs <- as.matrix()         # time-series of conditional correlations, 'rho'  (Size = T x N(N-1)/2)
  # STCC$stderr <- as.matrix()            # Residual Std Err matrix calculated from the returned Hessian matrix  
  # -- Elements returned from the Optim() function: -- #
  # STCC$value <- as.numeric()            # numeric value returned from the optim() function for the log-liklihood value
  # STCC$hess <- as.matrix()              # Hessian matrix returned from optim() (or optimHess())
  
  
  ### ============================================== ###
  ###         STEC list object definition            ###
  ### ============================================== ###
  
  # STEC <- list()
  # STEC$type <- as.numeric()             # 0:"CEC" (only 1 par = p1),  1:"STEC" (up to 5 pars)
  # STEC$pars <- as.vector()              # A vector containing the parameters, p1, p2, Speed, Loc1, <Loc2> (all scalars)
  # STEC$linpars <- as.vector() or NULL   # linpar1, linpar2,...  Note: NULL (or non-existing) => All functions will operate with linearised = FALSE
  # STEC$speedoption <- as.integer()      # speedoptions are: 1=gamma 2=gamma/sd(st) 3=exp(eta) 4=...
  # STEC$st <- as.vector()                # transition variable.  Default = seq(1,Tobs)/Tobs, i.e. linear
  # STEC$shape <- as.vector()             # vector of integers describing the transitions (G's)
  #                                       # e.g. c(NA) => No transitions.  c(1,2,3) => 3 transitions
  #                                       # 1=Single location, 2=Double location, 3=Single location, with Squared transition
  # STEC$condcorrs <- as.matrix()         # time-series of conditional correlations, 'rho'  (Size = T x N(N-1)/2)
  # STEC$stderr <- as.matrix()            # Residual Std Err matrix calculated from the returned Hessian matrix
  # -- Elements returned from the Optim() function: -- #
  # STEC$value <- as.numeric()            # numeric value returned from the optim() function for the log-liklihood value
  # STEC$hess <- as.matrix()              # Hessian matrix returned from optim() (or optimHess())
  
  
  
  ### ++++++++++++ MULTIVARIATE OBJECTS - START ++++++++++++ ###
  
  #nTV <- list(N)              # List of TV() list objects, where N = total number of objects
  #nTV[[n]]                    # Reference to one specific TV() list object
  #nTV$delta0 <- as.vector()   # Vector containing c(TV[[1]]$delta0, TV[[2]]$delta0,..., TV[[N]]$delta0)
  #nTV$pars <- as.vector()     # Vector containing c(TV[[1]]$pars, TV[[2]]$pars,..., TV[[N]]$pars)
  #nTV$hess <- as.matrix       # Hessian matrix of all TV$pars, returned from optim() (or optimHess())
  #nTV$stderr <- as.matrix     # Residual Std Err matrix calculated from the returned Hessian matrix
  
  #nGARCH <- list(N)              # List of GARCH() list objects, where N = total number of objects
  #nGARCH[[n]]                    # Reference to one specific GARCH() list object
  #nGARCH$pars <- as.vector()     # Vector containing c(GARCH[[1]]$pars, GARCH[[2]]$pars,..., GARCH[[N]]$pars)
  #nGARCH$hess <- as.matrix       # Hessian matrix of all GARCH$pars (with Omega removed), returned from optim() (or optimHess())
  #nGARCH$stderr <- as.matrix     # Residual Std Err matrix calculated from the returned Hessian matrix
  
  
  # Note: STCC is by definition a multivariate object.  We cannot have a univariate STCC
  
  ### ++++++++++++  MULTIVARIATE OBJECTS - END  ++++++++++++ ###
  
}

#----------------------------------------------------------#
#               Utility Functions 
#----------------------------------------------------------#
zUtilFn <- TRUE
if (zUtilFn) {
  
fit_acf <- function(pars,ac_data,rtn_est=FALSE) {
  # INPUT:  pars = vector of: [alpha,beta]
  #         ac_data = vector of 'actual' AutoCorr data

  maxLag <- length(ac_data)
  ac_est <- vector()  # The calculated estimate of Auto Corr
  a <- pars[1]  #alpha
  b <- pars[2]  #beta
  rtn <- 1e10  #large value to indicate error to optim

  # Parameter Boundary Checks:
  if (a < 0) return(rtn)
  if (b < 0) return(rtn)
  if (a+b > 1) return(rtn)

  ac_1 <- a*(1-b^2-a*b)/(1-b^2-2*a*b)

  for (t in 2:maxLag) {
    ac_est[t-1] <- ac_1*((a+b)^(t-1))
  }
  
  if (rtn_est) return(ac_est) else {
    act_est <- (ac_data[2:maxLag]-ac_est)^2
    rtn <- sum(act_est)
  }
  
  #return:
  rtn
}
#
fit_acf_c <- cmpfun(fit_acf)


fit_acf_1 <- function(pars,ac_data) {
  # INPUT:  pars = vector of: [decay_rate]
  #         ac_data = vector of 'actual' AutoCorr data
  
  maxLag <- length(ac_data)
  ac_est <- vector()   # The calculated estimate of Auto Corr
  ac_t1 <- ac_data[1]  #ac_t1
  decay <- pars[1]  #decay rate
  rtn <- 1e10  #large value to indicate error to optim
  
  # Parameter Boundary Checks:
  # Try fixing ac_t1 to be the first actual parameter?
  #ac_t1 <- ac_data[1]
  if (decay < 0) return(rtn)
  if (decay > 1) return(rtn)
  
  ac_est[1] <- ac_t1
  for (t in 2:maxLag) {
    ac_est[t-1] <- ac_t1*(decay^(t-1))
  }
  #return:
  act_est <- (ac_data-ac_est)^2
  rtn <- sum(act_est)
  
}

fit_acf_2 <- function(pars,ac_data) {
  # INPUT:  pars = vector of: [ac_t1,decay_rate]
  #         ac_data = vector of 'actual' AutoCorr data
  
  maxLag <- length(ac_data)
  ac_est <- vector()  # The calculated estimate of Auto Corr
  ac_t1 <- pars[1]  #ac_t1
  decay <- pars[2]  #decay rate
  rtn <- 1e10  #large value to indicate error to optim
  
  # Parameter Boundary Checks:
  # keep params positive & less than 1
  if (decay < 0) return(rtn)
  if (decay > 1) return(rtn)
  if (ac_t1 < 0) return(rtn)
  if (ac_t1 > 1) return(rtn)
  
  ac_est[1] <- ac_t1
  for (t in 2:maxLag) {
    ac_est[t-1] <- ac_t1*(decay^(t-1))
  }
  #return:
  act_est <- (ac_data-ac_est)^2
  rtn <- sum(act_est)
  
}
  
  
validateTV <- function(tv) {
  ### Test if the passed-in tv list object is valid ###
  
  errmsg <- vector("character")
  
  #1: Validate $delta0:
  if (is.null(tv$delta0)) {
    errmsg <- "You must provide $delta0! "
    cat("\n", errmsg, "\nYour object is:\n")
    return(tv)  
  } else { 
    if (length(tv$delta0)>1) errmsg <- c(errmsg, "$delta0 must be a single value vector! ")
    else if (is.na(tv$delta0[1])) errmsg <- c(errmsg, "$delta0 must be a valid number, cannot be NA")
  }
  
  #2: Validate $pars: when $pars are provided
  if (!is.null(tv$pars)) {
    errmsg <- c(errmsg, "\nYou have provided $pars...\n ")
    if (is.null(tv$shape)) errmsg <- c(errmsg,"but no $shape\n") else {
      if (is.na(tv$shape)) errmsg <- c(errmsg,"but $shape is NA\n")
      if (min(tv$shape)<1 || max(tv$shape)>3) errmsg <- c(errmsg,"but $shape contains invalid data - only 1, 2, or 3 supported\n")
      if (length(tv$pars)/length(tv$shape)!=4) errmsg <- c(errmsg,"but the number of pars doesn't match the $shape provided \n  (Did you forget to include a NA value for loc2?)\n")
    }
    if (is.null(tv$st)) errmsg <- c(errmsg,"but no $st\n") else {
      if (is.na(tv$shape)) errmsg <- c(errmsg,"but $st is NA\n")
    }
    if (is.null(tv$speedoption)) errmsg <- c(errmsg,"but no $speedoption\n") else {
      if (is.na(tv$speedoption)) errmsg <- c(errmsg,"but $speedoption is NA\n")
      if (min(tv$speedoption)<1 || max(tv$speedoption)>3) errmsg <- c(errmsg,"but $speedoption contains invalid data - only 1, 2, or 3 supported\n")
    }
    cat("\n", errmsg, "\nYour object is:\n")
    return(tv)
  }
  
  #3: Validate $pars: when $pars are Not provided
  if (is.null(tv$pars)) {
    errmsg <- c(errmsg, "\nYou have not provided $pars...\n ")
    if (!is.null(tv$shape)) errmsg <- c(errmsg,"but have provided $shape\n") 
    if (!is.null(tv$st)) errmsg <- c(errmsg,"but have provided $st\n") 
    if (!is.null(tv$speedoption)) errmsg <- c(errmsg,"but have provided $speedoption\n") 
    cat("\n", errmsg, "\nYour object is:\n")
    return(tv)
  }
  
  
  cat("\n", errmsg, "\nYour object is:\n")
  return(tv)
  
}  #End: Validate(TV)

matrix.sort <- function (mat, sort, decreasing = FALSE) {
  # This function will sort a matrix by a column
  m <- do.call("order", c(as.data.frame(mat[, sort]), decreasing = decreasing))
  mat[m, ]
}

unregisterDoParallel <- function(cluster=NULL) {
  # This function will release all variables & memory tied up by a call to foreach(..)
  # Very useful when testing, as the environment builds up quickly!
  if(!is.null(cluster)) stopCluster(cluster)
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

setParScale <- function(delta,TV) {
  ###
  ###  Set up the parameter scaling for optim.  This is Critical for accuracy!
  ###
    
  # handle variance targeting:
  if (is.null(TV$var_target)) TV$var_target <- TRUE  #Default value for TV
  
	if (is.null(TV$pars)) parScale <- ceiling(TV$delta0) else {
		# Divide: d0/TV$delta0, d/delta, sp/maxSpeed, l1/1, l2/1
		maxSpeed <- switch(TV$speedoption,1000,(1000/sd(TV$st)),7.0,0.30)
		parScale <- rep(c(delta,maxSpeed,1,1),length(TV$shape))
		#Now we need to remove the elements that are NA in TV$pars - to match optimpars
		parScale <- parScale[!is.na(TV$pars)]
		if (isTRUE(TV$var_target)) parScale <- c(ceiling(TV$delta0), parScale)
	}
  
  # Now handle linearised parameters:
  if(!is.null(TV$linpars)) {
    # Remove the NA elements
    linPars <- TV$linpars[!is.na(TV$linpars)]
    # Set parScale to match the length of optimpars
    parScale <- c(parScale,ceiling(linPars))
  }
  
  #Return:
  parScale
}

calculateOmega <- function(alpha,beta,delta=0) {
  # When VarienceTargetting is ON: omega is calculated as below:
  # Since the 'delta' param is only used for Garch Type=2 (GJR), we default it to 0
  #Return:
  omega <- 1-alpha-beta-0.5*delta
}

calculate_local_variance <- function(e,localvarwindow=300) {
  # Local Varience Targetting is ON: We can estimate Garch parameters, while minimising the effect of shock-persistence on alpha & beta
  # omega can then be calculated as: localVar * (1-alpha-beta-0.5*delta)
  
  Tobs <- length(e)
  localVar <- vector("numeric",Tobs)
  
  for (t in 1:Tobs) {
    if (t < (localvarwindow/2)) {
      localVar[t] <- var(e[1:localvarwindow])
    } else if (t > (Tobs - localvarwindow/2)) {
      localVar[t] <- var(e[(Tobs - localvarwindow):Tobs])
    } else localVar[t] <- var(e[(t - localvarwindow/2):(t + localvarwindow/2)])
  }
    
  # Return:
  localVar
}

calculate_local_kurtosis <- function(e,localkurtwindow=300) {
  Tobs <- length(e)
  localKurt <- vector("numeric",Tobs)
  
  for (t in 1:Tobs) {
    if (t < (localkurtwindow/2)) {
      localKurt[t] <- kurtosis(e[1:localkurtwindow])
    } else if (t > (Tobs - localkurtwindow/2)) {
      localKurt[t] <- kurtosis(e[(Tobs - localkurtwindow):Tobs])
    } else localKurt[t] <- kurtosis(e[(t - localkurtwindow/2):(t + localkurtwindow/2)])
  }
  
  # Return:
  localKurt
}


garchpars_to_optimpars <- function(garch) {
  ## Call this on a Garch object to get the required optimpars, which you can then pass to the optim() function.
  
  # Check if variance targetting is being used. If not, set default to optimise ALL parameters.
  if (is.null(garch$var_target)) garch$var_target <- FALSE
  
  optimpars <- vector("numeric")
  if (!garch$var_target) {
    # VarienceTargetting is OFF: All Garch parameters are optimised
    if (garch$type==1) {
      optimpars[1] <- garch$pars[1]
      optimpars[2] <- garch$pars[2]
      optimpars[3] <- garch$pars[3]
    } else if (garch$type==2) {
      optimpars[1] <- garch$pars[1]
      optimpars[2] <- garch$pars[2]
      optimpars[3] <- garch$pars[3]
      optimpars[4] <- garch$pars[4]
    } else if (garch$type==3) {
      optimpars[1] <- garch$pars[1]
      optimpars[2] <- garch$pars[3]
      optimpars[3] <- garch$pars[4]
    }
  } else {
    # VarienceTargetting is ON: omega is not passed in as a free param - it is calculated inside the Loglik function.
    if (garch$type==1) {
      optimpars[1] <- garch$pars[2]
      optimpars[2] <- garch$pars[3]
    } else if (garch$type==2) {
      optimpars[1] <- garch$pars[2]
      optimpars[2] <- garch$pars[3]
      optimpars[3] <- garch$pars[4]
    } else if (garch$type==3) {
      optimpars[1] <- garch$pars[3]
      optimpars[2] <- garch$pars[4]
    }
  }
  #Return:
  optimpars
  
}  # End: garchpars_to_optimpars

optimpars_to_garchpars <- function(garch,optimpars) {
  ## Call this on a Garch object, with the optimpars returned from the optim() function.
  ## This function will work out which parameter is which and put them back into Garch$pars.
  
  # Check if variance targetting is being used. If not, set default to optimise ALL parameters.
  if (is.null(garch$var_target)) garch$var_target <- FALSE
  
  pars <- vector("numeric",4L)
  if (!garch$var_target) {
    # VarienceTargetting is OFF: All Garch parameters are optimised
    if (garch$type==1) {
      pars[1] <- optimpars[1]
      pars[2] <- optimpars[2]
      pars[3] <- optimpars[3]
      pars[4] <- NA      
    } else if (garch$type==2) {
      pars[1] <- optimpars[1]
      pars[2] <- optimpars[2]
      pars[3] <- optimpars[3]
      pars[4] <- optimpars[4]
    } else if (garch$type==3) {
      pars[1] <- optimpars[1]
      pars[2] <- NA
      pars[3] <- optimpars[2]
      pars[4] <- optimpars[3]
    }
  } else {
    # VarienceTargetting is ON: omega is not passed in as a free param - it is calculated inside the Loglik function.
    if (garch$type==1) {
      pars[1] <- calculateOmega(optimpars[1],optimpars[2])
      pars[2] <- optimpars[1]
      pars[3] <- optimpars[2]
      pars[4] <- NA
    } else if (garch$type==2) {
      pars[1] <- calculateOmega(optimpars[1],optimpars[2])
      pars[2] <- optimpars[1]
      pars[3] <- optimpars[2]
      pars[4] <- optimpars[3]
    } else if (garch$type==3) {
      alpha <- 0
      pars[1] <- calculateOmega(alpha,optimpars[1],optimpars[2])
      pars[2] <- NA
      pars[3] <- optimpars[1]
      pars[4] <- optimpars[2]
    }
  }
  #Return:
  pars
} # End: optimpars_to_garchpars

displayProbabilities <- function(probDist,TestValues,kurt,pers,subtitle="") {
  ##  Just a quick-n-dirty function to calculate and display p-values & probability histograms,
  ##  from the data created by the Calculate
  
  require(graphics)
  
  cat("\n")
  print("Histograms of the H01, H02 & H03 Probability Distributions have been plotted.")
  
  testScenario <- paste0("K=",kurt," pers=0.",pers)
  histogramBreaks <- 40
  
  #H03 Case:
  Pval_TR2 <- mean(probDist[,"Pval_TR2.H03"],na.rm = TRUE)
  ValidTR2 <- probDist[,"Stat_TR2.H03"]
  ValidTR2 <- ValidTR2[!is.na(ValidTR2)]
  numValidTR2 <- length(ValidTR2)
  ptitle <- paste("LM3 for TR2",testScenario, sep = " - ")
  hist(ValidTR2, breaks = histogramBreaks, main="")
  title(main=ptitle, sub=subtitle)
  #
  Pval_Robust <- mean(probDist[,"Pval_Robust.H03"],na.rm = TRUE)
  ValidRobust <- probDist[,"Stat_Robust.H03"]
  ValidRobust <- ValidRobust[!is.na(ValidRobust)]
  numValidRobust <- length(ValidRobust)
  ptitle <- paste("LM3 for ROBUST",testScenario, sep = " - ")
  hist(ValidRobust, breaks = histogramBreaks, main="")
  title(main=ptitle, sub=subtitle)
  #
  cat("\n ", testScenario," (LM3) LM-TR2 =",TestValues$LMTR2.H03,"/ p =",Pval_TR2,"/ Valid TR2 Tests =", numValidTR2, " ||  LM-Rob =",TestValues$LMRobust.H03,"/ p =",Pval_Robust,"/ Valid Robust Tests =", numValidRobust)
  
  #H02 Case:
  Pval_TR2 <- mean(probDist[,"Pval_TR2.H02"],na.rm = TRUE)
  ValidTR2 <- probDist[,"Stat_TR2.H02"]
  ValidTR2 <- ValidTR2[!is.na(ValidTR2)]
  numValidTR2 <- length(ValidTR2)
  ptitle <- paste("LM2 for TR2",testScenario, sep = " - ")
  hist(ValidTR2, breaks = histogramBreaks, main="")
  title(main=ptitle, sub=subtitle)
  #
  Pval_Robust <- mean(probDist[,"Pval_Robust.H02"],na.rm = TRUE)
  ValidRobust <- probDist[,"Stat_Robust.H02"]
  ValidRobust <- ValidRobust[!is.na(ValidRobust)]
  numValidRobust <- length(ValidRobust)
  ptitle <- paste("LM2 for ROBUST",testScenario, sep = " - ")
  hist(ValidRobust, breaks = histogramBreaks, main="")
  title(main=ptitle, sub=subtitle)
  #
  cat("\n ", testScenario," (LM2) LM-TR2 =",TestValues$LMTR2.H02,"/ p =",Pval_TR2,"/ Valid TR2 Tests =", numValidTR2, " ||  LM-Rob =",TestValues$LMRobust.H02,"/ p =",Pval_Robust,"/ Valid Robust Tests =", numValidRobust)
  
  #H01 Case:
  Pval_TR2 <- mean(probDist[,"Pval_TR2.H01"],na.rm = TRUE)
  ValidTR2 <- probDist[,"Stat_TR2.H01"]
  ValidTR2 <- ValidTR2[!is.na(ValidTR2)]
  numValidTR2 <- length(ValidTR2)
  ptitle <- paste("LM1 for TR2",testScenario, sep = " - ")
  hist(ValidTR2, breaks = histogramBreaks, main="")
  title(main=ptitle, sub=subtitle)
  #
  Pval_Robust <- mean(probDist[,"Pval_Robust.H01"],na.rm = TRUE)
  ValidRobust <- probDist[,"Stat_Robust.H01"]
  ValidRobust <- ValidRobust[!is.na(ValidRobust)]
  numValidRobust <- length(ValidRobust)
  ptitle <- paste("LM1 for ROBUST",testScenario, sep = " - ")
  hist(ValidRobust, breaks = histogramBreaks, main="")
  title(main=ptitle, sub=subtitle)
  #
  cat("\n ", testScenario," (LM1) LM-TR2 =",TestValues$LMTR2.H01,"/ p =",Pval_TR2,"/ Valid TR2 Tests =", numValidTR2, " ||  LM-Rob =",TestValues$LMRobust.H01,"/ p =",Pval_Robust,"/ Valid Robust Tests =", numValidRobust)
  
  cat("\n")
  
} # End: displayProbabilities

displayProbabilities_1 <- function(probDist,TestValues,confidenceLevel,subtitle="") {
  ##  Just a quick-n-dirty function to calculate and display p-values & probability histograms,
  ##  from the data created by the Calculate
  
  require(graphics)
  
  cat("\n")
  print("Histograms Probability Distributions have been plotted.")
  
  confidenceLevel <- as.character(confidenceLevel)
  testScenario <- paste0("K=3 a+b=0.",confidenceLevel)
  histogramBreaks <- 50
  
  Pval_TR2 <- mean(probDist[,"Pval_TR2"],na.rm = TRUE)
  ValidTR2 <- probDist[,"Stat_TR2"]
  ValidTR2 <- ValidTR2[!is.na(ValidTR2)]
  numValidTR2 <- length(ValidTR2)
  ptitle <- paste("TR2",testScenario, sep = " - ")
  hist(ValidTR2, breaks = histogramBreaks, main="")
  title(main=ptitle, sub=subtitle)
  #
  Pval_Robust <- mean(probDist[,"Pval_Robust"],na.rm = TRUE)
  ValidRobust <- probDist[,"Stat_Robust"]
  ValidRobust <- ValidRobust[!is.na(ValidRobust)]
  numValidRobust <- length(ValidRobust)
  ptitle <- paste("ROBUST",testScenario, sep = " - ")
  hist(ValidRobust, breaks = histogramBreaks, main="")
  title(main=ptitle, sub=subtitle)
  #
  cat("\n", testScenario, ", ",subtitle)
  cat("\n LM-TR2 =",TestValues$LMTR2,"/ p =",Pval_TR2,"/ Valid TR2 Tests =", numValidTR2)
  cat("\n LM-Rob =",TestValues$LMRobust,"/ p =",Pval_Robust,"/ Valid Robust Tests =", numValidRobust)
  cat("\n")
  
}

runSpecificationTests <- function(e,tv,tv_startvalue,testorder,data_start,data_end,numloops) {
###
###  This sub routine encapsulates a section of common code that needs to be run many times during
###  the Function specification phase.
###
  
  TV <- tv
  
  functionOrder <- paste0("TV",as.character(length(TV$shape)))  #e.g. TV0 or TV1

  ptitle <- paste0("WIG20 - ",functionOrder, " (",data_start,"-",data_end,") - ",testorder)
  plot(TV$condvars,type='l',main=ptitle)
  
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\n\n")
  
  refTestValues <- list()
  refTestValues$LMTR2 <- CalculateTestStat_TV(e,TV,testorder,usetest="TR2")
  refTestValues$LMRobust <- CalculateTestStat_TV(e,TV,testorder,usetest="ROBUST")
  
  #Now setup which tests will be used in the CalcProbabilityDist() function
  useTest <- c("NONE")
  if (!is.na(refTestValues$LMTR2)) useTest <- c("TR2")
  if (!is.na(refTestValues$LMRobust)) useTest <- c(useTest,"ROBUST")
  if (length(useTest)==1) stop("Couldn't calculate any of the Reference Test Values")
  
  ## FOR a+b=0.95
  
  # ----- Set the parameters to pass to CalcProbabilityDist()  ------ #
  # Reset the TV parameters & TV value to the same starting values used to generate RefTest:
  TV$delta0 <- tv_startvalue$delta0
  TV$pars <- tv_startvalue$pars
  TV$linpars <- tv_startvalue$linpars
  TV$value <- NULL
  
  saveAs <- paste0("ProbDist_",functionOrder,"_GarchRealData_SimuK3-95_Data",data_start,"-",data_end,"_",testorder,".RDS")
  saveAs <- file.path("Output",saveAs)
  
  # # ----- Load the '95' generated data with Garch  ------ #
  RefData_WithGarch <- readRDS("GeneratedRefData_WithGarch_K3-95_T4777_Rows2000_Seed5000.RDS") # 4777 rows, 7000 cols =  observations of GARCH data
  RefData_WithGarch <- RefData_WithGarch[data_start:data_end,1:numloops]
  
  ### Generate the Probability Distribution and save to a file:
  if (file.exists(saveAs)) probDist95 <- readRDS(saveAs) else {
    timestamp(prefix = "##-- ",suffix = " --##")  
    probDist95 <- CalcProbabilityDist(TV,RefData_WithGarch,refTestValues,useTest,testorder,saveAs,numCores,numloops) 
    timestamp(prefix = "##-- ",suffix = " --##")
  }

  #     ----- Calculate p-value of refTestValues using TestValue Distribution from RefData      ------ #
  #       Note: The data matrix returned from CalcProbabilityDist is different in 'functions_tvgjr_v2'
  #       If a Test Failed, the TestStat value will be NA (and the equivalent P_Val = 0)
  #       Pval = Sum of all Pval_<Test> / Number of valid Tests, i.e. Mean(Pval, na.rm=TRUE)
  displayProbabilities_1(probDist95,refTestValues,95,ptitle)
  

  ## FOR a+b=0.97
  
  # ----- Set the parameters to pass to CalcProbabilityDist()  ------ #
  # Reset the TV parameters & TV value to the same starting values used to generate RefTest:
  TV$delta0 <- tv_startvalue$delta0
  TV$pars <- tv_startvalue$pars
  TV$linpars <- tv_startvalue$linpars
  TV$value <- NULL
  
  saveAs <- paste0("ProbDist_",functionOrder,"_GarchRealData_SimuK3-97_Data",data_start,"-",data_end,"_",testorder,".RDS")
  saveAs <- file.path("Output",saveAs)
  
  # # ----- Load the '97' generated data with Garch  ------ #
  RefData_WithGarch <- readRDS("GeneratedRefData_WithGarch_K3-97_T4777_Rows2000_Seed5000.RDS") # 4777 rows, 7000 cols =  observations of GARCH data
  RefData_WithGarch <- RefData_WithGarch[data_start:data_end,1:numloops]
  
  ### Generate the Probability Distribution and save to a file:
  if (file.exists(saveAs)) probDist97 <- readRDS(saveAs) else {
    timestamp(prefix = "##-- ",suffix = " --##")  
    probDist97 <- CalcProbabilityDist(TV,RefData_WithGarch,refTestValues,useTest,testorder,saveAs,numCores,numloops) 
    timestamp(prefix = "##-- ",suffix = " --##")
  }
 
  displayProbabilities_1(probDist97,refTestValues,97,ptitle)
  
  rm(RefData_WithGarch)
}

runLM123Tests <- function(e,tv,data_start,data_end,numloops,series_name,kurt,pers) {
  ###
  ###  This sub routine encapsulates a section of common code that needs to be run many times during
  ###  the Function specification phase, but using the comparative LM tests.
  ###
  
  TV <- tv
  
  functionOrder <- paste0("TV",as.character(length(TV$shape)))  #e.g. TV0 or TV1
  
  ptitle <- paste0(series_name," - ",functionOrder, " (",data_start,"-",data_end,") - LM123")
  plot(TV$condvars,type='l',main=ptitle)
  
  cat("\n", ptitle, "Logliklihood Value: ", TV$value, "\n\n")
  
  # Set up what test statistics are being used: TR2, Robust, order of approximation (3 -> H03, 2 -> H02, 1 -> H01)
  useTest <- "TR2"
  refTestValues <- list()
  refTestValues$LMTR2.H01 <- CalculateTestStat_TV(e,TV,"H01",useTest)
  refTestValues$LMTR2.H02 <- CalculateTestStat_TV(e,TV,"H02",useTest)
  refTestValues$LMTR2.H03 <- CalculateTestStat_TV(e,TV,"H03",useTest)
  useTest <- "ROBUST"
  refTestValues$LMRobust.H01 <- CalculateTestStat_TV(e,TV,"H01",useTest)
  refTestValues$LMRobust.H02 <- CalculateTestStat_TV(e,TV,"H02",useTest)
  refTestValues$LMRobust.H03 <- CalculateTestStat_TV(e,TV,"H03",useTest)
  
  
  # ----- Set the parameters to pass to CalcProbabilityDist()  ------ #
  # K= 'kurt' P= 'pers':  garch$pars = (o,a,b,NA)
  kurtPers <- paste0("K",kurt,"_P",pers)
  dataRange <- paste0("_T",data_start,"-",data_end)
  
  # Reset the TV parameters to the same starting values used to generate RefTest:
  TV$delta0 <- tv$delta0
  TV$pars <- tv$pars
  TV$value <- NULL
  TV$error <- NULL
  
  saveAs <- paste0(series_name,"-ProbDist_",functionOrder,"_GarchRealData_Simu",kurtPers,dataRange,".RDS")
  saveAs <- file.path("Output",saveAs)
  
  TestOrder <- NULL  # This will force all Tests to be run!
  useTest <- NULL    # This will force all Tests to be run!
  
  # # ----- Load the generated data with Garch  ------ #
  dataFile <- paste0("GeneratedRefData_WithGarch_",kurtPers,"_T6149_Rows2000_Seed7000.RDS")
  RefData_WithGarch <- readRDS(dataFile) # 6149 rows, 2000 cols =  observations of GARCH data
  RefData_WithGarch <- RefData_WithGarch[data_start:data_end,1:numloops]
  
  ### Generate the Probability Distribution and save to a file:
  timestamp(prefix = "##-- ",suffix = " --##")
  #probDist <- CalcProbabilityDist_LM(TV,RefData_WithGarch,refTestValues,useTest,TestOrder,saveAs,numCores,numLoops) 
  timestamp(prefix = "##-- ",suffix = " --##")
  
  ### Read it back from a file, if you've already generated it!
  probDist <- readRDS(saveAs)
  
  #     ----- Calculate p-value of refTestValues using TestValue Distribution from RefData      ------ #
  #       Note: The data matrix returned from CalcProbabilityDist is different in 'functions_tvgjr_v2'
  #       If a Test Failed, the TestStat value will be NA (and the equivalent P_Val = 0)
  #       Pval = Sum of all Pval_<Test> / Number of valid Tests, i.e. Mean(Pval, na.rm=TRUE)
  
  testScenario <- paste0(kurtPers,dataRange)
  displayProbabilities(probDist,refTestValues,testScenario,ptitle)
  
  
  
} # End: runLM123Tests

} #End: UtilFn
#==========================================================#

#----------------------------------------------------------#
#               Anna's Tricky Functions 
#----------------------------------------------------------#
zAnnaFn <- TRUE
if (zAnnaFn) {

myQ.EC <- function(N)
{
  # compute matrix of eigenvectors (columns of Q)
  # N = number of series
  Q <- matrix(0,N,N)
  for (i in 1:N)
  {
    if (i==1) Q[,i] <- (1/sqrt(N))*matrix(1,N,1)
    else
    {
      tmp <- 1/(sqrt((N-i+2)*(N-i+1)))
      for (j in seq((i-1),N))
      {
        if (j==(i-1)) Q[j,i] <- tmp*(N-i+1)
        else Q[j,i] <- tmp*(-1)
      }
    }
  }
  #Return
  Q # NXN
}  #End: myQ.EC(N)

myL.EC <- function(N,rho){
  L <- rep((1-rho),N)
  L[1] <- L[1]+rho*N
  #Return vector:
  L # Nx1
}

myVecl <- function(M){
  idx = 0
  N <- ncol(M)
  vM <- matrix(0,N*(N-1)/2,1)
  for (idxCol in seq(1,N-1)){
    for (idxRow in seq(idxCol+1,N)){
      idx <- idx+1
      vM[idx,1] <- M[idxRow,idxCol]
    }
  }
  #Return:
  vM
}

myUnVecl <- function(vM){
  k <- length(vM)
  N <- (1+sqrt(1+8*k))/2
  M <- matrix(0,N,N)
  M[lower.tri(M,FALSE)] <- vM
  M <- M + t(M) + diag(N)
  #Return:
  M
}

myFilter <- function(mX,vB){
  # mX -- T x s matrix
  # vB -- 1 x s vector of coefficients
  # does AR type filtering with lag 1 only
  # output -- mY -- T x s matrix
  # mY[1,] = 0...0
  # mY[t,s] = mX[t,s]+vB[s]*mY[t-1,s]
  s <- length(vB)
  Tobs <- NROW(mX)
  # GLEN: add error check, NCOL(mX)=s
  mY <- matrix(mX[1,],1,s)
  for (t in 2:Tobs){
    mY <- rbind(mY,mX[t,]+vB*mY[(t-1),])
  }
  # Return:
  mY
}

myVecd <- function(M){
  # vectorises main diagonal and off diagonal
  N <- NCOL(M)
  if (N==3){
    vM <- matrix(0,nrow=N*(N+1)/2,ncol=1)
    vM[1:3] <- diag(M)^2
    vM[4] <- sqrt(2)*M[3,2]
    vM[5] <- sqrt(2)*M[3,1]
    vM[6] <- sqrt(2)*M[2,1]
    return(vM)
  }
  else return(NULL)
  
}

} #End: AnnaFn
#==========================================================#

#----------------------------------------------------------#
#    TVGJR - Derivative Functions & 'g' Calculators
#----------------------------------------------------------#
zDerivFn <- TRUE
if (zDerivFn) {

myG <- function(speed,loc1,loc2,st,shape,speedoption){
  # input:	pars = speed,loc1,<loc2>
  #   			st   = transition variable
  #         shape =
  #         speedoption = 1: gamma, 2: gamma/sd(st), 3: exp(eta)
  # output:	Tx1 matrix
  
  st_c <- switch(shape,st-loc1,(st-loc1)*(st-loc2),(st-loc1)*(st-loc1))
  G_inv <- switch(speedoption, 1+exp(-speed*st_c), 1+exp(-speed/sd(st)*st_c), 1+exp(-exp(speed)*st_c))
  return(1/G_inv)
  
}

calculate_g.old <- function(tv) {
  ###  Calculate the Conditional Variances of a TV model
  ###  Input: TV object
  ###  Output: column-based vector 'g', of Length = number of observations

  g <- rep(tv$delta0,length(tv$st))

  if (!is.null(tv$shape)) {
    shape <- tv$shape
    pars <- tv$pars
    speedoption <- tv$speedoption
    st <- tv$st
    for (i in seq_along(shape)) {
      delta_pos <- 4*(i-1)+1
      speed_pos <- 4*(i-1)+2
      loc1_pos <- 4*(i-1)+3
      loc2_pos <- 4*(i-1)+4
      
      delta <- pars[delta_pos]
      speed <- pars[speed_pos]
      loc1 <- pars[loc1_pos]
      loc2 <- pars[loc2_pos]
      
      #Only pass the Speed & Location pars to myG()
      Gi <- myG(speed,loc1,loc2,st,shape[i],speedoption) # Tx1 each
      g <- g + delta*Gi
    }
  }  #End: If (!is.null(shape))
  
  if (!is.null(tv$linpars)) {
  # Note: Deliberately not using a for..loop below for performance reasons
    if (!is.na(tv$linpars[1])) g <- g + tv$linpars[1]
    if (!is.na(tv$linpars[2])) g <- g + tv$linpars[2]*tv$st
    if (!is.na(tv$linpars[3])) g <- g + tv$linpars[3]*tv$st*tv$st
    if (!is.na(tv$linpars[4])) g <- g + tv$linpars[4]*tv$st*tv$st*tv$st
  }
  
  #Return:
  g
}

dg.dt.old <- function(tv){
  # derivative of g(t) w.r.t. TV-parameters (d0,d1,eta,c1,<c2>,...,<lin1>,<lin2>,...)
  # output: T x (#of TV pars (incl. delta0) in Ho) matrix
  
  pars <- tv$pars
  shape <- tv$shape
  speedoption <- tv$speedoption
  st <- tv$st
  Tobs <- length(st)
  
  # if linparsMode == 2, then we need the derivative of TV$pars & linpars
  rtn <- matrix(nrow=Tobs,ncol=length(pars)+1)
  rtn[,1] <- 1
  
  if (!is.null(shape)) {
    for (i in seq_along(shape))
    {
      delta_pos <- 4*(i-1)+1
      speed_pos <- 4*(i-1)+2
      loc1_pos <- 4*(i-1)+3
      loc2_pos <- 4*(i-1)+4
      
      delta <- pars[delta_pos]
      speed <- pars[speed_pos]
      loc1 <- pars[loc1_pos]
      loc2 <- pars[loc2_pos]
      
      # Define: Speed Transform = 'speed_transf'
      speed_transf <- switch(speedoption,speed,speed/sd(st),exp(speed))
      # Define: Transition Variable-Location, st-c = 'st_c'
      st_c <- switch(shape[i],st-loc1,(st-loc1)*(st-loc2),(st-loc1)*(st-loc1))
      
      Gi <- myG(speed,loc1,loc2,st,shape[i],speedoption) # Tx1 each
      
      rtn[,(delta_pos+1)] <- Gi
      rtn[,(speed_pos+1)] <- delta*st_c*speed_transf*Gi*(1-Gi)
      
      if (shape[i]==1) {
        rtn[,(loc1_pos+1)] <- -delta*speed_transf*Gi*(1-Gi)
      } else if (shape[i]==2) {
        rtn[,(loc1_pos+1)] <- -delta*speed_transf*(st-loc2)*Gi*(1-Gi)
        rtn[,(loc2_pos+1)] <- -delta*speed_transf*(st-loc1)*Gi*(1-Gi)
      } else if (shape[i]==3) {
        rtn[,(loc1_pos+1)] <- -delta*speed_transf*2*(st-loc1)*Gi*(1-Gi)
      }
      
    } # End: for (i in 1:length(shape))
    
  } # End: If (!is.null(shape))
  
  # Now remove the NA columns (+ 1 because rtn includes derivative of delta0):
  if (!is.null(pars) && any(is.na(pars))){
    NA_idx <- which(is.na(pars)) + 1
    #Return:
    rtn <- rtn[,-NA_idx]  # T x (#of TV pars incl. delta0, but excluding any 'NA' parameters)
  } 

  # Now calculate the derivatives of the linearised parameters (if any):
  if (!is.null(tv$linparsMode)) {
    # Initialise the rtn_lin matrix with data = NA
    rtn_lin <- matrix(data=NA,nrow=Tobs,ncol=4)
    if (!is.na(tv$linpars[1])) rtn_lin[,1] <- 1 
    if (!is.na(tv$linpars[2])) rtn_lin[,2] <- tv$st
    if (!is.na(tv$linpars[3])) rtn_lin[,3] <- tv$st*tv$st
    if (!is.na(tv$linpars[4])) rtn_lin[,4] <- tv$st*tv$st*tv$st
    # Now remove all columns that still = NA:
    NA_idx <- which(is.na(rtn_lin))
    rtn_lin <- rtn_lin[,-NA_idx]
    rtn <- cbind(rtn,rtn_lin)
  }

  #Return:
  rtn  # T x (#of TV pars incl. delta0)
  
}


calculate_g <- function(tv) {
  ###  Calculate the Conditional Variances of a TV model
  ###  Input: TV object
  ###  Output: column-based vector 'g', of Length = number of observations
  
  g <- rep(tv$delta0,tv$Tobs)

  if (!is.null(tv$pars)) {
    myPars <- matrix(data=tv$pars, nrow=length(tv$shape), ncol=4, byrow=TRUE, dimnames=list(NULL,c("delta","speed","loc1","loc2")))
    # initialise some variables
    stdev_st <- sd(tv$st)
    st_c <- 0
    speed_transf <- 0
    Gi <- 0
    for (i in seq_along(tv$shape)) {
        st_c <- switch(tv$shape[i],tv$st-myPars[i,"loc1"],(tv$st-myPars[i,"loc1"])*(tv$st-myPars[i,"loc2"]),(tv$st-myPars[i,"loc1"])*(tv$st-myPars[i,"loc1"]))  
        switch (tv$speedoption,
          Gi <- 1/(1+exp(-myPars[i,"speed"]*st_c)),
          Gi <- 1/(1+exp(-myPars[i,"speed"]/stdev_st*st_c)),
          Gi <- 1/(1+exp(-exp(myPars[i,"speed"])*st_c))
        ) #end: switch()
        g <- g + myPars[i,"delta"]*Gi
    }
  }
  
  if (!is.null(tv$linpars)) {
    # Note: Deliberately not using a for..loop below for performance reasons
    if (!is.na(tv$linpars[1])) g <- g + tv$linpars[1]
    if (!is.na(tv$linpars[2])) g <- g + tv$linpars[2]*tv$st
    if (!is.na(tv$linpars[3])) g <- g + tv$linpars[3]*tv$st*tv$st
    if (!is.na(tv$linpars[4])) g <- g + tv$linpars[4]*tv$st*tv$st*tv$st
  }
  
  #Return:
  g
}

dg.dt <- function(tv){
  # derivative of g(t) w.r.t. TV-parameters (d0,d1,eta,c1,<c2>,...,<lin1>,<lin2>,...)
  # output: T x (#of TV pars (incl. delta0) in Ho) matrix
  
  # We need the derivative of TV$pars - We will process TV$linpars(if any) later
  rtn <- matrix(nrow=tv$Tobs,ncol=length(tv$pars)+1)  # Note: pars may still contain NA values here
  rtn[,1] <- 1
  
  if (!is.null(tv$pars)) {
    
    # get the tvpars as a matrix for convenience
    myPars <- matrix(data=tv$pars, nrow=length(tv$shape), ncol=4, byrow=TRUE, dimnames=list(NULL,c("delta","speed","loc1","loc2")))
    # initialise some variables
    stdev_st <- sd(tv$st)
    st_c <- 0
    speed_transf <- 0
    Gi <- 0
    
      for (i in seq_along(tv$shape)) {
    
        # Define: Transition Variable-Location, st-c = 'st_c'
        st_c <- switch(tv$shape[i],tv$st-myPars[i,"loc1"],(tv$st-myPars[i,"loc1"])*(tv$st-myPars[i,"loc2"]),(tv$st-myPars[i,"loc1"])*(tv$st-myPars[i,"loc1"]))
        
        switch (tv$speedoption,
          {speed_transf <- myPars[i,"speed"]
           Gi <- 1/(1+exp(-myPars[i,"speed"]*st_c))},
          {speed_transf <- myPars[i,"speed"]/stdev_st
           Gi <- 1/(1+exp(-myPars[i,"speed"]/stdev_st*st_c))},
          {speed_transf <- exp(myPars[i,"speed"])
           Gi <- 1/(1+exp(-exp(myPars[i,"speed"])*st_c))}
        ) #end: switch()
  
        deriv_const <- myPars[i,"delta"]*speed_transf*Gi*(1-Gi)
        
        startPos <- 4*(i-1)+2
        rtn[,startPos] <- Gi
        rtn[,startPos+1] <- deriv_const*st_c
        
        switch (tv$shape[i],
          rtn[,startPos+2] <- -deriv_const,
          {rtn[,startPos+2] <- -deriv_const*(tv$st-myPars[i,"loc2"])
           rtn[,startPos+3] <- -deriv_const*(tv$st-myPars[i,"loc1"])},
          rtn[,startPos+2] <- -deriv_const*2*(tv$st-myPars[i,"loc1"])
        )
      } # End: for (i in seq_along(shape))

      # Now remove any NA columns:
      if (any(is.na(tv$pars))){
        NA_idx <- which(is.na(tv$pars)) + 1  #(+ 1 because rtn includes derivative of delta0)
        #Return:
        rtn <- rtn[,-NA_idx]  
      }
  }
  
  # Now calculate the derivatives of the linearised parameters (if any):
  if (!is.null(tv$linpars)) {
    # Initialise the rtn_lin matrix with data = NA
    rtn_lin <- matrix(data=NA,nrow=tv$Tobs,ncol=4)
    if (!is.na(tv$linpars[1])) rtn_lin[,1] <- 1 # I think we need to drop this!
    if (!is.na(tv$linpars[2])) rtn_lin[,2] <- tv$st
    if (!is.na(tv$linpars[3])) rtn_lin[,3] <- tv$st*tv$st
    if (!is.na(tv$linpars[4])) rtn_lin[,4] <- tv$st*tv$st*tv$st
    # Now remove all columns that still = NA:
    NA_idx <- which(is.na(tv$linpars))
    rtn_lin <- rtn_lin[,-NA_idx]
    rtn <- cbind(rtn,rtn_lin)
  }
  
  #Return:
  dimnames(rtn) <- NULL
  rtn  # T x (#of TV pars incl. delta0)
  
}  # End dg.dt_new

dg.dt2.H0 <- function(st){
  # purpose:	recursively calculates the partial derivative dg/dTheta2
  #           in the linearized model,
  #			Theta2=parameters from the linearized TV component, d0,d1,d2,d3
  # output:	Tx4 matrix, each row = dg/dd0, dg/dd1, dg/dd2, dg/dd3
  rtn <- matrix(0,length(st),4)
  rtn[,1] <- 1
  rtn[,2] <- st
  rtn[,3] <- st*st
  rtn[,4] <- st*st*st
  
  #Return:
  rtn # Tx4
}

} #End: DerivFn
#==========================================================#

#----------------------------------------------------------#
#    TVGJR - TEST Functions 
#----------------------------------------------------------#
zTestFn <- TRUE
if (zTestFn) {

myTest.TV.noGARCH.TR2.old <- function(e,tv,TestOrder="H03"){
  # TR2-form of the LM test
  # No GARCH in the H0 model, no GARCH in the H1 model, i.e. no GARCH estimated, potential GARCH is ignored
  # H0: gt=delta0
  # H1: gt=delta0 + TV$pars
  # e    = raw data
  # shape = shape of TV under H0
  # speedoption = 1:gamma,2:gamma/sd(st),3:exp(eta)
  # st = transition variable in TV and used in test
  
  #H03 Test => TV$linpars must be c(NA,1,1,NA)
  #H02 Test => TV$linpars must be c(NA,1,NA,NA)
  #H01 Test => TV$linpars must be NULL
  
  TestStat <- NA
  
  Tobs <- length(tv$st)
  g <- calculate_g(tv)
  
  # 1. Get squared standardised residuals minus one, compute SSR0
  psi2_1 <- (e*e/g)-1
  SSR0 <- sum(psi2_1*psi2_1)    # Scalar
  
  # 2. Regress psi2_1 on 1/gt*(dgdt and dgdt2.H0) (here t=TV component in H0, t2=linearised (tested) TV component)
  #dgdt = Tx1 or Tx4 or Tx7 or... NCOL(dgdt) increases with the order of TV function. 
  dgdt <- dg.dt(tv)
  dgdt2 <- dg.dt2.H0(tv$st) # Tx4 but first column is ones, same as in dgdt, so remove that
  
  #dgdt2 is always a Tx1 matrix - need "drop=F" to prevent R reducing result to a Vector type
  if(TestOrder=="H03") dgdt2 <- dgdt2[,4,drop=FALSE]           
  else if(TestOrder=="H02") dgdt2 <- dgdt2[,3,drop=FALSE]
  else if(TestOrder=="H01") dgdt2 <- dgdt2[,2,drop=FALSE]
  else if(TestOrder=="H0") dgdt2 <- dgdt2[,2:4]
  
  X <- cbind(dgdt,dgdt2)
  X <- X/g
  
  # 3. Compute test statistic
  #SSR1 <- sum((psi2_1-X%*%solve(t(X)%*%X)%*%t(X)%*%psi2_1)^2)
  #Test$LMstat = T*(SSR0-SSR1)/SSR0
  Xinv <- NULL
  try(Xinv <- solve(t(X)%*%X))
  if(is.null(Xinv)) {
    warning("\n myTest.TV.noGARCH.TR2() threw an error doing solve() - returning NA")
    rm(psi2_1,g,X,dgdt,dgdt2)
    return(NA)
  }
  #psi2_1 <- as.matrix(psi2_1)
  SSR1 <- sum((psi2_1-X%*%Xinv%*%t(X)%*%psi2_1)^2)
  TestStat <- Tobs*(SSR0-SSR1)/SSR0

  #Tidy up & release memory before returning:
  rm(psi2_1,g,X,Xinv,dgdt,dgdt2)
 
  #Return:
  TestStat
  
}  #End: myTest.TV.noGARCH.TR2()

myTest.TV.noGARCH.TR2 <- function(e,tv,TestOrder="H0"){
  # TR2-form of the LM test
  # No GARCH in the H0 model, no GARCH in the H1 model, i.e. no GARCH estimated, potential GARCH is ignored
  # H0: gt=delta0
  # H1: gt=delta0 + TV$pars
  # e    = raw data
  # shape = shape of TV under H0
  # speedoption = 1:gamma,2:gamma/sd(st),3:exp(eta)
  # st = transition variable in TV and used in test
  
  #H03 Test => TV$linpars must be c(NA,1,1,NA)
  #H02 Test => TV$linpars must be c(NA,1,NA,NA)
  #H01 Test => TV$linpars must be NULL
  
  TestStat <- NA
  
  g <- calculate_g(tv)
  
  # 1. Get squared standardised residuals minus one, compute SSR0
  psi2_1 <- e*e/g-1
  SSR0 <- sum(psi2_1*psi2_1)    # Scalar
  
  # 2. Regress psi2_1 on 1/gt*(dgdt and dgdt2.H0) (here t=TV component in H0, t2=linearised (tested) TV component)
  #dgdt = Tx1 or Tx4 or Tx7 or... NCOL(dgdt) increases with the order of TV function. 
  dgdt <- dg.dt(tv)
  dgdt2 <- dg.dt2.H0(tv$st) # Tx4 but first column is ones, same as in dgdt, so remove that
  
  #dgdt2 is always a Tx1 matrix - need "drop=F" to prevent R reducing result to a Vector type
  if(TestOrder=="H0") dgdt2 <- dgdt2 <- dgdt2[,2:4]           
  else if(TestOrder=="H01") dgdt2 <- dgdt2[,2,drop=FALSE]
  else if(TestOrder=="H02") dgdt2 <- dgdt2[,3,drop=FALSE]
  else if(TestOrder=="H03") dgdt2 <- dgdt2[,4,drop=FALSE]
  
  X <- cbind(dgdt,dgdt2)/g
  
  # 3. Compute test statistic
  Xinv <- NULL
  try(Xinv <- solve(crossprod(X,X)))
  if(is.null(Xinv)) {
    warning("\n myTest.TV.noGARCH.TR2() threw an error doing solve() - returning NA")
    rm(psi2_1,g,X,dgdt,dgdt2)
    return(NA)
  }
  SSR1 <- sum((psi2_1-X%*%Xinv%*%t(X)%*%psi2_1)^2)
  # Return:
  TestStat <- tv$Tobs*(SSR0-SSR1)/SSR0
  
}  #End: myTest.TV.noGARCH.TR2_new()

myTest.TV.noGARCH.robust <- function(e,tv,TestOrder="H0"){
  # Robust form of the LM test
  # No GARCH in the H0 model, no GARCH in the H1 model, i.e. no GARCH estimated, potential GARCH is ignored
  # Ho: gt=delta0
  # H1: gt=delta0 + TV-terms
  # pars = estimate of delta0
  # e    = raw data
  # tv$shape = shape of TV under H0
  # tv$speedoption = 1:gamma,2:gamma/sd(st),3:exp(eta)
  # tv$st = transition variable in TV and used in test
  
  #H03 Test => TV$linpars <- c(NA,1,1,NA)
  #H02 Test => TV$linpars <- c(NA,1,NA,NA)
  #H01 Test => TV$linpars <- NULL
  
  TestStat <- NA
  
  # 1. Get squared standardised residuals minus one
  # 2. Regress r2 = 1/gt*dgdt2.Ho on r1 = 1/gt*dgdt, and compute residual vectors wt

  dgdt <- dg.dt(tv)
  dgdt2 <- dg.dt2.H0(tv$st) # Tx4 but first column is ones, same as in dgdt, so remove that

  #dgdt2 - need "drop=F" to prevent R reducing result to a Vector type when only 1 column returned
  if(TestOrder=="H03") dgdt2 <- dgdt2[,4,drop=FALSE]           
  else if(TestOrder=="H02") dgdt2 <- dgdt2[,3,drop=FALSE]
  else if(TestOrder=="H01") dgdt2 <- dgdt2[,2,drop=FALSE]
  else if(TestOrder=="H0") dgdt2 <- dgdt2[,2:4]

  g <- calculate_g(tv)
  X <- dgdt/g 
  
  # Note: Solve can fail, which means we can't calculate a TestStat:
  Xinv <- NULL
  #try(Xinv <- solve(t(X)%*%X))
  try(Xinv <- solve(crossprod(X,X)))
  if(is.null(Xinv)) {
    warning("\n myTest.TV.noGARCH.robust() threw an error doing solve() - returning NA")
    rm(g,X,dgdt,dgdt2)
    return(NA)
  }
  XXXX <- X%*%Xinv%*%t(X)
  Y <- as.matrix(dgdt2/g)
  w <- as.matrix(Y-XXXX%*%Y)

  #3. Regress 1 on (psi2-1)*w, and compute SSR
  psi2_1 <- as.vector(e*e/g-1)
  X <- psi2_1*w  #psi2_1 must be a vector for this!!
  
  #4. Compute test statistic:

  Xinv <- NULL
  #try(Xinv <- solve(t(X)%*%X))
  try(Xinv <- solve(crossprod(X,X)))
  if(is.null(Xinv)) {
    warning("\n myTest.TV.noGARCH.Robust() threw an error doing solve() - returning NA")
    rm(psi2_1,g,w,Y,X,dgdt,dgdt2)
    return(NA)
  }

  TestStat <- tv$Tobs-sum(diag(tv$Tobs)-(X%*%Xinv%*%t(X)))
  
  #Tidy up & release memory before returning:
  rm(tv,psi2_1,g,w,Y,X,Xinv,dgdt,dgdt2)
  
  #Return:
  TestStat
  
}  #End: myTest.TV.noGARCH.robust()



#TODO: Refactor or Do NOT use!
myTest.TV.noGARCH.F <- function(tv_r,tv_u,e){
  # F-form of the LM test
  # No GARCH in the H0 model, no GARCH in the H1 model, i.e. no GARCH estimated, potential GARCH is ignored
  # e    = raw data
  # shape = shape of TV under H0
  # speedoption = 1:gamma,2:gamma/sd(st),3:exp(eta)
  # st = transition variable in TV and used in test
  # tv_r = Restricted function, tv_u = Unrestricted function
  # Note: All the TV$... non-linear elements are common to both functions, only the linpars[] are different
  
  shape <- tv_u$shape
  speedoption <- tv_u$speedoption
  st <- tv_u$st
  sp_pos <- tv_u$speed_pos
  Tobs <- length(st)
  
  # Create the restricted & unrestricted TV params to be estimated
  tv_r$pars <- c(tv_r$pars,tv_r$linpars)
  tv_u$pars <- c(tv_u$pars,tv_u$linpars)
  
  #Estimate both & extract the pars vectors:
  myControl = list(fnscale = -1, ndeps = rep.int(1e-4, length(tv_r$pars)), maxit = 500, REPORT = 100)
  tv_r <- EstimateTV(tv_r,e,myControl,TRUE)
  cat("\ntv_r$pars=",tv_r$pars)
  cat("\ntv_r$value=",tv_r$value)
  
  myControl = list(fnscale = -1, ndeps = rep.int(1e-4, length(tv_u$pars)), maxit = 500, REPORT = 100)
  tv_u <- EstimateTV(tv_u,e,myControl,TRUE)
  cat("\ntv_u$pars=",tv_u$pars)
  cat("\ntv_u$value=",tv_u$value)
  
  tv_r$linpars <- tail(tv_r$pars,length(tv_r$linpars))
  tv_u$linpars <- tail(tv_u$pars,length(tv_u$linpars))
  
  g_r <- rep(tv_r$pars[1],Tobs)
  g_u <- rep(tv_u$pars[1],Tobs)
  
  
  if (!is.null(shape))
  {
    for (i in seq_along(shape))
    {
      startPos = sp_pos[i]
      endPos = if(shape[i]==2) startPos+2 else startPos+1
      #Only pass the Speed & Location pars to myG()
      Gi_r <- myG(tv_r$pars[startPos:endPos],st,shape[i],speedoption) # Tx1 each
      g_r <- g_r + tv_r$pars[startPos-1]*Gi_r
      Gi_u <- myG(tv_u$pars[startPos:endPos],st,shape[i],speedoption) # Tx1 each
      g_u <- g_u + tv_u$pars[startPos-1]*Gi_u
    }
  }  #End: if (!is.null(shape))
  
  for (n in seq_along(tv_r$linpars)) g_r <- g_r + tv_r$linpars[n]*(st^n)
  for (n in seq_along(tv_u$linpars)) g_u <- g_u + tv_u$linpars[n]*(st^n)
  
  SSRr <- sum(((e^2/g_r)-1)^2)
  SSRu <- sum(((e^2/g_u)-1)^2)
  
  cat("\nSSRr, SSRu",SSRr,SSRu)
  
  Test <- list()
  Test$Stat <- (SSRr-SSRu)/(SSRu/(Tobs-length(tv_u$pars)))
  Test$Pval <- pf(Test$Stat, df1 = 1, df2 = Tobs-length(tv_u$pars), lower.tail = FALSE)
  
  #Return:
  Test
}  #End: myTest.TV.noGARCH.F

#TODO: Refactor or Do NOT use!
myTest.TV.noGARCH.LR <- function(tv_r,tv_u,e){
  # Liklihood Ratio Test
  # No GARCH in the H0 model, no GARCH in the H1 model, i.e. no GARCH estimated, potential GARCH is ignored
  # e    = raw data
  # shape = shape of TV under H0
  # speedoption = 1:gamma,2:gamma/sd(st),3:exp(eta)
  # st = transition variable in TV and used in test
  # tv_r = Restricted function, tv_u = Unrestricted function
  # Note: All the TV$... non-linear elements are common to both functions, only the linpars[] are different
  
  
  shape <- tv_u$shape
  speedoption <- tv_u$speedoption
  st <- tv_u$st
  sp_pos <- tv_u$speed_pos
  Tobs <- length(st)
  
  # Create the restricted & unrestricted TV params to be estimated
  tv_r$pars <- c(tv_r$pars,tv_r$linpars)
  tv_u$pars <- c(tv_u$pars,tv_u$linpars)
  
  #Estimate both & extract the pars vectors:
  myControl = list(fnscale = -1, ndeps = rep.int(1e-4, length(tv_r$pars)), maxit = 500, REPORT = 100)
  tv_r <- EstimateTV(tv_r,e,myControl,TRUE)
  cat("\ntv_r$pars=",tv_r$pars)
  cat("\ntv_r$value=",tv_r$value)
  
  myControl = list(fnscale = -1, ndeps = rep.int(1e-4, length(tv_u$pars)), maxit = 500, REPORT = 100)
  tv_u <- EstimateTV(tv_u,e,myControl,TRUE)
  cat("\ntv_u$pars=",tv_u$pars)
  cat("\ntv_u$value=",tv_u$value)
  
  Test <- list()
  Test$Stat <- 2*(tv_u$value - tv_r$value)
  Test$Pval <- pchisq(Test$Stat, 1, lower.tail = FALSE)
  
  #Return:
  Test
}  #End: myTest.TV.noGARCH.LR()

myTest.CCCvSTCC.LM.old <- function(e,H0,H1,testorder=1) {
  
  CCC <- H0
  STCC <- H1
  
  N <- NCOL(e)  # Number of series
  Tobs <- NROW(e)
  P <- CCC$P
  st <- STCC$st
  testOrder <- testorder    # Can be 1 or 2, indicating single or double transition

  # 1. Create M_cor:
  if(T){
    Pinv <- try(solve(P))
    if (!is.matrix(Pinv)) stop("LM Test: Can't invert the correlation matrix P")
    # Create an Identity Matrix:
    I <- diag(N)
    
    # Construct the K matrix
    K <- NULL
    for (i in 1:N) {
      # block rows
      Krow <- NULL
      for (j in 1:N) {
        # block columns
        block <- matrix(0,N,N)
        block[j,i] <- 1
        Krow <- cbind(Krow,block)
      }
      K <- rbind(K,Krow)
    }
  }
  M_cor <- Pinv %x% Pinv + (Pinv %x% I) %*% K %*% (Pinv %x% I)
  
  # 2. Create the Information Matrix for the Correlation, and extract the SE corner of its inverse:
  if(T){
    # Construct the U matrix: Dimensions = N^2 x N(N-1)/2
    U <- NULL
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        block <- matrix(0,N,N)
        block[i,j] <- block[j,i] <- 1
        Ucol <- as.vector(block)
        U <- cbind(U,Ucol)
      }
    }
    
    IM_cor <- matrix(0,(testOrder+1)*N*(N-1)/2,(testOrder+1)*N*(N-1)/2)
    for (t in 1:Tobs) {
      # Calculate v_rho_t:
      if (testOrder==1) v_rho_t <- c(1,st[t]) else if (testOrder==2) v_rho_t <- c(1,st[t],st[t]^2) else stop("LM Test: testOrder must be 1 or 2!")
      # Calculate x_rho_t:
      x_rho_t <- (-0.5*v_rho_t) %x% t(U)
      IM_cor <- IM_cor + (x_rho_t %*% M_cor %*% t(x_rho_t))
    }
    IM_cor <- IM_cor/Tobs
    # inverse of the information matrix
    IM_cor_inv <- try(solve(IM_cor))
    if (!is.matrix(IM_cor_inv)) stop("LM Test: Can't invert the information matrix IM_cor")
  }
  # Extract testOrder(N(N-1))/2 x testOrder(N(N-1))/2 SE block,  # testOrder=1, N=3, => 3 x 3, testOrder=2, N=3, => 6 x 6
  block_start <- (N*(N-1)/2) + 1
  block_end <- NCOL(IM_cor_inv)
  IM_inv_SE <- IM_cor_inv[(block_start:block_end),(block_start:block_end)]
  
  # 3. Calculate the derivative of the log likelihood w.r.t correlation parameters under test:
  if(T){
    # z_t is the time series of volatility standardised returns or returns if no vol.model
    z_t <- e
    dll_d_rho <- matrix(0,testOrder*N*(N-1)/2, 1)
    for (t in 1:Tobs) {
      if (testOrder==1) v_rho_t <- c(1,st[t]) else if (testOrder==2) v_rho_t <- c(1,st[t],st[t]^2) else stop("LM Test: testOrder must be 1 or 2!")
      dll_d_rho <- dll_d_rho + (-0.5 * (v_rho_t[2:(testOrder+1)] %x% t(U)) %*% (as.vector(Pinv)-(Pinv %x% Pinv) %*% (z_t[t,] %x% z_t[t,])) )
    }
  }
  # dll_dv_rho, Dimension: testorder*N(N-1)/2 x 1
  
  # 4. Return: LM test statistic 
  LM <- (t(dll_d_rho)%*%IM_inv_SE%*%dll_d_rho)/Tobs
  # Return:
  LM
  
}  # End: myTest.CCCvSTCC.LM()

myTest.CCCvSTCC.LM <- function(e,H0,H1,testorder=1) {
  
  nGARCH <- H0$nGARCH
  nTV <- H0$nTV
  CCC <- H0$CCC 
  STCC <- H1
  
  if (is.null(CCC)) stop("LM Test: need CCC as input (H0)")
  if (is.null(STCC$st)) stop("LM Test: need transition variable for the test as input (H1)")
  
  N <- NCOL(e)  # Number of series
  Tobs <- NROW(e)
  P <- CCC$P
  st <- STCC$st
  testOrder <- testorder    # Can be 1 or 2, indicating single or double transition
  I <- diag(nrow = N,ncol = N) # identity matrix
  Pinv <- try(solve(P))
  if (!is.matrix(Pinv)) stop("LM Test: Can't invert the correlation matrix P")
  
  #Construct the K matrix: N^2 x N^2
  if(T){
    K <- NULL
    for (i in 1:N) {
      # block rows
      Krow <- NULL
      for (j in 1:N) {
        # block columns
        block <- matrix(0,N,N)
        block[j,i] <- 1
        Krow <- cbind(Krow,block)
      }
      K <- rbind(K,Krow)
    }
  }
  
  #Construct the U matrix:Dimensions = N^2 x N(N-1)/2
  if(T){
    U <- NULL
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        block <- matrix(0,N,N)
        block[i,j] <- block[j,i] <- 1
        Ucol <- as.vector(block)
        U <- cbind(U,Ucol)
      }
    }
  }
  
  # 2. Create the Information Matrix for the Correlation, and extract the SE corner of its inverse:
  if(T){
    
    # Initialise variables:
    IM_cor <- matrix(0,(testOrder+1)*N*(N-1)/2,(testOrder+1)*N*(N-1)/2)
    g <- matrix(1,Tobs,N) # for now, no TV estimated, so these are 1, T x N
    h <- matrix(1,Tobs,N) # place for estimated GARCH variances, ones at this stage, T x N
    if (!is.null(nGARCH)){
      # Extract the number of Garch pars: (Assume that Garch$type=0 has 1 par only, Garch$type=1 has 3, etc. )
      Num_garch_pars <- 0
      for (n in 1:N) Num_garch_pars <- Num_garch_pars + length(nGARCH[[n]]$pars)
      IM_garch <- matrix(0,Num_garch_pars,Num_garch_pars)
      IM_garch_cor <- matrix(0,Num_garch_pars,(testOrder+1)*N*(N-1)/2) 
    }
    
    # Calculate real values
    if (!is.null(nGARCH)){
      beta <- rep(0,N) # estimated beta coefficients from GARCH equations
      for (n in 1:N){
        h[,n] <- nGARCH[[n]]$condvars # univariate GARCH volatilities, T x N 
        if(nGARCH[[n]]$type==0) beta[n] <- 0 else beta[n] <- nGARCH[[n]]$pars[3]
      }
    }
    if (!is.null(nTV)){
      for (n in 1:N){
        g[,n] <- nTV[[n]]$condvars # univariate TV volatilities, T x N 
      }
    }
    w <- e/sqrt(g) # returns standardised by g(t), T x N
    z <- w/sqrt(h) # returns standardised by GARCH volatilities and TV g(t), T x N

    # partial derivatives of P w.r.t correlation parameters, T x 1 or T x 2
    # these are for alternative parameters, later add column of ones at the front
    if (testOrder==1) v_rho <- cbind(st) else if (testorder==2) v_rho <- cbind(st,st^2) else stop("LM Test: testOrder must be 1 or 2!")
    
    # score for rho
    zKRONz <- NULL
    for (t in 1:Tobs) zKRONz <- cbind(zKRONz, z[t,] %x% z[t,]) # (N^2 x T), each col = "z_t kron z_t"
    
    dlldrho_A <- -0.5*t(U)%*%( as.vector(Pinv)%x%t(rep(1,Tobs))-(Pinv%x%Pinv)%*%(zKRONz)  )%*%v_rho # N(N-1)/2 x testorder
    dlldrho_A <- vec(dlldrho_A) # testorder*N(N-1)/2 x 1, SUM OVER TIME
    v_rho <- cbind(rep(1,Tobs),v_rho) # T x 2 or T x 3, now add column of ones at the front
    
    # partial derivatives for GARCH and "x" matrices for GARCH parts    
    if (!is.null(nGARCH)){
      One_33 <- matrix(1,3,3)
      One_31 <- matrix(1,3,1)
      One_13 <- matrix(1,1,3)
      One_11 <- matrix(1,1,1)
      # partial derivatives of h1,h2,...,hN w.r.t garch_pars
      v_garch <- NULL
      h_scale <- NULL
      beta_scale <- NULL
      for (n in 1:N){
        if(nGARCH[[n]]$type==0) {
            v_garch <- cbind(v_garch,c(0,rep(1,(Tobs-1)))) 
            beta_scale <- c(beta_scale,beta[n])
            h_scale <- cbind(h_scale,h[,n])
          } else {
            v_garch <- cbind(v_garch,c(0,rep(1,(Tobs-1))),c(0,w[1:(Tobs-1),n]^2),c(0,h[1:(Tobs-1),n])) # T x Num_garch_pars, each row = "1~w(i,t-1)^2~h(i,t-1)", i=1,...,N    
            beta_scale <- beta_scale <- c(beta_scale,(beta[n] %x% rep(1,3)))
            h_scale <- cbind(h_scale,h[,n] %x% t(rep(1,3)))
          }
      }
      dhdt <- myFilter(v_garch,beta_scale) # T x Num_garch_pars, each row = "(dh(i,t)/dtheta(i))' ", i=1,...,N
      # matrix containing the x_t's (garch part)
      x_garch <- -0.5*dhdt/(h_scale) # T x Num_garch_pars, each row = "x_it'", i=1,...,N
    }
    
    
    # Information matrix:
    # Blocks involving GARCH only or GARCH and Correlation:
    if (!is.null(nGARCH)){
      # IM_garch (Num_garch_pars x Num_garch_pars), SUM OVER TIME
      I.P.Pinv <- I + P*Pinv
      I.P.Pinv_scale <- NULL
      scaleFactor <- NULL
      for (i in 1:N) {
        # i = row index
        I.P.Pinv_scale_row <- NULL
        for (j in 1:N) {
          # j = col index
          if(nGarch[[i]]$type==1 && nGarch[[j]]$type==1) scaleFactor <- One_33 else
            if(nGarch[[i]]$type==1 && nGarch[[j]]$type==0) scaleFactor <- One_31 else
              if(nGarch[[i]]$type==0 && nGarch[[j]]$type==1) scaleFactor <- One_13 else
                if(nGarch[[i]]$type==0 && nGarch[[j]]$type==0) scaleFactor <- One_11 else stop("Something went wrong!")
          I.P.Pinv_scale_row <- cbind(I.P.Pinv_scale_row, (I.P.Pinv[i,j] %x% scaleFactor))
        }
        I.P.Pinv_scale <- rbind(I.P.Pinv_scale,I.P.Pinv_scale_row)
      }
      
      IM_garch <- ((t(x_garch)%*%x_garch) * I.P.Pinv_scale) / Tobs
    
      # IM_garch_cor, 3N x (testorder+1)*N(N-1)/2
      Mhelp1 <- NULL
      for (n in 1:N){
        Mhelp1 <- rbind(Mhelp1, (Pinv[n,]%x%I[n,]+I[n,]%x%Pinv[n,])) # N x N^2
      }
      Mhelp1 <- -0.5*(Mhelp1%*%U)
      Mhelp1_scale <- NULL  # (Num_Garch_Pars x N(N-1)/2)
      for (n in 1:N) {
        if (nGarch[[n]]$type==0) Mhelp1_scale <- rbind(Mhelp1_scale,Mhelp1[n,]) else
          if (nGarch[[n]]$type==1) Mhelp1_scale <- rbind(Mhelp1_scale, (Mhelp1[n,] %x% One_31))
      }

      Mhelp2 <- t(t(v_rho)%*%x_garch)/Tobs # Num_garch_pars x 2 (or x3 if TestOrder=2), SUM OVER TIME
      IM_garch_cor <- NULL
      for (i in 1:NCOL(v_rho)){
        IM_garch_cor <- cbind(IM_garch_cor,((Mhelp2[,i]%x%t(rep(1,(N*(N-1))/2))) * Mhelp1_scale))  # (Num_Garch_Pars x (testOrder+1)*N(N-1)/2),  SUM OVER TIME
      }
    }
    
    # Blocks involving TV only or TV and Correlation:
    if (!is.null(nTV)){
      # to be written...
      IM_tv <- NULL
      IM_tv_cor <- NULL
    }
    # Blocks involving TV and GARCH:
    if (!is.null(nGARCH) && !is.null(nTV)){
      # to be written...
      IM_garch_tv <- NULL
    }
    # Block involving Correlation only:
    # IM_cor (testorder+1)*N(N-1)/2, SUM OVER TIME
    Mhelp3 <- t(U)%*%(Pinv%x%Pinv+(Pinv%x%I)%*%K%*%(Pinv%x%I))%*%U # N(N-1)/2 x N(N-1)/2
    Mhelp4 <- t(v_rho)%*%v_rho/Tobs # 2x2 or 3x3, SUM OVER TIME
    IM_cor <- 0.25*(Mhelp4%x%Mhelp3) # (testorder+1)*N(N-1)/2 x (testorder+1)*N(N-1)/2
    
    # IM
    #IM <- rbind(cbind(IM_garch,IM_garch_cor),cbind(t(IM_garch_cor),IM_cor)) # the whole IM, not really needed...
    # block corresponding to correlations of the inverse of the IM matrix:
    if (!is.null(nGARCH) && !is.null(nTV)){
      # to be written...
      
    }
    if (is.null(nGARCH) && !is.null(nTV)){
      # to be written
      IM_tv_inv <- try(solve(IM_tv))
      if (!is.matrix(IM_tv_inv)) stop("LM Test: Can't invert the information matrix IM_tv")
      IM_inv <- solve(IM_cor-t(IM_tv_cor)%*%IM_tv_inv%*%IM_tv_cor)
    }
    if (!is.null(nGARCH) && is.null(nTV)){
      IM_garch_inv <- try(solve(IM_garch))
      if (!is.matrix(IM_garch_inv)) stop("LM Test: Can't invert the information matrix IM_garch")
      IM_inv <- solve(IM_cor-t(IM_garch_cor)%*%IM_garch_inv%*%IM_garch_cor)
    }
    if (is.null(nGARCH) && is.null(nTV)){
      IM_inv <- solve(IM_cor)
    }
    if (!is.matrix(IM_inv)) stop("LM Test: Can't invert the information matrix IM_inv")
    # block corresponding to the corr.parameters that are set to zero under null
    SE_dim <- testOrder*(N*(N-1))/2
    block_start <- NCOL(IM_inv)-SE_dim+1
    block_end <- NCOL(IM_inv)
    IM_inv_SE <- IM_inv[(block_start:block_end),(block_start:block_end)]
    # LM statistic
    LM <- (1/Tobs)*t(dlldrho_A)%*%IM_inv_SE%*%dlldrho_A
    # Return:
    LM
  }
  

}  # End: myTest.CCCvSTCC.LM.new2()

myTest.CCCvSTCC.LM_c <- cmpfun(myTest.CCCvSTCC.LM)

myTest.STEC.LM <- function(e,tv,garch,stec){
  # test of H0:CEC vs H1:STEC
  # NOTE: myTest.STEC.LM only written for no-TV no-GARCH CEC-only case
  
  TestStat <- NA
  
  N <- NCOL(e)
  Tobs <- NROW(e)
  
  STEC <- stec
  st <- STEC$st
  Q <- myQ.EC(N)
  rho <- STEC$pars[1]
  psi0 <- 1+(N-1)*rho
  
  w <- e%*%Q         # TxN,  wt = Q'zt = Q' Sinv Dinv et -- -- -- for now S=D=I !!!!!!! CHANGE FOR TV and GARCH
  s_psi0 <- s_psi1 <- s_psi2 <- vector("numeric",length=Tobs)
  
  # s_psi0[t] <- -0.5*(1/psi0-(N-1)/(N-psi0))-0.5*((N-1)/(N-psi0)^2*w[t,]%*%t(w[t,])-w[t,1]^2*(1/psi0^2+(N-1)/(N-psi0)^2))  #Tx1
  
  a <- -0.5*(1/psi0-(N-1)/(N-psi0))
  b <- (N-1)/(N-psi0)^2
  c <- (1/psi0^2+(N-1)/(N-psi0)^2)
  
  #TODO: Can we optimise the loop below?
  for (t in seq(1,Tobs)){
    ww <- (w[t,,drop=FALSE])%*%t(w[t,,drop=FALSE])
    s_psi0[t] <- a -0.5*(b * ww - w[t,1]^2 * c)  #Tx1
    s_psi1[t] <- a*st[t]-0.5*(b * ww - w[t,1]^2 * c)*st[t]
    s_psi2[t] <- a*st[t]^2 -0.5*(b * ww - w[t,1]^2 * c)*st[t]^2
  }
  
  Jpsi0psi0 <- sum(s_psi0^2)
  Jpsi0psi1 <- sum(s_psi0*s_psi1)
  
  Jpsi0psi2 <- sum(s_psi0*s_psi2)
  
  Jpsi1psi1 <- sum(s_psi1^2)
  Jpsi1psi2 <- sum(s_psi1*s_psi2)
  
  Jpsi2psi2 <- sum(s_psi2^2)
  
  A <- matrix(Jpsi0psi0, nrow=1, ncol=1)
  B <- matrix(c(Jpsi0psi1,Jpsi0psi2), nrow=1, ncol=2)
  C <- t(B)
  D <- matrix(c(Jpsi1psi1,Jpsi1psi2,Jpsi1psi2,Jpsi2psi2),nrow=2, ncol=2, byrow=TRUE)
  
  Iblock <- (D-C%*%(A^(-1))%*%B)
  try(Iblockinv <- solve(Iblock)) # 2x2
  
  sblock <- matrix(c(sum(s_psi1),sum(s_psi2)), nrow=2, ncol=1)  # 2x1
  
  try(TestStat <- t(sblock)%*%Iblockinv%*%sblock)
  
  #Return:
  TestStat
  
}  #End: myTest.STEC.LM()

#NOTE: Do NOT use!
myTest.CEC.LM <- function(e,tv,garch,stec,opt){
  # test of Ho: CEC vs H1: CCC
  # NOTE: myTest.CEC.LM only written for no-TV no-GARCH CEC-only case
  
  TestStat <- NA
  
  N <- NCOL(e)
  Tobs <- NROW(e)
  
  STEC <- stec
  Q <- myQ.EC(N)
  rho <- STEC$pars[1]
  L <- myL.EC(N,rho)
  P <- Q%*%diag(L)%*%t(Q)
  
  # option 1:
  if (opt == 1){
    Z <- Zt <- matrix(0,nrow=N,ncol=N)
    #for (t in seq(1,Tobs)){ Zt <- Zt + e[t,]%*%t(e[t,]) }
    #TODO: Generalize the code below for all N
    #Faster performance, knowing N=3:
    a <- e[,1]
    b <- e[,2]
    c <- e[,3]
    aa <- sum(a*a)
    bb <- sum(b*b)
    cc <- sum(c*c)
    ab <- sum(a*b)
    bc <- sum(b*c)
    ac <- sum(a*c)
    
    Zt[1,] <- c(aa,ab,ac)
    Zt[2,] <- c(ab,bb,bc)
    Zt[3,] <- c(ac,bc,cc)
    
    Z <- Zt/Tobs
    
    II <- rep(1,N)%*%t(rep(1,N)) # NxN, 11' 
    I <- diag(1,nrow=N,ncol=N)   # NXN, I
    S <- 2*rbind( cbind((1-rho^2)*I+rho^2*II , rho*sqrt(2)*((rho-1)*I+II)), cbind(rho*sqrt(2)*((rho-1)*I+II) , P+rho^2*II) ) 
    Sinv <- solve(S)
    
    vZ <- myVecd(Z)
    vP <- myVecd(P)
    
    TestStat <- Tobs/2*(t(vZ-vP))%*%Sinv%*%(vZ-vP)
  }
  # option 2:
  if (opt == 2){
    Z <- Zt <- matrix(0,nrow=N,ncol=N)
    St <- matrix(0,ncol=N*(N+1)/2,nrow=N*(N+1)/2)
    #TODO: Optimize the code below to remove loop!
    for (t in seq(1,Tobs)){
      Zt <- Zt + e[t,]%*%t(e[t,])
      St <- St + myVecd(e[t,]%*%t(e[t,]))%*%t(myVecd(e[t,]%*%t(e[t,])))
    }
    
    Z <- Zt/Tobs
    S <- St/Tobs
    Sinv <- solve(S)
    
    vZ <- myVecd(Z)
    vP <- myVecd(P)
    
    TestStat <- Tobs/2*(t(vZ-vP))%*%Sinv%*%(vZ-vP)
  }
  # return:
  TestStat
  
}  #End: myTest.CEC.LM()

} # End: TestFn
#==========================================================#

#----------------------------------------------------------#
#    TVGJR - Log-Liklihood Calculator Functions 
#----------------------------------------------------------#
zLogLikFn <- TRUE
if (zLogLikFn) {

  
myLogLik.tv_univar <- function(optimpars,e,tv,return_ll=TRUE) {
  # model: TV (no GARCH), univariate
  # purpose: compute loglikelihood (and related components) value
  # input: optimpars   -- parameter vector including: (d0,delta(i),speed(i),loc1(i),<loc2(i)>) for all G(i)
  #        e           -- returns (not standardised)
  #        tv          -- TV object 
  # output:       return_ll   -- return value indicator:
  #               TRUE  : sum of loglikelihood(t) values over t=1...T
  #               FALSE : conditional variances, T rows
  # Note:  tv$var_target: TRUE => delta0 is free => include delta0 in optimpars
  
  #Create a new TV list within the scope of this function:
  TV <- tv
  
  # Check if variance targetting is being used. If not, set default to optimise ALL parameters.
  if (is.null(TV$var_target)) TV$var_target <- TRUE
  
  # Extract the parameters from optimpars back into the tv() object - then use tv() to calculate ll:
  if(!is.null(TV$linpars)) {  
      TV$delta0 <- optimpars[1]
      if(!is.null(TV$shape)) {
        lengthTVpars <- length(TV$pars[!is.na(TV$pars)])
        TV$pars[!is.na(TV$pars)] <- optimpars[2:(lengthTVpars+1)]    
      }
      TV$linpars[!is.na(TV$linpars)] <- tail(optimpars,length(TV$linpars[!is.na(TV$linpars)]))
  } else if (isTRUE(TV$var_target)) {
    # delta0 is a free parameter
    TV$delta0 <- optimpars[1]
    if(!is.null(TV$shape)) TV$pars[!is.na(TV$pars)] <- tail(optimpars,-1)
  } else {
    # TV$delta0 is used in calculation of LogLiklihood, but is not in optimpars
    if(!is.null(TV$shape)) TV$pars[!is.na(TV$pars)] <- optimpars else stop("You can't estimate TV with no $shape (just $delta0) and var_target = FALSE!")
  }
  
  # Check 1. Check that delta0 is positive
  if (TV$delta0 < 0) return(err_output)
  
  # Checks 2 - 6 (if we have any more parameters)
  if (!is.null(TV$shape)) {
    
    pars <- TV$pars
    speedoption <- TV$speedoption
    shape <- TV$shape
    st <- TV$st
    
    sp_pos <- l1_pos <- l2_pos <- vector("numeric")
    for (i in seq_along(shape))
    {
      #delta_pos <- c(delta_pos, 4*(i-1)+1)
      sp_pos <- c(sp_pos, 4*(i-1)+2)
      l1_pos <- c(l1_pos, 4*(i-1)+3)
      l2_pos <- c(l2_pos, 4*(i-1)+4)  #May be NA
    }
    
    # Check 2: Check the boundary values for speed params:
    #speedoption -- 1=gamma, 2=gamma/std(st), 3=exp(eta), 4=1/lambda^2
    maxSpeed <- switch(speedoption,1000,(1000/sd(st)),7.0,0.30)
    if (max(pars[sp_pos]) > maxSpeed) return(err_output)
    if (min(pars[sp_pos]) < 0) return(err_output)
    
    # Check 3: Check the loc1 locations fall within min-max values of st
    # We must have at least 1 loc1 to be inside this shape..loop, so no need to check if loc1 contains a valid value:
    if (min(pars[l1_pos]) < min(st)) return(err_output)
    if (max(pars[l1_pos]) > max(st)) return(err_output)
    
    # Switch OFF Warnings() for the loc2 checks as the NA element causes many warnings
    options(warn=-1)
    
    # Check 4: Check that loc1.1 < loc2.1 < loc3.1... for all G(i)
    # Method: Subtract loc1_pos vector from loc2_pos vector and ensure it is positive:
    tmp <- pars[l2_pos] - pars[l1_pos]
    # Note: tmp will contain NA wherever a loc2 element was NA - we can strip these out:
    tmp <- tmp[!is.na(tmp)]
    if (any(tmp)<0) return(err_output)
    
    # Check 5: Check the loc2 locations fall within min-max values of st
    # Confirm we have at least one valid numeric loc 2, before checking min & max:
    if (is.finite(min(pars[l2_pos],na.rm=TRUE))) {
      if (min(pars[l2_pos]) < min(st)) return(err_output)
      if (max(pars[l2_pos]) > max(st)) return(err_output)
    }
    
    # Switch ON Warnings() again
    options(warn=0)
    
    # Check 6: Check that loc1.1 < loc1.2 where 2 locations exist... for all G(i)
    # We do need to have at least 2 loc1.n locations for this error check
    tmp <- pars[l1_pos]
    if (length(tmp) > 1) {
      for (n in 2:length(tmp)) if (tmp[n-1] > tmp[n]) return(err_output)
    }
    
    
  }  # End: If (!is.null(shape))
  
  ###  ALL ERROR CHECKS PASSED  ###
  
  # recursively (over t=1...Tobs) compute g(t)=delta0 [+ delta1*G1(t) [+ delta2*G2(t) [+ ... ]]]
  g <- calculate_g(TV)
  
  #Check that g is positive! - No negative elements allowed:
  if (min(g) < 0) return(err_output)
  
  #Return g, if requested:
  if (!return_ll) return(g)
  
  # Calculate ll, if we need to return it:
  ll <- sum(-0.5*log(2*pi) - 0.5*log(g) - 0.5*(e*e)/g)
  
}  #End: myLogLik.tv_univar
  
myLogLik.garch_univar <- function(optimpars,w,garch,return_ll=TRUE){
  # model: GARCH
  # purpose: compute loglikelihood (and related components) value
  # input: optimpars = parameter vector: (<omega>,<alpha>,beta,<delta>)
  #                    Note: The actual optimpars parameters passed-in depend on the Garch$type and Garch$var_target
  #                    garch$pars will always contain 4 elements as follows:
  #                       omega <- optimpars[1]
  #                       alpha <- optimpars[2]
  #                       beta <- optimpars[3]
  #                       delta <- optimpars[4]
  #        w = returns/sqrt(g) = e/sqrt(g)
  #        type = 1:GARCH or 2:GJR or 3:GJR with alpha=0
  
  err_output <- -1e9
  
  # Check if variance targetting is being used. If not, set default to optimise ALL parameters.
  if (is.null(garch$var_target)) garch$var_target <- FALSE
  
  # Extract parameters from optimpars & do error-checking:
  
  if (!garch$var_target) {
    # VarienceTargetting is OFF: All Garch parameters are optimised
    if (garch$type==1) {
      omega <- optimpars[1]
      alpha <- optimpars[2]
      beta <- optimpars[3]
      delta <- 0
    } else if (garch$type==2) {
      omega <- optimpars[1]
      alpha <- optimpars[2]
      beta <- optimpars[3]
      delta <- optimpars[4]
      #Type2-Specific Error Checks:
      if (alpha+delta < 0) return(err_output)
      if (alpha+delta/2+beta >= 1) return(err_output)
    } else if (garch$type==3) {
      omega <- optimpars[1]
      alpha <- 0
      beta <- optimpars[2]
      delta <- optimpars[3]
    }
  } else {
    # VarienceTargetting is ON: omega is not passed in as a free param - it is calculated below.
    if (garch$type==1) {
      alpha <- optimpars[1]
      beta <- optimpars[2]
      delta <- 0
    } else if (garch$type==2) {
      alpha <- optimpars[1]
      beta <- optimpars[2]
      delta <- optimpars[3]
      #Type2-Specific Error Checks:
      if (alpha+delta < 0) return(err_output)
      if (alpha+delta/2+beta >= 1) return(err_output)
    } else if (garch$type==3) {
      alpha <- 0
      beta <- optimpars[1]
      delta <- optimpars[2]
    }
    omega <- calculateOmega(alpha,beta,delta)
  }
  
  ## Universal Error Checks:
  if (omega <= 0) return(err_output)
  if (alpha < 0) return(err_output)
  if (beta < 0) return(err_output)
  ## End of Error Checks 
  
  Tobs <- length(w)
  
  h <- rep(0,Tobs)
  h[1] <- sum(w*w)/Tobs
  # min(w[t-1],0) => Only return negative values or 0
  for(t in 2:Tobs) h[t] <- omega + alpha*(w[t-1])^2 + beta*h[t-1] + delta*(min(w[t-1],0))^2
  if (!return_ll) return(h)
  
  #Return
  ll <- sum(-0.5*log(2*pi) -0.5*log(h) -0.5*(w*w)/h)
  
}  #End: myLogLik.garch_univar()

myLogLik.garch_univar_targetomega <- function(optimpars,w,garch,localvar,return_ll=TRUE){
  # model: GARCH
  # purpose: compute loglikelihood (and related components) value
  # input: optimpars = parameter vector: (<omega>,<alpha>,beta,<delta>)
  #        w = returns/sqrt(g) = e/sqrt(g)
  #                    Note: The actual optimpars parameters passed-in depend on the Garch$type and Garch$var_target
  #                    garch$pars will always contain 4 elements as follows:
  #                       omega <- optimpars[1]
  #                       alpha <- optimpars[2]
  #                       beta <- optimpars[3]
  #                       delta <- optimpars[4]
  #        garch$type = 1:GARCH or 2:GJR or 3:GJR with alpha=0
  #        localvar = vector; The local variance around each observation 
  
  err_output <- -1e9
  
  # Extract parameters from optimpars & do error-checking:
    
  # Omega is not passed in as a free param - it is calculated below, to target the variance of a local window
  if (garch$type==1) {
    alpha <- optimpars[1]
    beta <- optimpars[2]
    delta <- 0
  } else if (garch$type==2) {
    alpha <- optimpars[1]
    beta <- optimpars[2]
    delta <- optimpars[3]
    #Type2-Specific Error Checks:
    if (alpha+delta < 0) return(err_output)
    if (alpha+delta/2+beta >= 1) return(err_output)
  } else if (garch$type==3) {
    alpha <- 0
    beta <- optimpars[1]
    delta <- optimpars[2]
  }
  
  ## Universal Error Checks:
  if (alpha < 0) return(err_output)
  if (beta < 0) return(err_output)
  if (alpha+beta+0.5*delta >=1) return(err_output)
  ## End of Error Checks 
  
  Tobs <- length(w)
  omega <- localvar*(1-alpha-beta-(delta/2))
  
  h <- rep(0,Tobs)
  h[1] <- sum(w*w)/Tobs
  # min(w[t-1],0) => Only return negative values or 0
  for(t in 2:Tobs) {
    h[t] <- omega[t] + alpha*(w[t-1])^2 + beta*h[t-1] + delta*(min(w[t-1],0))^2
  }
  if (!return_ll) return(h)
  
  #Return
  sum(-0.5*log(2*pi) -0.5*log(h) -0.5*(w*w)/h)
  
}  #End: myLogLik.garch_univar_targetomega()


myLogLik.arch_univar <- function(optimpars,w,arch,return_ll=TRUE){
  # model: ARCH
  # purpose: compute loglikelihood (and related components) value
  # input: optimpars = parameter vector: (<omega>,<alpha>,beta,<delta>)
  #                    Note: The actual optimpars parameters passed-in depend on the arch$type and arch$var_target
  #                    arch$pars will contain n elements
  #        w = returns/sqrt(g) = e/sqrt(g)
  #        type = 1:ARCH
  
  err_output <- -1e9
  
  ## Universal Error Checks:
  if (min(optimpars) <= 0) return(err_output)
  ## End of Error Checks 
  
  Tobs <- length(w)
  h <- rep(0,Tobs)

  #ARCH Type: 1
  if (arch$type==1){
    h[1:(length(optimpars)-1)] <- sum(w*w)/Tobs
   
    for(t in length(optimpars):Tobs){
      h[t] <- optimpars[1] + optimpars[2]*(w[t-1])^2 + optimpars[3]*(w[t-2])^2  + optimpars[4]*(w[t-3])^2 
      + optimpars[5]*(w[t-4])^2 + optimpars[6]*(w[t-5])^2  
    }
  }
  
  #ARCH Type: 2
  if (arch$type==2){
    h[1] = sum(w*w)/Tobs
    for(t in 2:Tobs){
      h[t] <- optimpars[1]
      for(i in 1:(t-1)){
        tmp <- optimpars[2]*optimpars[3]^(i-1)*(w[t-i])^2 
        if (tmp >= 1e-7) h[t] <- h[t] + tmp else break
      }
    }
  }
  
  
  if (!return_ll) return(h)
  
  #Return
  ll <- sum(-0.5*log(2*pi) -0.5*log(h) -0.5*(w*w)/h)
  
}  #End: myLogLik.arch_univar()


myLogLik.garch_univar_c <- cmpfun(myLogLik.garch_univar)

myLogLik.arch_univar_c <- cmpfun(myLogLik.arch_univar)


myLogLik.stec <- function(pars,z,stec,return_ll=TRUE){
  # model: STEC
  # input: pars        -- c(rho1,rho2,speed,loc1,[loc2])
  #        z           -- volatility standardised returns (matrix TxN)
  #        shape       -- scalar, 1: (st-c), 2: (st-c1)(st-c2), 3: (st-c)^2
  #        speedoption -- scalar, 1=gamma, 2=gamma/std(st), 3=exp(eta), 4=1/lambda^2
  #        st          -- transition variable
  #        output: sum of log-likelihoods over time (ll) = (Default)
  #        or correlations T x N(N-1)/2 (Pt)
  
  STEC <- stec
  
  err_output <- -1e9
  Tobs <- NROW(z)
  N <- NCOL(z)
  
  #Initialise return value to the Error:
  ll <- err_output
  
  # - - - CEC - - -
  if (STEC$type==0){
    rho <- pars
    Q <- myQ.EC(N)
    L <- myL.EC(N,rho)
    P <- Q%*%diag(L)%*%t(Q)
    if (min(L) <= 0) return(err_output)
    
    # - - - P(t) and loglik-value
    llt <- rep(0,Tobs)
    Pinv <- Q%*%diag(1/L)%*%t(Q)
    detP <- prod(L)
    for(t in seq(1,Tobs))
    {
      llt[t] <- -0.5*log(detP)-0.5*( t(z[t,])%*%(Pinv)%*%z[t,])
    }
  }
  
  # - - - STEC - - -
  if (STEC$type==1){
    st <- STEC$st
    shape <- STEC$shape
    speedoption <- STEC$speedoption
    
    rho <- pars[1]
    rhostar <- pars[2]
    rhot <- rep(0,Tobs)
    L1 <- myL.EC(N,(rho-rhostar))
    if (min(L1) <= 0) return(err_output)
    L2 <- myL.EC(N,(rho+rhostar))
    if (min(L2) <= 0) return(err_output)
    Q <- myQ.EC(N)
    
    speed <- pars[3]
    loc1 <- pars[4]
    loc2 <- pars[5]  #May be NA
    
    # - - - Gt - - -
    Gt <- myG(speed,loc1,loc2,st,shape,speedoption)
    
    # - - - P(t) and loglik-value
    llt  <- rep(0,Tobs)
    for(t in seq(1,Tobs))
    {
      rhot[t] <- (1-Gt[t])*(rho-rhostar) + Gt[t]*(rho+rhostar)
      Lt <- myL.EC(N,rhot[t])
      Ptinv <- Q%*%diag(1/Lt)%*%t(Q)
      detPt <- prod(Lt)
      llt[t] <- -0.5*log(detPt)-0.5*( t(z[t,])%*%(Ptinv)%*%z[t,])
    }
  }
  
  ll <- sum(llt)
  if (return_ll) return(ll)
  else {
    if (STEC$type==0) Pt <- matrix(rho,nrow=Tobs,ncol=(N*(N-1)/2))
    if (STEC$type==1) Pt <- matrix(rhot,nrow=Tobs,ncol=(N*(N-1)/2),byrow=FALSE)
    return(Pt)
  }
} #End: myLogLik.stec

myLogLik.stcc <- function(pars,z,stcc,return_ll=TRUE){
  # model: STCC
  # input: pars        -- c(rhovec1,rhovec2,speed,loc1,[loc2])
  #        z           -- volatility standardised returns (matrix TxN)
  #        shape       -- scalar, 1: (st-c), 2: (st-c1)(st-c2), 3: (st-c)^2
  #        speedoption -- scalar, 1=gamma, 2=gamma/std(st), 3=exp(eta), 4=1/lambda^2
  #        st          -- transition variable
  #        output: sum of log-likelihoods over time (ll) = (Default)
  #        or correlations T x N(N-1)/2 (Pt)
  
  
  # functionsPath <- file.path(dirname(getwd()),"Functions")
  # functionsFile <- file.path(functionsPath,"functions_tvgjr_v4.r")
  # source(functionsFile, local=TRUE)
  
  STCC <- stcc
  speedoption <- STCC$speedoption
  shape <- STCC$shape
  st <- STCC$st
    
  if(shape==2) numTRpars <- 3 else numTRpars <- 2
  lengthPars <- length(pars)
    if(shape==2) {
      l2_pos <- lengthPars
      l1_pos <- lengthPars -1
      sp_pos <- lengthPars -2
    } else {
      l2_pos <- NA
      l1_pos <- lengthPars 
      sp_pos <- lengthPars -1
    }

  Tobs <- NROW(z)
  N <- NCOL(z)
  
  #Initialise return value to the Error:
  ll <- err_output
  
  ## ERROR CHECKS & CONSTRAINTS  ##
  
  # Switch OFF Warnings() for the loc2 checks as the NA element causes many warnings
  options(warn=-1)
  
  # Checks 2 - 6 (if we have any more parameters)
  if (!is.null(shape)) {
      
    # Check 2: Check the boundary values for speed params:
    #speedoption -- 1=gamma, 2=gamma/std(st), 3=exp(eta), 4=1/lambda^2
    maxSpeed <- switch(speedoption,1000,(1000/sd(st)),7.0,0.30)
    if (max(pars[sp_pos]) > maxSpeed) return(err_output)
    if (min(pars[sp_pos]) < 0) return(err_output)
    
    # Check 3: Check the locations fall within min-max values of st
    # We must have at least 1 loc1 to be inside this shape..loop, so no need to check if loc1 contains a valid value:
    if (min(pars[l1_pos]) < min(st)) return(err_output)
    if (max(pars[l1_pos]) > max(st)) return(err_output)
    if (!is.na(l2_pos)) {
      if (min(pars[l2_pos]) < min(st)) return(err_output)
      if (max(pars[l2_pos]) > max(st)) return(err_output)  
    }
    
    # Switch ON Warnings() again
    options(warn=0)
    
  }  # End: If (!is.null(shape))
  
  ###  ALL ERROR CHECKS PASSED  ###
    
  # - - - CCC - - -
  if (STCC$type==0){
    vP <- c(0,0,0) #pars
    mP <- myUnVecl(vP)
    eig <- eigen(mP,symmetric=TRUE,only.values = TRUE)
    if (min(eig$values) <= 0) return(err_output)
    
    # - - - P(t) and loglik-value
    llt <- rep(0,Tobs)
    mPinv <- solve(mP)
    detmP <- det(mP)
    #TODO: Can we optimise the Loop below?
    const <- -0.5*log(detmP)
    for(t in seq(1,Tobs))
    {
      llt[t] <- const - 0.5*(t(z[t,])%*%(mPinv)%*%z[t,])
    }
    
  }
  
  # - - - STCC - - -
  if (STCC$type==1){
    vP1 <- pars[1:length(STCC$P1)]
    vP2 <- pars[(length(STCC$P1)+1):(2*length(STCC$P2))]
    TRpars <- tail(pars,numTRpars)
    
    Pt <- matrix(0,Tobs,N*(N-1)/2)
    speed <- TRpars[1]
    loc1 <- TRpars[2]
    if(length(TRpars)==3) loc2 <- TRpars[3] else loc2 <- NA
    
    mP1 <- myUnVecl(vP1)
    eig1 <- eigen(mP1,symmetric=TRUE,only.values=TRUE)
    # Check for SPD - positive-definite check:
    if (min(eig1$values) <= 0) return(err_output)
    mP2 <- myUnVecl(vP2)
    eig2 <- eigen(mP2,symmetric=TRUE,only.values=TRUE)
    # Check for SPD - positive-definite check:
    if (min(eig2$values) <= 0) return(err_output)
    
    # - - - Gt - - -
    Gt <- myG(speed,loc1,loc2,st,shape,speedoption)
    
    # - - - P(t) and loglik-value
    llt <- NULL
    
     for(t in seq(1,Tobs))
     {
       Pt[t,] <- (1-Gt[t])*vP1 + Gt[t]*vP2
       mPt <- myUnVecl(Pt[t,])
       #mPtinv <- solve(mPt)
       mPtinv <- chol2inv(chol(mPt))
       llt[t] <- -0.5*log(det(mPt)) -0.5*( t(z[t,])%*%(mPtinv)%*%z[t,])
     } # End: for loop
    
    #llt <- foreach (t = 1:Tobs, .combine=rbind, .verbose = FALSE) %dopar% { 
    #  Pt[t,] <- (1-Gt[t])*vP1 + Gt[t]*vP2
    #  mPt <- myUnVecl(Pt[t,])
    #  #mPtinv <- solve(mPt)
    #  mPtinv <- chol2inv(chol(mPt))
    #  -0.5*log(det(mPt)) -0.5*( t(z[t,])%*%(mPtinv)%*%z[t,])
    #} # End: foreach loop
    
    
  } # End: if(STCC$type==1)
  
   
  if (return_ll) {
    rm(Pt,STCC) 
    return(sum(llt))
    }
  else {
    if (STCC$type==0) Pt <- matrix(vP,nrow=Tobs,ncol=length(vP),byrow=TRUE)
    return(Pt)
  }

}  #End: myLogLik.stcc()


myLogLik.multivar.tv <- function(optimpars,e,ntv,P=NULL,return_ll=TRUE){
  ###  Inputs:  ###
  # optimpars       -- vector of params for all N TV objects (Pars will be the optimised pars from the Specification process)
  # e               -- matrix of N data-series x T obsevations (panel)
  # ntv             -- list of N x TV objects
  # P               -- correlation matrix

  Tobs <- NROW(e)   # time dimension
  N <- NCOL(e)      # number of variables
  
  h <- matrix(data = 1, nrow = Tobs, ncol = N)  # Not used.  Included to keep the log-liklihood calculation the same everywhere.
  g <- matrix(nrow = Tobs, ncol = N)
  if (is.null(P)) P <- matrix(data=1,nrow = Tobs, ncol = N)
  z <- matrix(nrow = Tobs, ncol = N)
  
  err_output <- -1e9
  nTV <- ntv
  
  endpos <- 0
  for (n in 1:N) {
    # Extract the TV$pars for each series & pass to optim:
    startpos <- endpos + 1
    endpos <- length(nTV[[n]]$pars[!is.na(nTV[[n]]$pars)])
    npars <- optimpars[startpos:endpos]
    # Get univariate estimate & update optimpars vector:
    tmp <- NULL
    try(tmp <- myLogLik.tv_univar(npars,e[,n],nTV[[n]],return_ll))
    
    if (!is.null(tmp)) if (tmp$convergence==0) { 
      # Update the optimpars vector and the g(t) matrix
      optimpars[startpos:endpos] <- tmp$par
      g[,n] <- myLogLik.tv_univar(tmp$par,e[,n],nTV[[n]],return_ll=FALSE)
      z[,n] <- e[,n]/sqrt(g[,n])
      #
    } else {
      # Failed to get an estimate for this series!!
      strWarning <- paste0("Data Series ",n," failed to converge.  We will continue multivariate optimisation using the starting params provided")
      warning(strWarning)
    }
    
  } #End: for(n in 1:N)
  

  # Calculate liklihood for full model
  ll <- err_output
  llt <- rep(0,Tobs)
  
  for(t in 1:Tobs) {
    mPt <- myUnVecl(P[t,])
    mPtinv <- solve(mPt)
    llt[t] <- -0.5*log(2*pi)*N -0.5*sum(log(h[t,])) -0.5*sum(log(g[t,])) -0.5*log(det(mPt))-0.5*(t(z[t,])%*%(mPtinv)%*%z[t,])
  }
  
  ll <- sum(llt)
  
  #Return:
  if (return_ll) rtn <- ll else rtn <- g
  rtn
  
}  #End: myLogLik.multivar.tv()


myLogLik.multivar.old <- function(npars,e,ntv,ngarch,stcc,target=NULL){
  ###  Inputs:  ###
  # npars                   -- vector of params for all N (matching 'target') that can be passed to optim?  TODO - Confirm this!!
  # target                  -- Used to...
  # var_target      -- FALSE means delta0 is Not free 
  
  ###  Common code - Start  ###
  if (is.null(target)) stop("You must supply a value for Target = one of 'TV', 'GARCH', 'STCC' or 'ALL'")
  
  Tobs <- NROW(e) # time dimension
  N <- NCOL(e)        # number of variables
  
  h <- matrix(nrow = Tobs, ncol = N)
  g <- matrix(nrow = Tobs, ncol = N)
  w <- matrix(nrow = Tobs, ncol = N)
  z <- matrix(nrow = Tobs, ncol = N)
  
  err_output <- -1e9
  return_ll <- FALSE  #This function only calculates cond.variances for now.
  
  nTV <- ntv
  nGARCH <- ngarch
  STCC <- stcc
  
  if(target=="TV") {
    
    for (n in 1:N) {
      #TV
      optimpars <- npars[4*(n-1)+1:4*(n-1)+4]
      optimpars <- optimpars[!is.na(optimpars)]
      #myLogLik.tv_univar <- function(optimpars,e,tv,return_ll=TRUE)
      tmp <- myLogLik.tv_univar(optimpars,e[,n],nTV[[n]],return_ll)
      #TODO: Ensure the loglik() functions return a 'known' value for error/failure to calc g
      if (length(tmp)==1) return(err_output) else
      {
        g[,n] <- tmp
        w[,n] <- e[,n]/sqrt(g[,n])
      }
      #GARCH
      optimpars <- nGARCH[[n]]$pars
      optimpars <- optimpars[!is.na(optimpars)]
      #myLogLik.garch_univar <- function(pars,w,garch,return_ll=TRUE){
      tmp <- myLogLik.garch_univar(optimpars,w[,n],nGARCH[[n]],return_ll)
      if (length(tmp)==1) return(err_output) else {
        h[,n] <- tmp
        z[,n] <- w[,n]/sqrt(h[,n])
      }
    } #End: for(n in 1:N)
    
    #STCC
    tmp <- myLogLik.stcc(STCC$pars,z,STCC,return_ll)
    if (length(tmp)==1) return(err_output) else P <- tmp
    
  } #End: if(target=="TV")
  
  if(target=="GARCH"){
    
    for (n in  1:N){
      g[,n] <- nTV[[n]]$condvars
      w[,n] <- e[,n]/sqrt(g[,n])
      
      optimpars <- npars[4*(n-1)+1:4*(n-1)+4]
      optimpars <- optimpars[!is.na(optimpars)]
      tmp <- myLogLik.garch_univar(optimpars,w[,n],nGARCH[[n]],return_ll)
      if (length(tmp)==1) return(err_output)
      else {
        h[,n]<-tmp
        z[,n] <- w[,n]/sqrt(h[,n])
      }
    }
    
    tmp <- myLogLik.stcc(STCC$pars,z,STCC,return_ll)
    if (length(tmp)==1) return(err_output)
    else P <- tmp
    
  } #End: if(target=="GARCH")
  
  if(target=="STCC"){
    
    for (n in  1:N) {
      g[,n] <- nTV[[n]]$condvars
      w[,n] <- e[,n]/sqrt(g[,n])
      h[,n] <- nGARCH[[n]]$condvars
      z[,n] <- w[,n]/sqrt(h[,n])
    }
    #TODO: What bit of npars goes in here??
    tmp <- myLogLik.stcc(npars,z,STCC,return_ll)
    if (length(tmp)==1) return(err_output) else P <- tmp
    
  }  #End: if(target=="STCC")
  
  if(target=="ALL"){
    
    ## Use tv_univar to get TV$condvars (g), then calculate w using this.
    ## Then use garch_univar to get GARCH$condvars (h), then calculate z using this.
    ## Use myLogLik.STCC to calc P using the STCC$pars & z
    ## Finally - calc the Optimum LogLikelihood using all these values
    ## !! Note: Ensure all values used exist in the pars vector passed into this function as the 1st argument
    ##          Format = (nTV[[1]]$pars,..,nTV[[N]]$pars,nGARCH[[1]]$pars,..,nGARCH[[N]]$pars,STCC$pars)
    
    for (n in 1:N) {
      #TV
      optimpars <- npars[4*(n-1)+1:4*(n-1)+4]
      optimpars <- optimpars[!is.na(optimpars)]
      #myLogLik.tv_univar <- function(optimpars,e,tv,return_ll=TRUE)
      tmp <- myLogLik.tv_univar(optimpars,e[,n],nTV[[n]],return_ll)
      #TODO: Ensure the loglik() functions return a 'known' value for error/failure to calc g
      if (length(tmp)==1) return(err_output) else
      {
        g[,n] <- tmp
        w[,n] <- e[,n]/sqrt(g[,n])
      }
      #GARCH
      start_pos <- N*4
      optimpars <- npars[start_pos+(4*(n-1))+start_pos+(1:4*(n-1)+4)]
      optimpars <- optimpars[!is.na(optimpars)]
      #myLogLik.garch_univar <- function(pars,w,garch,return_ll=TRUE){
      tmp <- myLogLik.garch_univar(optimpars,w[,n],nGARCH[[n]],return_ll)
      if (length(tmp)==1) return(err_output) else {
        h[,n] <- tmp
        z[,n] <- w[,n]/sqrt(h[,n])
      }
    } #End: for(n in 1:N)
    
    #STCC
    optimpars <- tail(npars,length(STCC$pars))
    optimpars <- optimpars[!is.na(optimpars)]
    tmp <- myLogLik.stcc(optimpars,z,STCC,return_ll)
    if (length(tmp)==1) return(err_output) else P <- tmp
    
  } #End: if(target=="ALL")
  
  
  # Calculate liklihood for full model
  ll <- err_output
  llt <- rep(0,Tobs)
  
  for(t in seq(1,Tobs))
  {
    mPt <- myUnVecl(P[t,])
    mPtinv <- solve(mPt)
    llt[t] <- -0.5*log(2*pi)*N -0.5*sum(log(h[t,])) -0.5*sum(log(g[t,])) -0.5*log(det(mPt))-0.5*(t(z[t,])%*%(mPtinv)%*%z[t,])
  }
  
  ll <- sum(llt)
  
  if (!return_ll) {
    #Return a list of g[],h[],P
    rtn <- list()
    rtn$g <- g
    rtn$h <- h
    rtn$P <- P
    return(rtn)
  }
  
  #Return:
  ll
  
}  #End: myLogLik.multivar()



} #End: LogLikFn
#==========================================================#

#===========================================================================#
#    Wrapper Functions below - These call the lower level functions above.
#   - Try to use these functions instead of call the ones above directly -
#===========================================================================#
zEstimateFn <- TRUE
if (zEstimateFn) {

EstimateTV <- function(e,tv,calcHess=FALSE) {
  TV <- tv
  if (is.null(TV$var_target)) TV$var_target <- TRUE  #Default value 
  if (is.null(TV$optimcontrol)) stop("TV$optimcontrol is mandatory.")
  
  ###
  ### ---  Set up 'optimpars' and call optim to calculate the estimate --- ###
  ###
  
  # First we check for the simple case of just delta0 provided, no TV$pars or TV$linpars:
  if (is.null(TV$shape)) {
    #TV Ordser 0 function - delta0 is the only param, so skip the optim()
    TV$delta0 <- var(e) #Quick estimate for delta0 only
    TV$condvars <- rep(TV$delta0,NROW(e))
    TV$value <- myLogLik.tv_univar(TV$delta0,e,TV,return_ll=TRUE)
    TV$error <- FALSE
    return(TV)
    #else run the optim() to get the estimate
  }

  # Next we need to check if we are handling linearized parameters (Note: var_target does not apply in this mode):
  if (!is.null(TV$linpars)) {
	  optimpars <- c(TV$delta0,TV$pars,TV$linpars)
    optimpars <- optimpars[!is.na(optimpars)]
    myControl <- TV$optimcontrol
  } else {
    # Handle variance targeting - if $var_target is False, strip out the delta0 from parameters & mycontrol
    if (TV$var_target) {
      optimpars <- c(TV$delta0,TV$pars) 
      myControl <- TV$optimcontrol
    } else {
      optimpars <- TV$pars
      myControl <- TV$optimcontrol
      myControl$parscale <- TV$optimcontrol$parscale[-1]
      myControl$ndeps <- TV$optimcontrol$ndeps[-1]
    }
    # Finally - strip out any NA values:
    optimpars <- optimpars[!is.na(optimpars)]
  }
 
  # Now call optim:
  tmp <- NULL
  #tmp <- optim(optimpars,myLogLik.tv_univar,e,TV,gr=NULL,method="BFGS",control=myControl,hessian=calcHess)
  try(tmp <- optim(optimpars,myLogLik.tv_univar,e,TV,gr=NULL,method="BFGS",control=myControl,hessian=calcHess))      
 
  
  ###    
  ### ---  Interpret the response from optim --- ###
  ###
  
  # An unhandled error could result in a NULL being returned by optim()
  if (is.null(tmp)) {
    TV$delta0 <- NA
    TV$pars <- NA
    if (!is.null(TV$linpars)) TV$linpars <- NA
    TV$value <- -1e9   #Set the value = the ErrorOutput value used by the myLogLik.tv_univar() function
    TV$error <- TRUE
    return(TV)
  }
  
  if (tmp$convergence==0) { 
    #Optim converged successfully => we expect tmp$par to have good estimates!
    TV$value <- tmp$value
    TV$error <- FALSE
    # Now get the conditional variances
    TV$condvars <- myLogLik.tv_univar(tmp$par,e,TV,return_ll=FALSE)
    
    #Update the TV object paramters using optimised pars
    if (!is.null(TV$linpars)) {
      TV$delta0 <- tmp$par[1]
      if(!is.null(TV$pars)) TV$pars[!is.na(TV$pars)] <- tmp$par[2:(length(TV$pars[!is.na(TV$pars)])+1)]
      TV$linpars[!is.na(TV$linpars)] <- tail(tmp$par,length(TV$linpars[!is.na(TV$linpars)]))
    } else {
      if (!TV$var_target) TV$pars[!is.na(TV$pars)] <- tmp$par else { 
      TV$delta0 <- tmp$par[1]
      TV$pars[!is.na(TV$pars)] <- tail(tmp$par,-1)  
      }
    }
    
    if (calcHess) {
      TV$hessian <- tmp$hessian
      try(TV$stderr <- sqrt(-diag(solve(tmp$hessian))))  #Use try(..) here as solve can fail!
      if (is.null(TV$stderr)) warning("Failed to calculate TV$stderr - probably a numerical solve() problem")
    }
  } else { 
    #Failed to converge 
    TV$delta0 <- NA
    TV$pars <- NA
    if (!is.null(TV$linpars)) TV$linpars <- NA
    TV$value <- err_output
    TV$error <- TRUE
    warning("EstimateTV() - failed to converge. Check the optim controls & starting params ")
    return(TV)
  }
  rm(tmp)
  
  #Return:
  TV
}  #End: EstimateTV()


EstimateGARCH <- function(e,garch,calcHess=FALSE,targetOmega=FALSE) {
  
  ###
  ### ---  Call optim to calculate the estimate --- ###
  ###
  GARCH <- garch
  if (is.null(GARCH$var_target)) GARCH$var_target <- FALSE  #Default value
  if (is.null(GARCH$optimcontrol)) stop("GARCH$optimcontrol is mandatory.")
  
  ## First, set the optimpars, based on Garch$type & Garch$var_target:
  optimpars <- garchpars_to_optimpars(GARCH)
  ## Second, set the optim Controls, based on Garch$var_target:
  if (!GARCH$var_target) myControl <- GARCH$optimcontrol else {
    myControl <- GARCH$optimcontrol
    myControl$parscale <- GARCH$optimcontrol$parscale[-1]
    myControl$ndeps <- GARCH$optimcontrol$ndeps[-1]
  }
  tmp <- NULL
  if (targetOmega) {
    try (tmp <- optim(optimpars,myLogLik.garch_univar_targetomega,e,GARCH,gr=NULL,method="BFGS",control=myControl,hessian=calcHess))
  } else try (tmp <- optim(optimpars,myLogLik.garch_univar_c,e,GARCH,gr=NULL,method="BFGS",control=myControl,hessian=calcHess))
  
  ###    
  ### ---  Interpret the response from optim --- ###
  ###
  
  # An unhandled error could result in a NULL being returned by optim()
  if (is.null(tmp)) {
    GARCH$value <- -111333 #Set the value = the ErrorOutput value used by the myLogLik.garch_univar() function
    GARCH$error <- TRUE
    return(GARCH)
  }
  
  if (tmp$convergence==0) {
    #Optim converged successfully => we expect tmp$par to have good estimates!
    GARCH$pars <- optimpars_to_garchpars(GARCH,tmp$par)
    GARCH$condvars <- myLogLik.garch_univar_c(tmp$par,e,GARCH,return_ll=FALSE)
    GARCH$value <- tmp$value
    GARCH$error <- FALSE
    if (calcHess) {
      GARCH$hess <-tmp$hessian
      GARCH$stderr <- NULL
      try(GARCH$stderr <- sqrt(-diag(solve(tmp$hessian))))
      if (is.null(GARCH$stderr)) warning("Failed to calculate GARCH$stderr - probably a numerical solve() problem")
    }
    
  } else {
    #Failed to converge - return specific error code in the $value
    GARCH$value <- -111000-tmp$convergence
    GARCH$error <- TRUE
    return(GARCH)
  }
  
  rm(tmp)
  
  #Return:
  GARCH
}  #End: EstimateGARCH()


EstimateARCH <- function(e,arch,calcHess=FALSE) {
  
  ###
  ### ---  Call optim to calculate the estimate --- ###
  ###
  ARCH <- arch
  #if (is.null(GARCH$var_target)) GARCH$var_target <- FALSE  #Default value
  if (is.null(ARCH$optimcontrol)) stop("ARCH$optimcontrol is mandatory.")
  
  ## First, set the optimpars, based on Garch$type & Garch$var_target:
  optimpars <- ARCH$pars
  ## Second, set the optim Controls, based on Garch$var_target:
  myControl <- ARCH$optimcontrol
  tmp <- NULL
  try (tmp <- optim(optimpars,myLogLik.arch_univar_c,e,ARCH,gr=NULL,method="BFGS",control=myControl,hessian=calcHess))
  
  ###    
  ### ---  Interpret the response from optim --- ###
  ###
  
  # An unhandled error could result in a NULL being returned by optim()
  if (is.null(tmp)) {
    ARCH$value <- -111333 #Set the value = the ErrorOutput value used by the myLogLik.garch_univar() function
    ARCH$error <- TRUE
    return(ARCH)
  }
  
  if (tmp$convergence==0) {
    #Optim converged successfully => we expect tmp$par to have good estimates!
    ARCH$pars <- tmp$par
    ARCH$condvars <- myLogLik.arch_univar_c(tmp$par,e,ARCH,return_ll=FALSE)
    ARCH$value <- tmp$value
    ARCH$error <- FALSE
    if (calcHess) {
      ARCH$hess <-tmp$hessian
      ARCH$stderr <- NULL
      try(ARCH$stderr <- sqrt(-diag(solve(tmp$hessian))))
      if (is.null(ARCH$stderr)) warning("Failed to calculate ARCH$stderr - probably a numerical solve() problem")
    }
    
  } else {
    #Failed to converge - return specific error code in the $value
    ARCH$value <- -111000-tmp$convergence
    ARCH$error <- TRUE
    return(ARCH)
  }
  
  rm(tmp)
  
  #Return:
  ARCH
}  #End: EstimateARCH()


EstimateTVGARCH <- function(e,series,maxLoops=20,verbose=FALSE) {
  ###
  ### ---  Call optim for TV & GARCH in a loop until Converged --- ###
  ###
  
  TV <- series$tv
  GARCH <- series$garch
  
  #Set the variance targetting:
  TV$var_target <- series$var_target
  GARCH$var_target <- series$var_target
  
  # Initialise Repeat..loop control values:
  loopcount <- 0
  ll_const <- -0.5*log(2*pi)
  e2 <- e*e
  TV_Converged <- FALSE
  GARCH_Converged <- FALSE
  
  repeat {
    
    # parameter values before optimization
    if (!series$var_target) G_prevPars <- GARCH$pars else G_prevPars <- GARCH$pars[-1]
    G_prevPars <- G_prevPars[!is.na(G_prevPars)]
    G_prevData <- e/sqrt(TV$condvars)
    
    GARCH <- EstimateGARCH(G_prevData,GARCH)
    
    if (!series$var_target) currentPars <- GARCH$pars else currentPars <- GARCH$pars[-1]
    currentPars <- currentPars[!is.na(currentPars)]
    # Check if the change in param value is small enough to exit:
    if (series$var_target) {
      if (all((G_prevPars-currentPars) < series$garchparLimit[-1])) GARCH_Converged <- TRUE } 
    else {
      # var_target = FALSE  
      if (all(abs(G_prevPars-currentPars) < series$garchparLimit)) GARCH_Converged <- TRUE 
    }
    if (verbose) cat("\nGARCH value:",GARCH$value, "Converged=",GARCH_Converged )
    
    if (verbose) cat("\n garchPars-GARCH$pars:",(G_prevPars-currentPars))
    if (verbose) cat("\n garchparLimit:",series$garchparLimit)
    
    
    # parameter values before optimization
    if (!series$var_target) TV_prevPars <- TV$pars else TV_prevPars <- c(TV$delta0,TV$pars)
    TV_prevPars <- TV_prevPars[!is.na(TV_prevPars)]
    
    TV <- EstimateTV(e/sqrt(GARCH$condvars),TV)
    
    if (!series$var_target) currentPars <- TV$pars else currentPars <- c(TV$delta0,TV$pars)
    currentPars <- currentPars[!is.na(currentPars)]
    # Check if the change in param value is small enough to exit:
    if (!series$var_target) {
      if (all(abs(TV_prevPars-currentPars) < series$tvparLimit[-1])) TV_Converged <- TRUE }
    else {
      # Var_target = TRUE
        if (all(abs(TV_prevPars-currentPars) < series$tvparLimit)) TV_Converged <- TRUE
    }
    if (verbose) cat("\nTV value:",TV$value, "Converged=", TV_Converged)

    if (verbose) cat("\n tvPars-TV$pars:",(TV_prevPars-currentPars))
    if (verbose) cat("\n tvparLimit:",series$tvparLimit)
    
    g <- TV$condvars
    h <- GARCH$condvars
    series$value <- sum(ll_const - 0.5*log(g) - 0.5*log(h) - 0.5*e2/(g*h))
    if (verbose) cat("\nSeries value:",series$value," in Loop:",loopcount)
    if (verbose) cat("\n")
    
    #Exit time:
    if ((loopcount==maxLoops) || (TV_Converged && GARCH_Converged)) {
     # Finally - Calculate the hessian matrix & std Errors
      
      GARCH$pars[!is.na(GARCH$pars)] <- G_prevPars
      GARCH <- EstimateGARCH(G_prevData,GARCH,calcHess = TRUE)
      
      if (!series$var_target) TV$pars[!is.na(TV$pars)] <- TV_prevPars else {
        TV$delta0 <- TV_prevPars[1]
        TV$pars[!is.na(TV$pars)] <- tail(TV_prevPars,-1)
      } 
      TV <- EstimateTV(e/sqrt(GARCH$condvars),TV,calcHess = TRUE)
      
      if (loopcount==maxLoops) {
        cat("\n\nMaximum Loops reached.  Convergence state is:")
        cat("\nTV_Converged:",TV_Converged,"\nGARCH_Converged:",GARCH_Converged)
      } else cat("\n\nTV & GARCH Convergence achieved in:",loopcount,"loops")
      
      break 
    }
    
    #Progress Bar:
    if (!verbose) cat(".")
    
    loopcount <- loopcount + 1
      
  } # End repeat..loop
  
  # Update the series with the optimised values
  series$tv <- TV
  series$garch <- GARCH
  
  #Return:
  series
  
} # End: EstimateTVGARCH()

EstimateSTCC <- function(stcc,z,mycontrol,calcHess=FALSE) {
  STCC <- stcc
  
  ###
  ### ---  Call optim to calculate the estimate --- ###
  ###
  

  if(STCC$shape==2) numTRpars <- 3 else numTRpars <- 2
  optimpars <- c(STCC$P1,STCC$P2,STCC$TRpars)
  
  tmp <- NULL
  #myLogLik.stcc <- function(pars,z,stcc,return_ll=TRUE)
  try(tmp <- optim(optimpars,myLogLik.stcc,z,STCC,gr=NULL,method="BFGS",control=mycontrol,hessian=calcHess))
  
  ### ---  Interpret the response from optim --- ###
  # An unhandled error could result in a NULL being returned by optim()
  if (is.null(tmp)) {
    STCC$value <- err_output 
    return(STCC)
  }
  
  #Optim converged successfully => we expect tmp$par to have good estimates!
  if (tmp$convergence==0) {
   
    if (!is.null(STCC$P1)) STCC$EST_P1 <- tmp$par[1:length(STCC$P1)]
    if (!is.null(STCC$P2)) STCC$EST_P2 <- tmp$par[(length(STCC$P1)+1):(2*length(STCC$P2))]
    if (!is.null(STCC$TRpars)) STCC$EST_TRpars <- tail(tmp$par,numTRpars)
    
    STCC$value <- tmp$value
    STCC$condcorrs <- myLogLik.stcc(tmp$par,z,STCC,return_ll=FALSE)
    if (calcHess) {
      STCC$hess <- tmp$hessian
      try (STCC$stderr <- sqrt(-diag(solve(tmp$hessian))))
      if (is.null(STCC$stderr)) warning("Failed to calculate STCC$stderr - probably a numerical solve() problem")
    }
  } else { 
    #Failed to converge
    STCC$value <- err_output
  }
  # For troubleshooting - attach the optim() ouput
  STCC$tmp <- tmp
  #Return:
  STCC

}  #End: EstimateSTCC()

EstimateSTEC <- function(stec,z,mycontrol,calcHess=FALSE) {
  STEC <- stec
  
  ###
  ### ---  Call optim to calculate the estimate --- ###
  ###
  
  optimpars <- STEC$pars
  optimpars <- optimpars[!is.na(optimpars)]
  
  #myLogLik.stec <- function(pars,z,stec,return_ll=TRUE){
  try(tmp <- optim(optimpars,myLogLik.stec,z,STEC,gr=NULL,method="BFGS",control=mycontrol,hessian=calcHess))
  
  ###    
  ### ---  Interpret the response from optim --- ###
  ###
  
  # An unhandled error could result in a NULL being returned by optim()
  if (is.null(tmp)) {
    STEC$value <- -111333  #Set the value = the Erroreturn_ll=FALSErOutput value used by the myLogLik.stcc() function
    return(STEC)
  }
  
  if (tmp$convergence==0) {
    #Optim converged successfully => we expect tmp$par to have good estimates!
    
    STEC$pars[!is.na(STEC$pars)] <- tmp$par
    STEC$value <- tmp$value
    STEC$condcorrs <- myLogLik.stec(tmp$par,z,STEC,return_ll=FALSE)
    if (calcHess) {
      STEC$hess <- tmp$hessian
      try(STEC$stderr <- sqrt(-diag(solve(tmp$hessian))))
      if (is.null(STEC$stderr)) warning("Failed to calculate STEC$stderr - probably a numerical solve() problem")
    }
  } else { 
    #Failed to converge
    STEC$value <- -111000-tmp$convergence
  }
  rm(tmp)
  
  #Return:
  STEC
}  #End: EstimateSTEC()

} #End EstimateFn
#==========================================================#

zCalcFn <- TRUE
if (zCalcFn) {
  
CalculateTestStat_TV <- function(e,tv,testorder,usetest) {
  # Inputs:  
  #         tv          - TV list object
  #         e           - vector of data
  #         testorder   - Character string containing "H01", "H02", "H03"
  #         useTest     - Character string = "TR2" or "ROBUST" (Not Both!!)
  # Outputs:
  #         TestStat    - Numeric value.  Note: This function can only return one value at a time!!
  #
  
  TestStat <- NA

  if (usetest=="TR2") TestStat <- myTest.TV.noGARCH.TR2(e,tv,testorder)
  if (usetest=="ROBUST") TestStat <- myTest.TV.noGARCH.robust(e,tv,testorder)

  #Return:
  TestStat
  
}  #End CalculateTestStat_TV()


CalcProbabilityDist <- function(tv,refdata_withgarch,reftest,usetest="TR2",testorder=NULL,
                                saveas=NULL,numcores=1,numLoops=10,useRcpp=FALSE) {
  
  ###==================================================================================###
  ###                             SETUP FOR B-LOOP
  ###==================================================================================###
  
  ## Setup the parallel backend environment ##
  Sys.setenv("MC_CORES" = numcores)
  cl <- makeCluster(numcores)
  registerDoParallel(cl, cores = numcores)

  TV <- tv          # Initialise the local TV object
  B <- numLoops
  
  LMRef <- reftest  
  
  # Which Tests do we want to run?
  # Note: The usetest parameter must be in uppercase!  Don't want to write code to force it! 'toupper(usetest)'
  runTR2 <- any(grepl("TR2",usetest, useBytes = TRUE))
  runRobust <- any(grepl("ROBUST",usetest, useBytes = TRUE))
  
  ## Setup the matrix to store the simulation results - depends on the Order of TV function
  numTestResults <- as.integer(runTR2)*3 + as.integer(runRobust)*3   # ([TestStat,Pval,RefStat] x 2 Tests) 
  LogLiklihood_value <- 1
  if (is.null(tv$pars)) numpars <- 1 else numpars <- 1 + length(tv$pars)
  if (!is.null(tv$linpars)) numpars <- numpars + length(tv$linpars)
  numResultCols <- 1 + numTestResults + LogLiklihood_value + numpars  # The first column is used for the b-loop index 'b'
  runSim <- matrix(NA,nrow=B,ncol=numResultCols)
  
  ###==================================================================================###
  ###                             START B-LOOP
  ###==================================================================================###
  runSim <- foreach(b = 1:B, .export = c("Test_TV_noGARCH"),  .inorder=FALSE, .combine=rbind, .verbose = FALSE) %dopar% {
    #for (b in 1:B) {
    functionsPath <- file.path(dirname(getwd()),"Functions")
    functionsFile <- file.path(functionsPath,"functions_tvgjr_v4.r")
    source(functionsFile, local=TRUE)

    # cppfunctionsFile <- file.path(functionsPath,"cFunctions_tvgjr.cpp")
    # Rcpp::sourceCpp(cppfunctionsFile, env=environment())
    
    runSimrow <- c(b,rep(NA,numTestResults))  
    
    sim_e <- as.vector(refdata_withgarch[,b])
    
    TV <- EstimateTV(sim_e,tv)  # Note: The tv must always be the same params that are passed in!  Only the sim_e changes!!
    if (!TV$error) {
      # Now run the requested Tests:
      if(useRcpp) {
        
        simTEST <- Test_TV_noGARCH(e, TV$Tobs,TV$delta0,TV$pars,TV$linpars,TV$shape,TV$speedoption,TV$st, testorder) 
        runSimrow[2:4] <- c(simTEST["LMTR2"],as.integer(simTEST["LMTR2"] > LMRef$LMTR2),LMRef$LMTR2) 
        runSimrow[5:7] <- c(simTEST["LMRobust"],as.integer(simTEST["LMRobust"] > LMRef$LMRobust),LMRef$LMRobust)
        
      } else {
        
        if (runTR2) simTEST <- CalculateTestStat_TV(sim_e,TV,testorder,"TR2")        # else simTEST <- NA from initialisation above
        runSimrow[2:4] <- c(simTEST,as.integer(simTEST > LMRef$LMTR2),LMRef$LMTR2) 
        if (runRobust) simTEST <- CalculateTestStat_TV(sim_e,TV,testorder,"ROBUST")  # else simTEST <- NA from initialisation above
        runSimrow[5:7] <- c(simTEST,as.integer(simTEST > LMRef$LMRobust),LMRef$LMRobust)
      } 
        
    } else {
      
      # There was an error estimating TV
      TV$delta0 <- tv$delta0
      TV$pars <- tv$pars
      TV$linpars <- tv$linpars
      TV$value <- NA
    }  


    # Return when using foreach(...):
    result <- c(runSimrow,TV$value,TV$delta0) 
    if (!is.null(TV$pars)) result <- c(result,TV$pars)
    if (!is.null(TV$linpars)) result <- c(result,TV$linpars)
    #Return:
    result
    
    # Return when using for(...)  -  Use this format for troubleshooting!
    # runSim[b,] <- c(runSimrow,TV$value,TV$delta0)
    # if (!is.null(TV$pars)) runSim[b,] <- c(runSim[b,],TV$pars)
    # if (!is.null(TV$linpars)) runSim[b,] <- c(runSim[b,],TV$linpars)
    
  } # End: foreach(b = 1:B,...
  
    # Matrix structured as follows:
    loopDataNames <- c("b")
    if (runTR2) loopDataNames <- c(loopDataNames,"Stat_TR2","Pval_TR2","Ref$LMTR2")
    if (runRobust) loopDataNames <- c(loopDataNames,"Stat_Robust","Pval_Robust","Ref$LMRobust")
    loopDataNames <- c(loopDataNames,"Estimated_LL")
    TVparNames <- c("TV$delta0")
    if(!is.null(TV$pars)) TVparNames <- c(TVparNames,paste0("TV$par",as.character(seq(1,length(TV$pars)))))
    if(!is.null(TV$linpars)) TVparNames <- c(TVparNames,paste0("TV$linpar",as.character(seq(1,length(TV$linpars)))))
    loopDataNames <- c(loopDataNames,TVparNames)
    try(colnames(runSim) <- loopDataNames)
  
  if (!is.null(saveas)) saveRDS(runSim,saveas)
  
  unregisterDoParallel(cl)
  rm(cl,refdata_withgarch,TV)
  
  #Return:
  runSim
  
} # End of CalcProbabilityDist()

CalcProbabilityDist_c <- cmpfun(CalcProbabilityDist)

CalcProbabilityDist_LM <- function(tv,refdata_withgarch,reftest,usetest=NULL,testorder=NULL,saveas=NULL,numcores=1,numLoops=10)
{
  ###==================================================================================###
  ###                             SETUP FOR B-LOOP
  ###==================================================================================###
  
  Sys.setenv("MC_CORES" = numcores)
  cl <- makeCluster(numcores)
  registerDoParallel(cl, cores = numcores)
  
  TV <- tv          # Initialise the local TV object
  B <- numLoops
  
  LMRef <- reftest  #Note: LMRef$... must contain values for all Tests, H01, H02 & H03 
  ## Setup the matrix to store the simulation results - depends on the Order of TV function
  numTestResults <- 13  # b_idx + (2 Results x 2 Tests x 3 testOrders)
  numRefTestAndLL <- 7  # (TR2 RefStat, Robust RefStat) x 3, LogLiklihood_value
  if (is.null(tv$pars)) numpars <- 1 else numpars <- 1 + length(tv$pars)
  numResultCols <- numTestResults+numRefTestAndLL+numpars
  runSim <- matrix(NA,nrow=B,ncol=numResultCols)
  
  ###==================================================================================###
  ###                             START B-LOOP
  ###==================================================================================###
  runSim <- foreach(b = 1:B, .inorder=FALSE, .combine=rbind, .verbose = FALSE) %dopar% {
    #for (b in 1:B) {
    functionsPath <- file.path(dirname(getwd()),"Functions")
    functionsFile <- file.path(functionsPath,"functions_tvgjr_v4.r")
    source(functionsFile, local=TRUE)
    
    runSimrow <- c(b,rep(NA,numTestResults-1))
    
    sim_e <- as.vector(refdata_withgarch[,b])
    
    TV <- EstimateTV(sim_e,tv)  # Note: The tv must always be the same params that are passed in!  Only the sim_e changes!!
    if (!TV$error) {
      # Now run the requested Tests:  (Note: set testorder=Null to run all tests at once)
      #H01
      if (is.null(testorder) || testorder=="H01"){
        simTEST <- CalculateTestStat_TV(sim_e,TV,"H01","TR2")
        if (!is.na(simTEST)) runSimrow[2:3] <- c(simTEST,as.integer(simTEST > LMRef$LMTR2.H01)) 
        #
        simTEST <- CalculateTestStat_TV(sim_e,TV,"H01","ROBUST")
        if (!is.na(simTEST)) runSimrow[4:5] <- c(simTEST,as.integer(simTEST > LMRef$LMRobust.H01))
      }
      #H02
      if (is.null(testorder) || testorder=="H02"){
        simTEST <- CalculateTestStat_TV(sim_e,TV,"H02","TR2")
        if (!is.na(simTEST)) runSimrow[6:7] <- c(simTEST,as.integer(simTEST > LMRef$LMTR2.H02))
        #
        simTEST <- CalculateTestStat_TV(sim_e,TV,"H02","ROBUST")
        if (!is.na(simTEST)) runSimrow[8:9] <- c(simTEST,as.integer(simTEST > LMRef$LMRobust.H02))
      }
      #H03
      if (is.null(testorder) || testorder=="H03"){
        simTEST <- CalculateTestStat_TV(sim_e,TV,"H03","TR2")
        if (!is.na(simTEST)) runSimrow[10:11] <- c(simTEST,as.integer(simTEST > LMRef$LMTR2.H03))
        #
        simTEST <- CalculateTestStat_TV(sim_e,TV,"H03","ROBUST")
        if (!is.na(simTEST)) runSimrow[12:13] <- c(simTEST,as.integer(simTEST > LMRef$LMRobust.H03))
      }
      rm(sim_e,simTEST)
    }
    # Return when using foreach(...):
    result <- c(runSimrow,LMRef$LMTR2.H01,LMRef$LMTR2.H02,LMRef$LMTR2.H03,LMRef$LMRobust.H01,LMRef$LMRobust.H02,LMRef$LMRobust.H03,TV$value,TV$delta0,TV$pars)
    
    # Return when using for(...)  -  Use this format for troubleshooting!
    #runSim[b,] <- c(runSimrow,LMRef$LMTR2.H01,LMRef$LMTR2.H02,LMRef$LMTR2.H03,LMRef$LMRobust.H01,LMRef$LMRobust.H02,LMRef$LMRobust.H03,TV$value,TV$delta0,TV$pars)
    
  } # end of runSim <- foreach(b = 1:B,...
  
  
  if (!is.null(saveas))
  {
    # Matrix structured as follows:
    loopDataNames <- c("b")
    loopDataNames <- c(loopDataNames,"Stat_TR2.H01","Pval_TR2.H01","Stat_Robust.H01","Pval_Robust.H01")
    loopDataNames <- c(loopDataNames,"Stat_TR2.H02","Pval_TR2.H02","Stat_Robust.H02","Pval_Robust.H02")
    loopDataNames <- c(loopDataNames,"Stat_TR2.H03","Pval_TR2.H03","Stat_Robust.H03","Pval_Robust.H03")
    loopDataNames <- c(loopDataNames,"Ref$LMTR2.H01","Ref$LMTR2.H02","Ref$LMTR2.H03","Ref$LMRobust.H01","Ref$LMRobust.H02","Ref$LMRobust.H03","Estimated_LL")
    TVparNames <- c("TV$delta0")
    if(!is.null(TV$pars)) TVparNames <- c(TVparNames,paste0("TV$par",as.character(seq(1,length(TV$pars)))))
    loopDataNames <- c(loopDataNames,TVparNames)
    colnames(runSim) <- loopDataNames
    saveRDS(runSim,saveas)
  }
  
  unregisterDoParallel(cl)
  rm(refdata_withgarch,TV)
  
  #Return:
  runSim
  
} # End of CalcProbabilityDist_LM()



} #End: CalcFn
#==========================================================#

zDataGenFn <- TRUE
if(zDataGenFn) {
 
  
GenerateRefData_WithGarch <- function(garch,seedstart=1,numseries=1000,tsamples=5000,
                                      saveas="RefData_WithGarch.RDS") {
  
  ###
  #   Note: R uses column-major as the default order for matrices
  #   Therefore, it makes programming/life easier if we build our reference data this way!
  ###
  
  par_o <- garch$pars[1]
  if (garch$type==3) par_a <- 0 else par_a <- garch$pars[2]
  par_b <- garch$pars[3]
  if (garch$type==2) par_d <- garch$pars[4] else par_d <- 0
  
  refData <- matrix(NA,tsamples,numseries)
  
  for (b in seq(1,numseries))
  {
    set.seed(seedstart+b)
    z <- rnorm(tsamples)
    discard <- 2000    #Discard the first 1000 iterations of h(t) to get a good starting value
    z_discard <- rnorm(discard)
    # recursive ht and et
    sim_e_garch <- rep(0,tsamples)
    # Set the first data points, then loop through the remainder
    ht_1 <- ht <- 1
    et_1 <- et <- sqrt(ht)*z_discard[1]
    for (t in 2:discard) {
      ht_1 <- ht <- par_o + par_a*(et_1)^2 + par_b*ht_1 + par_d*(min(et_1,0)^2)
      et_1 <- et <- sqrt(ht)*z_discard[t]
    }
    #
    ht <- par_o + par_a*(et_1)^2 + par_b*ht_1 + par_d*(min(et_1,0)^2)
    sim_e_garch[1] <- sqrt(ht)*z[1]
    for (t in 2:tsamples) {
      ht_1 <- ht <- par_o + par_a*(sim_e_garch[t-1])^2 + par_b*ht_1 + par_d*(min(sim_e_garch[t-1],0)^2)
      sim_e_garch[t] <- sqrt(ht)*z[t]
    }
    refData[,b] <- sim_e_garch
  } # End of for (b in seq(1,numloops))
  
  saveRDS(refData,saveas)
  #Return:
  rm(refData)
  
} # End: GenerateRefData_WithGarch()
  
  
GenerateRefData_WithGarch_c <- cmpfun(GenerateRefData_WithGarch)

GenerateRefData_WithTVGarch <- function(garch,tv,seedstart=0,numloops=5000,tsamples=5000,saveas="GeneratedRefData_WithTVGarch.RDS")
{
  
  ###
  #   Note: R uses column-major as the default order for matrices
  #   Therefore, it makes programming/life easier if we build our reference data this way!
  ###

  if (garch$type==1){
    par_o <- garch$pars[1]
    par_a <- garch$pars[2]
    par_b <- garch$pars[3]
    par_d <- 0
  }else if (garch$type==2){
    par_o <- garch$pars[1]
    par_a <- garch$pars[2]
    par_b <- garch$pars[3]
    par_d <- garch$pars[4]
  }else if (garch$type==3){
    par_o <- calculateOmega(garch$pars[2],garch$pars[3],garch$pars[4])
    par_a <- garch$pars[2]
    par_b <- garch$pars[3]
    par_d <- garch$pars[4]
  }
  
  TV <- tv
  refData <- matrix(NA,tsamples,numloops)
  
  #z    # this will store iid data
  #e    # this will store data with GARCH
  #eps  # this will store correlated data with GARCH and TV
  
  for (b in seq(1,numloops))
  {
    set.seed(seedstart+b)
    z <- rnorm(tsamples)
    if (is.null(GARCH)) e<-z
    else{
      discard <- round(tsamples*0.25)    #Discard the first 25% of the data set
      z_discard <- rnorm(discard)
      # recursive ht and et
      e <- rep(0,tsamples)
      # Set the first data points, then loop through the remainder
      ht_1 <- ht <- 1
      et_1 <- et <- sqrt(ht)*z_discard[1]
      for (t in 2:discard) {
        ht_1 <- ht <- par_o + par_a*(et_1)^2 + par_b*ht_1 + par_d*(min(et_1,0)^2)
        et_1 <- et <- sqrt(ht)*z_discard[t]
      }
      #
      ht <- par_o + par_a*(et_1)^2 + par_b*ht_1 + par_d*(min(et_1,0)^2)
      e[1] <- sqrt(ht)*z[1]
      for (t in 2:tsamples) {
        ht_1 <- ht <- par_o + par_a*(e[t-1])^2 + par_b*ht_1 + par_d*(min(e[t-1],0)^2)
        e[t] <- sqrt(ht)*z[t]
      }
    } # GARCH data created => e complete
    
    if (is.null(TV)) eps <- e
    else{
      # - - - NO TV at all (no delta0)- - -
      if (is.null(TV$delta0)){
        gt <- 1
        eps <- e
      } else if (is.null(TV$pars)) {
        # - - - NO TV (only constant delta0)- - -
        gt <- TV$delta0
        eps <- e*sqrt(gt)
      } else {
        # - - - Full TV - - -
        gt <- calculate_g(TV)
        eps <- e*sqrt(gt)
      }
    } # TV data created -> eps complete
    
    refData[,b] <- eps
  } # End of for (b in seq(1,numloops))
  
  saveRDS(refData,saveas)
  #Return:
  refData
  
}  # End: GenerateRefData_WithTVGarch

GenerateRefData_WithTVGarch_c <- cmpfun(GenerateRefData_WithTVGarch)


GenerateRefData <- function(corr,garch,tv,e=NULL,e_discard=NULL,seedstart=0,numseries=10,tobs=1000,saveas="GenRefData.RDS")
{
  
## ----- Comments to explain how this function works ----- ##
##  This function is designed to generate a matrix of N data series, to be used as reference data with known qualities
##  Garch, TV, and Correlation information can be passed to the function as a parameter. None or all can be included.
##  If e is provided we will use it as the data - matrix of noise expected
##  If e_discard is provided we will use it as the garch discard data - matrix of noise expected
##  If e is not provided, then e & e_discard will be generated using the other parameters passed in.
  
# v   : this will store iid data
# z   : this will store correlated iid data 
# w   : this will store correlated data + GARCH 
# e   : this will store correlated data + GARCH + TV 

  
  Tobs <- tobs
  N <- numseries
  refData <- matrix(NA,Tobs,N)  # this will store the Tobs x N final generated data
  
  if (!is.null(e)) {
    v <- e
    v_discard <- e_discard
    discardObs <- NROW(v_discard)
  } else {
  
    # We need to have some data to discard, when generating Garch:
    discardObs <- 1500
    
    # Make some noise to start with!!
    set.seed(seedstart)
    dataVec <- rnorm((Tobs+discardObs)*N)
    v_vec <- dataVec[1:(Tobs*N)]
    # Use the tail of this vector for our discard:
    discard_vec <- tail(dataVec,(discardObs*N))
    
    v <- matrix(v_vec,nrow=Tobs,ncol=N)
    v_discard <- matrix(discard_vec,nrow=discardObs,ncol=N)
  }
  
  # Generate Correlated iid Data
  z <- v
  if (!is.null(corr)) {
    
      # - - - CCC - - -
      if (corr$type=="CCC"){
        CCC <- corr$CCC
        P <- CCC$P
        eig <- eigen(P)       # eigenvalues and eigenvectors
        Psqrt <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
        z <- t(Psqrt%*%t(v))
      }
      
      # - - - STCC - - -
      else if (corr$type=="STCC"){
        STCC <- corr$STCC
        rhovec1 <- as.vector(STCC$P1)
        rhovec2 <- as.vector(STCC$P2[(N*(N-1)/2+1):(N*(N-1))])
        speed <- STCC$pars[1]
        loc1 <- STCC$pars[2]
        loc2 <- STCC$pars[3]
        
        Gt <- myG(speed,loc1,loc2,STCC$st,STCC$shape,STCC$speedoption)
        P <- matrix(0,Tobs,N)
        for (t in 1:Tobs){
          P[t,] <- (1-Gt[t])*rhovec1 + Gt[t]*rhovec2
          mP <- myUnVecl(P[t,])
          eig <- eigen(mP)           # eigenvalues and eigenvectors
          mPtsqrt <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors) 
          z[t,] <- mPtsqrt%*%v[t,]
        }
      }

  } #End: Generate Correlated Data
  
  
  # Now add Garch and/or TV to the data (loop over N-series)
  for (n in 1:N) {
    
    w <- z
    # Generate GARCH in the data:
    if (!is.null(ngarch)) {
      garch <- ngarch[[n]]
      # Set Garch params based on type:
      if (garch$type==0){
        par_o <- garch$pars[1]
        par_a <- 0
        par_b <- 0
        par_d <- 0
      }else if (garch$type==1){
        par_o <- garch$pars[1]
        par_a <- garch$pars[2]
        par_b <- garch$pars[3]
        par_d <- 0
      }else if (garch$type==2){
        par_o <- garch$pars[1]
        par_a <- garch$pars[2]
        par_b <- garch$pars[3]
        par_d <- garch$pars[4]
      }else if (garch$type==3){
        par_o <- calculateOmega(garch$pars[2],garch$pars[3],garch$pars[4])
        par_a <- garch$pars[2]
        par_b <- garch$pars[3]
        par_d <- garch$pars[4]
      } else stop("Garch$type is missing or invalid!")
      
      # Initialise the first data points: e(t-1), h(t-1) using the v_discard[<discard>,N]
      ht_1 <- ht <- 1
      et_1 <- et <- v_discard[1,n]
      for (t in 2:discardObs) {
        ht_1 <- ht <- par_o + par_a*et_1*et_1 + par_b*ht_1 + par_d*(min(et_1,0)^2)
        et_1 <- et <- sqrt(ht)*v_discard[t,n]
      }
      
      # Now generate the actual Garch data:  
      ht_1 <- ht <- par_o + par_a*et_1*et_1 + par_b*ht_1 + par_d*(min(et_1,0)^2)
      w[1,n] <- sqrt(ht)*z[1,n]
      for (t in 2:Tobs) {
        ht_1 <- ht <- par_o + par_a*w[t-1,n]*w[t-1,n] + par_b*ht_1 + par_d*(min(w[t-1,n],0)^2)
        w[t,n] <- sqrt(ht)*z[t,n]
      }
    } #End: Generate Garch Data
    
    e <- w
    # Generate TV in the data:
    if (!is.null(ntv)) {
      tv <- ntv[[n]]
      # Note: calculation of 'g' depends on tv$st to get the length of data series
      tv$Tobs <- NROW(e)
      tv$st <- seq(1,tv$Tobs)/tv$Tobs
      if (is.null(tv$pars)) {
        # - - - NO TV (only constant delta0) - - -
        #g <- tv$delta0
        g <- rep(tv$delta0,Tobs)
        e[,n] <- w[,n]*sqrt(g)
      } else {
        # - - - Full TV - - -
        g <- calculate_g(tv)
        e[,n] <- w[,n]*sqrt(g)
      } 
    } #End: Generate TV Data 
    
    # Update the return matrix with this data series:
    refData[,n] <- e[,n]
    
    
  } #End: for (n in 1:N) 
  
  
  if(!is.null(saveas)) saveRDS(refData,saveas)
  #Return:
  refData
  
}  # End: GenerateRefData()

GenerateRefData_c <- cmpfun(GenerateRefData)

simPD <- function(n,tobs,corr,loops=100) {
  # This function requires the N, Tobs, & Correlation parameters to be set correctly, BEFORE calling it!
  # No params are passed & no error-checking is done.
  # This is just a dummy sub-routine to save me from 'copy-pasting' the same code too many times!
  # Setup parameters in generic format for the GenerateRefData() function: 
  
  functionsPath <- file.path(dirname(getwd()),"Functions")
  functionsFile <- file.path(functionsPath,"functions_tvgjr_v4.r")
  source(functionsFile, local=TRUE)
  
  N <- n
  Tobs <- tobs
  CORR <- corr
  simLoops <- loops
  
  corr_case <- CORR$corr_case
  
  # Parallel Loop Index = b
  b <- 0
  ccDist <- foreach(b = 1:simLoops, .inorder=FALSE, .combine=rbind, .verbose = FALSE) %dopar% {
    
    # GenerateRefData() uses N unique seeds: 
    # Therefore we need our seedstart to jump in blocks of N
    # e.g. N=2: SeedStart = c(1,3,5,7...)
    saveAs <- paste0("RefData_N",N,"_T",Tobs,"_C",corr_case,".RDS")
    seedStart <- ((b-1)*N) + 1
    refData <- GenerateRefData(corr=CORR,garch=NULL,tv=NULL,seedstart=seedStart,numseries=N,tobs=Tobs,saveas=saveAs)
    # Estimate the constant correlation:
    CCC <- list()
    CCC$P <- cor(refData)
    # Provide an STCC object for H1
    STCC <- list()
    STCC$st <- (1:Tobs)/Tobs  # Linear time-based transition
    
    # Initialise the Test Results:
    K1 <- K2 <- NaN
    K1 <- try(myTest.CCCvSTCC.LM(e=refData,H0=CCC,H1=STCC,testorder=1))
    K2 <- try(myTest.CCCvSTCC.LM(e=refData,H0=CCC,H1=STCC,testorder=2))
    
    # Return:
    c(K1,K2)
  }
  
  # Save Distribution Results to CSV:
  csvName1 <- paste0("SimulatedDist_N",N,"_T",Tobs,"_Corr1","_K1.csv")
  write.csv(ccDist[,1],csvName1,row.names = FALSE)
  csvName2 <- paste0("SimulatedDist_N",N,"_T",Tobs,"_Corr1","_K2.csv")
  write.csv(ccDist[,2],csvName2,row.names = FALSE)
  
  # Unregister the Parallel Backend Process:
  unregisterDoParallel()
  
  #Return
  TRUE
}

simPD_c <- cmpfun(simPD)

##TODO: Fix this up!:
myGenerateData <- function(Tobs,N,tv=list(),garch=list(),corr=list())
{
  # input: Tobs      -- time dimension
  #        N             -- number of series
  #        TV.pars_all    -- N-list of TV parameters: delta0, (delta1,speed,loc) for G1, for G2, etc...
  #        TV.shape       -- N-list of vectors of shape indicators for each G in each g, NA: No shape, 1: (st-c), 2: (st-c1)(st-c2), 3: (st-c)^2
  #        TV.speedoption -- N- list/vector, indicator for each series, NA:No speedoption, 1=gamma, 2=gamma/std(st), 3=exp(eta), 4=1/lambda^2
  #        TV.st          -- N-list/matrix of transition variables for each series for their TV-component
  #        GARCH.type -- N-list of indicators for each series of their GARCH specification, NA: no GARCH, 1: GARCH, 2: GJR-GARCH
  #        // CORR.type  -- N-list of indicator for correlation model, NA: No type, 0:EC, 1:CCC, 2:STEC, 3:STCC, ...
  #        //CORR.pars  -- parameters for the correlation model
  #
  #       Replace CORR with STEC, or STCC or ...  as needed
  #
  #
  
  TV <- tv
  GARCH <- garch
  CORR <- corr
  
  set.seed(1)
  wstar = matrix(rnorm(Tobs*N),nrow=Tobs,ncol=N) # iid errors, unocorrelated
  z <- wstar # this will store correlated iid data
  e <- wstar # this will store correlated data with GARCH
  eps <- wstar # this will store correlated data with GARCH and TV
  
  #-----------------------------------#
  # * * * correlation structure * * * #
  #-----------------------------------#
  
  # - NO CORRELATION (i.e. CCC with P=I)                              <- - - changed 2016_04_15
  if (is.na(CORR$type)) z <- wstar
  
  # - - - CCC - - -
  else if (CORR$type==1){
    rho <- CORR$pars
    P <- myUnVecl(rho)            # N x N
    eig <- eigen(P)               # eigenvalues and eigenvectors
    Psqrt <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors) # N x N
    z <- t(Psqrt%*%t(wstar))
  }
  
  # - - - STCC - - -
  else if (CORR$type==3){
    rhovec1 <- CORR$pars[1:(N*(N-1)/2)]
    rhovec2 <- CORR$pars[(N*(N-1)/2+1):(N*(N-1))]
    speed <- CORR$pars[(N*(N-1)+1)]
    loc1 <- CORR$pars[(N*(N-1)+2):length(CORR$pars)]
    loc2 <- NA
    Gt <- myG(speed,loc1,loc2,CORR$st,CORR$shape,CORR$speedoption)
    for (t in seq(1,Tobs)){
      Pt[t,] <- (1-Gt[t])*rhovec1 + Gt[t]*rhovec2
      mPt <- myUnVecl(Pt[t,])
      eig <- eigen(mPt)           # eigenvalues and eigenvectors
      mPtsqrt <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors) # N x N
      z[t,] <- mPtsqrt%*%wstar[t,]
    }
  }
  
  #--------------------------------#
  #   * * * GARCH structure * * *  #
  #--------------------------------#
  for (n in 1:N){
    # - - - NO GARCH - - -
    if (is.na(GARCH$type)) e[,n] <- z[,n]
    # - - - GARCH(1) or GJR(2) - - -
    else if (GARCH$type==1 || GARCH$type==2) {
      discard <- 1000
      z_discard <- matrix(rnorm(discard),nrow=discard,ncol=1)
      # Dt
      o <- GARCH$pars[1]
      a <- GARCH$pars[2]
      b <- GARCH$pars[3]
      if (GARCH$type==2) d <- GARCH$pars[4] else d <- 0
      # recursive ht and et
      for (t in 1:discard){
        if (t==1) ht <- 1
        else ht <- o + a*(et_1)^2 + b*ht_1 + d*(min(et_1,0)^2)
        et <- sqrt(ht)*z_discard[t]
        ht_1 <- ht
        et_1 <- et
      }
      for (t in 1:Tobs){
        if (t==1) ht <- o + a*(et_1)^2 + b*ht_1 + d*(min(et_1,0)^2)
        else ht <- o + a*(e[t-1,n])^2 + b*ht_1 + d*(min(e[t-1,n],0))^2
        e[t,n] <- sqrt(ht)*z[t,n]
        ht_1 <- ht
      }
    }  #End: if (GARCH$type==1 || GARCH$type==2) {
    
  }  #End: for (n in 1:N) (GARCH Structure)
  
  
  #----------------------------#
  #  * * * TV structure * * *  #
  #----------------------------#
  
  for (n in 1:N){
    # - - - NO TV at all (no delta0)- - -
    if (is.null(TV$delta0)){
      gt <- 1
      eps[,n] <- e[,n]
    }
    # - - - NO TV (only constant deta0)- - -
    else if (is.null(TV$pars) && is.null(TV$linpars)){
      gt <- TV$delta0
      eps[,n] <- e[,n]*sqrt(gt)
    }
    # - - - NO TV (only constant delta0 & linearised params)- - -
    else if (is.null(TV$pars) && length(TV$linpars)>0){
      gt <- TV$delta0
      for (i in seq_along(TV$linpars)) gt <- gt + TV$linpars[i]*(TV$st^i)
      eps[,n] <- e[,n]*sqrt(gt)
    }
    # - - - TV - - -
    else {
      gt <- calculate_g(TV)
      eps[,n] <- e[,n]*sqrt(gt)
    }
    
  }
  
  #return:
  eps
  
}  #End: My generate data

##TODO: Fix this up!:
myGenSTECData <- function(Tobs,N,tv=list(),garch=list(),stec=list(),seed = 1)
{
  # input: Tobs      -- time dimension
  #        N             -- number of series
  #        TV.pars_all    -- N-list of TV parameters: delta0, (delta1,speed,loc) for G1, for G2, etc...
  #        TV.shape       -- N-list of vectors of shape indicators for each G in each g, NA: No shape, 1: (st-c), 2: (st-c1)(st-c2), 3: (st-c)^2
  #        TV.speedoption -- N- list/vector, indicator for each series, NA:No speedoption, 1=gamma, 2=gamma/std(st), 3=exp(eta), 4=1/lambda^2
  #        TV.st          -- N-list/matrix of transition variables for each series for their TV-component
  #        GARCH.type -- N-list of indicators for each series of their GARCH specification, NA: no GARCH, 1: GARCH, 2: GJR-GARCH
  #        // CORR.type  -- N-list of indicator for correlation model, NA: No type, 0:EC, 1:CCC, 2:STEC, 3:STCC, ...
  #        //CORR.pars  -- parameters for the correlation model
  #
  #       Replace CORR with STEC, or STCC or ...  as needed
  #
  #
  
  TV <- tv
  GARCH <- garch
  STEC <- stec
  
  #------------------------------------------#
  #  Handle NULL inputs:
  #------------------------------------------#
  if (is.null(STEC$type)) STEC$type <- NA
  if (is.null(GARCH$type)) GARCH$type <- NA
  if (is.null(TV$shape)) TV$shape <- NA
  #------------------------------------------#
  
  set.seed(seed)
  wstar = matrix(rnorm(Tobs*N),nrow=Tobs,ncol=N) # iid errors, unocorrelated
  z <- wstar # this will store correlated iid data
  e <- wstar # this will store correlated data with GARCH
  eps <- wstar # this will store correlated data with GARCH and TV
  
  #-----------------------------------#
  # * * * correlation structure * * * #
  #-----------------------------------#
  
  # - NO CORRELATION (i.e. CEC with P=I)                              <- - - changed 2016_04_15
  if (is.na(STEC$type)) z <- wstar
  # - - - CEC - - -
  else if (STEC$type==0){
    rho <- STEC$pars[1]
    Q <- myQ.EC(N)                # N x N eigenvectors
    L <- myL.EC(N,rho)            # N x 1 eigenvalues
    Lsqrt = diag(sqrt(L))         # N x N
    mPsqrt <- Q%*%Lsqrt%*%t(Q)    # square root of P is Q*sqrt(L)*t(Q)
    z <- t(mPsqrt%*%t(wstar))     # T x N
  }
  # - - - STEC - - -
  else if (STEC$type==1){
    rho <- STEC$pars[1]
    rhostar <- STEC$pars[2]
    speed <- STEC$pars[3]
    loc1 <- STEC$pars[4]
    loc2 <- STEC$pars[5]
    Q <- myQ.EC(N)                   # N x N eigenvectors
    Gt <- myG(speed,loc1,loc2,STEC$st,STEC$shape,STEC$speedoption) # T x 1
    rhot <- (1-Gt)*(rho-rhostar) + Gt*(rho+rhostar) # T x 1
    for (t in seq(1,Tobs)){
      Lt <- myL.EC(N,rhot[t])        # N x 1 eigenvalues
      Ltsqrt <- diag(sqrt(Lt))        # N x N
      mPtsqrt <- Q%*%Ltsqrt%*%t(Q)   # square root of Pt is Q*sqrt(Lt)*t(Q)
      z[t,] <- mPtsqrt%*%wstar[t,]
    }
  }
  
  
  #--------------------------------#
  #   * * * GARCH structure * * *  #
  #--------------------------------#
  for (n in 1:N){
    # - - - NO GARCH - - -
    if (is.na(GARCH$type)) e[,n] <- z[,n]
    # - - - GARCH(1) or GJR(2) - - -
    else if (GARCH$type==1 || GARCH$type==2) {
      discard <- 1000
      z_discard <- matrix(rnorm(discard),nrow=discard,ncol=1)
      # Dt
      o <- GARCH$pars[1]
      a <- GARCH$pars[2]
      b <- GARCH$pars[3]
      if (GARCH$type==2) d <- GARCH$pars[4] else d <- 0
      # recursive ht and et
      for (t in 1:discard){
        if (t==1) ht <- 1
        else ht <- o + a*(et_1)^2 + b*ht_1 + d*(min(et_1,0)^2)
        et <- sqrt(ht)*z_discard[t]
        ht_1 <- ht
        et_1 <- et
      }
      for (t in 1:Tobs){
        if (t==1) ht <- o + a*(et_1)^2 + b*ht_1 + d*(min(et_1,0)^2)
        else ht <- o + a*(e[t-1,n])^2 + b*ht_1 + d*(min(e[t-1,n],0))^2
        e[t,n] <- sqrt(ht)*z[t,n]
        ht_1 <- ht
      }
    }  #End: if (GARCH$type==1 || GARCH$type==2) {
    
  }  #End: for (n in 1:N) (GARCH Structure)
  
  
  #----------------------------#
  #  * * * TV structure * * *  #
  #----------------------------#
  
  for (n in 1:N){
    # - - - NO TV at all (no delta0)- - -
    if (is.null(TV$delta0)){
      gt <- 1
      eps[,n] <- e[,n]
    }
    # - - - NO TV (only constant deta0)- - -
    else if (is.null(TV$pars) && is.null(TV$linpars)){
      gt <- TV$delta0
      eps[,n] <- e[,n]*sqrt(gt)
    }
    # - - - NO TV (only constant deta0 & linearised params)- - -
    else if (is.null(TV$pars) && length(TV$linpars)>0){
      gt <- TV$delta0
      for (i in seq_along(TV$linpars)) gt <- gt + TV$linpars[i]*(TV$st^i)
      eps[,n] <- e[,n]*sqrt(gt)
    }
    # - - - TV - - -
    else {
      gt <- calculate_g(TV)
      eps[,n] <- e[,n]*sqrt(gt)
    }
  }
  
  #return:
  eps
  
}  #End: My generate data

} #End: zDataGenFn
#===========================================================================#


#==============================================================================================#
#    Note: Functions below are broken!  Need to fix, including remove linearised param         #
#==============================================================================================#
zGridSearchFn <- TRUE
if(zGridSearchFn) {

myGridSearchTV_StartPars <- function(tv=list(),e,parBoundary,gridsize=10,linearised=FALSE)
{
  #Inputs:
  #tv = A fully populated TV list object, tv$pars[] will be put into a grid & estimated
  #     Note: A maximum of 7 parameters is supported - restriction imposed by nested For..Loop structrure
  #e  = dataset we want to use for param search
  #parBoundary = list containing upper & lower bounds for each parameter
  #gridsize = Dimension along 1-side of grid, e.g. gridsize = 10 => 10x10 grid for 2 pars, i.e. 100 estimations
  #Outputs:
  #myGridSearchTV_StartPars = matrix[(No. of estimations), TV$pars] => Returns the starting parameters in a matrix for each grid point
  
  TV <- tv
  
  if (linearised) TV$pars <- c(TV$pars,TV$linpars)
  numTVpars <- length(TV$pars)
  if (numTVpars > 7) stop("Sorry mate, I can only handle 7 parameters at the moment.  :(")
  
  gridsteps <- gridsize - 1
  #Sometimes we may want to run a more granular search across 1 or more params:
  fineGridSize <- 0
  fineGrid <- (gridsize * fineGridSize) - 1
  
  fineGrid <- gridsteps
  
  #The paramater boundaries depend on the TV_Order (determined by TV$shape) & the linearised parameter
  #Set up the boundaries & step size for each loop:
  par1_by <- ((parBoundary$par1[2]-parBoundary$par1[1])/gridsteps)
  par2_by <- ((parBoundary$par2[2]-parBoundary$par2[1])/gridsteps)
  try(par3_by <- ((parBoundary$par3[2]-parBoundary$par3[1])/fineGrid))
  try(par4_by <- ((parBoundary$par4[2]-parBoundary$par4[1])/gridsteps))
  try(par5_by <- ((parBoundary$par5[2]-parBoundary$par5[1])/gridsteps))
  try(par6_by <- ((parBoundary$par6[2]-parBoundary$par6[1])/fineGrid))  #Speed of second transition
  try(par7_by <- ((parBoundary$par7[2]-parBoundary$par7[1])/gridsteps))
  
  #Define the Output matrix:
  numEstimations <- (gridsize + fineGridSize)^numTVpars
  myStartPars <- matrix(0,numEstimations,numTVpars)  #(Number-of-Estimations) x (Number-of-Params)
  estimateCount <- 0
  
  cat("\nCreating Grid Param Matrix...")
  #Note: We will always have a minimum of 2 params: d0+lin1
  #      So we only need to check if params 3-7 have been provided
  for (par1 in seq(parBoundary$par1[1],parBoundary$par1[2],par1_by))
  {
    for (par2 in seq(parBoundary$par2[1],parBoundary$par2[2],par2_by))
    {
      if (numTVpars >= 3)
      {
        for (par3 in seq(parBoundary$par3[1],parBoundary$par3[2],par3_by))
        {
          if (numTVpars >= 4)
          {
            for (par4 in seq(parBoundary$par4[1],parBoundary$par4[2],par4_by))
            {
              if (numTVpars >= 5)
              {
                for (par5 in seq(parBoundary$par5[1],parBoundary$par5[2],par5_by))
                {
                  if (numTVpars >= 6)
                  {
                    for (par6 in seq(parBoundary$par6[1],parBoundary$par6[2],par6_by))
                    {
                      if (numTVpars == 7)
                      {
                        for (par7 in seq(parBoundary$par7[1],parBoundary$par7[2],par7_by))
                        {
                          estimateCount <- estimateCount + 1
                          myStartPars[estimateCount,] <- c(par1,par2,par3,par4,par5,par6,par7)
                          cat(".")
                        }
                      } else
                      {
                        estimateCount <- estimateCount + 1
                        myStartPars[estimateCount,] <- c(par1,par2,par3,par4,par5,par6)
                        cat(".")
                      }  #End: if (numTVpars === 7) .. else
                    }  #End: for (par6 in seq...)
                    
                  } else
                  {
                    estimateCount <- estimateCount + 1
                    myStartPars[estimateCount,] <- c(par1,par2,par3,par4,par5)
                    cat(".")
                  }  #End: if (numTVpars >= 6) .. else
                }  #End: for (par5 in seq...)
              } else
              {
                estimateCount <- estimateCount + 1
                myStartPars[estimateCount,] <- c(par1,par2,par3,par4)
                cat(".")
              }  #End: if (numTVpars >= 5) .. else
            } #End: for (par4 in seq...)
          } else
          {
            estimateCount <- estimateCount + 1
            myStartPars[estimateCount,] <- c(par1,par2,par3)
            cat(".")
          }  #End: if (numTVpars >= 4) .. else
        }  #End: for (par3 in seq...)
      } else
      {
        estimateCount <- estimateCount + 1
        myStartPars[estimateCount,] <- c(par1,par2)
        cat(".")
      }  #End: if (numTVpars >= 3) .. else
    }  #End: for (par2 in seq...)
  }  #End: for (par1 in seq...)
  
  #Return:
  myStartPars
  
}  # End: myGridSearchTV_StartPars()

myGridSearchTV_FixedPars <- function(tv=list(),e,parBoundary,gridsize=10,linearised=FALSE)
{
  #Inputs:
  #tv = A fully populated TV list object, tv$pars[] will be put into a grid & estimated
  #parBoundary = list containing upper & lower bounds for each parameter
  #gridsize = Dimension along 1-side of grid, e.g. gridsize = 10 => 10x10 grid for 2 pars, i.e. 100 estimations
  #Outputs:
  #myStartPars = matrix[(No. of estimations), TV$pars] => Returns the starting parameters in a matrix for each grid point
  
  TV <- tv
  
  numTVpars <- length(TV$linpars)
  
  gridsteps <- gridsize - 1
  
  #The paramater boundaries depend on the TV_Order (determined by TV$shape) & the linearised parameter
  #Set up the boundaries & step size for each loop:
  linpar1_by <- ((parBoundary$linpar1[2]-parBoundary$linpar1[1])/gridsteps)
  linpar2_by <- ((parBoundary$linpar2[2]-parBoundary$linpar2[1])/gridsteps)
  
  #Define the Output matrix:
  numEstimations <- gridsize^numTVpars
  myStartPars <- matrix(0,numEstimations,numTVpars+length(TV$pars))  #(Number-of-Estimations) x (Number-of-Params)
  estimateCount <- 0
  
  cat("\nCreating Grid Param Matrix...")
  #Note: We are fixing the non-linear params and just 'gridding' the 2 linear params
  for (linpar1 in seq(parBoundary$linpar1[1],parBoundary$linpar1[2],linpar1_by))
  {
    
    for (linpar2 in seq(parBoundary$linpar2[1],parBoundary$linpar2[2],linpar2_by))
    {
      estimateCount <- estimateCount + 1
      myStartPars[estimateCount,] <- c(TV$pars[1],TV$pars[2],TV$pars[3],TV$pars[4],TV$pars[5],TV$pars[6],TV$pars[7],linpar1,linpar2)
      cat(".")
    }
    
  }  #End: for (linpar1 in seq...)
  
  
  #Return:
  myStartPars
  
}  # End: myGridSearchTV_FixedPars()

myGridSearchTV_Estimates <- function(tv,e,startPars,linearised=FALSE)
{
  TV <- tv
  numTVpars <- length(TV$pars)
  
  cat("\nEstimating Grid...")
  myControl = list(fnscale = -1, ndeps = rep.int(1e-4, numTVpars), maxit = 500, REPORT = 100)
  
  listItem <- vector("numeric")
  rtnForeach <- list()
  rtnForeach <- foreach(n = 1:NROW(startPars), .combine = rbind) %dopar%
  {
    TV$pars <- startPars[n,]
    tmpTV <- TV
    try(tmpTV <- EstimateTV(TV,e,myControl,linearised))
    if (is.null(tmpTV)) Optimvalue <- -111222 else Optimvalue <- tmpTV$value
    listItem <- c(startPars[n,],tmpTV$par, Optimvalue)
    rm(tmpTV)
    #Internal Return from Foreach:
    listItem
  }
  
  #Final return - list:
  rtnForeach
  
}  #End: myGridSearchTV_Estimates()

} #End: GridSearchFn


