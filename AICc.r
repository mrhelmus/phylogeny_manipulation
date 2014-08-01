##generic
AICc <- function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
  UseMethod("AICc", mod)
}

##phylolm objects
AICc.phylolm<-function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){

    if(identical(nobs, NULL)) {n <- length(mod$residuals)} else {n <- nobs}
    LL <- logLik(mod)$logLik
    K <- logLik(mod)$df  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    return(AICc)
  }

##phyloglm objects
AICc.phyloglm<-function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){

    if(identical(nobs, NULL)) {n <- length(mod$residuals)} else {n <- nobs}
    LL <- logLik(mod)$logLik
    K <- logLik(mod)$df  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    return(AICc)
  }
