######AIC functions ripped from AICcmodavg should do model averaging and AICc tables for phyloglm objects,
#but does not account for differences in the phylogeny fit (alpha/lambda values)
##generic
AICc <- function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){
  UseMethod("AICc", mod)
}

AICc.phyloglm<-function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){

    if(identical(nobs, NULL)) {n <- length(mod$residuals)} else {n <- nobs}
    LL <- logLik(mod)$logLik
    K <- logLik(mod)$df  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    return(AICc)
  }

AICc.phylolm<-function(mod, return.K = FALSE, second.ord = TRUE, nobs = NULL, ...){

    if(identical(nobs, NULL)) {n <- length(mod$residuals)} else {n <- nobs}
    LL <- logLik(mod)$logLik
    K <- logLik(mod)$df  #extract correct number of parameters included in model - this includes LM
    if(second.ord == TRUE) {AICc <- -2*LL+2*K*(n/(n-K-1))}  else{AICc <- -2*LL+2*K}
    if(return.K == TRUE) AICc[1] <- K #attributes the first element of AICc to K
    return(AICc)
  }


### AIC table to compare models

aictab <- function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...) {
  ##format list according to model class
  cand.set <- formatCands(cand.set)
  UseMethod("aictab", cand.set)
}


  ##phyloglm
aictab.AICphyloglm <-
  function(cand.set, modnames = NULL, second.ord = TRUE, nobs = NULL, sort = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }


    ##add check to see whether response variable is the same for all models
    check.resp <- lapply(X = cand.set, FUN = function(b) formula(b)[2])
    if(length(unique(check.resp)) > 1) stop("\nYou must use the same response variable for all models\n")

    Results <- data.frame(Modnames = modnames)                    #assign model names to first column
    Results$K <- unlist(lapply(cand.set, AICc, return.K = TRUE, second.ord = second.ord, nobs = nobs))     #extract number of parameters
    Results$AICc <- unlist(lapply(cand.set, AICc, return.K = FALSE, second.ord = second.ord, nobs = nobs))  #extract AICc                                      #
    Results$Delta_AICc <- Results$AICc - min(Results$AICc)            #compute delta AICc
    Results$ModelLik <- exp(-0.5*Results$Delta_AICc)                #compute model likelihood required to compute Akaike weights
    Results$AICcWt <- Results$ModelLik/sum(Results$ModelLik)        #compute Akaike weights

   ##check if some models are redundant
    if(length(unique(Results$AICc)) != length(cand.set)) warning("\nCheck model structure carefully as some models may be redundant\n")

    Results$LL <- unlist(lapply(X= cand.set, FUN = function(i) logLik(i)[1]))

    ##rename correctly to AIC
    if(second.ord == FALSE) {
      colnames(Results)[1:6] <- c("Modnames", "K", "AIC", "Delta_AIC", "ModelLik", "AICWt")
    }

    if(sort)  {
      Results <- Results[rev(order(Results[, 6])),] 	  #if sort=TRUE, models are ranked based on Akaike weights
      Results$Cum.Wt <- cumsum(Results[, 6])                        #display cumulative sum of Akaike weights
    } else {Results$Cum.Wt <- NULL}


    class(Results) <- c("aictab", "data.frame")
    return(Results)
  }





##generic
modavg <- function(cand.set, parm, modnames = NULL, second.ord = TRUE,
                   nobs = NULL, uncond.se = "revised", conf.level = 0.95,
                   exclude = NULL, warn = TRUE, ...){
  cand.set <- formatCands(cand.set)
  UseMethod("modavg", cand.set)
}


##utility function to format candidate list of models to new class
##extract class from models in list and create new class
formatCands <- function(cand.set) {

  ##extract model class
  if(!is.list(cand.set)) stop("\n\'cand.set\' needs to be a list of candidate models\n")
  n.mods <- length(cand.set)
  all.mods <- lapply(cand.set, class)
  check.class <- unique(all.mods)
  out.class <- NULL
  ##for "coxph", c("coxph.null", "coxph"), c("clogit", "coxph")
  if(all(regexpr("coxph", check.class) != -1)) {
    out.class <- "coxph"
  }

  ##if NULL
  if(is.null(out.class)) {
    if(length(check.class) > 1) stop("\nFunctions do not support mixture of model classes\n")
    out.class <- unlist(check.class)
  }

  ##rename class
  mod.class.new <- c(paste("AIC", paste(out.class, collapse = "."), sep =""))

  ##add to list
  new.cand.set <- cand.set

  ##new S3 class
  class(new.cand.set) <- mod.class.new
  return(new.cand.set)
}





##utility functions used with modavg( ) to accomodate different specifications of interaction terms (e.g., A:B, B:A, A*B, B*A)
##in models of same set

####################################################
##function to reverse terms in interaction
reverse.parm <- function(parm) {

  ##check if ":" appears in term
  val <- grep(pattern = ":", x = parm)

  ##set value to NULL
  parm.alt <- NULL

  ##if ":" appears, then reverse interaction term
  if(length(val) > 0) {

    ##additional check if interaction involves more than 2 terms
    check.terms <- unlist(strsplit(x = parm, split = ":"))

    ##number of terms in interaction
    n.check.terms <- length(check.terms)

    ##issue warning if more than 2 terms are involved
    if(n.check.terms > 2) warning("\nThis function only supports two terms in an interaction:\n",
                                  "for more complex interactions, either create terms manually before analysis \n",
                                  "or double-check that models have been correctly included in model-averaging table\n")

    ##reverse order of interaction
    parm.alt.tmp <- rep(NA, n.check.terms)

    for (b in 1:n.check.terms) {
      parm.alt.tmp[b] <- check.terms[n.check.terms - b + 1]
    }

    ##paste terms together
    parm.alt <- paste(parm.alt.tmp, collapse = ":")

    return(parm.alt)
  }
}

#reverse.parm(parm = "BARE:AGE")




####################################################
##function to reverse order of exclude terms with colon or asterisk
reverse.exclude <- function(exclude) {

  ##remove all leading and trailing white space and within parm
  exclude <- lapply(exclude, FUN = function(i) gsub('[[:space:]]+', "", i))

  ##determine which terms are interactions with colons
  which.inter <- grep(pattern = ":", x = exclude)
  n.inter <- length(which.inter)

  ##list to hold reverse terms
  excl.list.alt <- list( )
  excl.list.alt2 <- list( )
  inter.star <- list( )

  ##if there are interaction terms with colons
  if (n.inter > 0) {

    ##create list for interaction
    rev.inter <- exclude[which.inter]

    ##create list to hold results
    excl.rev.list <- list( )
    for (b in 1:length(rev.inter)) {
      excl.rev.list[b] <- strsplit(x = rev.inter[b][[1]], split = ":")
    }

    ##add interaction with asterisk
    inter.star <- lapply(excl.rev.list, FUN = function(i) paste(i, collapse = " * "))

    ##additional check if interaction involves more than 2 terms
    n.check.terms <- unlist(lapply(excl.rev.list, length))

    ##issue warning if more than 2 terms are involved
    if(any(n.check.terms > 2)) warning("\nThis function only supports two terms in an interaction:\n",
                                       "for more complex interactions, either create terms manually before analysis \n",
                                       "or double-check that models have been correctly excluded in model-averaging table\n")

    ##iterate over each item in excl.rev.list
    for(k in 1:n.inter) {

      inter.id <- excl.rev.list[k][[1]]
      n.elements <- length(inter.id)

      ##reverse order of interaction
      parm.alt.tmp <- rep(NA, n.elements)

      for (b in 1:n.elements) {
        parm.alt.tmp[b] <- inter.id[n.elements - b + 1]
      }

      ##paste terms together
      excl.list.alt[k] <- paste(parm.alt.tmp, collapse = ":")
      excl.list.alt2[k] <- paste(parm.alt.tmp, collapse = " * ")

    }
  }


  ##determine which terms are interactions with asterisk
  which.inter.star <- grep(pattern = "\\*", x = exclude)
  n.inter.star <- length(which.inter.star)

  ##set lists to hold values
  inter.space <- list( )
  inter.nospace <- list( )
  ##list to hold reverse terms
  excl.list.alt.star <- list( )
  excl.list.alt.star2 <- list( )


  ##if there are interaction terms with asterisks
  if (n.inter.star > 0) {

    ##create list for interaction
    rev.inter <- exclude[which.inter.star]

    ##create vector to hold results
    excl.rev.list <- list( )
    for (b in 1:length(rev.inter)) {
      excl.rev.list[b] <- strsplit(x = rev.inter[b][[1]], split = "\\*")
    }

    ##paste interaction term with space
    inter.space <- lapply(excl.rev.list, FUN = function(i) paste(i, collapse = " * "))
    inter.nospace <- lapply(excl.rev.list, FUN = function(i) paste(i, collapse = ":"))

    ##additional check if interaction involves more than 2 terms
    n.check.terms <- unlist(lapply(excl.rev.list, length))

    ##issue warning if more than 2 terms are involved
    if(any(n.check.terms > 2)) warning("\nThis function only supports two terms in an interaction:\n",
                                       "for more complex interactions, either create terms manually before analysis \n",
                                       "or double-check that models have been correctly excluded in model-averaging table\n")


    ##iterate over each item in excl.rev.list
    for(k in 1:n.inter.star) {

      inter.id <- excl.rev.list[k][[1]]
      n.elements <- length(inter.id)

      ##reverse order of interaction
      parm.alt.tmp <- rep(NA, n.elements)

      for (b in 1:n.elements) {
        parm.alt.tmp[b] <- inter.id[n.elements - b + 1]
      }

      ##paste terms together
      excl.list.alt.star[k] <- paste(parm.alt.tmp, collapse = " * ")
      excl.list.alt.star2[k] <- paste(parm.alt.tmp, collapse = ":")

    }
  }


  ##add step to replicate each term with colon and asterisk

  ##combine into exclude
  exclude.out <- unique(c(exclude, excl.list.alt, excl.list.alt2, inter.space, inter.nospace, excl.list.alt.star, excl.list.alt.star2, inter.star))
  if(length(exclude.out) == 0) {exclude.out <- NULL}
  return(exclude.out)
}



##phyloglm
modavg.AICphyloglm <-
  function(cand.set, parm, modnames = NULL, second.ord = TRUE,
           nobs = NULL, uncond.se = "revised", conf.level = 0.95,
           exclude = NULL, warn = TRUE, ...){

    ##check if named list if modnames are not supplied
    if(is.null(modnames)) {
      if(is.null(names(cand.set))) {
        modnames <- paste("Mod", 1:length(cand.set), sep = "")
        warning("\nModel names have been supplied automatically in the table\n")
      } else {
        modnames <- names(cand.set)
      }
    }

       #MAY NEED TO EDIT THIS IF PHYLOGLM MODELS ARE GENERALIZED FROM LOGISTIC
    ##check that link function is the same for all models
#    check.link <- unlist(lapply(X = cand.set, FUN=function(i) i$link))
#    unique.link <- unique(x = check.link)
#    if(length(unique.link) > 1) stop("\nIt is not appropriate to compute a model-averaged beta estimate\n",
#                                         "from models using different link functions\n")


#####MODIFICATIONS BEGIN#######
    ##remove all leading and trailing white space and within parm
    parm <- gsub('[[:space:]]+', "", parm)

    ##reverse parm
    reversed.parm <- reverse.parm(parm)
    exclude <- reverse.exclude(exclude = exclude)
#####MODIFICATIONS END######


    ##extract model formula for each model in cand.set
    mod_formula <- lapply(cand.set, FUN = function(i) rownames(summary(i)$coefficients))

    nmods <- length(cand.set)

    ##setup matrix to indicate presence of parms in the model
    include <- matrix(NA, nrow = nmods, ncol = 1)
    ##add a check for multiple instances of same variable in given model (i.e., interactions)
    include.check <- matrix(NA, nrow = nmods, ncol = 1)

    ##iterate over each formula in mod_formula list
    for (i in 1:nmods) {
      idents <- NULL
      idents.check <- NULL
      form <- mod_formula[[i]]

######################################################################################################
######################################################################################################
###MODIFICATIONS BEGIN
      ##iterate over each element of formula[[i]] in list
      if(is.null(reversed.parm)) {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")== "-1" , 0, 1)
        }
      } else {
        for (j in 1:length(form)) {
          idents[j] <- identical(parm, form[j]) | identical(reversed.parm, form[j])
          idents.check[j] <- ifelse(attr(regexpr(parm, form[j], fixed=TRUE), "match.length")=="-1" & attr(regexpr(reversed.parm, form[j],
                                                                  fixed=TRUE), "match.length")=="-1" , 0, 1)
        }
      }
###MODIFICATIONS END
######################################################################################################
######################################################################################################


      include[i] <- ifelse(any(idents==1), 1, 0)
      include.check[i] <- ifelse(sum(idents.check)>1, "duplicates", "OK")
    }

  #####################################################
    ##exclude == NULL; warn=TRUE:  warn that duplicates occur and stop
    if(is.null(exclude) && identical(warn, TRUE)) {
      ##check for duplicates in same model
      if(any(include.check == "duplicates")) {
        stop("\nSome models include more than one instance of the parameter of interest. \n",
             "This may be due to the presence of interaction/polynomial terms, or variables\n",
             "with similar names:\n",
             "\tsee \"?modavg\" for details on variable specification and \"exclude\" argument\n")
      }

    }

    ##exclude == NULL; warn=FALSE:  compute model-averaged beta estimate from models including variable of interest,
    ##assuming that the variable is not involved in interaction or higher order polynomial (x^2, x^3, etc...),
    ##warn that models were not excluded
    if(is.null(exclude) && identical(warn, FALSE)) {
      if(any(include.check == "duplicates")) {
        warning("\nMultiple instances of parameter of interest in given model is presumably\n",
                "not due to interaction or polynomial terms - these models will not be\n",
                "excluded from the computation of model-averaged estimate\n")
      }

    }

    ##warn if exclude is neither a list nor NULL
    if(!is.null(exclude)) {
      if(!is.list(exclude)) {stop("\nItems in \"exclude\" must be specified as a list\n")}
    }


    ##if exclude is list
    if(is.list(exclude)) {

      ##determine number of elements in exclude
      nexcl <- length(exclude)

      ##check each formula for presence of exclude variable extracted with formula( )
      not.include <- lapply(cand.set, FUN = formula)

      ##set up a new list with model formula
      forms <- list()
      for (i in 1:nmods) {
        form.tmp <- as.character(not.include[[i]])[3]
        if(attr(regexpr("\\+", form.tmp), "match.length")==-1) {
          forms[i] <- form.tmp
        } else {forms[i] <- strsplit(form.tmp, split=" \\+ ")}
      }

      ##additional check to see whether some variable names include "+"
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\+", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nPlease avoid \"+\" in variable names\n")

      ##additional check to determine if intercept was removed from models
      check.forms <- unlist(lapply(forms, FUN=function(i) any(attr(regexpr("\\- 1", i), "match.length")>0)[[1]]))
      if (any(check.forms==TRUE)) stop("\nModels without intercept are not supported in this version, please use alternative parameterization\n")


      ##search within formula for variables to exclude
      mod.exclude <- matrix(NA, nrow=nmods, ncol=nexcl)

      ##iterate over each element in exclude list
      for (var in 1:nexcl) {

        ##iterate over each formula in mod_formula list
        for (i in 1:nmods) {
          idents <- NULL
          form.excl <- forms[[i]]

          ##iterate over each element of forms[[i]]
          for (j in 1:length(form.excl)) {
            idents[j] <- identical(exclude[var][[1]], gsub("(^ +)|( +$)", "", form.excl[j]))
          }
          mod.exclude[i,var] <- ifelse(any(idents==1), 1, 0)
        }

      }

      ##determine outcome across all variables to exclude
      to.exclude <- rowSums(mod.exclude)


      ##exclude models following models from model averaging
      include[which(to.exclude>=1)] <- 0


    }



    ##add a check to determine if include always == 0
    if (sum(include)==0) {stop("\nParameter not found in any of the candidate models\n") }

    new.cand.set <- cand.set[which(include==1)] #select models including a given parameter
    new.mod.name <- modnames[which(include==1)]    #update model names

    new_table <- aictab(cand.set=new.cand.set, modnames=new.mod.name,
                        second.ord=second.ord, nobs=nobs, sort=FALSE)  #recompute AIC table and associated measures
    new_table$Beta_est <- unlist(lapply(new.cand.set, FUN=function(i) coef(i)[paste(parm)])) #extract beta estimate for parm
    new_table$SE <- unlist(lapply(new.cand.set, FUN=function(i) sqrt(diag(vcov(i)))[paste(parm)]))

    ##compute model-averaged estimates, unconditional SE, and 95% CL
    if(second.ord == TRUE) {
      Modavg_beta <- sum(new_table$AICcWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICcWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICcWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }


    } else {
      Modavg_beta <- sum(new_table$AICWt*new_table$Beta_est)

      ##unconditional SE based on equation 4.9 of Burnham and Anderson 2002
      if(identical(uncond.se, "old")) {
        Uncond_SE <- sum(new_table$AICWt*sqrt(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2))
      }

      ##revised computation of unconditional SE based on equation 6.12 of Burnham and Anderson 2002; Anderson 2008, p. 111
      if(identical(uncond.se, "revised")) {
        Uncond_SE <- sqrt(sum(new_table$AICWt*(new_table$SE^2 + (new_table$Beta_est- Modavg_beta)^2)))
      }


    }

    zcrit <- qnorm(p=(1-conf.level)/2, lower.tail=FALSE)
    Lower_CL <- Modavg_beta - zcrit*Uncond_SE
    Upper_CL <- Modavg_beta + zcrit*Uncond_SE
    out.modavg <- list("Parameter"=paste(parm), "Mod.avg.table" = new_table, "Mod.avg.beta" = Modavg_beta,
                       "Uncond.SE" = Uncond_SE, "Conf.level" = conf.level, "Lower.CL" = Lower_CL,
                       "Upper.CL" = Upper_CL)

    class(out.modavg) <- c("modavg", "list")
    return(out.modavg)
  }


`dredge` <-
function(global.model, beta = FALSE, evaluate = TRUE, rank = "AICc",
		 fixed = NULL, m.max = NA, m.min = 0, subset, marg.ex = NULL,
		 trace = FALSE, varying, extra, ct.args = NULL, ...) {

	gmEnv <- parent.frame()
	gmCall <- .getCall(global.model)
	gmNobs <- nobs(global.model)

	if (is.null(gmCall)) {
		gmCall <- substitute(global.model)
		if(!is.call(gmCall)) {
			stop("need a 'global.model' with call component. Consider using ",
				if(inherits(global.model, c("gamm", "gamm4")))
					"'uGamm'" else "'updateable'")
		}
		#"For objects without a 'call' component the call to the fitting function \n",
		#" must be used directly as an argument to 'dredge'.")
		# NB: this is unlikely to happen:
		if(!is.function(eval(gmCall[[1L]], parent.frame())))
			gettext('could not find function "%s"', deparse(gmCall[[1L]],
				control = NULL), domain = "R")
	} else {
		# if 'update' method does not expand dots, we have a problem
		# with expressions like ..1, ..2 in the call.
		# So, try to replace them with respective arguments in the original call
		is.dotted <- grep("^\\.\\.", sapply(as.list(gmCall), deparse))
		if(length(is.dotted) > 0L) {
			substGmCall <- substitute(global.model)
			if(is.name(substGmCall)) {
				.cry(NA, "call to 'global.model' contains '...' arguments and cannot be updated: %s", deparse(gmCall, control = NULL))
			} else gmCall[is.dotted] <-
				substitute(global.model)[names(gmCall[is.dotted])]
		}

		## object from 'run.mark.model' has $call of 'make.mark.model' - fixing it here:
		if(inherits(global.model, "mark") && gmCall[[1L]] == "make.mark.model") {
			gmCall <- call("run.mark.model", model = gmCall, invisible = TRUE)
		}

	}

	LL <- .getLik(global.model)
	logLik <- LL$logLik
	lLName <- LL$name

	# *** Rank ***
	rank.custom <- !missing(rank)

	if(!rank.custom && lLName == "qLik") {
		rank <- "QIC"
		.cry(NA, "using 'QIC' instead of 'AICc'", warn = TRUE)
	}

	rankArgs <- list(...)
	IC <- .getRank(rank, rankArgs)
	ICName <- as.character(attr(IC, "call")[[1L]])

	if(length(tryCatch(IC(global.model), error = function(e) {
		e$call <- do.call(substitute, list(attr(IC, "call"), list(x = as.name("global.model"))))
		stop(e)
	})) != 1L) {
		.cry(NA, "result of '%s' is not of length 1", deparse(attr(IC,
			"call"), control = NULL)[1L])
	}

	allTerms <- allTerms0 <- getAllTerms(global.model, intercept = TRUE,
		data = eval(gmCall$data, envir = gmEnv))

	# Intercept(s)
	interceptLabel <- attr(allTerms, "interceptLabel")
	if(is.null(interceptLabel)) interceptLabel <- "(Intercept)"
	nInts <- sum(attr(allTerms, "intercept"))

	#XXX: use.ranef <- FALSE
	#if(use.ranef && inherits(global.model, "mer")) {
		#allTerms <- c(allTerms, paste("(", attr(allTerms0, "random.terms"), ")",
			#sep = ""))
	#}


	# Check for na.omit
	if(!(gmNA.action <- .checkNaAction(cl = gmCall, what = "'global.model'")))
		.cry(NA, attr(gmNA.action, "message"))


	if(names(gmCall)[2L] == "") gmCall <-
		match.call(gmCall, definition = eval(gmCall[[1L]], envir = parent.frame()),
				   expand.dots = TRUE)


	# TODO: other classes: model, fixed, etc...
	gmCoefNames <- fixCoefNames(names(coeffs(global.model)))

	n.vars <- length(allTerms)

	if(isTRUE(rankArgs$REML) || (isTRUE(.isREMLFit(global.model)) && is.null(rankArgs$REML)))
		.cry(NA, "comparing models fitted by REML", warn = TRUE)

	if (beta && is.null(tryCatch(beta.weights(global.model), error = function(e) NULL,
		warning = function(e) NULL))) {
		.cry(NA, "do not know how to calculate beta weights for '%s', argument 'beta' ignored",
			 class(global.model)[1L], warn = TRUE)
		beta <- FALSE
	}

	m.max <- if (missing(m.max)) (n.vars - nInts) else min(n.vars - nInts, m.max)

	# fixed variables:
	if (!is.null(fixed)) {
		if (inherits(fixed, "formula")) {
			if (fixed[[1L]] != "~" || length(fixed) != 2L)
				.cry(NA, "'fixed' should be a one-sided formula", warn = TRUE)
			fixed <- as.vector(getAllTerms(fixed))
		} else if (identical(fixed, TRUE)) {
			fixed <- as.vector(allTerms[!(allTerms %in% interceptLabel)])
		} else if (!is.character(fixed)) {
			.cry(NA, paste("'fixed' should be either a character vector with",
						   " names of variables or a one-sided formula"))
		}
		if (!all(fixed %in% allTerms)) {
			.cry(NA, "not all terms in 'fixed' exist in 'global.model'", warn = TRUE)
			fixed <- fixed[fixed %in% allTerms]
		}
	}
	fixed <- c(fixed, allTerms[allTerms %in% interceptLabel])
	n.fixed <- length(fixed)
	termsOrder <- order(allTerms %in% fixed)
	allTerms <- allTerms[termsOrder]
	gmFormulaEnv <- environment(as.formula(formula(global.model), env = gmEnv))
	# TODO: gmEnv <- gmFormulaEnv ???

	### BEGIN Manage 'varying'
	## @param:	varying
	## @value:	varying, varying.names, variants, nvariants, nvarying
	if(!missing(varying) && !is.null(varying)) {
		nvarying <- length(varying)
		varying.names <- names(varying)
		fvarying <- unlist(varying, recursive = FALSE, use.names = FALSE)
		vlen <- vapply(varying, length, 1L)
		nvariants <- prod(vlen)
		variants <- as.matrix(expand.grid(split(seq_len(sum(vlen)),
			rep(seq_len(nvarying), vlen))))

		flat.variant.Vvals <- unlist(lapply(varying, .makeListNames),
			recursive = FALSE, use.names = FALSE)

	} else {
		variants <- varying.names <- NULL
		nvariants <- 1L
		nvarying <- 0L
	}
	## END: varying

	## BEGIN Manage 'extra'
	## @param:	extra, global.model, gmFormulaEnv,
	## @value:	extra, nextra, extraNames, nullfit_
	if(!missing(extra) && length(extra) != 0L) {
		# a cumbersome way of evaluating a non-exported function in a parent frame:
		extra <- eval(as.call(list(call("get", ".get.extras", envir = call("asNamespace",
															 .packageName), inherits = FALSE),
					 substitute(extra), r2nullfit = TRUE)), parent.frame())

		#extra <- eval(call(".get.extras", substitute(extra), r2nullfit = TRUE), parent.frame())
		if(any(c("adjR^2", "R^2") %in% names(extra))) {
			nullfit_ <- null.fit(global.model, evaluate = TRUE, envir = gmFormulaEnv)
		}
		applyExtras <- function(x) unlist(lapply(extra, function(f) f(x)))
		extraResult <- applyExtras(global.model)
		if(!is.numeric(extraResult))
			.cry(NA, "function in 'extra' returned non-numeric result")

		nextra <- length(extraResult)
		extraNames <- names(extraResult)
	} else {
		nextra <- 0L
		extraNames <- character(0L)
	}
	## END: manage 'extra'

	nov <- as.integer(n.vars - n.fixed)
	ncomb <- (2L ^ nov) * nvariants

	if(nov > 31L) .cry(NA, "number of predictors (%d) exceeds allowed maximum of 31", nov)
	#if(nov > 10L) warning(gettextf("%d predictors will generate up to %.0f combinations", nov, ncomb))
	nmax <- ncomb * nvariants
	if(evaluate) {
		ret.nchunk <- 25L
		ret.ncol <- n.vars + nvarying + 3L + nextra
		ret <- matrix(NA_real_, ncol = ret.ncol, nrow = ret.nchunk)
		coefTables <- vector(ret.nchunk, mode = "list")
	} else {
		ret.nchunk <- nmax
	}

	calls <- vector(mode = "list", length = ret.nchunk)

	## BEGIN: Manage 'subset'
	## @param:	hasSubset, subset, allTerms, [interceptLabel],
	## @value:	hasSubset, subset
	if(missing(subset))  {
		hasSubset <- 1L
	} else {
		if(!tryCatch(is.language(subset) || is.matrix(subset), error = function(e) FALSE))
			subset <- substitute(subset)

		if(is.matrix(subset)) {
			dn <- dimnames(subset)
			#at <- allTerms[!(allTerms %in% interceptLabel)]
			n <- length(allTerms)
			if(is.null(dn) || any(sapply(dn, is.null))) {
				di <- dim(subset)
				if(any(di != n)) stop("unnamed 'subset' matrix does not have both dimensions",
					" equal to number of terms in 'global.model': %d", n)

				dimnames(subset) <- list(allTerms, allTerms)
			} else {
				if(!all(unique(unlist(dn)) %in% allTerms))
					warning("at least some dimnames of 'subset' matrix do not ",
					"match term names in 'global.model'")

				subset0 <- subset
				subset <- matrix(subset[
					match(allTerms, rownames(subset)),
					match(allTerms, colnames(subset))],
					dimnames = list(allTerms, allTerms),
					nrow = n, ncol = n)
				tsubset <- t(subset)
				nas <- is.na(subset)
				i <- lower.tri(subset) & is.na(subset) & !t(nas)
				ti <- t(i)
				subset[i] <- subset[ti]
				subset[ti] <- NA
			}
			if(any(!is.na(subset[!lower.tri(subset)]))) {
				warning("non-missing values exist outside the lower triangle of 'subset'")
				subset[!lower.tri(subset)] <- NA
			}
			mode(subset) <- "logical"
			hasSubset <- 2L # subset as matrix
		} else {
			if(inherits(subset, "formula")) {
				if (subset[[1L]] != "~" || length(subset) != 2L)
					stop("'subset' formula should be one-sided")
				subset <- subset[[2L]]
			}
			subset <- as.expression(subset)
			ssValidNames <- c("comb", "*nvar*")

			subsetExpr <- subset[[1L]]

			## subset X
			#gloFactorTable <- t(attr(terms(global.model), "factors")[-1L, ] != 0)
			gloFactorTable <- t(attr(terms(reformulate(allTerms0[!(allTerms0
				%in% interceptLabel)])), "factors") != 0)

			rownames(gloFactorTable) <- allTerms0[!(allTerms0 %in% interceptLabel)]


			subsetExpr <- .substFun4Fun(subsetExpr, ".", function(x, fac, at, vName) {
				if(length(x) != 2L) .cry(x, "exactly one argument needed, %d given.", length(x) - 1L)
				if(length(x[[2L]]) == 2L && x[[2L]][[1L]] == "+") {
					fun <- "all"
					sx <- as.character(x[[2L]][[2L]])
				} else {
					fun <- "any"
					sx <- as.character(x[[2L]])
				}
				#print(sx)
				dn <- dimnames(fac)
				#print(dn)
				#browser()
				if(!(sx %in% dn[[2L]])) .cry(x, "unknown variable name '%s'", sx)
				as.call(c(as.name(fun), call("[", vName, as.call(c(as.name("c"),
					match(dn[[1L]][fac[, sx]], at))))))
			}, gloFactorTable, allTerms, as.name("comb"))

			subsetExpr <- .subst4Vec(subsetExpr, allTerms, as.name("comb"))


			if(nvarying) {
				ssValidNames <- c("cVar", "comb", "*nvar*")
				subsetExpr <- .substFun4Fun(subsetExpr, "V", function(x, cVar, fn) {
					if(length(x) > 2L) .cry(x, "discarding extra arguments", warn = TRUE)
					i <- which(fn == x[[2L]])[1L]
					if(is.na(i)) .cry(x, "'%s' is not a valid name of 'varying' element",
									  as.character(x[[2L]]), warn = TRUE)
					call("[[", cVar, i)
				}, as.name("cVar"), varying.names)
				if(!all(all.vars(subsetExpr) %in% ssValidNames))
					subsetExpr <- .subst4Vec(subsetExpr, varying.names,
											 as.name("cVar"), fun = "[[")
			}
			ssVars <- all.vars(subsetExpr)
			okVars <- ssVars %in% ssValidNames
			if(!all(okVars)) stop("unrecognized names in 'subset' expression: ",
				prettyEnumStr(ssVars[!okVars]))

			ssEnv <- new.env(parent = .GlobalEnv)
			ssFunc <- setdiff(all.vars(subsetExpr, functions = TRUE), ssVars)
			if("dc" %in% ssFunc) assign("dc", .subset_dc, ssEnv)

			hasSubset <- if(any(ssVars == "cVar")) 4L else # subset as expression
				3L # subset as expression using 'varying' variables

		}
	} # END: manage 'subset'

	comb.sfx <- rep(TRUE, n.fixed)
	comb.seq <- if(nov != 0L) seq_len(nov) else 0L
	k <- 0L
	extraResult1 <- integer(0L)
	ord <- integer(ret.nchunk)

	argsOptions <- list(
		response = attr(allTerms0, "response"),
		intercept = nInts,
		interceptLabel = interceptLabel,
		random = attr(allTerms0, "random"),
		gmCall = gmCall,
		gmEnv = gmEnv,
		allTerms = allTerms0,
		gmCoefNames = gmCoefNames,
		gmDataHead = if(!is.null(gmCall$data)) {
			if(eval(call("is.data.frame", gmCall$data), gmEnv))
				eval(call("head", gmCall$data, 1L), gmEnv) else gmCall$data
			} else NULL,
		gmFormulaEnv = gmFormulaEnv
		)

	matchCoefCall <- as.call(c(alist(matchCoef, fit1, all.terms = allTerms,
		  beta = beta, allCoef = TRUE), ct.args))

	# TODO: allow for 'marg.ex' per formula in multi-formula models
	if(missing(marg.ex) || (!is.null(marg.ex) && is.na(marg.ex))) {
		newArgs <- makeArgs(global.model, allTerms, rep(TRUE, length(allTerms)),
							argsOptions)
		formulaList <- if(is.null(attr(newArgs, "formulaList"))) newArgs else
			attr(newArgs, "formulaList")

		marg.ex <- unique(unlist(lapply(sapply(formulaList, formulaMargChk,
			simplify = FALSE), attr, "marg.ex")))
		if(!length(marg.ex)) marg.ex <- NULL else
			cat("Marginality exceptions:", sQuote(marg.ex), "\n")
	}
	###

	retColIdx <- if(nvarying) -n.vars - seq_len(nvarying) else TRUE

	prevJComb <- 0L
	for(iComb in seq.int(ncomb)) {
		jComb <- ceiling(iComb / nvariants)
		if(jComb != prevJComb) {
			isok <- TRUE
			prevJComb <- jComb

			comb <- c(as.logical(intToBits(jComb - 1L)[comb.seq]), comb.sfx)
			nvar <- sum(comb) - nInts

			if(nvar > m.max || nvar < m.min ||
			   switch(hasSubset,
					FALSE,
					!all(subset[comb, comb], na.rm = TRUE),
					!.evalExprIn(subsetExpr, env = ssEnv, enclos = parent.frame(),
						comb = comb, `*nvar*` = nvar),
					FALSE
			   )) {
				isok <- FALSE
				next;
			}
			newArgs <- makeArgs(global.model, allTerms[comb], comb, argsOptions)
			formulaList <- if(is.null(attr(newArgs, "formulaList"))) newArgs else
				attr(newArgs, "formulaList")


			if(!all(vapply(formulaList, formulaMargChk, logical(1L), marg.ex)))  {
				isok <- FALSE; next;
			}
			if(!is.null(attr(newArgs, "problems"))) {
				print.warnings(structure(vector(mode = "list",
					length = length(attr(newArgs, "problems"))),
						names = attr(newArgs, "problems")))
			} # end if <problems>

			cl <- gmCall
			cl[names(newArgs)] <- newArgs
		} #  end if(jComb != prevJComb)

		if(!isok) next;
		## --- Variants ---------------------------
		clVariant <- cl
		if (nvarying) {
			cvi <- variants[(iComb - 1L) %% nvariants + 1L, ]

			if(hasSubset == 4L &&
				!.evalExprIn(subsetExpr, env = ssEnv, enclos = parent.frame(),
					comb = comb, `*nvar*` = nvar, cVar = flat.variant.Vvals[cvi]))
						next;
			clVariant[varying.names] <- fvarying[cvi]
		}

		if(trace) {
			cat(iComb, ": "); print(clVariant)
			utils::flush.console()
		}

		if(evaluate) {
			# begin row1: (clVariant, gmEnv, modelId, IC(), applyExtras(),
			#              nextra, allTerms, beta,
			#              if(nvarying) variantsIdx[v] else NULL
			fit1 <- tryCatch(eval(clVariant, gmEnv), error = function(err) {
				err$message <- paste(conditionMessage(err), "(model",
					iComb, "skipped)", collapse = "")
				class(err) <- c("simpleError", "warning", "condition")
				warning(err)
				return(NULL)
			})

			if (is.null(fit1)) next;

			if(nextra != 0L) {
				extraResult1 <- applyExtras(fit1)
				if(length(extraResult1) < nextra) {
					tmp <- rep(NA_real_, nextra)
					tmp[match(names(extraResult1), names(extraResult))] <- extraResult1
					extraResult1 <- tmp
				}
				#row1 <- c(row1, extraResult1)
			}

			#mcoef1 <- matchCoef(fit1, all.terms = allTerms, beta = beta,
			#	allCoef = TRUE)
			mcoef1 <- eval(matchCoefCall)

			ll <- logLik(fit1)
			nobs1 <- nobs(fit1)
			if(nobs1 != gmNobs) warning(gettextf(
				"number of observations in model #%d (%d) differs from that in the global model (%d)",
				iComb, nobs1, gmNobs))

			row1 <- c(mcoef1[allTerms], extraResult1,
				df = attr(ll, "df"), ll = ll, ic = IC(fit1)
			)
			## end -> row1

			k <- k + 1L # all OK, add model to table

			ret.nrow <- nrow(ret)
			if(k > ret.nrow) { # append if necesarry
				nadd <- min(ret.nchunk, nmax - ret.nrow)
				coefTables <- c(coefTables, vector(nadd, mode = "list"))
				ret <- rbind(ret, matrix(NA, ncol = ret.ncol, nrow = nadd),
					deparse.level = 0L)
					calls <- c(calls, vector("list", nadd))
				ord <- c(ord, integer(nadd))
			}
			ret[k, retColIdx] <- row1
			coefTables[[k]] <- attr(mcoef1, "coefTable")
		} else { # if !evaluate
			k <- k + 1L # all OK, add model to table
		}
		ord[k] <- iComb
		calls[[k]] <- clVariant
	} ### for (iComb ...)

	if(k == 0L) stop("result is empty")
	names(calls) <- ord
	if(!evaluate) return(calls[seq_len(k)])

	if(k < nrow(ret)) {
		i <- seq_len(k)
		ret <- ret[i, , drop = FALSE]
		ord <- ord[i]
		calls <- calls[i]
		coefTables <- coefTables[i]
	}

	if(nvarying) {
		varlev <- ord %% nvariants; varlev[varlev == 0L] <- nvariants
		ret[, n.vars + seq_len(nvarying)] <- variants[varlev, ]
	}

	ret <- as.data.frame(ret)
	row.names(ret) <- ord

	# Convert columns with presence/absence of terms to factors
	tfac <- which(!(allTerms %in% gmCoefNames))

	ret[tfac] <- lapply(ret[tfac], factor, levels = NaN, labels = "+")

	i <- seq_along(allTerms)
	v <- order(termsOrder)
	ret[, i] <- ret[, v]
	allTerms <- allTerms[v]
	colnames(ret) <- c(allTerms, varying.names, extraNames, "df", lLName, ICName)

	if(nvarying) {
		variant.names <- vapply(flat.variant.Vvals, function(x) if(is.character(x)) x else
			deparse(x, control = NULL, width.cutoff = 20L)[1L], character(1L))
		vnum <- split(seq_len(sum(vlen)), rep(seq_len(nvarying), vlen))
		names(vnum) <- varying.names
		for (i in varying.names) ret[, i] <-
			factor(ret[, i], levels = vnum[[i]], labels = variant.names[vnum[[i]]])

	}

	o <- order(ret[, ICName], decreasing = FALSE)
	ret <- ret[o, ]
	coefTables <- coefTables[o]

	ret$delta <- ret[, ICName] - min(ret[, ICName])
	ret$weight <- exp(-ret$delta / 2) / sum(exp(-ret$delta / 2))

	structure(ret,
		class = c("model.selection", "data.frame"),
		model.calls = calls[o],
		global = global.model,
		global.call = gmCall,
		terms = structure(allTerms, interceptLabel = interceptLabel),
		rank = IC,
		rank.call = attr(IC, "call"),
		beta = beta,
		call = match.call(expand.dots = TRUE),
		coefTables = coefTables,
		nobs = gmNobs,
		vCols = varying.names
	)
} ######

# Generalized Variance-Inflation Factors (Henric Nilsson and John Fox)

vif<-function(mod, ...){
	UseMethod("vif")
}


vif.phyloglm <- function(mod, ...) {
	if (any(is.na(coef(mod))))
		stop ("there are aliased coefficients in the model")
	v <- vcov(mod)
	assign <- attributes(mod$X)$assign
	if (names(coefficients(mod)[1]) == "(Intercept)") {
		v <- v[-1, -1]
		assign <- assign[-1]
	} else {warning("No intercept: vifs may not be sensible.")}
	terms <- colnames(mod$X)[colnames(mod$X)!="(Intercept)"]
	n.terms <- length(terms)
	if (n.terms < 2) stop("model contains fewer than 2 terms")
	R <- cov2cor(v)
	detR <- det(R)
	result <- matrix(0, n.terms, 3)
	rownames(result) <- terms
	colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
	for (term in 1:n.terms) {
		subs <- which(assign == term)
		result[term, 1] <- det(as.matrix(R[subs, subs])) *
			det(as.matrix(R[-subs, -subs])) / detR
		result[term, 2] <- length(subs)
	}
	if (all(result[, 2] == 1)) {result <- result[, 1]
  } else {
  result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
	return(result)
}





