###------------- Code for meta-analyzing cohorts with confounder imbalance -----------------
#
message("====================================================================")
message("                     CIMBAL v0.7 is loaded")
message("====================================================================")
message("If you use this software, please cite:")
message("Ray et al.(2022) Meta-analysis under Imbalance in Measurement ")
message("    of Confounders in Cohort Studies Using Only Summary-level Data.")
message("--------------------------------------------------------------------")
message("For updated citation, check https://github.com/RayDebashree/CIMBAL")
message("--------------------------------------------------------------------")
message("")



###---------------------Begin: Intermediate functions to implement imputed adjusted estimates---------------------

# calculate meta-analyzed-imputed beta estimate using meta-analyzed full cohort info
.beta.var.tilde.meta <- function(est.unadj.2.vec, var.unadj.2.vec, est.adj.2.vec, var.adj.2.vec, est.unadj.1.vec, var.unadj.1.vec, cov.unadj.adj){
		# est.unadj.2.vec: vector of unadjusted beta=log(OR) estimates from cohorts with full info
		# var.unadj.2.vec: vector of unadjusted variance estimates from cohorts with full info
		# est.adj.2.vec  : vector of adjusted beta=log(OR) estimates from cohorts with full info
		# var.adj.2.vec	 : vector of adjusted variance estimates from cohorts with full info
		# est.unadj.1.vec: vector of unadjusted beta=log(OR) estimates from cohorts with incomplete info
		# var.unadj.1.vec: vector of unadjusted variance estimates from cohorts with incomplete info
	M <- length(est.unadj.2.vec)
	if(length(var.unadj.2.vec)!=M | length(est.adj.2.vec)!=M | length(var.adj.2.vec)!=M) 
		stop("Vectors of adjusted and unadjusted estimates and their variances from cohorts with complete confounder info should be of same length (same as the no. of complete cohorts)")
	if(sum(is.na(est.unadj.2.vec))>0 | sum(is.na(var.unadj.2.vec))>0 | sum(is.na(est.adj.2.vec))>0 | sum(is.na(var.adj.2.vec))>0 | sum(is.na(est.unadj.1.vec))>0 | sum(is.na(var.unadj.1.vec))>0)
		stop('There seems to be missing information for log odds estimates or their variances from cohorts needed for the correction approach.')
	# meta-analyzed info for the cohorts with full info
	unadj.meta.2 <- meta.fixed.invvar(beta.vec=est.unadj.2.vec, var.vec=var.unadj.2.vec)
		est.unadj.meta.2 <- unadj.meta.2$beta.meta
		var.unadj.meta.2 <- unadj.meta.2$var.meta
	adj.meta.2 <- meta.fixed.invvar(beta.vec=est.adj.2.vec, var.vec=var.adj.2.vec)
		est.adj.meta.2 <- adj.meta.2$beta.meta
		var.adj.meta.2 <- adj.meta.2$var.meta
	# meta-analyzed info for the cohorts with incomplete info
	unadj.meta.1 <- meta.fixed.invvar(beta.vec=est.unadj.1.vec, var.vec=var.unadj.1.vec)
		est.unadj.meta.1 <- unadj.meta.1$beta.meta
		var.unadj.meta.1 <- unadj.meta.1$var.meta
	# imputed meta-analyzed estimates
	est.tilde.meta <- est.adj.meta.2 - est.unadj.meta.2 + est.unadj.meta.1
	var.tilde.meta <- var.adj.meta.2 + var.unadj.meta.2 + var.unadj.meta.1 - 2*cov.unadj.adj
	if(var.tilde.meta<=0){
		warning("imputed adjusted meta-analyzed variance estimate for cohorts with incomplete confounder info came out non-positive - returning NA.")
		var.tilde.meta <- NA
	}
	return(list(est.tilde.meta.1=est.tilde.meta, var.tilde.meta.1=var.tilde.meta, cov.unadj.adj=cov.unadj.adj))
}

# function for different checks needed to implement cimbal function
.checks <- function(dat){
    # check colnames
    if(sum(colnames(dat) %in% c("cohort", "samplesize", "b.unadj", "se.unadj", "b.adj", "se.adj"))!=ncol(dat))
        stop("Please make sure that the input matrix or dataframe has 6 columns with names 'cohort', 'samplesize', 'b.unadj', 'se.unadj', 'b.adj', 'se.adj'. Cohorts without b.adj or se.adj should report NA. All cohorts with non-missing b.unadj, se.unadj, b.adj and se.adj will be used to obtain imputed adjusted estimates for those without complete information.")
    # check unique cohort names
    if(length(unique(dat$cohort))!=nrow(dat))
        stop("Please ensure that cohort names/numbers are unique and are not duplicated.")
    #
    if(sum(is.na(dat$b.unadj))>0 | sum(is.na(dat$se.unadj))>0){
        remove.rows <- which(is.na(dat$b.unadj) | is.na(dat$se.unadj))
        dat <- dat[-remove.rows,]
        warning(paste0("All cohorts must have the unadjusted estimates and SEs. Some cohorts seem to have missing information about unadjusted estimates. Removed ",length(remove.rows)," such cohort(s)..."))
    }
    K <- nrow(dat)
    if(K<=1)
        stop("Need estimates from at least 2 cohorts to continue.")
    if(sum(is.na(dat$b.adj))==K | sum(is.na(dat$se.adj))==K){
        message("Current summary data from cohorts:"); print(dat)
        stop("There is no cohort with full confounder info that can be used to obtain imputed adjusted estimates for those without complete information.")
    }
    if(sum(is.na(dat$b.adj))==0 | sum(is.na(dat$se.adj))==0){
        message("Current summary data from cohorts:"); print(dat)
        stop("There is no cohort with incomplete information that needs correction using CIMBAL. Use function meta.fixed.invvar() for meta-analysis.")
    }
    return(dat)
}

# calculate covariance of b.uandj and b.adj from estimated correlation from full cohorts
.cov.unadj.adj <- function(est.unadj.2.vec, est.adj.2.vec, var.unadj.2.vec, var.adj.2.vec, var.meta.unadj.2, var.meta.adj.2, ncohort.thresh, mute.msgs, limitat0=TRUE){
	# est.adj.2.vec: vector of adjusted beta estimates corresponding to M cohorts with complete confounder info
	# est.unadj.2.vec: vector of unadjusted beta estimates corresponding to M cohorts with complete confounder info
	M <- length(est.unadj.2.vec)
	covar <- 0
	# recommended no. of cohorts to use for calculating this covariance
	ncohort.thresh.rec <- 25	
	if(ncohort.thresh < ncohort.thresh.rec) 
		if(!mute.msgs) warning(paste0("Parameter 'ncohort.thresh'=",ncohort.thresh," was specified. Note we recommend using at least ",ncohort.thresh.rec," cohorts to estimate the covariance between unadjusted and adjusted estimates."))
		
	if(M >= ncohort.thresh){
	    covar <- cor(est.adj.2.vec,est.unadj.2.vec) * var.meta.unadj.2*var.meta.adj.2 * sum(sqrt(1/var.unadj.2.vec * 1/var.adj.2.vec))
	}else{ 
		if(!mute.msgs) message("Not enough cohorts to estimate covariance between b.adj and b.unadj. Theoretically (for linear models) and empirically (for logistic model), we found this covariance should be ≥0. So, 0 is assigned as covariance.")
	}  
	
	if(covar < 0 & limitat0){
		warning(paste0("Estimated covariance between b.adj and b.unadj using ",M," complete cohorts is <0 that may result in underestimated variance of the final meta-analyzed beta estimate and possibly inflated hypothesis tests. Theoretically (for linear models) and empirically (for logistic model), we found this covariance should be ≥0. So, 0 is assigned as covariance."))
		covar <- 0
	}

	return(list(covar=covar, corr=cor(est.adj.2.vec,est.unadj.2.vec)))
}


# calculate meta-analysis of beta-estimates and their variance estimates which includes imputed estimates using CIMBAL
# (note: here we need to account for covariances between imputed and existing fully adjusted estimates)
# (note: inverse variances are no longer optimal weights due to dependence between beta estimates)
.meta.fixed.invwt.cimbal <- function(dat, dat.corr, fullconf, noconf, cov.unadj.adj, ncohort.thresh, mute.msgs, returnSE=FALSE){
	# meta-analyzed estimates from the two categories of cohorts
        est.adj.2.vec <- dat$b.adj[fullconf]
        var.adj.2.vec <- (dat$se.adj[fullconf])^2
	meta.adj.2 <- unlist(meta.fixed.invvar(est.adj.2.vec, var.adj.2.vec), use.names=F)
	meta.corradj.1 <- c(dat.corr$b.adj[nrow(dat.corr)], (dat.corr$se.adj[nrow(dat.corr)])^2)
		
	# some quantities needed for meta-analysis
	    est.unadj.1.vec <- dat$b.unadj[noconf]
	    var.unadj.1.vec <- (dat$se.unadj[noconf])^2
	meta.unadj.1 <- unlist(meta.fixed.invvar(est.unadj.1.vec, var.unadj.1.vec), use.names=F)
	    est.unadj.2.vec <- dat$b.unadj[fullconf]
	    var.unadj.2.vec <- (dat$se.unadj[fullconf])^2
	meta.unadj.2 <- unlist(meta.fixed.invvar(est.unadj.2.vec, var.unadj.2.vec), use.names=F)
	# meta-analysis weights
	wt.meta.corradj.1 <- cov.unadj.adj/(meta.unadj.1[2] + meta.unadj.2[2])
	wt.meta.adj.2 <- 1 - wt.meta.corradj.1
	if(wt.meta.corradj.1<0|wt.meta.adj.2<0) message("Warning: one of the fixed-effect meta-analysis weights is negative!")
	
	# meta-analyzed estimates
	beta.meta <- wt.meta.corradj.1*meta.corradj.1[1] + wt.meta.adj.2*meta.adj.2[1]
	var.meta <- meta.adj.2[2] - (cov.unadj.adj^2)/(meta.unadj.1[2] + meta.unadj.2[2])
		
	# return meta-analyzed estimates
	if(returnSE){
		return(list(beta.meta=beta.meta, se.meta=sqrt(var.meta)))
	}else{
		return(list(beta.meta=beta.meta, var.meta=var.meta))
	}
}


###---------------------End: Intermediate functions to implement imputed adjusted estimates---------------------



#################################################################################################################
### Main function to implement CIMBAL: correction for confounder imbalance and subsequent meta-analysis
#################################################################################################################
### Inputs:
###		dat: input data must be a dataframe  with column names "cohort", "samplesize", "b.unadj", "se.unadj", "b.adj", "se.adj"
###			 (use NA's for the missing values, e.g., NA in b.adj and se.adj columns for cohorts without full confounder info)
### Outputs:
###		dat.imputed: the imputed version of 'dat' data frame
###			(i.e., NA's in b.adj and se.adj columns for cohorts without full confounder info will be filled)
###		beta.meta.cimbal: the beta=log(OR) estimate after meta-analyzing CIMBAL-imputed adjusted estimates and adjusted estimates from
###			cohorts with full confounder info
###		SE.meta.cimbal: the SE estimate after meta-analyzing ...

# impute individual-cohort-level beta and variance estimates using meta-analyzed full cohort info
impute.summstat.single <- function(est.unadj.2.vec, var.unadj.2.vec, est.adj.2.vec, var.adj.2.vec, est.unadj.1, var.unadj.1, mute.msgs=FALSE, ncohort.thresh=25){
    # est.unadj.2.vec: vector of unadjusted beta=log(OR) estimates from cohorts with full info
    # var.unadj.2.vec: vector of unadjusted variance estimates from cohorts with full info
    # est.adj.2.vec  : vector of adjusted beta=log(OR) estimates from cohorts with full info
    # var.adj.2.vec	 : vector of adjusted variance estimates from cohorts with full info
    # est.unadj.1	 : the unadjusted beta estimate from the cohort with incomplete info (the one cohort to be "imputed" for confounder imbalance)
    # var.unadj.1	 : the unadjusted variance estimate from the cohort with incomplete info (the one cohort to be "imputed" for confounder imbalance)
    # input checks
    M <- length(est.unadj.2.vec)
    if(length(var.unadj.2.vec)!=M | length(est.adj.2.vec)!=M | length(var.adj.2.vec)!=M) 
        stop("Vectors of adjusted and unadjusted log odds estimates and their variances from cohorts with complete confounder info should be of same length (same as the no. of complete cohorts)")
    if(length(est.unadj.1)>1 | length(var.unadj.1)>1)
        stop("There should be log odds estimate or variance estimate for exactly 1 cohort with incomplete confounder info.")
    if(sum(is.na(est.unadj.2.vec))>0 | sum(is.na(var.unadj.2.vec))>0 | sum(is.na(est.adj.2.vec))>0 | sum(is.na(var.adj.2.vec))>0 | is.na(est.unadj.1) | is.na(var.unadj.1))
        stop('There seems to be missing information for log odds estimates and/or their variances from cohorts needed for the correction approach.')
    
    # meta-analyzed info for the cohorts with full info
    unadj.meta.2 <- meta.fixed.invvar(beta.vec=est.unadj.2.vec, var.vec=var.unadj.2.vec)
    est.unadj.meta.2 <- unadj.meta.2$beta.meta
    var.unadj.meta.2 <- unadj.meta.2$var.meta
    adj.meta.2 <- meta.fixed.invvar(beta.vec=est.adj.2.vec, var.vec=var.adj.2.vec)
    est.adj.meta.2 <- adj.meta.2$beta.meta
    var.adj.meta.2 <- adj.meta.2$var.meta
    
    # obtain covariance between adjusted and unadjusted estimates from complete cohorts
    cov.unadj.adj <- .cov.unadj.adj(est.unadj.2.vec=est.unadj.2.vec, est.adj.2.vec=est.adj.2.vec,
                                    var.unadj.2.vec=var.unadj.2.vec, var.adj.2.vec=var.adj.2.vec, 
                                    var.meta.unadj.2=var.unadj.meta.2, var.meta.adj.2=var.adj.meta.2, ncohort.thresh, mute.msgs)$covar
    
    
    # imputed estimates
    est.tilde.1 <- est.adj.meta.2 - est.unadj.meta.2 + est.unadj.1
    var.tilde.1 <- var.adj.meta.2 + var.unadj.meta.2 + var.unadj.1 - 2*cov.unadj.adj
    if(var.tilde.1<=0){
        warning(paste("Imputed adjusted variance estimate came out non-positive for Cohort",i,"- returning NA."))
        var.tilde.1 <- NA
    }
    return(list(beta.adj.imputed=est.tilde.1, var.adj.imputed=var.tilde.1))
}

# calculate meta-analysis of beta-estimates and their variance estimates
meta.fixed.invvar <- function(beta.vec, var.vec, returnSE=FALSE){
	weights <- 1/var.vec
	beta.meta <- sum(weights*beta.vec)/sum(weights)
	var.meta <- 1/sum(weights)
	if(returnSE){
		return(list(beta.meta=beta.meta, se.meta=sqrt(var.meta)))
	}else{
		return(list(beta.meta=beta.meta, var.meta=var.meta))
	}
}

# the main cimbal function
cimbal <- function(dat, mute.msgs=FALSE, ncohort.thresh=25){
    ### Checks
    # check input type 
    dat <- .checks(dat)
    K <- nrow(dat)
    # which cohorts do not have confounder information
    noconf <- which((is.na(dat$b.adj) | is.na(dat$se.adj)) & !is.na(dat$b.unadj) & !is.na(dat$se.unadj))
    # which cohorts have full confounder information
    fullconf <- which(!is.na(dat$b.adj) & !is.na(dat$se.adj) & !is.na(dat$b.unadj) & !is.na(dat$se.unadj))
    if(length(noconf)+length(fullconf)!=K)
        stop("Some problem with the current data format. Total no. of cohorts with and without complete information does not match with total no. of cohorts available.")
    
    ### obtain imputed adjusted beta, SE estimates
    dat.corr <- as.data.frame(dat)
    # total samplesize of the incomplete cohorts
    n.noconf <- sum(dat$samplesize[noconf])
    # cohort name of the meta-analyzed incomplete cohort
    name.noconf <- paste(dat$cohort[noconf],collapse="+")
    # initialize the imputed dataframe to return
    # remove the individual incomplete cohorts (since meta-analyzed imputed estimate of incomplete cohorts will be used)
    dat.corr <- dat.corr[-noconf,]
    # add a row in the end that will correspond to estimates from the meta-analyzed imputed cohort
    dat.corr[nrow(dat.corr)+1,] <- rep(NA,6)
    dat.corr$cohort[nrow(dat.corr)] <- name.noconf		
    dat.corr$samplesize[nrow(dat.corr)] <- n.noconf	
    # fill in the imputed beta and variance estimates for the meta-analyzed cohort
    # meta-analyzed unadjusted variances from incomplete cohort
    est.unadj.1.vec <- dat$b.unadj[noconf]
    var.unadj.1.vec <- (dat$se.unadj[noconf])^2
    meta.unadj.1 <- unlist(meta.fixed.invvar(est.unadj.1.vec, var.unadj.1.vec), use.names=F)
    # meta-analyzed unadjusted variances from complete cohort to be used to get cov.unadj.adj
    est.unadj.2.vec <- dat$b.unadj[fullconf]
    var.unadj.2.vec <- (dat$se.unadj[fullconf])^2
    meta.unadj.2 <- unlist(meta.fixed.invvar(est.unadj.2.vec, var.unadj.2.vec), use.names=F)
    # meta-analyzed unadjusted variances from complete cohort to be used to get cov.unadj.adj
    est.adj.2.vec <- dat$b.adj[fullconf]
    var.adj.2.vec <- (dat$se.adj[fullconf])^2
    meta.adj.2 <- unlist(meta.fixed.invvar(est.adj.2.vec, var.adj.2.vec), use.names=F)
    cov.unadj.adj <- .cov.unadj.adj(est.unadj.2.vec=dat$b.unadj[fullconf], est.adj.2.vec=dat$b.adj[fullconf],
                                    var.unadj.2.vec=var.unadj.2.vec, var.adj.2.vec=var.adj.2.vec, 
                                    var.meta.unadj.2=meta.unadj.2[2], var.meta.adj.2=meta.adj.2[2], ncohort.thresh, mute.msgs)$covar
    
    tilde.meta.1 <- unlist(.beta.var.tilde.meta(est.unadj.2.vec=dat$b.unadj[fullconf], 
                                                var.unadj.2.vec=(dat$se.unadj[fullconf])^2, 
                                                est.adj.2.vec=dat$b.adj[fullconf], 
                                                var.adj.2.vec=(dat$se.adj[fullconf])^2, 
                                                est.unadj.1.vec=dat$b.unadj[noconf], 
                                                var.unadj.1.vec=(dat$se.unadj[noconf])^2, cov.unadj.adj=cov.unadj.adj), use.names=FALSE)
    
    dat.corr[nrow(dat.corr),c('b.adj','se.adj')] <- c(tilde.meta.1[1], sqrt(tilde.meta.1[2]))
    dat.corr[nrow(dat.corr),c('b.unadj','se.unadj')] <- c(meta.unadj.1[1], sqrt(meta.unadj.1[2]))
    
    # return outputs final meta-analyzed estimates
    meta.cimbal <- .meta.fixed.invwt.cimbal(dat, dat.corr, fullconf, noconf, cov.unadj.adj, ncohort.thresh, mute.msgs, returnSE=TRUE)
    
    return(list(dat.imputed=dat.corr, beta.meta.cimbal=meta.cimbal$beta.meta, se.meta.cimbal=meta.cimbal$se.meta))
}


