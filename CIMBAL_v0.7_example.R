#-------- Download to local directory & source from there OR directly source CIMBAL from github
# source("CIMBAL_v0.7.R")
require(devtools)
source_url("https://github.com/RayDebashree/CIMBAL/blob/main/CIMBAL_v0.7.R?raw=TRUE")
#--------

### For an example, let's first simulate a toy set of summary 
### statistics from 60 cohorts on a binary outcome (Y), binary exposure (X) 
### and 2 binary confounders (C1, C2)

	# function to simulate data with Y, X, C1 and C2 for sample size n
	getdata <- function(n, a, b, bx, model="glm"){	
	    # simulate confounders
    	c1 <- rbinom(n,1,0.1)
    	c2 <- rbinom(n,1,0.6) 
    	# get the binary exposure with prob px
	    px <- 1/(1+exp(-(a[1]+a[2]*c1+a[3]*c2)))        
	    x <- rbinom(n,1,px)   
	    # simulate binary response Y if "glm" else normal response
	    if(model=="glm"){
	        py <- 1/(1+exp(-(b[1]+b[2]*c1+b[3]*c2+bx*x)))         
	        y <- rbinom(n,1,py)        	
	    }else{
	        e <- rnorm(n,0,1)
	        y <- b[1]+b[2]*c1+b[3]*c2+bx*x+e
	    }
	    return(data.frame(y,x,c1,c2))
	}
	
	# function to get summary statistics from cohorts
	getcohorts <- function(samplesizes, a, b, bx, model, type){
        # no. of cohorts
        K <- length(samplesizes)
        # initialize the data frame of summary stats
        mydat <- as.data.frame(cbind(1:K, samplesizes))
        colnames(mydat) <- c("cohort","samplesize")
        #cohortnames <- sapply(1:K, function(i) paste0("Cohort",i))
        mydat$cohort <- 1:K
        mydat$se.adj <- mydat$b.adj <- mydat$se.unadj <- mydat$b.unadj <- NA
        # initialize the full data
        set.seed(2022)
        for (cohort in 1:K){
            ### simulate the data on cohort of sample size n
            dat1 <- getdata(samplesizes[cohort], a, b, bx)
            
            ### analysis output of dat1
            if(type=="unadj") outu <- glm(y~x, data=dat1, family="binomial", 
control=list(maxit=1e4))
            if(type=="padj") outu <- glm(y~x+c1, data=dat1, family="binomial", 
control=list(maxit=1e4))
            outa <- glm(y~x+c1+c2, data=dat1, family="binomial", 
control=list(maxit=1e4))
			rm(dat1)
            mydat[cohort,3:6] <- c(coef(summary(outu))['x',c('Estimate',
'Std. Error')], coef(summary(outa))['x',c('Estimate','Std. Error')])
        }
        return(mydat)
    }

### Simulating 60 cohorts and obtaining their unadjusted and adjusted
### estimates; data generation models used are:
###     model for X: logit(P(X=1)) = a0 + a1*C1 + a2*C2
###     model for Y: logit(P(Y=1)) = b0 + b1*C1 + b2*C2 + bx*X
samplesizes <- rep(150, 60)
a <- c(log(0.5/0.5), 0.5, 0.5) 
b <- c(log(0.3/0.7), 0.5, 0.5)
bx <- log(1)
mydat <- getcohorts(samplesizes, a, b, bx, model="glm", type="unadj")
### randomly assign 30 cohorts to have only unadjusted estimates
set.seed(1)
noconfcohorts <- sort(sample(1:60, size=30, replace=F))
mydat[which(mydat$cohort %in% noconfcohorts), c('b.adj','se.adj')] <- NA
mydat

### Implementing CIMBAL to impute adjusted estimates for the 
### combined no-confounder cohort and meta-analyze adjusted 
### estimates from all 60 cohorts
out <- cimbal(dat=mydat)
# imputed adjusted estimate for combined 30 no-confounder cohorts
out$dat.imputed[nrow(out$dat.imputed),]
# final exposure-outcome adjusted effect estimate meta-analyzing all
# 60 cohorts
out$beta.meta.cimbal

# corresponding adjusted SE estimate 
out$se.meta.cimbal

### One can also implement CIMBAL if some cohorts report partially 
### adjusted estimates and others report fully adjusted estimates
mydat <- getcohorts(samplesizes, a, b, bx, model="glm", type="padj")
# randomly assign 30 cohorts to have only C1-adjusted estimates
set.seed(1)
noconfcohorts <- sort(sample(1:60, size=30, replace=F))
mydat[which(mydat$cohort %in% noconfcohorts), c('b.adj','se.adj')] <- NA
# implement CIMBAL
out <- cimbal(dat=mydat)
# imputed fully adjusted estimate for combined 30 cohorts with only 
# C1 confounder (note, although output says 'b.unadj' or 'se.unadj',
# they really show partially adjusted or C1-adjusted estimates
# and not unadjusted estimates)
out$dat.imputed[nrow(out$dat.imputed),]
# final exposure-outcome fully adjusted effect estimate meta-analyzing 
# all 60 cohorts
out$beta.meta.cimbal
# corresponding adjusted SE estimate 
out$se.meta.cimbal 
