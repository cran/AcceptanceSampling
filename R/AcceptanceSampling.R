## AcceptanceSampling.R --- 
##
## Author: Andreas Kiermeier
##
## Created: 08 Mar 2007 
##
## Purpose: A package to provide functionality for creating and
##          evaluating acceptance sampling plans.
## ----------------------------------------------------------------------

## ----------------------------------------------------------------------
## Class definitions
## ----------------------------------------------------------------------

setClass("OC2c", representation(n="numeric", ## A vector of sample sizes at each
                              ## stage of sampling - NOT comulative sample size
                              c="numeric", ## vector of acceptance numbers for
                              ## each stage of sampling. Accept if actual number
                              ## of defectives/defects is <= c
                              r="numeric", ## vector of rejection numbers for
                              ## each stage of sampling. Reject if actual number
                              ## of defectives/defects is >= r
                              type="character",
                              paccept="numeric",
                              "VIRTUAL"),
         validity=function(object){
           if(any(is.na(object@n)) | any(is.na(object@c)) |
              any(is.na(object@r)))
             return("Missing values in 'n', 'c', or 'r'")
           ## Check that n, c and r are of the same length
           l <- length(object@n)
           if (l != length(object@c) | l != length(object@r))
             return("n, c and r must be of same length.")
           ## Check that the acceptance numbers make sense
           if (any(object@c < 0) | any(object@c > object@n))
             return("Acceptance number(s) 'c' must be in the range [0,n]")
           ## Check that the rejection numbers make sense
           if (any(object@r < 0) | any(object@r > object@n))
             return("Rejection number(s) 'r' must be in the range [0,n].")
           if (any(object@r <= object@c))
             return("Rejection number(s) 'r' must be greater than acceptance number(s) 'c'.")
           ## For double sampling (or more) make sure that acceptance and
           ## rejection number are non-decreasing and non-increasing, respectively.
           if (l > 1) {
             if (any(diff(object@c)<0) )
               return("'c' must be non-decreasing")
             if (any(diff(object@r)<0) )
               return("'r' must be non-decreasing")
           }
           ## Check that a decision is made on the last sample
           if (object@r[l] != object@c[l] + 1)
             return("Decision from last sample cannot be made: r != c+1")
           ## Otherwise things seem fine.
           return(TRUE)
         })

setClass("OCbinomial",
         representation("OC2c",
                        pd="numeric"),
         contains="OC2c",
         prototype=list("OC2c", type="binomial", pd=seq(0,1,by=0.01)),
         validity=function(object){
           ## Check that the proportion of defectives make sense
           if (any(is.na(object@pd)))
             return("Missing values in 'pd'")
           if (any(object@pd < 0.) | any(object@pd > 1.) )
             return("Proportion defectives must be in the range [0,1]")
         })

setClass("OChypergeom",
         representation("OC2c",
                        N="numeric",
                        pd="numeric"),
         contains="OC2c",
         prototype=list("OC2c", type="hypergeom",N=100,pd=(0:100)/100),
         validity=function(object){
           ## Check that the population size of of length 1
           if (length(object@N) > 1)
             return("Length of population size 'N' != 1")
           if (is.na(object@N))
             return("Missing value in 'N'")
           ## Check that the population size is non-negative
           if (object@N < 0.)
             return("Population size 'N' must be non-negative")
           ## Check that the proportion of defectives make sense
           if (any(is.na(object@pd)))
             return("Missing value in 'pd'")
           if (any(object@pd < 0.) | any(object@pd > 1) )
             return("Proportion defectives 'pd' must be in the range [0,1]")
         })

setClass("OCpoisson",
         representation("OC2c",
                        pd="numeric"),
         contains="OC2c",
         prototype=list("OC2c", type="poisson",pd=seq(0,1,0.01)),
         validity=function(object){
           ## Check that the proportion of defectives make sense
           if (any(is.na(object@pd)))
             return("Missing values in 'pd'")
           if (any(object@pd < 0.))
             return("Rate of defects 'pd' must be non-negative")
         })


## ----------------------------------------------------------------------
## Methods to create new object and calculate P(accept)
## Only OC2c to be exported
## other functions are helpers only
## ----------------------------------------------------------------------

OC2c <- function(n,c,r=if (length(c)==1) c+1 else NULL,
               type=c("binomial","hypergeom", "poisson"), ...){
  ## Decide on what 'type' to use
  type <- match.arg(type)
  OCtype <- paste("OC",type,sep="")

  ## Create a new object of that type
  obj <- new(OCtype, n=n, c=c, r=r, type=type, ...)
  
  ## Evaluate the probability of acceptance for this type and given
  ## pd.
  ## First get the generic calculation function
  OCtype <- getFromNamespace(paste("calc.",OCtype,sep=""),
                ns="AcceptanceSampling")

  ## now, based on the type, decide on what to pass to the function
  ## Only need to check for existing type since new() would have stuffed up
  ## if we don't have a class for the type.
  if (type =="binomial")
    obj@paccept <- OCtype(n=obj@n, c=obj@c, r=obj@r, pd=obj@pd) 
  if (type =="hypergeom")
    obj@paccept <- OCtype(n=obj@n, c=obj@c, r=obj@r, N=obj@N, D=obj@pd*obj@N) 
  if (type =="poisson")
    obj@paccept <- OCtype(n=obj@n, c=obj@c, r=obj@r, pd=obj@pd) 

  obj
}



calc.OCbinomial <- function(n,c,r,pd)
{
  ## n needs to be cumulative since c and r are specified that way too.
  n <- cumsum(n)
  ## Get a list with a vector for each pd.
  ## Length of vector equals number of samples, e.g. double = length 2.
  p.accept <- lapply(pd, FUN=function(el) pbinom(q=c, size=n, prob=el))
  p.unsure <- lapply(pd, FUN=function(el) {
    pbinom(q=(r-1), size=n, prob=el) - pbinom(q=c, size=n, prob=el)})

  ## Now combine the sampling stages via helper function
  pa <- mapply(FUN=calc.paccept, p.accept=p.accept, p.unsure=p.unsure)
  pa
}


calc.OChypergeom <- function(n,c,r,N,D)
{
  ## n needs to be cumulative since c and r are specified that way too.
  n <- cumsum(n)
  ## Get a list with a vector for each pd.
  ## Length of vector equals number of samples, e.g. double = length 2.
  ## Note that c, n, D, and N are not required to be integers. the
  ## Hypergeometric function deals with them "appropriately" and we
  ## adopt the same treatment here.
  ## white ball is defect; black is non-defect
  p.accept <- lapply(D, FUN=function(el) phyper(q=c, m=el, n=N-el, k=n))
  p.unsure <- lapply(D, FUN=function(el) {
    phyper(q=(r-1), m=el, n=N-el, k=n) - phyper(q=c, m=el, n=N-el, k=n)})

  ## Now combine the sampling stages via helper function
  pa <- mapply(FUN=calc.paccept, p.accept=p.accept, p.unsure=p.unsure)
  pa
}

calc.OCpoisson <- function(n,c,r,pd)
{
  ## n needs to be cumulative since c and r are specified that way too.
  n <- cumsum(n)
  ## Get a list with a vector for each pd.
  ## Length of vector equals number of samples, e.g. double = length 2.
  ## The rate of defects is given per item.  Need to convert to
  ## rate per sample size (multiply by n)
  p.accept <- lapply(pd, FUN=function(el) ppois(q=c, lambda=el*n) )
  p.unsure <- lapply(pd, FUN=function(el) {
    ppois(q=(r-1), lambda=el*n) - ppois(q=c, lambda=el*n)})

  ## Now combine the sampling stages via helper function
  pa <- mapply(FUN=calc.paccept, p.accept=p.accept, p.unsure=p.unsure)
  pa
}

calc.paccept <- function(p.accept, p.unsure)
{
  ## Purpose: From to matrices (one column per sample stage) calculate the
  ##          overall probability of acceptance
  ##          Not to be used directly by the user.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## p.accept: P(accept) at each stage, i, of sampling, i.e. P(X_i <= c_i)
  ## p.unsure: P(unsure) at each stage, i, of sampling, i.e. P(c_i < X_i < r_i)
  ## ----------------------------------------------------------------------
  ## Author: Andreas Kiermeier, Date: 14 May 2007, 15:56
  pu <- c(1, cumprod(p.unsure))
  pa <- sum(p.accept*pu[-length(pu)])
}




## ----------------------------------------------------------------------
## Printing methods and functions
## ----------------------------------------------------------------------

OC2c.show.default <-
  function(object){
    if(length(object@n)==0){
      x <- matrix(rep(NA,3), ncol=1)
    }
    else
      x <- rbind(object@n, object@c, object@r)
    dimnames(x) <- list(c("Sample size(s)", "Acc. Number(s)",
                          "Rej. Number(s)"),
                        paste("Sample", 1:ncol(x)))
    show(x)
  }

OC2c.show.prob <-
  function(object) {
    if (object@type=="binomial") {
      x <- cbind(object@pd, object@paccept)
      colnames(x) <- c("Prop. defective","P(accept)")
    }
    else if (object@type=="hypergeom"){
      x <- cbind(object@pd*object@N, object@pd, object@paccept)
      colnames(x) <- c("Pop. Defectives", "Pop. Prop. defective","P(accept)")
    }
    else if (object@type=="poisson"){
      x <- cbind(object@pd, object@paccept)
      colnames(x) <- c("Rate of defects","P(accept)")
    }
    else
      stop("No full print method defined for this type")
    
    rownames(x) <- rep("", length(object@paccept))
    show(x)
  }


setMethod("show", "OC2c",
          function(object){
            cat(paste("Acceptance Sampling Plan (",object@type,")\n\n",sep=""))
            OC2c.show.default(object)
          })

setMethod("show", "OChypergeom",
          function(object){
            cat(paste("Acceptance Sampling Plan (",
                      object@type," with N=",object@N,")\n\n",sep=""))
            OC2c.show.default(object)
          })

setMethod("summary", "OC2c",
          function(object, full=FALSE){
            cat(paste("Acceptance Sampling Plan (",object@type,")\n\n",sep=""))
            OC2c.show.default(object)
            if (full){
              cat("\nDetailed acceptance probabilities:\n\n")
              OC2c.show.prob(object)
            }
          })

setMethod("summary", "OChypergeom",
          function(object, full=FALSE){
            cat(paste("Acceptance Sampling Plan (",
                      object@type," with N=",object@N,")\n\n",sep=""))
            OC2c.show.default(object)
            if (full){
              cat("\nDetailed acceptance probabilities:\n\n")
              OC2c.show.prob(object)
            }
          })


## ----------------------------------------------------------------------
## Plotting methods
## ----------------------------------------------------------------------

setMethod("plot", signature(x="OCbinomial", y="missing"),
          function(x, y, type="o", ylim=c(0,1),...){
            plot(x@pd, x@paccept, type=type,
                 xlab="Proportion defective", ylab="P(accept)",
                 ylim=ylim, ...)
          })

setMethod("plot", signature(x="numeric", y="OCbinomial"),
          function(x, y, type="o", ylim=c(0,1),...){
            plot(x, y@paccept, type=type,
                 ylab="P(accept)", ylim=ylim, ...)
          })


setMethod("plot", signature(x="OChypergeom", y="missing"),
          function(x, type="p", ylim=c(0,1), axis=c("pd","D","both"), ...){
            xs <- match.arg(axis)

            if (xs=="pd")
              plot(x@pd, x@paccept, type=type,
                   xlab=paste("Proportion of population defectives (N=",x@N,")",sep=""),
                   ylab="P(accept)", ylim=ylim, ...)
            else if (xs=="D")
              plot(x@pd*x@N, x@paccept, type=type,
                   xlab=paste("Population defectives, D (N=",x@N,")",sep=""),
                   ylab="P(accept)", ylim=ylim, ...)
            else if (xs=="both") {
              plot(x@pd, x@paccept, type=type,
                   xlab=paste("Proportion of population defectives, (N=",x@N,")",sep=""),
                   ylab="P(accept)", ylim=ylim, mar=c(5,4,5,2)+0.1,...)
              ax <- axis(1)
              axis(3, at=ax, labels=ax*x@N)
              mtext(paste("Population defectives, D (N=",x@N,")",sep=""),
                    side=3, line=3)
            }
          })

setMethod("plot", signature(x="numeric", y="OChypergeom"),
          function(x, y, type="p", ylim=c(0,1), ...){
            plot(x, y@paccept, type=type,
                 ylab="P(accept)", ylim=ylim, ...)
          })


setMethod("plot", signature(x="OCpoisson", y="missing"),
          function(x, y, type="o", ylim=c(0,1),...){
            plot(x@pd, x@paccept, type=type,
                 xlab="Rate of defects", ylab="P(accept)",
                 ylim=ylim, ...)
          })

setMethod("plot", signature(x="numeric", y="OCpoisson"),
          function(x, y, type="o", ylim=c(0,1),...){
            plot(x, y@paccept, type=type,
                 ylab="P(accept)", ylim=ylim, ...)
          })


## ----------------------------------------------------------------------
## Methods to evaluation risk points
## All these functions are helpers only and should not be exported
## "assess" methods are exported
## ----------------------------------------------------------------------

check.paccept <-
  function(pa){
    ## Purpose: Utility function to check that supplied P(accept) values
    ##          fall within [0,1]
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## pa: a vector of P(accept) values
    ## ----------------------------------------------------------------------
    ## Author: Andreas Kiermeier, Date: 16 May 2007, 10:19
    if (any(pa < 0) | any(pa > 1))
      return(FALSE)
    return(TRUE)
  }

check.quality <-
  function(pd, type){
    ## Purpose: Utility function to check that supplied Proportion defective
    ##          values fall within
    ##          [0,1] for the binomial or hypergeometric
    ##          [0,inf] for the poisson
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## pd: a vector of proportion defective values
    ## ----------------------------------------------------------------------
    ## Author: Andreas Kiermeier, Date: 16 May 2007, 10:19
    if (any(pd < 0))
      return(FALSE)
    if ((type=="binomial" | type=="hypergeom") & any(pd > 1))
      return(FALSE)
    return(TRUE)
  }

assess.OC2c <-
  function(object, PRP, CRP){
    ## Purpose: This is the function that does the work.
    ##          Evaluate whether a particular sampling plan can meet
    ##          specified producer and/or consumer risk points
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## object: An object of class OC2c
    ## PRP   : Producer risk point in the form c(pdefect, paccept)
    ## CRP   : Consumer risk point in the form c(pdefect, paccept)
    ## print : Print the result
    ## ----------------------------------------------------------------------
    ## Author: Andreas Kiermeier, Date: 16 May 2007, 10:19
    
    planOK <- TRUE
    ## Check that what we are given is OK
    if (missing(PRP))
      PRP <- rep(NA,3)
    else if (!missing(PRP)){
      if( !check.quality(PRP[1], type=object@type) |
         !check.paccept(PRP[2]) )
        stop("Quality and/or desired P(accept) out of bound")

      ## Get the appropriate function for the distribution and
      ## calculate the P(accept)
      calc.pa <- getFromNamespace(paste("calc.OC",object@type,sep=""),
                     ns="AcceptanceSampling")

      pa <- switch(object@type,
                   binomial=calc.pa(object@n, object@c, object@r, PRP[1]),
                   hypergeom=calc.pa(object@n, object@c, object@r, object@N, PRP[1]*object@N),
                   poisson=calc.pa(object@n, object@c, object@r, PRP[1]))
      
      PRP <- c(PRP, pa)

      ## Check that the plan meets the desired point
      ## For PRP have to have P(accept) greater than desired prob.
      if (pa >= PRP[2])
        planOK <- TRUE
      else
        planOK <- FALSE
    }

    
    if (missing(CRP))
      CRP <- rep(NA,3)
    else if (!missing(CRP)){
      if( !check.quality(CRP[1], type=object@type) |
         !check.paccept(CRP[2]) )
        stop("Quality and/or desired P(accept) out of bound")
      ## Get the appropriate function for the distribution and
      ## calculate the P(accept)
      calc.pa <- getFromNamespace(paste("calc.OC",object@type,sep=""),
                     ns="AcceptanceSampling")
      pa <- switch(object@type,
                   binomial=calc.pa(object@n, object@c, object@r, CRP[1]),
                   hypergeom=calc.pa(object@n, object@c, object@r, object@N, CRP[1]*object@N),
                   poisson=calc.pa(object@n, object@c, object@r, CRP[1]))

      CRP <- c(CRP, pa)
      ## Check that the plan meets the desired point
      ## For CRP have to have P(accept) less than desired prob.
      if (pa <= CRP[2])
        planOK <- planOK & TRUE
      else
        planOK <- planOK & FALSE
    }
    return(list(OK=planOK, PRP=PRP, CRP=CRP))
  }

setGeneric("assess", function(object, PRP, CRP, print=TRUE)
           standardGeneric("assess"))

setMethod("assess", signature(object="OC2c"),
          function(object, PRP, CRP, print)
          {
            ## Purpose: Evaluate whether a particular sampling plan can meet
            ##          specified producer and/or consumer risk points
            ## ----------------------------------------------------------------------
            ## Arguments:
            ## object: An object of class OCbinomial
            ## PRP   : Producer risk point in the form c(pdefect, paccept)
            ## CRP   : Consumer risk point in the form c(pdefect, paccept)
            ## print : Print the result
            ## ----------------------------------------------------------------------
            ## Author: Andreas Kiermeier, Date: 16 May 2007, 10:19

            if(!hasArg(PRP) & !hasArg(CRP))
              stop("At least one risk point, PRP or CRP, must be specified")

            plan <- assess.OC2c(object, PRP, CRP)
            if (print) {
              cat(paste("Acceptance Sampling Plan (",object@type,")\n\n",sep=""))
              OC2c.show.default(object)
              cat(paste("\nPlan", ifelse(plan$OK, "CAN","CANNOT"),
                        "meet desired risk point(s):\n\n"))

              ## Both PRP and CRP
              if(hasArg(PRP) & hasArg(CRP))
                RP <- cbind(PRP=plan$PRP, CRP=plan$CRP)
              ## Only PRP
              else if (hasArg(PRP))
                RP <- cbind(PRP=plan$PRP)
              ## Only CRP
              else if (hasArg(CRP))
                RP <- cbind(CRP=plan$CRP)

              rownames(RP) <- c("       Quality", "  RP P(accept)", "Plan P(accept)")
              show(t(RP))
            }

            if(object@type=="hypergeom")
              return(invisible(c(list(n=object@n, c=object@c, r=object@r,
                                      n=object@N), plan)))
            else
              return(invisible(c(list(n=object@n, c=object@c, r=object@r), plan)))
          })


### Local Variables:
### comment-start: "## "
### mode: auto-fill
### fill-column: 80
### End:
