##' Risk predictions to work with riskRegression package
##'
##' Risk predictions to work with riskRegression package
##'
##' @param object phreg/binreg/cifreg object
##' @param newdata newdata
##' @param times times for predictions
##' @param ... additional arguments to lower level functions
##' @author Thomas Scheike
##' @aliases predictRisk predictRisk.phreg predictRisk.binreg predictRisk.cifreg
##'   predictRisk.cifregFG predictRisk.recreg
##' @export
predictRisk <- function(object, ...) UseMethod("predictRisk")

##' @export
predictRisk.phreg <- function(object,newdata,times=NULL,...)
{# {{{
	pcif <- predict(object,newdata,times=times,se=0)
	return(c(1-pcif$surv))
}# }}}

##' @export
predictRisk.binreg <- function(object,newdata,cause,times=NULL,...)
{# {{{
	pcif <- predict(object,newdata,se=0)
	return(c(pcif))
}# }}}

##' @export
predictRisk.cifreg <- function(object,newdata,cause,times=NULL,...)
{# {{{b
	pcif <- predict(object,newdata,times=times,se=0)$cif
	return(c(pcif))
}# }}}

##' @export
predictRisk.cifregFG <- function(object,newdata,cause,times=NULL,...)
{# {{{
	pcif <- predict(object,newdata,times=times,se=0)$cif
	return(c(pcif))
}# }}}

##' @export
predictRisk.recreg <- function(object,newdata,cause,times=NULL,...)
{# {{{
	pcif <- predict(object,newdata,times=times,se=0)$cumhaz
	return(c(pcif))
}# }}}
