

gsub2 <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
      pattern <- rep(pattern, length.out=length(replacement))
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- sub(pattern[i], replacement[i], result, ...)
  }
  result
}

##' @export
procform <- function(formula=NULL, sep="\\|", nsep=1, return.formula=FALSE, data=NULL,
             no.match=TRUE, regex=FALSE, return.list=TRUE, specials=NULL,...) {# {{{
    res <- NULL
    if (is.null(formula)) {
        res <- colnames(data)
    } else if (is.character(formula)) {
        if (is.null(data)) {
            res <- unique(formula)
        } else {
            yy <-c()
            for (y0 in formula) {
                y0orig <- y0
                if (!regex) y0 <- glob2rx(y0)
                npos <- grep(y0,names(data),perl=mets.options()$regex.perl)
                if (no.match && length(npos)==0) {
                    yy <- union(yy, y0orig)
                } else {
                    yy <- union(yy,names(data)[npos])
                }
            }
            res <- unique(yy)
        }
    }
    if (is.numeric(formula)) res <- colnames(data)[formula]
    if (is.character(res)) {
        if (!return.list) return(res)
        if (return.formula) return(as.formula(paste("~",paste(res,collapse="+"))))
        return(list(response=res,predictor=NULL,filter=NULL))
    }


    ## Add parantheses around quotes if it is not a function call
    if (inherits(formula,"formula")) {

        st <- Reduce(paste,deparse(formula))
        strsplit(st,"\"")
        quotepos <- gregexpr("[\"']",st)[[1]]
        if (quotepos[1]>0) {
            sts <- strsplit(st,"[\"']")[[1]]
            foundsep <- any(grepl("|", sts, fixed=TRUE))
            p <- length(quotepos)
            ##repl <- rep(c("(\"","\")"),p)
            for (i in seq(p/2)*2-1) {
                sts[i] <- paste0(sts[i],"(\"")
                sts[i+1] <- paste0(sts[i+1],"\")")
            }
            ## To handle regular expression entered as strings in the formula, we add a 'filter' expression at the end of the formula
            if (!foundsep) sts <- c(sts,"|1")
            formula <- as.formula(paste(sts,collapse=""))
        }
    }

    aa <- attributes(terms(formula,data=data,specials="regex"))
    if (aa$response == 0) {
        res <- NULL
    } else {
        res <- paste(deparse(formula[[2]]), collapse = "")
    }
    filter.expression <- NULL
    foundsep <- FALSE
    pred <- filter <- c()
    if (!is.null(sep) && length(aa$term.labels) > 0) {
        foundsep <- any(grepl(sep,aa$term.labels))
        if (foundsep) {
            if (nsep>1) {
                xc <- gsub(" ","",unlist(lapply(aa$term.labels, function(z) strsplit(z,sep)[[1]])))
                pred <- xc[1]
                filter <- xc[-1]
            } else {
                xc <- gsub(" ","",unlist(lapply(aa$term.labels, function(z) {
                    spl <- regexpr(sep,z) ## first appearance
                    pred <- substr(z,1,spl-1)
                    filter <- substr(z,spl+1,nchar(z))
                    return(c(pred,filter))
                })))
                pred <- xc[1]
                filter <- xc[2]
            }
            if (any(pred==".")) {
                f <- as.formula(paste0(paste0(c(res,filter),collapse="+"),"~."))
                x <- attributes(terms(f,data=data))$term.labels
                pred <- x
            }
            if (filter%in%c("0","-1")) {
                filter <- list()
                filter.expression <- NULL
            } else {
                filter.expression <- parse(text=filter)
                filter <- as.list(filter)
            }
        }
    }
    if (!foundsep) pred <- aa$term.labels

    expandst <- function(st) {
        st <- res <- unlist(strsplit(gsub(" ","",st),"\\+"))
        if (any(unlist(lapply(st, function(x) grepl("^\\(",x))))) {
            res <- c()
            for (x in st) {
                if (grepl("^\\(",x)) {
                    x <- gsub('\\"',"",x)
                    x <- gsub('^\\(',"",x)
                    x <- gsub('\\)$',"",x)
                    res <- c(res,unlist(procform(x,data=data,regex=regex, no.match=FALSE)$response))
                } else {
                    res <- c(res,x)
                }
                res <- unique(res)
            }
        }
        return(res)
    }
    res <- expandst(res)
    pred <- expandst(pred)
    if (any(res==".")) {
        diffset <- c(".",setdiff(pred,res))
        res <- setdiff(union(res,colnames(data)),diffset)
    }
    filter <- lapply(filter, expandst)
    if (!is.null(specials)) {
        foundspec <- replicate(length(specials),c())
        names(foundspec) <- specials
        rmidx <- c()
        spec <- paste0("^",specials,"\\(")
        val <- lapply(spec, function(x) which(grepl(x,pred)))
        for (i in seq_along(val)) {
            if (length(val[[i]])>0) { # special function found
                rmidx <- c(rmidx,val[[i]])
                cleanpred <- gsub("\\)$","",gsub(spec[i],"",pred[val[[i]]]))
                foundspec[[i]] <- c(foundspec[[i]],cleanpred)
            }
        }
        if (length(rmidx)>0)
            pred <- pred[-rmidx]
        if (length(pred)==0) pred <- NULL
        specials <- foundspec
        for (i in seq_along(specials)) if (is.null(specials[[i]])) specials[i] <- NULL
        if (length(specials)==0) specials <- NULL
    }
    if (return.formula) {
        if (foundsep && !is.null(filter)) {
            filter <- lapply(filter,
                             function(z) as.formula(paste0("~", paste0(z,collapse="+"))))
        }
        if (length(pred)>0)
            pred <- as.formula(paste0("~", paste0(pred,collapse="+"), collapse=""))
        if (length(res)>0)
            res <- as.formula(paste0("~", paste0(res,collapse="+"), collapse=""))
        if (!is.null(specials)) {
            specials <- lapply(specials,function(x)
                              as.formula(paste0("~", paste0(x,collapse="+"), collapse="")))
        }
    }
    res <- list(response=res,
                predictor=pred,
                filter=filter,
                filter.expression=filter.expression,
                specials=specials)
    if (!return.list) return(unlist(unique(res)))
    return(res)
}

##' @export
procformdata <- function(formula,data,sep="\\|", na.action=na.pass, do.filter=TRUE, ...) {
    res <- procform(formula,sep=sep,data=data,return.formula=TRUE,...)
    if (inherits(res,"formula")) {
        res <- list(response=res)
    }
    y <- x <- NULL
    filter <- res$filter.expression
    if (!do.filter) {
        filter <- NULL
    }
    ### when filter.expression is expression(1) then also no filter, ts
    if ((!missing(filter))) if (!is.null(filter)) if (as.character(filter)=="1") filter <- NULL

    if (length(res$response)>0) {
        if (is.null(filter)) y <- model.frame(res$response,data=data,na.action=na.action)
        else y <- model.frame(res$response,data=subset(data,eval(filter)),na.action=na.action)
    }
    if (length(res$predictor)>0) {
        if (is.null(filter)) x <- model.frame(res$predictor,data=data,na.action=na.action)
        else x <- model.frame(res$predictor,data=subset(data,eval(filter)),na.action=na.action)

    }

    specials <- NULL
    if (!is.null(res$specials)) {
        specials <- lapply(res$specials,
                      function(x) {
                          if (is.null(filter)) model.frame(x,data=data,na.action=na.action)
                          else model.frame(x,data=subset(data,eval(filter)),na.action=na.action)
                      })
    }

    if (!do.filter) {
        group <- lapply(res$filter, function(x) model.frame(x,data=data,na.action=na.action))
        return(list(response=y,predictor=x,group=group,specials=specials))
    }
    return(list(response=y,predictor=x,specials=specials))
}# }}}

procform2 <- function(y,x=NULL,z=NULL,...) {# {{{
    yx <- procform(y,return.formula=FALSE,...)
    y <- yx$response
    x0 <- yx$predictor
    z0 <- NULL
    if (length(yx$filter)>0) z0 <- yx$filter[[1]]
    if (is.null(x) && length(y)>0) x <- x0
    if (NCOL(x)==0) x <- NULL
    if (length(y)==0) y <- x0
    if (!is.null(x)) {
        x <- unlist(procform(x,return.formula=FALSE,...))
    }
    if (is.null(z)) z <- z0
    if (!is.null(z)) {
        zz <- unlist(procform(z,return.formula=FALSE,...))
    }
    if (is.null(z)) z <- z0
    return(list(y=y,x=x,z=z))
}# }}}

##' @export
procform3 <- function(y,x=NULL,z=NULL,...) {# {{{
    yx <- procform(y,return.formula=FALSE,...)
    x0 <- yx$predictor
    y  <- yx$response
    if (is.null(yx$predictor)) { x0 <- yx$response ; y <- NULL}

    if (is.null(y))
    if (!is.null(x)) {
        x <- procform(x,return.formula=FALSE,...)
        y <- c(x$predictor,x$response)
    }
    return(list(y=y,x=x0,z=z))
}# }}}


Specials <- function(f,spec,split1=",",split2=NULL,...) {# {{{
  tt <- terms(f,spec)
  pos <- attributes(tt)$specials[[spec]]
  if (is.null(pos)) return(NULL)
  x <- rownames(attributes(tt)$factors)[pos]
  st <- gsub(" ","",x) ## trim
  spec <- unlist(strsplit(st,"[()]"))[[1]]
  res <- substr(st,nchar(spec)+2,nchar(st)-1)
  if (!is.null(split1))
      res <- unlist(strsplit(res,split1))
  res <- as.list(res)
  for (i in seq_along(res)) {
      if (length(grep("~",res[[i]]))>0) {
          res[[i]] <- as.formula(res[[i]])
      }
  }
  return(res)
}# }}}

decomp.specials <- function (x, pattern = "[()]", sep = ",", ...)
  {
    st <- gsub(" ", "", x)
    if (!is.null(pattern))
      st <- rev(unlist(strsplit(st, pattern, ...)))[1]
    unlist(strsplit(st, sep, ...))
  }


model.extract2 <- function(frame, component) {
  component <- as.character(substitute(component))
  if (component %in% c("response", "offset")) {
    return(do.call(
      model.extract,
      list(frame = frame, component = component)
    ))
  }
  vname <- paste0("(", component, ")")
  if (!(vname %in% names(frame))) {
    regex <- paste0("^", component, "\\(.*\\)$")
    if (any(grepl(regex, names(frame)))) {
      vname <- grep(regex, names(frame))
      if (length(vname) > 1) stop("model.extract2: non-unique component")
    }
  }
  rval <- frame[[vname]]

  if (!is.null(rval)) {
    if (length(rval) == nrow(frame)) {
      names(rval) <- attr(frame, "row.names")
    } else if (is.matrix(rval) && nrow(rval) == nrow(frame)) {
      t1 <- dimnames(rval)
      dimnames(rval) <- list(
        attr(frame, "row.names"),
        t1[[2L]]
      )
    }
  }
  return(rval)
}

# #' Extract design matrix from data.frame and formula
# #' @title Extract design matrix
# #' @param formula formula
# #' @param data data.frame
# #' @param intercept (logical) If FALSE an intercept is not included in the
# #'   design matrix
# #' @param response (logical) if FALSE the response variable is dropped
# #' @param rm_envir Remove environment
# #' @param ... additional arguments (e.g, specials such weights, offsets, ...)
# #' @param specials character vector specifying functions in the formula that
# #'   should be marked as special in the [terms] object
# #' @param specials.call (call) specials optionally defined as a call-type
# #' @param levels a named list of character vectors giving the full set of
# #'   levels
# #'   to be assumed for each factor
# #' @param design.matrix (logical) if FALSE then only response and specials are
# #'   returned. Otherwise, the design.matrix `x` is als part of the returned
# #'   object.
# #' @return An object of class 'mets.design'
# #' @examples
# #' n <- 1e3
# #' a <- rbinom(n, 1, 0.5)
# #' t <- rweibull(n, shape=0.5, exp(3 + a))
# #' d <- data.frame(time=t, status=TRUE, a=a, id=seq(n), x=rnorm(n))
# #' des <- proc_design(
# #' Surv(time, status) ~ x + strata(a) + cluster(id),
# #'   specials = c("strata", "cluster"),
# #'   data=d)
# #' new <- update_design(des, head(d))
# #' new$strata
# #' new$cluster
# #' new$y
# #' new$x
proc_design <- function(formula, data, ..., # nolint
                   intercept = FALSE,
                   response = TRUE,
                   rm_envir = FALSE,
                   specials = NULL,
                   specials.call = NULL,
                   levels = NULL,
                   design.matrix = TRUE) {
  dots <- substitute(list(...))
  if ("subset" %in% names(dots)) stop(
    "subset is not an allowed specials argument for `design`"
  )
  tt <- terms(formula, data = data, specials = specials)
  term.labels <- attr(tt, "term.labels") # predictors

  if (response && inherits(
      try(model.frame(update(tt, ~1), data = data), silent = TRUE),
      "try-error"
  )) {
      response <- FALSE
  }
  # delete response to generate design matrix when making predictions
  if (!response) tt <- delete.response(tt)

  sterm.list <- c()
  if (length(specials) > 0) {
    des <- attr(tt, "factors")
    for (s in specials) {
      sterm <- rownames(des)[attr(tt, "specials")[[s]]]
      sterm.list <- c(sterm.list, sterm)
    }
    if (length(sterm.list) > 0) {
      # create formula without specials
      if ((nrow(attr(tt, "factors")) - attr(tt, "response")) ==
          length(sterm.list)) {
        # only specials on the rhs, remove everything
        print(formula)
          formula <- update(formula, ~1)
      } else {
          # predictors without the specials
          term.labels <- setdiff(
              term.labels,
              unlist(sterm.list)
          )
          ## formula <- update(tt, reformulate(term.labels))
      }
      upd <- paste(" ~ . - ", paste(sterm.list, collapse = " - "))
      formula <- update(formula, upd)
    }
  }

  xlev <- levels
  xlev[["response_"]] <- NULL
  if (!design.matrix) { # only extract specials, response
    des <- attr(tt, "factors")
    fs <- update(formula, ~1)
    if (length(sterm.list) > 0) {
      # formula with only special-terms
      fs <- reformulate(paste(sterm.list, collapse = " + "))
      fs <- update(formula, fs)
    }
    mf <- model.frame(fs, data=data, ...)
  } else { # also extract design matrix
    mf <- model.frame(tt,
        data = data, ...,
        xlev = xlev,
        drop.unused.levels = FALSE
        )
    if (is.null(xlev)) {
      xlev <- .getXlevels(tt, mf)
    }
    xlev0 <- xlev
  }

  y <- NULL
  if (response) {
    y <- tryCatch(
      model.response(mf, type = "any"),
          error = function(...) NULL
      )
      if (is.factor(y) || is.character(y)) {
          ylev <- levels[["response_"]]
          if (!is.null(ylev)) {
              factor(y, levels = ylev)
          } else {
              ylev <- if (is.factor(y)) {
                  levels(y)
              } else if (is.character(y)) {
                  levels(as.factor(y))
              }
              levels[["response_"]] <- ylev
          }
      }
  }

  has_intercept <- attr(tt, "intercept") == 1L
  specials <- union(
    specials,
    names(dots)[-1] # removing "" at first position when calling dots, which
  ) # is a call object

  specials.list <- c()
  specials.var <- c()
  if (length(specials) > 0) {
    for (s in specials) {
        w <- eval(substitute(model.extract2(mf, s), list(s = s)))
        specials.list <- c(specials.list, list(w))
        specials.var <- c(
            specials.var,
            list(unlist(Specials(tt, spec = s)))
        )
    }
    names(specials.var) <- specials
    names(specials.list) <- specials
    if (length(sterm.list) > 0) {
      if (design.matrix) {
        xlev0[sterm.list] <- NULL
        mf <- model.frame(formula(delete.response(terms(formula))),
                          data = data, ...,
                          xlev = xlev0,
                          drop.unused.levels = FALSE
                          )
      }
    }
  }

  if (!is.null(specials.call)) {
    specials.list2 <- eval(specials.call, data)
    for (n in names(specials.list2)) {
      if (is.null(specials.list[[n]])) {
        specials.list[[n]] <- specials.list2[[n]]
      }
    }
  }

  if (design.matrix) {
    x <- model.matrix(mf, data = data, xlev = xlev0)
    if (!intercept && has_intercept) {
      has_intercept <- FALSE
      x <- x[, -1, drop = FALSE]
    }
  } else {
    term.labels <- NULL
    x <- NULL
  }

  if (rm_envir) attr(tt, ".Environment") <- NULL
  if (is.null(specials.call)) specials.call <- dots

  xlev[["response_"]] <- levels[["response_"]]
  res <- c(
    list(
      formula = formula, # formula without specials
      terms = tt,
      term.labels = term.labels,
      levels = xlev,
      x = x, y = y,
      design.matrix = design.matrix,
      intercept = has_intercept,
      data = data[0, ], ## Empty data.frame to capture structure of data
      specials = specials,
      specials.var = specials.var,
      specials.call = specials.call
    ),
    specials.list
  )
  return(structure(res, class="mets.design"))
}

update_design <- function(object, data = NULL, response=FALSE,  ...) {
  if (is.null(data)) data <- object$data
  return(
    proc_design(object$terms,
      data = data,
      design.matrix = object$design.matrix,
      levels = object$levels,
      response = response,
      intercept = object$intercept,
      specials = object$specials,
      specials.call = object$specials.call
    )
  )
}

