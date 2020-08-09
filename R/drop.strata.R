##' @export
drop.strata <- function(x) {
    drop.specials(x, c('strata', 'strataC'))
}

##' @export
drop.specials <- function(x, components, ...) {
    variables <- c()
    for (comp in components) {
        mm <- lava::Specials(x, comp, ...)
        vars <- unlist(lapply(mm, function(x) strsplit(x,',')))
        variables <- c(variables, list(vars))
        for (i in mm) {
            newf <- as.formula(paste0('.~.-', comp, '(',i,')'))
            x <- update(x, newf)
        }
    }
    names(variables) <- components
    return( structure(x, variables=variables) )
}
