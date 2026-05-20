### drop.strata <- function(x) {
###    drop.specials(x, c('strata', 'strataC'))
### }

##' Remove Special Terms from a Formula
##'
##' Removes terms such as \code{strata()} or \code{cluster()} from a
##' formula/terms object.
##'
##' @param x a formula object.
##' @param components character vector of special term names to remove.
##' @param ... additional arguments passed to \code{lava::Specials}.
##' @return The updated formula with an attribute \code{"variables"} containing
##'   the extracted variable names from removed specials.
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
