.onAttach <- function(libname, pkgname) {
    # Suppress the warning about overwriting methods
    suppressWarnings({
        registerS3method("[", "Event", `[.Event`,
            envir = asNamespace(pkgname)
        )
        registerS3method("as.character", "Event", as.character.Event,
            envir = asNamespace(pkgname)
        )
        registerS3method("as.matrix", "Event", as.matrix.Event,
            envir = asNamespace(pkgname)
        )
        registerS3method("format", "Event", format.Event,
            envir = asNamespace(pkgname)
        )
        registerS3method("print", "Event", print.Event,
            envir = asNamespace(pkgname)
        )
        registerS3method("rbind", "Event", rbind.Event,
            envir = asNamespace(pkgname)
        )
        registerS3method("as.data.frame", "Event", as.data.frame.Event,
            envir = asNamespace(pkgname)
        )
        registerS3method("is.na", "Event", is.na.frame.Event,
            envir = asNamespace(pkgname)
        )
        registerS3method("length", "Event", length.Event,
            envir = asNamespace(pkgname)
        )
        registerS3method("c", "Event", c.Event,
            envir = asNamespace(pkgname)
        )
    })
}


.onLoad <- function(libname, pkgname) {

}
