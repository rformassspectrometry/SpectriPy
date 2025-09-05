#' @rdname MsBackendPy
setReplaceMethod("spectraVariableMapping", "Spectra", function(object, value) {
    spectraVariableMapping(object@backend) <- value
    object
})

#' @rdname MsBackendPy
setMethod("spectraVariableMapping", "Spectra", function(object) {
    spectraVariableMapping(object@backend)
})
