#' @title Low level functions to convert between Spectra and matchms Spectrum
#'
#' @name rspec_to_pyspec
#'
#' @description
#'
#' The `rspec_to_pyspec()` and `pyspec_to_rspec()` functions allow to convert
#' R [Spectra::Spectra()] objects into
#' [matchms](https://github.com/matchms/matchms) Python `matchms.Spectrum`
#' objects. These functions are designed for
#' **advanced users or developers** who want/need to integrate Python/matchms
#' functionality into R using *reticulate*. All other users should use the
#' dedicated R functions within this package that take care of running the
#' Python code in the correct Python environment.
#'
#' Parameter `mapping` allows to define which spectra variables (metadata)
#' should be copied between the R and Python spectra. Only provided spectra
#' variables will be copied to R respectively Python. `mapping` also defines
#' the mapping between the `Spectra`'s spectra variables and the Spectrum
#' metadata. The names of the character vector `mapping` are the R spectra
#' variables and the values the corresponding names in the Python's Spectrum
#' metadata. See the output of the `spectraVariableMapping()` function for the
#' default variables and the mapping of the names.
#'
#' The `spectraVariableMapping()` function provides a default mapping of some
#' core `Spectra` variables based on this [definition in matchms](https://github.com/matchms/matchms/blob/master/matchms/data/known_key_conversions.csv).
#' The function returns a named vector that can be directly used as parameter
#' `mapping` in the `rspec_to_pyspec()` and `pyspec_to_rspec()` functions.
#'
#' @param .check Optionally disable input parameter checking. Input parameter
#'     checking should only disabled for very good reasons.
#'
#' @param BPPARAM Optional parallel processing setup.
#' 
#' @param pythonLibrary character(1) specifying the Python library to be used
#'
#' @param x For `rspec_to_pyspec()`: `Spectra` object. For `pyspec_to_rspec()`:
#'     a Python list of matchms Spectrum objects.
#'
#' @param ... ignored.
#'
#' @return For `rspec_to_pyspec()`: Python array of Spectrum objects, same
#'     length than `x`. For `pyspec_to_rspec()`: [Spectra::Spectra()] with the
#'     converted spectra. For `spectraVariableMapping()`: named `character`
#'     vector with names being `Spectra` variable names and values the
#'     corresponding names in `matchms`.
#'
#' @author Michael Witting, Johannes Rainer
#'
#' @export
#'
#' @importFrom reticulate r_to_py import py_to_r
#'
#' @importFrom BiocParallel SerialParam bplapply
#'
#' @importMethodsFrom Spectra spectrapply
#'
#'
NULL

#' @importMethodsFrom Spectra spectraVariableMapping
#'
#' @exportMethod spectraVariableMapping
#'
#' @rdname rspec_to_pyspec
#' 
#' @examples
#' ## List the default spectra variables and their mapping.
#' spectraVariableMapping(pythonLibrary = "matchms")
#' spectraVariableMapping(pythonLibrary = "spectrum_utils")
setGeneric("spectraVariableMapping", function(pythonLibrary, ...)
    standardGeneric("spectraVariableMapping"))
setMethod("spectraVariableMapping", "character", 
    function(pythonLibrary = c("matchms", "spectrum_utils"), ...) {
    
    pythonLibrary <- match.arg(pythonLibrary)
    
    if (pythonLibrary == "matchms")
        mapping <- .SPECTRA_2_MATCHMS
    if (pythonLibrary == "spectrum_utils")
        mapping <- .SPECTRA_2_SPECTRUM_UTILS
    
    ## return the mapping
    mapping
})

.SPECTRA_2_MATCHMS <- c(
    precursorMz = "precursor_mz",
    precursorIntensity = "precursor_intensity",
    precursorCharge = "charge",
    rtime = "retention_time",
    collisionEnergy = "collision_energy",
    isolationWindowTargetMz = "isolation_window_target_mz",
    ## polarity = "ionmode", # Disabling since matchms does not support int.
    msLevel = "ms_level"
)

.SPECTRA_2_SPECTRUM_UTILS <- c(
    precursorMz = "precursor_mz",
    ##precursorIntensity = "",
    precursorCharge = "precursor_charge",
    rtime = "retention_time"
    ##collisionEnergy = "",
    ##isolationWindowTargetMz = "",
    ##polarity = ""
    ##msLevel = "",
)

#' @rdname rspec_to_pyspec
#'
#' @importFrom methods is
#'
#' @export
#' 
#' @examples
#' library(Spectra)
#' 
#' DF <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     name = c("Caffeine", "Caffeine", "1-Methylhistidine"),
#'     precursorMz = c(195.0877, 195.0877, 170.0924)
#' )
#' DF$intensity <- list(
#'     c(340.0, 416, 2580, 412),
#'     c(388.0, 3270, 85, 54, 10111),
#'     c(3.407, 47.494, 3.094, 100.0, 13.240))
#' DF$mz <- list(
#'     c(135.0432, 138.0632, 163.0375, 195.0880),
#'     c(110.0710, 138.0655, 138.1057, 138.1742, 195.0864),
#'     c(109.2, 124.2, 124.5, 170.16, 170.52))
#' sps <- Spectra(DF)
#'
#' ## apply the function rspec_to_pyspec on sps
#' rspec_to_pyspec(x = sps, pythonLibrary = "matchms")
#' ## gives [Spectrum(precursor m/z=195.09, 4 fragments between 135.0 and 195.1), Spectrum(precursor m/z=195.09, 5 fragments between 110.1 and 195.1), Spectrum(precursor m/z=170.09, 5 fragments between 109.2 and 170.5)]
#' 
#' rspec_to_pyspec(x = sps, pythonLibrary = "spectrum_utils")
#' ## gives [<spectrum_utils.spectrum.MsmsSpectrum>, <spectrum_utils.spectrum.MsmsSpectrum>, <spectrum_utils.spectrum.MsmsSpectrum>]
#'
rspec_to_pyspec <- function(x, pythonLibrary = c("matchms", "spectrum_utils"), ##mapping = spectraVariableMapping(),
    ##reference = import("matchms"), ## set import(, convert = TRUE/FALSE)???
    BPPARAM = SerialParam(), .check = TRUE, ...) {
    
    ## check arguments
    if (.check && !is(x, "Spectra"))
        stop("'x' should be a Spectra object.")
    
    pythonLibrary <- match.arg(pythonLibrary)
    
    ## apply on each spectra of x the conversion from R to Python
    plist <- spectrapply(object = x, 
        FUN = function(x_i) .single_rspec_to_pyspec(x = x_i, 
            pythonLibrary = pythonLibrary, ...), 
        BPPARAM = BPPARAM)
    r_to_py(unname(plist))
}

#' @rdname rspec_to_pyspec
#'
#' @importFrom Spectra concatenateSpectra
#'
#' @export
#' 
#' @examples 
#' library(Spectra)
#' 
#' DF <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     name = c("Caffeine", "Caffeine", "1-Methylhistidine"),
#'     precursorMz = c(195.0877, 195.0877, 170.0924)
#' )
#' DF$intensity <- list(
#'     c(340.0, 416, 2580, 412),
#'     c(388.0, 3270, 85, 54, 10111),
#'     c(3.407, 47.494, 3.094, 100.0, 13.240))
#' DF$mz <- list(
#'     c(135.0432, 138.0632, 163.0375, 195.0880),
#'     c(110.0710, 138.0655, 138.1057, 138.1742, 195.0864),
#'     c(109.2, 124.2, 124.5, 170.16, 170.52))
#' sps <- Spectra(DF)
#'
#' ## python library: matchms
#' ## apply the function rspec_to_pyspec on sps
#' py_obj <- rspec_to_pyspec(x = sps, pythonLibrary = "matchms")
#' pyspec_to_rspec(py_obj, pythonLibrary = "matchms")
#' 
#' ## python library: spectrum_utils
#' py_obj <- rspec_to_pyspec(x = sps, pythonLibrary = "spectrum_utils")
#' pyspec_to_rspec(py_obj, pythonLibrary = "spectrum_utils")
#' 
pyspec_to_rspec <- function(x, 
        BPPARAM = SerialParam(), .check = TRUE, 
        pythonLibrary = c("matchms", "spectrum_utils"), ...) {
    
    ## check the arguments
    pythonLibrary <- match.arg(pythonLibrary)
    
    ## obtain the mapping
    args <- list(...)
    if (!("spectraVariables" %in% names(args)))
        spectraVars <- spectraVariableMapping(pythonLibrary = pythonLibrary)
    if ("spectraVariables" %in% names(args))
        spectraVars <- args[["spectraVariables"]]
    
    if (!(is(x, "list") | is(x, "python.builtin.list")))
      stop("'x' is expected to be a Python list.")
    x <- py_to_r(x)
    
    ## check for matchms
    if (pythonLibrary == "matchms" && .check  && !all(vapply(x, function(z)
        is(z, "matchms.Spectrum.Spectrum"), logical(1))))
        stop("'x' is expected to be a Python list of matchms Spectrum objects.")
    
    if (pythonLibrary == "spectrum_utils" && .check && !all(vapply(x, function(z)
        is(z, "spectrum_utils.spectrum.MsmsSpectrum"), logical(1))))
        stop("'x' is expected to be a Python list of spectrum_utils MsmsSpectrum objects.")
    
    ## apply the Python to R Spectra conversion for each element of x and 
    ## combine the list entries
    
    spectra_list <- bplapply(X = x, 
        FUN = function(x_i) .single_pyspec_to_rspec(x = x_i,
            pythonLibrary = pythonLibrary, spectraVariables = spectraVars), 
        BPPARAM = BPPARAM)
    do.call(concatenateSpectra, spectra_list)    
}

#' @description
#'
#' Function to convert a **single** R Spectra object (of length 1) into a
#' Python matchms Spectrum using the `reticulate` package.
#'
#' @details
#' Converting to a `spectrum_utils`' `MsmsSpectrum` object might change the 
#' `m/z` values (approximately by +- 1e-06).
#' 
#' @param x `Spectra` object **of length 1!**.
#'
#' @param spectraVariables named `character` vector defining the spectra
#'     varibles that should be stored as metadata in `matchms`' metadata. Names
#'     are expected to be the spectra variable names and values the
#'     corresponding metadata fields in `matchms`. Defaults to
#'     [spectraVariableMapping()]. If `spectraVariables = character()` no
#'     metadata will be stored.
#'
#' @param reference Reference to Python environment matchms
#'
#' @return `Spectrum` Single Python Spectrum
#'
#' @author Michael Witting, Johannes Rainer
#'
#' @importMethodsFrom Spectra spectraData
#'
#' @importMethodsFrom Spectra peaksData
#'
#' @importMethodsFrom Spectra mz
#'
#' @importMethodsFrom Spectra intensity
#'
#' @importFrom reticulate np_array r_to_py
#' @importFrom Spectra precursorMz precursorCharge
#' 
#' @examples
#' library(Spectra)
#' 
#' DF <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     name = c("Caffeine", "Caffeine", "1-Methylhistidine"),
#'     precursorMz = c(195.0877, 195.0877, 170.0924)
#' )
#' DF$intensity <- list(
#'     c(340.0, 416, 2580, 412),
#'     c(388.0, 3270, 85, 54, 10111),
#'     c(3.407, 47.494, 3.094, 100.0, 13.240))
#' DF$mz <- list(
#'     c(135.0432, 138.0632, 163.0375, 195.0880),
#'     c(110.0710, 138.0655, 138.1057, 138.1742, 195.0864),
#'     c(109.2, 124.2, 124.5, 170.16, 170.52))
#' sps <- Spectra(DF)
#'
#' .single_rspec_to_pyspec(x = sps[1, ], pythonLibrary = "matchms")
#' ## gives Spectrum(precursor m/z=195.09, 4 fragments between 135.0 and 195.1)
#' .single_rspec_to_pyspec(x = sps[1, ], pythonLibrary = "spectrum_utils")
#' ## gives <spectrum_utils.spectrum.MsmsSpectrum object at ...>
#' 
#' @noRd
.single_rspec_to_pyspec <- function(x, pythonLibrary = c("matchms", "spectrum_utils"), ...) { 
    
    ## match argument pythonLibrary
    pythonLibrary <- match.arg(pythonLibrary)
    
    ## obtain the mapping
    args <- list(...)
    if (!("spectraVariables" %in% names(args)))
        spectraVars <- spectraVariableMapping(pythonLibrary = pythonLibrary)
    if ("spectraVariables" %in% names(args))
        spectraVars <- args[["spectraVariables"]]
    
    ## run import() to import the specified Python library, making it available
    ## for use from R
    reference <- import(module = pythonLibrary, convert = TRUE) ## set import(, convert = TRUE/FALSE)???
    
    pks <- unname(peaksData(x, c("mz", "intensity")))[[1L]]
    
    if (pythonLibrary == "matchms") {
        if (length(spectraVars)) {
            slist <- as.list(spectraData(x, columns = names(spectraVars)))
            ## ## Seems matchms.Spectrum does not support NA retention times?
            ## if (any(names(slist) == "rtime") && is.na(slist$rtime))
            ##     slist$rtime <- 0
            names(slist) <- spectraVars
            res <- reference$Spectrum(mz = np_array(pks[, 1L]),
                intensities = np_array(pks[, 2L]),
                metadata = r_to_py(slist))
        } else 
            res <- reference$Spectrum(mz = np_array(pks[, 1L]),
                intensities = np_array(pks[, 2L]))    
    }
    
    if (pythonLibrary == "spectrum_utils") {
        res <- reference$spectrum$MsmsSpectrum(mz = np_array(pks[, 1L]),
            intensity = np_array(pks[, 2L]), 
            identifier = "Spectra_to_MsmsSpectrum", 
            precursor_mz = precursorMz(x), 
            precursor_charge = precursorCharge(x))
        
        ## add spectraVariables
        slist <- as.list(spectraData(x, columns = names(spectraVars)))
        names(slist) <- spectraVars
        for (i in seq_along(slist))
            res[[names(slist)[i]]] <- slist[i]
    }
    
    ## return the Python object
    res
    
}

#' @description
#'
#' Function to convert a single Python Spectrum object into an R Spectra using
#' the `reticulate` package.
#' 
#' @details
#' Converting to a `spectrum_utils`' `MsmsSpectrum` object might change the 
#' `m/z` values (approximately by +- 1e-06). 
#'
#' @param x `Spectrum` Single Python Spectrum.
#'
#' @param spectraVariables named `character` vector with the names of the
#'     spectra variables that should be extracted. Names are expected to be
#'     the spectra variable names (in `Spectra`) and values the corresponding
#'     names of the variables within the Python Spectrum.
#'
#' @param reference Reference to Python environment matchms
#'
#' @return `Spectra` single R Spectra
#'
#' @author Michael Witting, Johannes Rainer
#'
#' @importFrom IRanges NumericList
#'
#' @importFrom S4Vectors DataFrame
#'
#' @importMethodsFrom Spectra Spectra
#'
#' @noRd
#' 
#' @examples
#' library(Spectra)
#' 
#' DF <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     name = c("Caffeine", "Caffeine", "1-Methylhistidine"),
#'     precursorMz = c(195.0877, 195.0877, 170.0924)
#' )
#' DF$intensity <- list(
#'     c(340.0, 416, 2580, 412),
#'     c(388.0, 3270, 85, 54, 10111),
#'     c(3.407, 47.494, 3.094, 100.0, 13.240))
#' DF$mz <- list(
#'     c(135.0432, 138.0632, 163.0375, 195.0880),
#'     c(110.0710, 138.0655, 138.1057, 138.1742, 195.0864),
#'     c(109.2, 124.2, 124.5, 170.16, 170.52))
#' sps <- Spectra(DF)
#'
#' ## python library: matchms
#' ## apply the function rspec_to_pyspec on sps
#' py_obj <- rspec_to_pyspec(x = sps, pythonLibrary = "matchms")
#' .single_pyspec_to_rspec(py_to_r(py_obj[0]), pythonLibrary = "matchms")
#' 
#' ## python library: spectrum_utils
#' py_obj <- rspec_to_pyspec(x = sps, pythonLibrary = "spectrum_utils")
#' .single_pyspec_to_rspec(py_obj[0], pythonLibrary = "spectrum_utils")
#' 
.single_pyspec_to_rspec <- function(x, 
    ##spectraVariables = spectraVariableMapping(), 
    pythonLibrary = c("matchms", "spectrum_utils"), ...) {
        
    ## match argument pythonLibrary
    pythonLibrary <- match.arg(pythonLibrary)
    
    ## obtain the mapping
    args <- list(...)
    if (!("spectraVariables" %in% names(args)))
        spectraVars <- spectraVariableMapping(pythonLibrary = pythonLibrary)
    if ("spectraVariables" %in% names(args))
        spectraVars <- args[["spectraVariables"]]
    
    if (pythonLibrary == "matchms") {
        
        plist <- x$metadata
        vars <- spectraVars[spectraVars %in% names(plist)]
        if (length(vars)) {
            rlist <- lapply(vars, function(z) plist[z])
            ## Drop NULL variables.
            spd <- DataFrame(rlist[lengths(rlist) > 0])
            if (!nrow(spd))
                spd <- DataFrame(msLevel = NA_integer_)
        } else
            spd <- DataFrame(msLevel = NA_integer_)
        spd$mz <- NumericList(as.numeric(x$peaks$mz), compress = FALSE)
        spd$intensity <- NumericList(as.numeric(x$peaks$intensities),
                                     compress = FALSE)
    }
    
    if (pythonLibrary == "spectrum_utils") {
        plist <- names(x)
        vars <- spectraVars[spectraVars %in% plist]## vars <- plist[plist %in% spectraVars]
        if (length(vars)) {
            rlist <- lapply(vars, function(z) unlist(unname(x[[z]])))
            ## Drop NULL variables.
            spd <- DataFrame(rlist[lengths(rlist) > 0])
            if (!nrow(spd))
                spd <- DataFrame(msLevel = NA_integer_)
        } else
            spd <- DataFrame(msLevel = NA_integer_)
        spd$mz <- NumericList(as.numeric(x$mz), compress = FALSE)
        spd$intensity <- NumericList(as.numeric(x$intensity), compress = FALSE)
    }
    
    Spectra(spd)
}

#' Extract all spectraData and all mz and intensity values, give them to
#' Python to create an array of Spectrum. Could be faster because loop is
#' performed in Python rather than in R.
#'
#' @noRd
.multi_rspec_to_pyspec <- function() {
}

#' Extract data in Python for all elements: python function should return a
#' list with the metadata (convert to DataFrame) and the m/z and intensity
#' values. Would avoid loop in R alltogether.
#'
#' @noRd
.multi_pyspec_to_rspec <- function() {
}

