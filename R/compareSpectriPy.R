#' @title Spectra similarity calculations using Python's matchms library
#'
#' @name compareSpectriPy
#'
#' @description
#'
#' The `compareSpectriPy()` function allows to calculate spectral similarity
#' scores using the `calculate_scores()` function of the Python
#' [matchms.similarity](https://matchms.readthedocs.io/en/latest/api/matchms.similarity.html).
#' module.
#'
#' Selection and configuration of the algorithm can be performed with one of the
#' *parameter* objects/functions:
#'
#' - `CosineGreedy()`: calculate the *cosine similarity score* between
#'   spectra. The score is calculated by finding the best possible matches
#'   between peaks of two spectra. Two peaks are considered a potential match if
#'   their m/z ratios lie within the given `tolerance`. The underlying peak
#'   assignment problem is here solved in a *greedy* way. This can perform
#'   notably faster, but does occasionally deviate slightly from a fully correct
#'   solution (as with the `CosineHungarian` algorithm). In practice this
#'   will rarely affect similarity scores notably, in particular for smaller
#'   tolerances. The algorithm can be configured with parameters `tolerance`,
#'   `mz_power` and `intensity_power` (see parameter description for more
#'   details). See also
#'   [matchms CosineGreedy](https://matchms.readthedocs.io/en/latest/api/matchms.similarity.CosineGreedy.html) for more information.
#'
#' - `CosineHungarian()`: calculate the *cosine similarity score* as with
#'   `CosineGreedy`, but using the Hungarian algorithm to find the best
#'   matching peaks between the compared spectra. The algorithm can be
#'   configured with parameters `tolerance`, `mz_power` and `intensity_power`
#'   (see parameter description for more details). See also
#'   [matchms CosingHungarian](https://matchms.readthedocs.io/en/latest/api/matchms.similarity.CosineHungarian.html) for more information.
#'
#' - `ModifiedCosine()`: The modified cosine score aims at quantifying the
#'   similarity between two mass spectra. The score is calculated by finding
#'   the best possible matches between peaks of two spectra. Two peaks are
#'   considered a potential match if their m/z ratios lie within the given
#'   `tolerance`, or if their m/z ratios lie within the tolerance once a
#'   mass shift is applied. The mass shift is simply the difference in
#'   precursor-m/z between the two spectra.
#'   See also [matchms ModifiedCosine](https://matchms.readthedocs.io/en/latest/api/matchms.similarity.ModifiedCosine.html) for more information.
#'
#' - `NeutralLossesCosine()`: The neutral losses cosine score aims at
#'   quantifying the similarity between two mass spectra. The score is
#'   calculated by finding the best possible matches between peaks of two
#'   spectra. Two peaks are considered a potential match if their m/z ratios lie
#'   within the given `tolerance` once a mass shift is applied. The mass shift
#'   is the difference in precursor-m/z between the two spectra. See also
#'   [matchms NeutralLossesCosine](https://matchms.readthedocs.io/en/latest/api/matchms.similarity.NeutralLossesCosine.html) for more information.
#'
#' - `FingerprintSimilarity()`: Calculate similarity between molecules based
#'   on their fingerprints. For this similarity measure to work, fingerprints
#'   are expected to be derived by running *add_fingerprint()*. See also
#'   [matchms FingerprintSimilarity](https://matchms.readthedocs.io/en/latest/api/matchms.similarity.FingerprintSimilarity.html) for more information.
#'
#'
#' @note
#'
#' Parameters and algorithms are named as originally defined in the *matchms*
#' library (i.e. all parameters in *snake_case* while *CamelCase* is used
#' for the algorithms.
#'
#' @param x A [Spectra::Spectra()] object.
#'
#' @param y A [Spectra::Spectra()] object to compare against. If missing,
#'   spectra similarities are calculated between all spectra in `x`.
#'
#' @param param One of the parameter classes listed above (such as
#'   `CosineGreedy`) defining the similarity scoring function in Python
#'   and its parameters.
#'
#' @param tolerance `numeric(1)`: tolerated differences in the peaks' m/z. Peaks
#'   with m/z differences `<= tolerance` are considered matching.
#'
#' @param mz_power `numeric(1)`: the power to raise m/z to in the cosine
#'   function. The default is 0, in which case the peak intensity products will
#'   not depend on the m/z ratios.
#'
#' @param ignore_peaks_above_precursor For `NeutralLossesCosine()`:
#'   `logical(1)`: if `TRUE` (the default), peaks with m/z values larger than
#'   the precursor m/z are ignored.
#'
#' @param intensity_power `numeric(1)`: the power to raise intensity to in the
#'   cosine function. The default is 1.
#'
#' @param ... ignored.
#'
#' @return `compareSpectriPy()` Returns a `numeric` matrix with the scores,
#'   with the number of rows equal to `length(x)` and the number of columns
#'   equal to `length(y)`.
#'
#' @author Carolin Huber, Michael Witting, Johannes Rainer, Helge Hecht,
#'   Marilyn De Graeve
#'
#' @seealso [Spectra::compareSpectra()] in the *Spectra* package for pure R
#'   implementations of spectra similarity calculations.
#'
#' @export
#'
#' @examples
#'
#' library(Spectra)
#' ## Create some example Spectra.
#' DF <- DataFrame(
#'     msLevel = c(2L, 2L, 2L),
#'     name = c("Caffeine", "Caffeine", "1-Methylhistidine"),
#'     precursorMz = c(195.0877, 195.0877, 170.0924)
#' )
#' DF$intensity <- list(
#'     c(340.0, 416, 2580, 412),
#'     c(388.0, 3270, 85, 54, 10111),
#'     c(3.407, 47.494, 3.094, 100.0, 13.240)
#' )
#' DF$mz <- list(
#'     c(135.0432, 138.0632, 163.0375, 195.0880),
#'     c(110.0710, 138.0655, 138.1057, 138.1742, 195.0864),
#'     c(109.2, 124.2, 124.5, 170.16, 170.52)
#' )
#' sps <- Spectra(DF)
#'
#' ## Calculate pairwise similarity beween all spectra within sps with
#' ## matchms' CosineGreedy algorithm
#' ## Note: the first compareSpectriPy will take longer because the Python
#' ## environment needs to be set up.
#' res <- compareSpectriPy(sps, param = CosineGreedy())
#' res
#'
#' ## Next we calculate similarities for all spectra against the first one
#' res <- compareSpectriPy(sps, sps[1], param = CosineGreedy())
#'
#' ## Calculate pairwise similarity of all spectra in sps with matchms'
#' ## ModifiedCosine algorithm
#' res <- compareSpectriPy(sps, param = ModifiedCosine())
#' res
#'
#' ## Note that the ModifiedCosine method requires the precursor m/z to be
#' ## known for all input spectra. Thus, it is advisable to remove spectra
#' ## without precursor m/z before using this algorithm.
#' sps <- sps[!is.na(precursorMz(sps))]
#' compareSpectriPy(sps, param = ModifiedCosine())
NULL

setGeneric("compareSpectriPy", function(x, y, param, ...) {
    standardGeneric("compareSpectriPy")
})

setClass(
    "CosineGreedy",
    slots = c(
        tolerance = "numeric", mzPower = "numeric", intensityPower = "numeric"
    ),
    prototype = prototype(tolerance = 0.1, mzPower = 0.0, intensityPower = 1.0),
    validity = function(object) {
        msg <- NULL
        if (length(object@tolerance) != 1 || object@tolerance < 0) {
            msg <- c("'tolerance' has to be a positive number of length 1")
        }
        if (length(object@mzPower) != 1) {
            msg <- c(msg, "'mz_power' has to be a number of length 1")
        }
        if (length(object@intensityPower) != 1) {
            msg <- c(msg, "'intensity_power' has to be a number of length 1")
        }
        msg
    }
)
setClass(
    "CosineHungarian",
    contains = "CosineGreedy")
setClass(
    "ModifiedCosine",
    contains = "CosineGreedy")
setClass(
    "NeutralLossesCosine",
    contains = "CosineGreedy",
    slots = c(ignorePeaksAbovePrecursor = "logical"),
    prototype = prototype(ignorePeaksAbovePrecursor = TRUE),
    validity = function(object) {
        msg <- NULL
        if (length(object@ignorePeaksAbovePrecursor) != 1) {
            msg <- paste0("'ignorePeaksAbovePrecursor' has to be a ",
                          "positive number of length 1")
        }
        msg
    }
)
setClass(
    "FingerprintSimilarity",
    contains = "CosineGreedy",
    slots = c(similarityMeasure = "character"),
    prototype = prototype(similarityMeasure = "jaccard"))

#' @rdname compareSpectriPy
#'
#' @importFrom methods new
#'
#' @export
CosineGreedy <- function(tolerance = 0.1, mz_power = 0.0,
                         intensity_power = 1.0) {
    new("CosineGreedy",
        tolerance = as.numeric(tolerance),
        mzPower = as.numeric(mz_power),
        intensityPower = as.numeric(intensity_power)
        )
}

#' @rdname compareSpectriPy
#'
#' @export
CosineHungarian <- function(tolerance = 0.1, mz_power = 0.0,
                            intensity_power = 1.0) {
    new("CosineHungarian",
        tolerance = as.numeric(tolerance),
        mzPower = as.numeric(mz_power),
        intensityPower = as.numeric(intensity_power)
    )
}

#' @rdname compareSpectriPy
#'
#' @export
ModifiedCosine <- function(tolerance = 0.1, mz_power = 0.0,
                                intensity_power = 1.0) {
    new("ModifiedCosine",
        tolerance = as.numeric(tolerance),
        mzPower = as.numeric(mz_power),
        intensityPower = as.numeric(intensity_power)
    )
}

#' @rdname compareSpectriPy
#'
#' @export
NeutralLossesCosine <- function(tolerance = 0.1, mz_power = 0.0,
                                intensity_power = 1.0,
                                ignore_peaks_above_precursor = TRUE) {
    new("NeutralLossesCosine",
        tolerance = as.numeric(tolerance),
        mzPower = as.numeric(mz_power),
        intensityPower = as.numeric(intensity_power),
        ignorePeaksAbovePrecursor = as.logical(ignore_peaks_above_precursor)
    )
}

## #' @rdname compareSpectriPy
## #'
## #' @export
## FingerprintSimilarity <- function(similarity_measure = "jaccard") {
##     new("FingerprintSimilarity",
##         similarityMeasure = similarity_measure
##         )
## }

#' @rdname compareSpectriPy
#'
#' @exportMethod compareSpectriPy
setMethod(
    "compareSpectriPy",
    signature = c(x = "Spectra", y = "Spectra", param = "CosineGreedy"),
    function(x, y, param, ...) {
        .compare_spectra_python(x, y, param)
    }
)

#' @rdname compareSpectriPy
setMethod(
    "compareSpectriPy",
    signature = c(x = "Spectra", y = "missing", param = "CosineGreedy"),
    function(x, y, param, ...) {
        .compare_spectra_python(x, y = NULL, param)
    }
)

#' Method to construct a python function name for a param object.
#'
#' @noRd
setGeneric("py_fun", function(object, ...) {
    standardGeneric("py_fun")
})
setMethod("py_fun", "CosineGreedy", function(object) {
    matchms_similarity$CosineGreedy(object@tolerance, object@mzPower,
                                    object@intensityPower)
})
setMethod("py_fun", "CosineHungarian", function(object) {
    matchms_similarity$CosineHungarian(object@tolerance, object@mzPower,
                                       object@intensityPower)
})
setMethod("py_fun", "ModifiedCosine", function(object) {
    matchms_similarity$ModifiedCosine(object@tolerance, object@mzPower,
                                      object@intensityPower)
})
setMethod("py_fun", "NeutralLossesCosine", function(object) {
    matchms_similarity$NeutralLossesCosine(object@tolerance, object@mzPower,
                                           object@intensityPower,
                                           object@ignorePeaksAbovePrecursor)
})

.fun_name <- function(x) {
    class(x)[1L]
}

#' Internal function to calculate similarities with Python's matchms. `Spectra`
#' will be converted to Python.
#'
#' @param x `Spectra`
#'
#' @param y `Spectra`
#'
#' @param param Parameter object.
#'
#' @param value `character(1)` defining which value should be returned from the
#'     Python call.
#'
#' @return a `numeric` `matrix` nrow being length of `x`, nrow length `y`.
#'
#' @noRd
#'
#' @author Carolin Huber, Johannes Rainer, Wout Bittremieux
#'
#' @importFrom reticulate r_to_py py_to_r
.compare_spectra_python <- function(x, y = NULL, param,
                                    mapping = spectraVariableMapping()) {
    ## Handle empty input.
    if (!length(x) || (!length(y) & !is.null(y))) {
        return(matrix(NA_real_, ncol = length(y), nrow = length(x)))
    }

    ## Convert R spectra to Python.
    py_x <- rspec_to_pyspec(x, mapping = mapping)
    py_y <- if (!is.null(y)) rspec_to_pyspec(y) else py_x
    is_symmetric <- is.null(y)

    ## Compute the similarity scores with matchms.
    scores <- matchms$calculate_scores(
        py_x, py_y, py_fun(param),
        is_symmetric = is_symmetric
        )
    ## Collect results
    py_to_r(scores$to_array()[paste0(.fun_name(param), "_score")])
}
