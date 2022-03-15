#' @title Spectra similarity calculations using matchms
#'
#' @name compareSpectriPy
#'
#' @description
#'
#' The `compareSpectriPy` function allows to calculate spectral similarity
#' scores using the `calculate_scores module` of the python
#' [matchms.similarity package](https://matchms.readthedocs.io/en/latest/api/matchms.similarity.html)
#' package.
#'
#' Selection and configuration of the algorithm can be performed with one of the
#' parameter objects:
#'
#' - `CosineGreedyParam`: calculate the *cosine similarity score* between
#'   spectra. The score is calculated by finding best possible matches between
#'   peaks of two spectra. Two peaks are considered a potential match if their
#'   m/z ratios lie within the given `tolerance`. The underlying peak assignment
#'   problem is here solved in a *greedy* way. This can perform notably faster,
#'   but does occasionally deviate slightly from a fully correct solution (as
#'   with the `CosineHungarianParam` algorithm). In practice this will rarely
#'   affect similarity scores notably, in particular for smaller tolerances. The
#'   algorithm can be configured with parameters `tolerance`, `mzPower` and
#'   `intensityPower` (see parameter description for more details).
#'
#' - `CosineHungarianParam`: calculate the *cosine similarity score* as with
#'   `CosineGreedyParam`, but using the Hungarian algorithm to find the best
#'   matching peaks between the compared spectra. The algorithm can be
#'   configured with parameters `tolerance`, `mzPower` and `intensityPower`
#'   (see parameter description for more details).
#'
#' - `ModifiedCosineParam`: The modified cosine score aims at quantifying the
#'   similarity between two mass spectra. The score is calculated by finding
#'   best possible matches between peaks of two spectra. Two peaks are
#'   considered a potential match if their m/z ratios lie within the given
#'   `tolerance`, or if their m/z ratios lie within the tolerance once a
#'   mass-shift is applied. The mass shift is simply the difference in
#'   precursor-m/z between the two spectra.
#'
#' - `NeutralLossesCosineParam`: The neutral losses cosine score aims at
#'   quantifying the similarity between two mass spectra. The score is
#'   calculated by finding best possible matches between peaks of two spectra.
#'   Two peaks are considered a potential match if their m/z ratios lie within
#'   the given `tolerance` once a mass-shift is applied. The mass shift is the
#'   difference in precursor-m/z between the two spectra.
#'
#' @param x A [Spectra()] object.
#'
#' @param y A [Spectra()] object to compare against. If missing spectra
#'   similarities are calculated between all spectr in `x`.
#'
#' @param tolerance `numeric(1)`: tolerated differences in peaks' m/z. Peaks
#'   with m/z differences `<= tolerance` are considered matching.
#'
#' @param mzPower `numeric(1)`: the power to raise m/z to in the cosine
#'   function. The default is 0, in which case the peak intensity products will
#'   not depend on the m/z ratios.
#'
#' @param intensityPower `numeric(1)`: the power to raise intensity to in the
#'   cosine function. The default is 1.
#'
#' @param ignorePeaksAbovePrecursor `logical(1)`: if `TRUE` (the default), peaks
#'   with m/z values larger then the precursor m/z are ignored.
#'
#' @return `compareSpectriPy` returns a `numeric` matrix with the scores,
#'   number of rows being equal to `length(x)` and number of columns equal to
#'   `length(y)`.
#'
#' @author Carolin Huber, Michael Witting, Johannes Rainer, Helge Hecht
#'
#' @seealso [compareSpectra()] in the `Spectra` package for pure R
#'     implementations of spectra similarity calculations.
#'
#' @export
#'
#' @importFrom reticulate py_run_string
#'
#' @examples
NULL

setGeneric("compareSpectriPy", function(x, y, param, ...)
    standardGeneric("compareSpectriPy"))

#' @importClassesFrom ProtGenerics Param
#'
#' @noRd
setClass("CosineGreedyParam",
         slots = c(
             tolerance = "numeric",
             mzPower = "numeric",
             intensityPower = "numeric"
         ),
         prototype = prototype(
             tolerance = 0.1,
             mzPower = 0.0,
             intensityPower = 1.0
         ),
         validity = function(object) {
             msg <- NULL
             if (length(object@tolerance) != 1 || object@tolerance < 0)
                 msg <- c("'tolerance' has to be a positive number of length 1")
             if (length(object@mzPower) != 1)
                 msg <- c(msg, "'mzPower' has to be a number of length 1")
             if (length(object@intensityPower) != 1)
                 msg <- c(msg,
                          "'intensityPower' has to be a number of length 1")
             msg
         })
setClass("CosineHungarianParam",
         contains = "CosineGreedyParam")
setClass("ModifiedCosineParam",
         contains = "CosineGreedyParam")
setClass("NeutralLossesCosineParam",
         contains = "CosineGreedyParam",
         slots = c(ignorePeaksAbovePrecursor = "logical"),
         prototype = prototype(ignorePeaksAbovePrecursor = TRUE),
         validity = function(object) {
             msg <- NULL
             if (length(object@ignorePeaksAbovePrecursor) != 1)
                 msg <- paste0("'ignorePeaksAbovePrecursor' has to be a ",
                               "positive number of length 1")
             msg
         })

#' @rdname compareSpectriPy
#'
#' @importFrom methods new
#'
#' @export
CosineGreedyParam <- function(tolerance = 0.1, mzPower = 0.0,
                              intensityPower = 1.0) {
    new("CosineGreedyParam", tolerance = as.numeric(tolerance),
        mzPower = as.numeric(mzPower),
        intensityPower = as.numeric(intensityPower))
}

#' @rdname compareSpectriPy
#'
#' @export
CosineHungarianParam <- function(tolerance = 0.1, mzPower = 0.0,
                                 intensityPower = 1.0) {
    new("CosineHungarianParam", tolerance = as.numeric(tolerance),
        mzPower = as.numeric(mzPower),
        intensityPower = as.numeric(intensityPower))
}

#' @rdname compareSpectriPy
#'
#' @export
ModifiedCosineParam <- function(tolerance = 0.1, mzPower = 0.0,
                                 intensityPower = 1.0) {
    new("ModifiedCosineParam", tolerance = as.numeric(tolerance),
        mzPower = as.numeric(mzPower),
        intensityPower = as.numeric(intensityPower))
}

#' @rdname compareSpectriPy
#'
#' @export
NeutralLossesCosineParam <- function(tolerance = 0.1, mzPower = 0.0,
                                     intensityPower = 1.0,
                                     ignorePeaksAbovePrecursor = TRUE) {
    new("NeutralLossesCosineParam", tolerance = as.numeric(tolerance),
        mzPower = as.numeric(mzPower),
        intensityPower = as.numeric(intensityPower),
        ignorePeaksAbovePrecursor = as.logical(ignorePeaksAbovePrecursor))
}

## #' @rdname compareSpectriPy
## #'
## #' @exportMethod compareSpectriPy
## setMethod(
##     "compareSpectriPy",
##     signature = c(x = "Spectra", y = "Spectra", param = "CosineGreedyParam"),
##     function(x, y, param, ...) {
##         ## Simply pass things along.
##     })
## setMethod(
##     "compareSpectriPy",
##     signature = c(x = "Spectra", y = "missing", param = "CosineGreedyParam"),
##     function(x, y, param, ...) {
##         ## Simply pass things along setting y = NULL.
##     })

#' helper function to extract cosine parameter settings for similarity functions
#' based on cosine.
#'
#' @noRd
.cosine_param_string <- function(x) {
    paste0("tolerance=", x@tolerance, ", mz_power=", x@mzPower,
           ", intensity_power=", x@intensityPower)
}
#' Could also define a method, but I guess that's overkill in this case.
#'
#' @noRd
.fun_name <- function(x) {
    sub("Param$", "", class(x)[1L])
}

#' (internal) helper method to build the python command to perform the
#' comparison. Each parameter class could (if needed) it's own implementation
#' to create the string. This methods will be called in the
#' `compare_spectra_python` function.
#'
#' @noRd
setGeneric("python_command", function(object, ...)
    standardGeneric("python_command"))
setMethod(
    "python_command",
    "CosineGreedyParam",
    function(object, input_param = "py_x, py_y", is_symmetric = "False") {
        FUN <- .fun_name(object)
        paste0("import matchms\n",
               "from matchms.similarity import ", FUN, "\n",
               "res = matchms.calculate_scores(", input_param,
               ", ", FUN, "(", .cosine_param_string(object), "), ",
               "is_symmetric=", is_symmetric, ")")
    })
setMethod(
    "python_command",
    "NeutralLossesCosineParam",
    function(object, input_param = "py_x, py_y", is_symmetric = "False") {
        FUN <- .fun_name(object)
        paste0("import matchms\n",
               "from matchms.similarity import ", FUN, "\n",
               "res = matchms.calculate_scores(", input_param,
               ", ", FUN, "(", .cosine_param_string(object), ", ",
               "ignore_peaks_above_precursor=",
               ifelse(object@ignorePeaksAbovePrecursor,
                      yes = "True", no = "False"),
               "), is_symmetric=", is_symmetric, ")")
    })

#' internal function to calculate similarities with python's matchms. `Spectra`
#' will be converted to python.
#'
#' @param x `Spectra`
#'
#' @param y `Spectra`
#'
#' @param param Parameter object.
#'
#' @return a `numeric` `matrix` nrow being length of `x`, nrow length `y`.
#'
#' @importFrom basilisk basiliskStart basiliskRun
#'
#' @noRd
#'
#' @author Carolin Huber, Johannes Rainer
.compare_spectra_python <- function(x, y = NULL, param, value = "score") {
    ## handle empty input
    if (!length(x) || (!length(y) & !is.null(y)))
        return(matrix(NA_real_, ncol = length(y), nrow = length(x)))

    cl <- basiliskStart(matchms_env)
    on.exit(basiliskStop(cl))

    basiliskRun(cl, function(x, y, param) {
        ref <- import("matchms")
        is_symmetric <- "False"
        py$py_x <- rspec_to_pyspec(x, reference = ref)
        if (is.null(y)) {
            py$py_y <- py$py_x
            is_symmetric <- "True"
        } else
            py$py_y <- rspec_to_pyspec(y, reference = ref)
        com <- python_command(param, is_symmetric = is_symmetric)
        ## Run the command. Result is in py$res
        py_run_string(com)
        ## Collect the results.
        py_run_string(
            paste0("sim = []\n",
                   "for x,y,z in res:\n  sim.append(z['", value, "'])"))
        matrix(unlist(py$sim, use.names = FALSE),
               nrow = length(x), byrow = TRUE)
    }, x = x, y = y, param = param)
}
