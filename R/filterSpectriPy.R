## #' @title Filter Spectra using matchms
## #'
## #' @name filterSpectriPy
## #'
## #' @description
## #' The `filterSpectriPy()` function allows to filter/process a `Spectra` object
## #' using the `select_by_intensity`, `select_by_mz`,
## #' `remove_peaks_around_precursor_mz`, and `normalize_intensities` of the python
## #' [matchms.filtering](https://matchms.readthedocs.io/en/latest/api/matchms.filtering.html)
## #' module.
## #'
## #' Selection and configuration of the algorithm can be performed with one of the
## #' parameter objects (equivalent to `matchms`' function names):
## #'
## #' - `select_by_intensity`: Keeps only the peaks within defined intensity range
## #' (keep if `intensity_from` >= intensity >= `intensity_to`).
## #'
## #' - `select_by_mz`: Keeps only the peaks between `mz_from` and `mz_to`
## #' (keep if `mz_from` >= m/z >= `mz_to`).
## #'
## #' - `remove_peaks_around_precursor_mz`: Removes the peaks that are within
## #' `mz_tolerance` (in Da) of the precursor mz, exlcuding the precursor peak.
## #'
## #' - `normalize_intensities`: Normalizes the intensities of peaks
## #' (and losses) to unit height.
## #'
## #' @param sps A [Spectra::Spectra()] object.
## #'
## #' @param param one of parameter classes listed above (such as
## #'   `select_by_intensity`) defining the filter/processing function in python
## #'   and its parameters.
## #'
## #' @param ... ignored.
## #'
## #' @return `filterSpectriPy()` returns a `Spectra` object on which the
## #' filtering/processing function has been applied
## #'
## #' @author Thomas Naake
## #'
## #' @seealso [Spectra::filterIntensity()], [Spectra::filterMzRange()],
## #' [Spectra::scalePeaks()] in the `Spectra` package for pure R
## #' implementations of filtering/processing calculations.
## #'
## #' @export
## #'
## #' @importFrom reticulate py_run_string
## #'
## #' @examples
## #' library(Spectra)
## #' ## create some example Spectra
## #' DF <- DataFrame(
## #'     msLevel = c(2L, 2L, 2L),
## #'     name = c("Caffeine", "Caffeine", "1-Methylhistidine"),
## #'     precursorMz = c(195.0877, 195.0877, 170.0924)
## #' )
## #' DF$intensity <- list(
## #'     c(340.0, 416, 2580, 412),
## #'     c(388.0, 3270, 85, 54, 10111),
## #'     c(3.407, 47.494, 3.094, 100.0, 13.240))
## #' DF$mz <- list(
## #'     c(135.0432, 138.0632, 163.0375, 195.0880),
## #'     c(110.0710, 138.0655, 138.1057, 138.1742, 195.0864),
## #'     c(109.2, 124.2, 124.5, 170.16, 170.52))
## #' sps <- Spectra(DF)
## #'
## #' ## process Spectra with matchms' select_by_intensity algorithm
## #' ## note: the first filterSpectriPy will take longer because the Python
## #' ## environment needs to be set up.
## #' filterSpectriPy(sps, param = select_by_intensity(intensity_from=50, intensity_to=400))
## #'
## #' ## Process Spectra with matchms' select_by_mz algorithm
## #' filterSpectriPy(sps, param = select_by_mz(mz_from=150, mz_to=450))
## #'
## #' ## Calculate pairwise similarity of all spectra in sps with matchms'
## #' ## remove_peaks_around_precursor_mz algorithm
## #' filterSpectriPy(sps, param = remove_peaks_around_precursor_mz(mz_tolerance=20))
## #'
## #' ## Calculate pairwise similarity of all spectra in sps with matchms'
## #' ## normalize_intensities algorithm
## #' filterSpectriPy(sps, normalize_intensities())
## NULL

## setGeneric("filterSpectriPy", function(sps, param, ...)
##     standardGeneric("filterSpectriPy"))

## #' @importClassesFrom ProtGenerics Param
## #'
## #' @noRd
## setClass("select_by_intensity",
##     slots = c(
##         intensity_from = "numeric", intensity_to = "numeric"),
##     prototype = prototype(
##              intensity_from = 20, intensity_to = 200),
##     validity = function(object) {
##         msg <- NULL
##         if (length(object@intensity_from) != 1 || object@intensity_from < 0)
##             msg <- c("'intensity_from' has to be a positive number of length 1")
##         if (length(object@intensity_to) != 1 || object@intensity_to < 0)
##             msg <- c("'intensity_to' has to be a positive number of length 1")
##         msg
##     }
## )
## setClass("select_by_mz",
##     slots = c(
##         mz_from = "numeric",
##         mz_to = "numeric"),
##     prototype = prototype(
##         mz_from = 150,
##         mz_to = 450),
##     validity = function(object) {
##         msg <- NULL
##         if (length(object@mz_from) != 1 || object@mz_from < 0)
##             msg <- c("'mz_from' has to be a positive number of length 1")
##         if (length(object@mz_to) != 1 || object@mz_to < 0)
##             msg <- c("'mz_to' has to be a positive number of length 1")
##         msg
##     }
## )
## setClass("remove_peaks_around_precursor_mz",
##     slots = c(
##         mz_tolerance = "numeric"),
##     prototype = prototype(
##             mz_tolerance = 20),
##     validity = function(object) {
##         msg <- NULL
##         if (length(object@mz_tolerance) != 1 || object@mz_tolerance < 0)
##             msg <- c("'mz_tolerance' has to be a positive number of length 1")
##         msg
##     }
## )
## setClass("normalize_intensities",
##     prototype = prototype(),
##     validity = function(object) {
##         msg <- NULL
##         msg
##     }
## )

## #' @rdname filterSpectriPy
## #'
## #' @param intensity_from `numeric(1)`: Set lower threshold for peak intensity.
## #' Default is 10.
## #'
## #' @param intensity_to `numeric(1)`: Set upper threshold for peak intensity.
## #' Default is 200.
## #'
## #' @importFrom methods new
## #'
## #' @export
## select_by_intensity <- function(intensity_from = 10, intensity_to = 200) {
##     new("select_by_intensity", intensity_from = as.numeric(intensity_from),
##         intensity_to = as.numeric(intensity_to))
## }

## #' @rdname filterSpectriPy
## #'
## #' @param mz_from `numeric(1)`: Set lower threshold for m/z peak positions.
## #' Default is 0.
## #'
## #' @param mz_to `numeric(1)`: Set upper threshold for m/z peak positions.
## #' Default is 1000.
## #'
## #' @export
## select_by_mz <- function(mz_from = 0, mz_to = 1000) {
##     new("select_by_mz", mz_from = as.numeric(mz_from),
##         mz_to = as.numeric(mz_to))
## }

## #' @rdname filterSpectriPy
## #'
## #' @param mz_tolerance `numeric(1)`: Tolerance of m/z values that are not
## #' allowed to lie within the precursor mz. Default is 17 Da.
## #'
## #' @export
## remove_peaks_around_precursor_mz <- function(mz_tolerance = 17) {
##     new("remove_peaks_around_precursor_mz",
##         mz_tolerance = as.numeric(mz_tolerance))
## }

## #' @rdname filterSpectriPy
## #'
## #' @export
## normalize_intensities <- function() {
##     new("normalize_intensities")
## }

## #' @rdname filterSpectriPy
## #'
## #' @exportMethod filterSpectriPy
## setMethod(
##     "filterSpectriPy",
##     signature = c(sps = "Spectra", param = "select_by_intensity"),
##     function(sps, param, ...) {
##         .filter_spectra_python(sps, param)
##     })

## #' @rdname filterSpectriPy
## #'
## #' @exportMethod filterSpectriPy
## setMethod(
##     "filterSpectriPy",
##     signature = c(sps = "Spectra", param = "select_by_mz"),
##     function(sps, param, ...) {
##         .filter_spectra_python(sps, param)
##     })

## #' @rdname filterSpectriPy
## #'
## #' @exportMethod filterSpectriPy
## setMethod(
##     "filterSpectriPy",
##     signature = c(sps = "Spectra", param = "remove_peaks_around_precursor_mz"),
##     function(sps, param, ...) {
##         .filter_spectra_python(sps, param)
##     })

## #' @rdname filterSpectriPy
## #'
## #' @exportMethod filterSpectriPy
## setMethod(
##     "filterSpectriPy",
##     signature = c(sps = "Spectra", param = "normalize_intensities"),
##     function(sps, param, ...) {
##         .filter_spectra_python(sps, param)
##     })


## #' helper function to extract parameter settings for filtering/processing
## #' functions.
## #'
## #' @noRd
## .select_by_intensity_param_string <- function(x) {
##     paste0("intensity_from=", x@intensity_from, ", intensity_to=", x@intensity_to)
## }
## .select_by_mz_param_string <- function(x) {
##     paste0("mz_from=", x@mz_from, ", mz_to=", x@mz_to)
## }
## .remove_peaks_around_precursor_mz_param_string <- function(x) {
##     paste0("mz_tolerance=", x@mz_tolerance)
## }
## .normalize_intensities_param_string <- function(x) {
##     paste0()
## }

## #' (internal) helper method to build the python command to perform the
## #' filtering/processing. Each parameter class could (if needed) it's own implementation
## #' to create the string. This methods will be called in the
## #' `filter_spectra_python` function.
## #'
## #' Generic "python_command" defined in `compareSpectriPy.R`
## #' @noRd
## setGeneric("python_command", function(object, ...)
##     standardGeneric("python_command"))
## setMethod(
##     "python_command",
##     "select_by_intensity",
##     function(object, input_param = "py_spectrum_in") {
##         FUN <- .fun_name(object)
##         paste0("import matchms\n",
##             "from matchms.filtering import ", FUN, "\n",
##             "res = [", FUN, "(s, ", .select_by_intensity_param_string(object), ") for s in ", input_param, "]\n")
##     })
## setMethod(
##     "python_command",
##     "select_by_mz",
##     function(object, input_param = "py_spectrum_in") {
##         FUN <- .fun_name(object)
##         paste0("import matchms\n",
##             "from matchms.filtering import ", FUN, "\n",
##             "res = [", FUN, "(s, ", .select_by_mz_param_string(object), ") for s in ", input_param, "]\n")
##     })
## setMethod(
##     "python_command",
##     "remove_peaks_around_precursor_mz",
##     function(object, input_param = "py_spectrum_in") {
##         FUN <- .fun_name(object)
##         paste0("import matchms\n",
##             "from matchms.filtering import ", FUN, "\n",
##             "res = [", FUN, "(s, ", .remove_peaks_around_precursor_mz_param_string(object), ") for s in ", input_param, "]\n")
##     })
## setMethod(
##     "python_command",
##     "normalize_intensities",
##     function(object, input_param = "py_spectrum_in") {
##         FUN <- .fun_name(object)
##         paste0("import matchms\n",
##                "from matchms.filtering import ", FUN, "\n",
##                "res = [", FUN, "(s, ", .normalize_intensities_param_string(object), ") for s in ", input_param, "]\n")
##     })

## #' internal function to filter/processing with python's matchms. `Spectra`
## #' will be converted to python `Spectrum` class and matchms' processing
## #' functions will be applied on the `Spectrum` objects. After processing, the
## #' matchms' `Spectrum` objects will be converted back to `Spectra` objects.
## #'
## #' @param sps `Spectra` object
## #'
## #' @param param Parameter object.
## #'
## #' @return a `Spectra` object
## #'
## #' @importFrom basilisk basiliskStart basiliskRun basiliskStop
## #'
## #' @importFrom reticulate py
## #'
## #' @noRd
## #'
## #' @author Thomas Naake, Johannes Rainer
## .filter_spectra_python <- function(sps, param) {
##     ## handle empty input
##     if (!length(sps))
##         return(Spectra())

##     cl <- basiliskStart(matchms_env)
##     on.exit(basiliskStop(cl))

##     basiliskRun(cl, function(sps, param) {
##         ref <- import("matchms")
##         vars <- c(precursorMz = "precursor_mz")
##         py$py_spectrum_in <- rspec_to_pyspec(sps,
##             reference = ref, mapping = vars)

##         ## run the command. Result is in py$res
##         com <- python_command(param)
##         py_run_string(com)

##         ## convert from Python Spectrum to R Spectra and return
##         pyspec_to_rspec(py$res, mapping = vars)

##     }, sps = sps, param = param)
## }
