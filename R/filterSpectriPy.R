#' @title Filter Spectra using Python's matchms library
#'
#' @name filterSpectriPy
#'
#' @description
#'
#' The `filterSpectriPy()` function allows to filter/process a `Spectra` object
#' using the `select_by_intensity()`, `select_by_mz()`,
#' `remove_peaks_around_precursor_mz()`, and `normalize_intensities()`
#' functions of the Python
#' [matchms.filtering](https://matchms.readthedocs.io/en/latest/api/matchms.filtering.html)
#' module.
#'
#' Selection and configuration of the algorithm can be performed with one of
#' the parameter objects (equivalent to *matchms*' function names):
#'
#' - `select_by_intensity()`: Keeps only the peaks within defined intensity
#'   range (keep if `intensity_from` >= intensity >= `intensity_to`). See also
#'   the respective [documentation in *matchms*](https://matchms.readthedocs.io/en/latest/api/matchms.filtering.peak_processing.select_by_intensity.html).
#'
#' - `select_by_mz()`: Keeps only the peaks between `mz_from` and `mz_to`
#'   (keep if `mz_from` >= m/z >= `mz_to`). See also the respective
#'   [documentation in *matchms*](https://matchms.readthedocs.io/en/latest/api/matchms.filtering.peak_processing.select_by_mz.html).
#'
#' - `remove_peaks_around_precursor_mz()`: Removes the peaks that are within
#'   `mz_tolerance` (in Da) of the precursor mz, excluding the precursor peak.
#'
#' - `normalize_intensities()`: Normalizes the intensities of peaks
#'   (and losses) to unit height.
#'
#' @note
#'
#' The first call to the `filterSpectriPy()` function can take longer because
#' the Python environment needs to be first set up.
#'
#' `filterSpectriPy()` first translates the `Spectra` to Python, applies the
#' filter functions from the *matchms* Python libraries and then translates
#' the filtered data back to a `Spectra` object. Thus, any spectra variables
#' other than those that are translated between R and Python will be lost
#' during the processing. Use [setSpectraVariableMapping()] to define which
#' spectra variables should be transferred/converted between R and Python.
#' See also examples below for more information.
#'
#' The [Spectra::Spectra()] object returned by `filterSpectriPy()` will
#' **always** use an in-memory backend (i.e. the [Spectra::MsBackendMemory()]),
#' independently of the backend used by the backend used by the input
#' `Spectra`.
#'
#' @param object A [Spectra::Spectra()] object.
#'
#' @param param one of parameter classes listed above (such as
#'     `select_by_intensity()`) defining the filter/processing function in
#'     Python and its parameters.
#'
#' @param mapping named `character()` defining which spectra variables/metadata
#'     should be converted between R and Python and how they should be renamed.
#'     Defaults to `spectraVariableMapping()`. See [setSpectraVariableMapping()]
#'     for more information.
#'
#' @param intensity_from `numeric(1)`: Set lower threshold for peak intensity.
#'     Default is 10.
#'
#' @param intensity_to `numeric(1)`: Set upper threshold for peak intensity.
#'     Default is 200.
#'
#' @param mz_from `numeric(1)`: Set lower threshold for m/z peak positions.
#'     Default is 0.
#'
#' @param mz_to `numeric(1)`: Set upper threshold for m/z peak positions.
#'     Default is 1000.
#'
#' @param mz_tolerance `numeric(1)`: Tolerance of m/z values that are not
#'     allowed to lie within the precursor mz. Default is 17 Da.
#'
#' @param ... ignored.
#'
#' @return `filterSpectriPy()` returns a `Spectra` object on which the
#'     filtering/processing function has been applied
#'
#' @author Thomas Naake
#'
#' @seealso
#'
#' - [Spectra::filterIntensity()], [Spectra::filterMzRange()],
#'   [Spectra::scalePeaks()] in the `Spectra` package for pure R
#'   implementations of filtering/processing calculations.
#'
#' - [rspec_to_pyspec()] or [pyspec_to_rspec()] for the functions used to
#'   translated the MS data between R and Python.
#'
#' @export
#'
#' @importFrom reticulate py_run_string
#'
#' @examples
#'
#' library(Spectra)
#'
#' ## create some example Spectra
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
#' ## Filter: select_by_intensity
#' res <- filterSpectriPy(
#'     sps, select_by_intensity(intensity_from = 15, intensity_to = 300))
#' ## Only mass peaks with intensities between the specified limits are
#' ## retained
#' intensity(res)
#' ## Compared to the original intensities
#' intensity(sps)
#'
#' ## Note that the spectra variable `"name"` was lost during conversion of
#' ## the MS data between R and Python:
#' sps$name
#' any(spectraVariables(res) == "name")
#'
#' ## Only spectra variables defined by `spectraVariableMapping()` are
#' ## converted and thus retained:
#' spectraVariableMapping()
#'
#' ## We can also pass a custom *spectra variable mapping* with the `mapping`
#' ## parameter to the `filterSpectriPy()` function. Below we create such
#' ## a mapping by adding the translation of a spectra variable `"name"` to
#' ## a metadata name `"compound_name"` to the default spectra variable
#' ## mapping `defaultSpectraVariableMapping()`.
#' map <- c(defaultSpectraVariableMapping(), name = "compound_name")
#' map
#'
#' ## Repeat the filtering operation passing this mapping information:
#' res <- filterSpectriPy(
#'     sps, select_by_intensity(intensity_from = 15, intensity_to = 300),
#'     mapping = map)
#' res$name
#'
NULL

setGeneric("filterSpectriPy", function(object, param, ...)
    standardGeneric("filterSpectriPy"))

setClass("filter_param")
setClass("select_by_intensity",
    contains = "filter_param",
    slots = c(
        intensity_from = "numeric", intensity_to = "numeric"),
    prototype = prototype(
             intensity_from = 20, intensity_to = 200),
    validity = function(object) {
        msg <- NULL
        if (length(object@intensity_from) != 1 || object@intensity_from < 0)
            msg <- c("'intensity_from' has to be a positive number of length 1")
        if (length(object@intensity_to) != 1 || object@intensity_to < 0)
            msg <- c("'intensity_to' has to be a positive number of length 1")
        msg
    }
)
setClass("select_by_mz",
    contains = "filter_param",
    slots = c(
        mz_from = "numeric",
        mz_to = "numeric"),
    prototype = prototype(
        mz_from = 150,
        mz_to = 450),
    validity = function(object) {
        msg <- NULL
        if (length(object@mz_from) != 1 || object@mz_from < 0)
            msg <- c("'mz_from' has to be a positive number of length 1")
        if (length(object@mz_to) != 1 || object@mz_to < 0)
            msg <- c("'mz_to' has to be a positive number of length 1")
        msg
    }
)
setClass("remove_peaks_around_precursor_mz",
    contains = "filter_param",
    slots = c(
        mz_tolerance = "numeric"),
    prototype = prototype(
            mz_tolerance = 20),
    validity = function(object) {
        msg <- NULL
        if (length(object@mz_tolerance) != 1 || object@mz_tolerance < 0)
            msg <- c("'mz_tolerance' has to be a positive number of length 1")
        msg
    }
)
setClass("normalize_intensities",
    contains = "filter_param",
    prototype = prototype(),
    validity = function(object) {
        msg <- NULL
        msg
    }
)

#' @rdname filterSpectriPy
#'
#' @importFrom methods new
#'
#' @export
select_by_intensity <- function(intensity_from = 10, intensity_to = 200) {
    new("select_by_intensity", intensity_from = as.numeric(intensity_from),
        intensity_to = as.numeric(intensity_to))
}

#' @rdname filterSpectriPy
#'
#' @export
select_by_mz <- function(mz_from = 0, mz_to = 1000) {
    new("select_by_mz", mz_from = as.numeric(mz_from),
        mz_to = as.numeric(mz_to))
}

#' @rdname filterSpectriPy
#'
#' @export
remove_peaks_around_precursor_mz <- function(mz_tolerance = 17) {
    new("remove_peaks_around_precursor_mz",
        mz_tolerance = as.numeric(mz_tolerance))
}

#' @rdname filterSpectriPy
#'
#' @export
normalize_intensities <- function() {
    new("normalize_intensities")
}

#' @rdname filterSpectriPy
#'
#' @exportMethod filterSpectriPy
setMethod(
    "filterSpectriPy",
    signature = c(object = "Spectra", param = "filter_param"),
    function(object, param, mapping = spectraVariableMapping(), ...) {
        .filter_spectra_python(object, param, mapping = mapping)
    })

#' Helper method to return a SpectrumProcessor that can be applied to the
#' list of spectra in Python
#'
#' @importFrom reticulate py_dict
#'
#' @noRd
setMethod("py_fun", "select_by_intensity", function(x) {
    matchms_filtering$SpectrumProcessor$create_partial_function(
        matchms_filtering$select_by_intensity,
        py_dict(c("intensity_from", "intensity_to"),
                c(x@intensity_from, x@intensity_to)))
})

setMethod("py_fun", "select_by_mz", function(x) {
    matchms_filtering$SpectrumProcessor$create_partial_function(
        matchms_filtering$select_by_mz,
        py_dict(c("mz_from", "mz_to"),
                c(x@mz_from, x@mz_to)))
})

setMethod("py_fun", "remove_peaks_around_precursor_mz", function(x) {
    matchms_filtering$SpectrumProcessor$create_partial_function(
        matchms_filtering$remove_peaks_around_precursor_mz,
        py_dict(c("mz_tolerance"),
                c(x@mz_tolerance)))
})

setMethod("py_fun", "normalize_intensities", function(x) {
    matchms_filtering$SpectrumProcessor$create_partial_function(
        matchms_filtering$normalize_intensities)
})

#' internal function to filter/processing with python's matchms. `Spectra`
#' will be converted to python `Spectrum` class and matchms' processing
#' functions will be applied on the `Spectrum` objects. After processing, the
#' matchms' `Spectrum` objects will be converted back to `Spectra` objects.
#'
#' @param sps `Spectra` object
#'
#' @param param Parameter object.
#'
#' @return a `Spectra` object
#'
#' @importFrom reticulate py
#'
#' @noRd
#'
#' @author Thomas Naake, Johannes Rainer
.filter_spectra_python <- function(x, param,
                                   mapping = spectraVariableMapping()) {
    ## handle empty input
    if (!length(x))
        return(Spectra())

    proc <- matchms_filtering$SpectrumProcessor$SpectrumProcessor(
        r_to_py(list(py_fun(param))))
    pyspec_to_rspec(
        proc$process_spectra(rspec_to_pyspec(x, mapping = mapping))[0],
        mapping = mapping)
}
