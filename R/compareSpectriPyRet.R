#' @title Spectra similarity calculations using matchms
#'
#' @name compareSpectriPyRet
#'
#' @description
#'
#' The `compareSpectriPyRet()` function allows to calculate spectral similarity
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
#' - `FingerprintSimilarityParam`: Calculate similarity between molecules based 
#'   on their fingerprints. For this similarity measure to work, fingerprints 
#'   are expected to be derived by running *add_fingerprint()*.
#'
#' @param x A [Spectra::Spectra()] object.
#'
#' @param y A [Spectra::Spectra()] object to compare against. If missing,
#'   spectra similarities are calculated between all spectra in `x`.
#'
#' @param param one of parameter classes listed above (such as
#'   `CosineGreedyParam`) defining the similarity scoring function in python
#'   and its parameters.
#'   
#' @param use_existing_env provide path to conda env. If nothing is provided
#'   a basilisk env is set up and used.
#'
#' @param tolerance `numeric(1)`: tolerated differences in peaks' m/z. Peaks
#'   with m/z differences `<= tolerance` are considered matching.
#'
#' @param mzPower `numeric(1)`: the power to raise m/z to in the cosine
#'   function. The default is 0, in which case the peak intensity products will
#'   not depend on the m/z ratios.
#'
#' @param ignorePeaksAbovePrecursor For `NeutralLossesCosineParam()`:
#'   `logical(1)`: if `TRUE` (the default), peaks with m/z values larger than
#'   the precursor m/z are ignored.
#'
#' @param intensityPower `numeric(1)`: the power to raise intensity to in the
#'   cosine function. The default is 1.
#'
#' @param ... ignored.
#'
#' @return `compareSpectriPyRet()` returns a `numeric` matrix with the scores,
#'   number of rows being equal to `length(x)` and number of columns equal to
#'   `length(y)`.
#'
#' @author Carolin Huber, Michael Witting, Johannes Rainer, Helge Hecht,
#'   Marilyn De Graeve
#'
#' @seealso [Spectra::compareSpectra()] in the *Spectra* package for pure R
#'   implementations of spectra similarity calculations.
#'
#' @export
#'
#' @importFrom reticulate py_run_string 
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
#'     c(3.407, 47.494, 3.094, 100.0, 13.240))
#' DF$mz <- list(
#'     c(135.0432, 138.0632, 163.0375, 195.0880),
#'     c(110.0710, 138.0655, 138.1057, 138.1742, 195.0864),
#'     c(109.2, 124.2, 124.5, 170.16, 170.52))
#' sps <- Spectra(DF)
#'
#' ## Calculate pairwise similarity beween all spectra within sps with
#' ## matchms' CosineGreedy algorithm
#' ## Note: the first compareSpectriPy will take longer because the Python
#' ## environment needs to be set up.
#' res <- compareSpectriPyRet(sps, param = CosineGreedyParam())
#' res
#'
#' ## Next we calculate similarities for all spectra against the first one
#' res <- compareSpectriPyRet(sps, sps[1], param = CosineGreedyParam())
#'
#' ## Calculate pairwise similarity of all spectra in sps with matchms'
#' ## ModifiedCosine algorithm
#' res <- compareSpectriPyRet(sps, param = ModifiedCosineParam(), use_existing_env = "my_env")
#' res
#'
#' ## Note that the ModifiedCosine method requires the precursor m/z to be
#' ## known for all input spectra. Thus, it is advisable to remove spectra
#' ## without precursor m/z before using this algorithm.
#' sps <- sps[!is.na(precursorMz(sps))]
#' compareSpectriPyRet(sps, param = ModifiedCosineParam())
NULL

setGeneric("compareSpectriPyRet", function(x, y, param, use_existing_env, ...)
  standardGeneric("compareSpectriPyRet"))

#' @rdname compareSpectriPyRet
#'
#' @exportMethod compareSpectriPyRet
setMethod(
  "compareSpectriPyRet",
  signature = c(x = "Spectra", y = "Spectra", param = "CosineGreedyParam", use_existing_env = "character"),
  function(x, y, param, use_existing_env, ...) {
    .compare_spectra_reticulate(x, y, param, use_existing_env)
  })

#' @rdname compareSpectriPyRet
setMethod(
  "compareSpectriPyRet",
  signature = c(x = "Spectra", y = "missing", param = "CosineGreedyParam", use_existing_env = "missing"),
  function(x, y, param, use_existing_env, ...) {
    .compare_spectra_reticulate(x, y = NULL, param, use_existing_env = NULL)
  })

#' @rdname compareSpectriPyRet
setMethod(
  "compareSpectriPyRet",
  signature = c(x = "Spectra", y = "Spectra", param = "CosineGreedyParam", use_existing_env = "missing"),
  function(x, y, param, use_existing_env, ...) {
    .compare_spectra_reticulate(x, y, param, use_existing_env = NULL)
  })

#' @rdname compareSpectriPyRet
setMethod(
  "compareSpectriPyRet",
  signature = c(x = "Spectra", y = "missing", param = "CosineGreedyParam", use_existing_env = "character"),
  function(x, y, param, use_existing_env, ...) {
    .compare_spectra_reticulate(x, y = NULL, param, use_existing_env)
  })

#' internal helper function which sets up a basilisk conda env unless a local 
#' conda env is provided
#'
#' @param use_existing_env `character`
#' 
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
#'
#' @importFrom reticulate py import use_condaenv
#'
#' @noRd
#'
#' @author Victor Chrone
.setup_conda_env <- function( use_existing_env = NULL) {
  # If the user provides an existing Conda environment, use reticulate
  if (!is.null(use_existing_env)) {
    message("Using existing Conda environment: ", use_existing_env)
    
    # Activate the provided Conda environment
    reticulate::use_condaenv(use_existing_env, required = TRUE)
    
    # Import required Python packages
    matchms <- reticulate::import("matchms", convert = TRUE)
    numpy <- reticulate::import("numpy", convert = TRUE)
    
    return(list(matchms = matchms, numpy = numpy))
  }
  
  # Otherwise, use basilisk to create an isolated environment
  message("No existing Conda environment provided. Using basilisk to set up an isolated Conda environment...")
  
  # Start the environment and install packages if needed
  cl <- basiliskStart(matchms_env)
  on.exit(basiliskStop(cl), add = TRUE)
  
  basiliskRun(cl, function() {
    reticulate::import("matchms")
    reticulate::import("numpy")
  })
  
  # Load the Python packages from the basilisk-managed environment
  matchms <- basiliskRun(cl, function() reticulate::import("matchms", convert = TRUE))
  numpy <- basiliskRun(cl, function() reticulate::import("numpy", convert = TRUE))
  
  return(list(matchms = matchms, numpy = numpy))
}

#' internal helper function to convert R Spectra object to a Python-compatible 
#' matchms.Spectrum format object using reticulate.
#'
#' @param spectra `Spectra`
#'
#' @return A Python list of matchms Spectrum objects
#'
#' @importFrom reticulate py import
#'
#' @noRd
#'
#' @author Victor Chrone
# Convert R Spectra object to Python-compatible format
.convert_to_python <- function(spectra) {
  
  # Import matchms and numpy package 
  matchms <- reticulate::import("matchms")
  similarity <- reticulate::import("matchms.similarity")
  np <- reticulate::import("numpy")  # Needed for array conversion
  
  py_list <- list()
  
  for (i in seq_along(spectra)) {
    py_list[[i]] <- matchms$Spectrum(
      mz = np$array(mz(spectra)[[i]], dtype = "float64"),  # Convert to NumPy array
      intensities = np$array(intensity(spectra)[[i]], dtype = "float64"),  # Convert to NumPy array
      metadata = list(precursor_mz = as.numeric(precursorMz(spectra)[i]))  # Ensure metadata is numeric
    )
  }
  
  return(reticulate::r_to_py(py_list))
}

#' internal function to calculate similarities with python's matchms. `Spectra`
#' using reticulate
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
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
#'
#' @importFrom reticulate py import
#'
#' @noRd
#'
#' @author Victor Chrone
.compare_spectra_reticulate <- function(x, y = NULL, param, use_existing_env = NULL) {
  
  #Set up Conda environment using basilisk if no local Conda environment is provided
  env_packages <- .setup_conda_env(use_existing_env = use_existing_env)
  
  # Import matchms and numpy package 
  matchms <- reticulate::import("matchms")
  similarity <- reticulate::import("matchms.similarity")
  np <- reticulate::import("numpy")  # Needed for array conversion
  
  # Convert R Spectra objects to Python-compatible format
  py_x <- .convert_to_python(x)
  py_y <- if (is.null(y)) py_x else .convert_to_python(y)
  
  # Determine if the comparison is symmetric
  is_symmetric <- is.null(y) || (length(x) == length(y))
  
  # Determine which similarity function to use
  similarity_function <- switch(
    class(param)[1],
    "CosineGreedyParam" = similarity$CosineGreedy,
    "CosineHungarianParam" = similarity$CosineHungarian,
    "ModifiedCosineParam" = similarity$ModifiedCosine,
    "NeutralLossesCosineParam" = similarity$NeutralLossesCosine,
    "FingerprintSimilarityParam" = similarity$FingerprintSimilarity,
    stop("Unknown similarity function")
  )
  
  # Initialize similarity function
  sim_func <- similarity_function(
    tolerance = param@tolerance,
    mz_power = param@mzPower,
    intensity_power = param@intensityPower
  )
  
  # Call calculate_scores from matchms, including `is_symmetric`
  scores_obj <- matchms$calculate_scores(py_x, py_y, similarity_function = sim_func, is_symmetric = is_symmetric)
  
  # Convert structured NumPy array to an R list
  scores_np <- reticulate::py_to_r(scores_obj$scores$to_array())
  
  # Determine correct matrix dimensions
  ncol_scores <- length(scores_np[[1]])
  nrow_scores <- length(scores_np) / ncol_scores
  
  # Initialize similarity matrix
  scores_matrix <- matrix(NA, nrow = nrow_scores, ncol = ncol_scores)
  
  # Extract only similarity scores (first value in tuple)
  for (i in 0:(nrow_scores - 1)) {
    for (j in 0:(ncol_scores - 1)) {
      scores_matrix[(i + 1), (j + 1)] <- as.numeric(scores_np[[i]][j][0])  # Extract first element (CosineGreedy_score)
    }
  }
  
  return(scores_matrix)
}

## from matchms import calculate_scores
## from matchms.importing import load_from_msp
## from ms2deepscore import MS2DeepScore
## from ms2deepscore.models import load_model

## model = load_model("load model file")
## similarity_measure = MS2DeepScore(model)
## scores_mat = similarity_measure.matrix(
## matchms_spectrum1,
## matchms_spectrum2
## )

