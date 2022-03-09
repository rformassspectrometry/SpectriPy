#' @title Low level functions to use compare functions of matchms in R
#'
#' @name compareSpectriPy
#'
#' @description
#'
#' The `compareSpectriPy` function allow to calculate spectral similarity values 
#' based on `calculate_scores module` of the python [matchms.similarity package](https://matchms.readthedocs.io/en/latest/api/matchms.similarity.html) package.
#' The function returns matrix that contains all scores calculated by the chosen
#' matchms. Per default, `CosineGreedy` is used.
#'
#' @param x A spectra object
#' 
#' @param y A spectra object to compare against, if left empty, x is compared against itself
#'
#' @param FUN Which comparison function of `matchms` shall be used. Default `CosineGreedy`, other possibilities `CosineHungarian`, `ModifiedCosine`, `NeutralLossesCosine`
#' 
#' @param tolerance Peaks will be considered a match when <= tolerance apart. Default is 0.1.
#'
#' @param mz_power The power to raise m/z to in the cosine function. The default is 0, in which case the peak intensity products will not depend on the m/z ratios.
#' 
#' @param intensity_power The power to raise intensity to in the cosine function. The default is 1.
#' 
#' @param ignore_peaks_above_precursor only used for `NeutralLossesCosine`
#' 
#' @param output Shall score or matches given in output matrix
#'
#' @return A matrix containing the calculated scores/ number of matched fragments
#'
#' @author Michael Witting, Johannes Rainer, Helge Hecht, Carolin Huber
#'
#' @export
#'
#' @importFrom reticulate py_run_string
#'
#' @examples
#'
compareSpectriPy <- function(x, 
                             y=NULL, 
                             FUN="CosineGreedy",
                             tolerance=0.2,
                             mz_power=0,
                             intensity_power=1,
                             ignore_peaks_above_precursor= FALSE,
                             output="scores"){
    
    # calculate half of the matrix if only one spectra object applied
    if(is.null(y)){
        is_symetric=TRUE
    }else{
        is_symetric=FALSE
    }
    
    # determine function
    if(!FUN %in% c("CosineGreedy", "CosineHungarian", "ModifiedCosine", "NeutralLossesCosine")){
        stop("Unknown spectral similarity function applied.")
    }
    
    # convert spectra objects to spectrum
    x_py <- rspec_to_pyspec(x) 
    
    if(is.null(y)){
        include_y_py <- ""
    }else{
        include_y_py <- ",y_py"
        y_py <- rspec_to_pyspec(y)
    }
    
    # modify function command depending on variables
    if(FUN=="NeutralLossesCosine"){
        FUN_defined <- paste0(FUN, "(tolerance=",tolerance, " mz_power=", mz_power, ", intensity_power=",intensity_power, ", ignore_peaks_above_precursor=", ignore_peaks_above_precursor,")")
    }else{
        FUN_defined <- paste0(FUN, "(tolerance=",tolerance, " mz_power=", mz_power, ", intensity_power=",intensity_power, ")")
    }
    
    pycommand <- paste0("result = matchms.calculate_scores( x_py", include_y_py, ", ", FUN_defined,", ", is_symetric,")")
    
    # run python command
    py_run_string(pycommand)
    
    # export scores
    py_run_string("res = []")
    if(output='score'){
        py_run_string("for x,y,z in result:
                        res.append(z['score'])")
    }
    # export matches
    if(output='matches'){
        py_run_string("for x,y,z in result:
                        res.append(z['matches'])")
    }
    
    # generate r matrix 
    score <- matrix(unlist(py$res), nrow=py$result$n_rows, ncol=py$result$n_cols)
    
    return(scores)
}
