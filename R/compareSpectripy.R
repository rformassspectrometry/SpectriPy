compareSpectripy <- function(queries, references=queries, tol=0.005, similarity_function="CosineGreedy", symmetric=F){
    # load conda environment with matchms and import
    use_condaenv("matchms")
    matchms <- import("matchms")
    
    # convert spectra to python spectrum
    spectrum_q <- unname(spectrapply(queries, rspect_to_pyspec))
    spectrum_l <- unname(spectrapply(references, rspect_to_pyspec))
    
    # import similarity functions
    similarity <- import("matchms.similarity")

    # perform comparison and copy results to python for read out
    if(similarity_function == "CosineGreedy"){
        message("Using msmatch CosineGreedy function")
        py$result <- matchms$calculate_scores(references = r_to_py(spectrum_q),
                                          queries = r_to_py(spectrum_l),
                                          similarity_function = similarity$CosineGreedy(tolerance = tol),
                                          is_symmetric = symmetric)
        
        # converting scores from msmatch
        #method <- py_run_file("python/format_scores.py")
        #res <- method$to_data_frame(py$result)
        
        py_run_string("score = []")
        py_run_string("for x,y,z in result:
                    score.append(z['score'])")
        res <- matrix(unlist(py$score), nrow=py$result$n_rows, ncol=py$result$n_cols)
        
        return(res)
        }

    }



