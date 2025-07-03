test_that(".initialize_libraries2 works", {
    a <- new.env()
    import_mock <- function(module, as = NULL, convert = TRUE,
                            delay_load = FALSE) module
    with_mocked_bindings(
        "import" = import_mock,
        .package = "reticulate",
        code = .initialize_libraries2(envir = a)
    )
    expect_equal(sort(ls(a)), sort(c("matchms", "matchms_similarity",
                                     "matchms_filtering", "spectrum_utils")))
})

test_that(".onLoad works", {
    ## Does not work because the `asNamespace()` call.
})
