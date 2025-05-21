test_that(".spectripy_env works", {
    res <- SpectriPy:::.spectripy_env()
    expect_equal(res, "r-spectripy")
    orig <- Sys.getenv("SPECTRIPY_ENV")
    Sys.setenv(SPECTRIPY_ENV="test")
    res <- SpectriPy:::.spectripy_env()
    expect_equal(res, "test")
    if (orig == "") Sys.unsetenv("SPECTRIPY_ENV")
    else Sys.setenv(SPECTRIPY_ENV=orig)
    res <- SpectriPy:::.spectripy_env()
    expect_equal(res, "r-spectripy")
})

test_that(".spectripy_use_conda works", {
    res <- SpectriPy:::.spectripy_use_conda()
    expect_true(res)
    orig <- Sys.getenv("SPECTRIPY_USE_CONDA")
    Sys.setenv(SPECTRIPY_USE_CONDA=FALSE)
    res <- SpectriPy:::.spectripy_use_conda()
    expect_false(res)
    Sys.setenv(SPECTRIPY_USE_CONDA=TRUE)
    res <- SpectriPy:::.spectripy_use_conda()
    expect_true(res)
    if (orig == "") Sys.unsetenv("SPECTRIPY_USE_CONDA")
    else Sys.setenv(SPECTRIPY_USE_CONDA=orig)
    res <- SpectriPy:::.spectripy_use_conda()
    expect_true(res)
})

test_that(".spectripy_use_system works", {
    res <- SpectriPy:::.spectripy_use_system()
    expect_false(res)
    orig <- Sys.getenv("SPECTRIPY_USE_SYSTEM")
    Sys.setenv(SPECTRIPY_USE_SYSTEM=TRUE)
    res <- SpectriPy:::.spectripy_use_system()
    expect_true(res)
    if (orig == "") Sys.unsetenv("SPECTRIPY_USE_SYSTEM")
    else Sys.setenv(SPECTRIPY_USE_SYSTEM=orig)
    res <- SpectriPy:::.spectripy_use_system()
    expect_false(res)
})

test_that(".py_check_install works", {
    virtualenv_install_mock <- function(envname = NULL, packages = NULL,
                                        ignore_installed = FALSE,
                                        pip_options = character(),
                                        requirements = NULL, ...,
                                        python_version = NULL) {
        message("Calling virtualenv_install_mock")
    }
    conda_install_mock <- function(envname = NULL, packages, forge = TRUE,
                                   channel = character(), pip = FALSE,
                                   pip_options = character(),
                                   pip_ignore_installed = FALSE,
                                   conda = "auto", python_version = NULL,
                                   additional_create_args = character(),
                                   additional_install_args = character(), ...) {
        message("Calling conda_install_mock")
    }
    ## Virtualenv
    expect_message(
        with_mocked_bindings(
            "virtualenv_install" = virtualenv_install_mock,
            .package = "reticulate",
            code = .py_check_install("matchmsddd", envname = "aaaa", FALSE)
        ),
        "Installing required library 'matchmsddd'"
    )
    expect_message(
        with_mocked_bindings(
            "virtualenv_install" = virtualenv_install_mock,
            .package = "reticulate",
            code = .py_check_install("matchmsddd", envname = "aaaa", FALSE)
        ),
        "Calling virtualenv_install_mock"
    )
    ## Conda
    expect_message(
        with_mocked_bindings(
            "conda_install" = conda_install_mock,
            .package = "reticulate",
            code = .py_check_install("matchmsddd", envname = "aaaa", TRUE)
        ),
        "Installing required library 'matchmsddd'"
    )
    expect_message(
        with_mocked_bindings(
            "conda_install" = conda_install_mock,
            .package = "reticulate",
            code = .py_check_install("matchmsddd", envname = "aaaa", TRUE)
        ),
        "Calling conda_install_mock"
    )
})

## test_that(".install_python_packages works", {
##     py_install_mock <- function(packages, envname = NULL,
##                                 method = c("auto", "virtualenv", "conda"),
##                                 conda = "auto", python_version = NULL,
##                                 pip = FALSE, ..., pip_ignore_installed = TRUE,
##                                 ignore_installed = FALSE) {
##         list(packages = packages, envname = envname, method = method,
##              pip = pip)
##     }
##     .is_matchms_available_mock <- function() FALSE
##     py_module_available_mock <- function(module) {
##         return(FALSE)
##     }

##     res <- with_mocked_bindings(
##         "py_module_available" = py_module_available_mock,
##         .package = "reticulate",
##         code = py_module_available("numpy")
##     )
##     expect_false(res)

##     res <- with_mocked_bindings(
##         ".is_matchms_available" = function() FALSE,
##         "py_install" = py_install_mock,
##         .package = "reticulate",
##         code = SpectriPy:::.install_python_packages("mytest", use_conda = TRUE)
##     )
##     expect_equal(res, list(packages = "matchms==0.28.2",
##                            envname = "mytest",
##                            method = "conda",
##                            pip = FALSE))

##     res <- with_mocked_bindings(
##         "py_module_available" = py_module_available_mock,
##         "py_install" = py_install_mock,
##         .package = "reticulate",
##         code = SpectriPy:::.install_python_packages("mytest", use_conda = FALSE)
##     )
##     expect_equal(res, list(packages = c("matchms==0.28.2", "numpy==2.0.2"),
##                            envname = "mytest",
##                            method = "virtualenv",
##                            pip = FALSE))

##     ## spectrum_utils missing
##     py_module_available_mock <- function(module) {
##         if (module == "matchms")
##             TRUE
##         else FALSE
##     }
##     res <- with_mocked_bindings(
##         "py_module_available" = py_module_available_mock,
##         "py_install" = py_install_mock,
##         .package = "reticulate",
##         code = SpectriPy:::.install_python_packages("mytest2", use_conda = TRUE)
##     )
##     expect_equal(res, list(packages = "matchms==0.28.2",
##                            envname = "mytest2",
##                            method = "conda",
##                            pip = FALSE))
##     res <- with_mocked_bindings(
##         "py_module_available" = py_module_available_mock,
##         "py_install" = py_install_mock,
##         .package = "reticulate",
##         code = SpectriPy:::.install_python_packages("mytest2", use_conda = FALSE)
##     )
##     expect_equal(res,
##                  list(packages = c("spectrum_utils==0.3.2", "numpy==2.0.2"),
##                       envname = "mytest2",
##                       method = "virtualenv",
##                       pip = FALSE))
## })

## test_that(".onLoad works", {
##     py_install_mock <- function(packages, envname = NULL,
##                                 method = c("auto", "virtualenv", "conda"),
##                                 conda = "auto", python_version = NULL,
##                                 pip = FALSE, ..., pip_ignore_installed = TRUE,
##                                 ignore_installed = FALSE) {
##         list(packages = packages, envname = envname, method = method,
##              pip = pip)
##     }
##     expect_message(
##         with_mocked_bindings(
##             "py_install" = py_install_mock,
##             ".spectripy_env" = function() "a",
##             ".spectripy_use_conda" = function() FALSE,
##             ".spectripy_use_system" = function() TRUE,
##             code = .onLoad(libname = "SpectriPy", pkgname = "SpectriPy")
##         ), "Using system Python")
##     expect_message(
##         with_mocked_bindings(
##             ".spectripy_env" = function() "a",
##             ".spectripy_use_conda" = function() TRUE,
##             ".spectripy_use_system" = function() FALSE,
##             .package = "SpectriPy",
##             code = .onLoad(libname = "SpectriPy", pkgname = "SpectriPy")
##         ), "Using conda environment 'a'")

## })
