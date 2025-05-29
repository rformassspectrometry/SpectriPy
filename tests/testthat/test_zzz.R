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
    expect_false(res)
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
    expect_false(res)
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

## test_that(".onLoad works", {
##     matchms <- NULL
##     res <- .onLoad("SpectriPy", "SpectriPy")
## })

test_that(".initialize_conda works", {
    ## Define mock functions for *reticulate*
    conda_create_mock <- function(envname = NULL, packages = NULL, ...,
                                  forge = TRUE, channel = character(),
                                  environment = NULL, conda = "auto",
                                  python_version = miniconda_python_version(),
                                  additional_create_args = character()) {
        TRUE
    }
    use_condaenv_mock <- function(condaenv = NULL, conda = "auto",
                                  required = NULL) {
        TRUE
    }
    ## Run test
    expect_message(
        res <- with_mocked_bindings(
            "conda_list" = function(conda = "auto") list(),
            "conda_create" = conda_create_mock,
            "use_condaenv" = use_condaenv_mock,
            ".py_check_install" = function(pkgs, envname = .spectripy_env(),
                                           use_conda = .spectripy_use_conda()) {
                c(pkgs, envname, use_conda)
            },
            code = SpectriPy:::.initialize_conda("aa")
        ), "Creating conda environment 'aa'")
    expect_equal(res, c(matchms = "matchms==0.30.0",
                        spectrum_utils = "spectrum_utils==0.3.2",
                        "aa", TRUE))
})

test_that(".initialize_virtualenv works", {
    ## Define mock functions for *reticulate*
    virtualenv_create_mock <-
        function(envname = NULL, python = virtualenv_starter(version),
                 ..., version = NULL, packages = "numpy", requirements = NULL,
                 force = FALSE,
                 module = getOption("reticulate.virtualenv.module"),
                 system_site_packages = getOption(
                     "reticulate.virtualenv.system_site_packages",
                     default = FALSE),
                 pip_version = getOption("reticulate.virtualenv.pip_version",
                                         default = NULL),
                 setuptools_version = getOption(
                     "reticulate.virtualenv.setuptools_version",
                     default = NULL),
                 extra = getOption("reticulate.virtualenv.extra",
                                   default = NULL)) TRUE
    virtualenv_exists_mock <- function(envname = NULL) FALSE
    use_virtualenv_mock <- function(virtualenv = NULL, required = NULL) TRUE
    ## Run test
    expect_message(
        res <- with_mocked_bindings(
            "virtualenv_create" = virtualenv_create_mock,
            "virtualenv_exists" = virtualenv_exists_mock,
            "use_virtualenv" = use_virtualenv_mock,
            ".py_check_install" = function(pkgs, envname = .spectripy_env(),
                                           use_conda = .spectripy_use_conda()) {
                c(pkgs, envname, use_conda)
            },
            code = SpectriPy:::.initialize_virtualenv("bb")
        ), "Creating virtual environment 'bb'")
    expect_equal(res, c(matchms = "matchms==0.30.0",
                        spectrum_utils = "spectrum_utils==0.3.2",
                        "bb", FALSE))
})

test_that(".initialize_libraries works", {
    matchms <- NULL
    matchms_similarity <- NULL
    matchms_filtering <- NULL
    spectrum_utils <- NULL
    ## import_mock <- function(module, as = NULL, convert = TRUE,
    ##                         delay_load = FALSE) module
    ## with_mocked_bindings(
    ##     "import" = import_mock,
    ##     .package = "reticulate",
    ##     code = .initialize_libraries()
    ## )
    ## expect_equal(matchms, "matchms")
})
