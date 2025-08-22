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
    orig <- getOption("spectripy.use.system")
    options(spectripy.use.system = TRUE)
})

test_that(".onLoad works", {
    ## Does not work because the `asNamespace()` call.
})

test_that(".is_spectripy_use_system works", {
    orig_opt <- getOption("spectripy.use_system")
    orig_env <- Sys.getenv("SPECTRIPY_USE_SYSTEM")

    options(spectripy.use.system = NULL)
    Sys.unsetenv("SPECTRIPY_USE_SYSTEM")
    expect_false(.is_spectripy_use_system())
    options(spectripy.use_system = "FALSE")
    expect_false(.is_spectripy_use_system())
    options(spectripy.use_system = "True")
    expect_true(.is_spectripy_use_system())
    Sys.setenv(SPECTRIPY_USE_SYSTEM="FALSE")
    expect_true(.is_spectripy_use_system())
    options(spectripy.use_system = NULL)
    expect_false(.is_spectripy_use_system())
    Sys.setenv(SPECTRIPY_USE_SYSTEM="1")
    expect_true(.is_spectripy_use_system())

    if (orig_env == "")
        Sys.unsetenv("SPECTRIPY_USE_SYSTEM")
    else Sys.setenv(SPECTRIPY_USE_SYSTEM=orig_env)
    options(spectripy.use_system = orig_opt)
})

test_that(".is_spectripy_use_conda works", {
    orig_opt <- getOption("spectripy.use_conda")
    orig_env <- Sys.getenv("SPECTRIPY_USE_CONDA")

    options(spectripy.use.conda = NULL)
    Sys.unsetenv("SPECTRIPY_USE_CONDA")
    expect_false(.is_spectripy_use_conda())
    options(spectripy.use_conda = "FALSE")
    expect_false(.is_spectripy_use_conda())
    options(spectripy.use_conda = "True")
    expect_true(.is_spectripy_use_conda())
    Sys.setenv(SPECTRIPY_USE_CONDA="FALSE")
    expect_true(.is_spectripy_use_conda())
    options(spectripy.use_conda = NULL)
    expect_false(.is_spectripy_use_conda())
    Sys.setenv(SPECTRIPY_USE_CONDA="1")
    expect_true(.is_spectripy_use_conda())

    if (orig_env == "")
        Sys.unsetenv("SPECTRIPY_USE_CONDA")
    else Sys.setenv(SPECTRIPY_USE_CONDA=orig_env)
    options(spectripy.use_conda = orig_opt)
})

test_that(".spectripy_env works", {
    orig_opt <- getOption("spectripy.env")
    orig_env <- Sys.getenv("SPECTRIPY_ENV")

    options(spectripy.env = NULL)
    Sys.unsetenv("SPECTRIPY_ENV")
    expect_equal(.spectripy_env(), "")

    options(spectripy.env = "myenvi")
    expect_equal(.spectripy_env(), "myenvi")
    Sys.setenv(SPECTRIPY_ENV="myotherenvi")
    expect_equal(.spectripy_env(), "myenvi")

    options(spectripy.env = NULL)
    expect_equal(.spectripy_env(), "myotherenvi")

    if (orig_env == "")
        Sys.unsetenv("SPECTRIPY_ENV")
    else Sys.setenv(SPECTRIPY_ENV=orig_env)
    options(spectripy.env = orig_opt)
})

test_that(".check_environment works", {
    sys_opt <- getOption("spectripy.use_system")
    sys_env <- Sys.getenv("SPECTRIPY_USE_SYSTEM")
    con_opt <- getOption("spectripy.use_conda")
    con_env <- Sys.getenv("SPECTRIPY_USE_CONDA")
    env_opt <- getOption("spectripy.env")
    env_env <- Sys.getenv("SPECTRIPY_ENV")
    ret_py <- Sys.getenv("RETICULATE_PYTHON")
    ret_env <- Sys.getenv("RETICULATE_PYTHON_ENV")
    ## cleaning options and environment
    options(spectripy.use_system=NULL)
    options(spectripy.use_conda=NULL)
    options(spectripy.env=NULL)
    Sys.unsetenv("SPECTRIPY_USE_SYSTEM")
    Sys.unsetenv("SPECTRIPY_USE_CONDA")
    Sys.unsetenv("SPECTRIPY_ENV")
    Sys.unsetenv("RETICULATE_PYTHON")
    Sys.unsetenv("RETICULATE_PYTHON_ENV")

    expect_silent(.check_environment())
    options(spectripy.use_system=TRUE)
    expect_warning(.check_environment(), "variable 'SPECTRIPY_USE_SYSTEM'")
    options(spectripy.use_system=NULL)
    options(spectripy.use_conda=1)
    expect_warning(.check_environment(), "variable 'SPECTRIPY_USE_CONDA'")
    options(spectripy.use_conda=NULL)
    options(spectripy.env="something")
    expect_warning(.check_environment(), "variable 'SPECTRIPY_ENV'")
    options(spectripy.env=NULL)
    Sys.setenv(RETICULATE_PYTHON="something")
    expect_message(.check_environment(), "Using Python defined by")
    Sys.unsetenv("RETICULATE_PYTHON")
    Sys.setenv(RETICULATE_PYTHON_ENV="something")
    expect_message(.check_environment(), "Using Python environment defined by")
    Sys.unsetenv("RETICULATE_PYTHON_ENV")

    ## Restoring options and environment
    options(spectripy.use_system = sys_opt)
    options(spectripy.use_conda = con_opt)
    options(spectripy.env = env_opt)
    if (sys_env == "")
        Sys.unsetenv("SPECTRIPY_USE_SYSTEM")
    else Sys.setenv(SPECTRIPY_USE_SYSTEM=sys_env)
    if (con_env == "")
        Sys.unsetenv("SPECTRIPY_USE_CONDA")
    else Sys.setenv(SPECTRIPY_USE_CONDA=con_env)
    if (env_env == "")
        Sys.unsetenv("SPECTRIPY_ENV")
    else Sys.setenv(SPECTRIPY_ENV=env_env)
    if (ret_py == "")
        Sys.unsetenv("RETICULATE_PYTHON")
    else Sys.setenv(RETICULATE_PYTHON=ret_py)
    if (ret_env == "")
        Sys.unsetenv("RETICULATE_PYTHON_ENV")
    else Sys.setenv(RETICULATE_PYTHON_ENV=ret_env)
})
