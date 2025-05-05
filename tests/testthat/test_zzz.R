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

test_that(".install_python_packages works", {
    #' mock for `py_module_available)_`, `py_install()`
    res <- with_mocked_bindings(
        "py_module_available" = function() {},
        "py_install" = function() {## return what happened.},
        code = alalalal
    )
})
