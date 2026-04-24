# Copy Python MS data structure on MS data replacement operations

Assigning a `MsBackendPy` to another variable name in R with e.g.
`a <- b` with `b` being a `MsBackendPy` object results in two different
`MsBackendPy` instances in R that point however to the **same** MS data
structure in Python. Changing thus one of the two `MsBackendPy` either
by subset or data replacement operations changes the shared MS data in
Python and thus inadvertedly also the data from any other `MsBackendPy`
instance in R pointing to the same variable in Python.

Setting `pyspec_copy_on_replace(TRUE)` **before** any data replacement
operation on a `MsBackendPy` instance causes the full MS data in Python
to be copied to a **new** Phython variable before replacing any values.
Thus, in situations described above, `a` and `b` would then point to two
different Python variables. After the subset or replacement operation
has been performed, the *copy-on-replace* setting can again be disabled
with `pyspec_copy_on_replace(FALSE)`.

In general, if there are multiple `MsBackendPy` instances (or `Spectra`
objects) pointing to the same Python variable, then
`pyspec_copy_on_replace(TRUE)` should be called before one of the
replacement methods `$<-`, `intensity<-`, `mz<-`, `peaksData<-`,
`spectraData<-` or
[`applyProcessing()`](https://rdrr.io/pkg/ProtGenerics/man/processingQueue.html)
is called on one of them. After the operation,
`pyspec_copy_on_replace(FALSE)` should be called.

See examples below for details.

Calling `pyspec_copy_on_replace()` without `TRUE` or `FALSE` returns a
`logical(1)` whether *copy-on-replace* is enabled or not.

The *copy-on-replace* approach avoids data corruptions by accidental
modification of common Python MS data shared by different `Spectra`
objects/ `MsBackendPy` instances. On the downside, it comes with a
slighlty lower performance (as the copy has to be cloned) and can result
in potentially multiple copies of the (same) MS data in Python.
*Copy-on-replace* should thus be used with care.

## Usage

``` r
pyspec_copy_on_replace(x = logical())
```

## Arguments

- x:

  [`logical()`](https://rdrr.io/r/base/logical.html)

## Value

If `pyspec_copy_on_replace()` is called without providing an input
parameter, a `logical(1)` is returned whether the *copy-on-replace*
option is enabled or not.

## Details

For the name of the Python MS data variable a number is appended to the
original variable name, ensuring that no other Python variable with the
same name exists (to avoid replacing any other variable in Python). For
example, if the original Python variable's name is `"my_data"`,
`"my_data_1"` is used if there is no other variable with that name. If
there is already another Python variable `"my_data_1"`, then
`"my_data_2"` is used instead.

*copy-on-replace* can also be enabled using an R option (e.g.
`options(pyspec_copy_on_replace = TRUE)`) or System environment variable
(e.g. `Sys.setenv(pyspec_copy_on_replace = TRUE)`).

## Author

Johannes Rainer

## Examples

``` r

## Is copy-on-replace enabled?
pyspec_copy_on_replace()
#> [1] FALSE

## Setting *copy-on-replace* to FALSE (the default)
pyspec_copy_on_replace(FALSE)

## Load an example Spectra object and change to use a `MsBackendPy`
library(Spectra)
library(MsBackendMgf)
fl <- system.file("extdata", "mgf", "test.mgf", package = "SpectriPy")
s <- Spectra(fl, source = MsBackendMgf())
#> Start data import from 1 files ... 
#> done
s <- setBackend(s, MsBackendPy(), pythonVariableName = "mgf_data")

## `s` is now a `Spectra` object with the data in Python, in a variable
## `"mgf_data"`.
s
#> MSn data (Spectra) with 100 spectra in a MsBackendPy backend:
#> Data stored in the "mgf_data" variable in Python
#> Processing:
#>  Switch backend from MsBackendMgf to MsBackendPy [Fri Apr 24 07:23:22 2026] 

## DATA CORRUPTION

## To show how data corruption can happen, we next create a subset of
## `s` and assign that to new variable `s_sub`. Both point to the same
## data in Python
s_sub <- s[1:10]
s_sub
#> MSn data (Spectra) with 10 spectra in a MsBackendPy backend:
#> Data stored in the "mgf_data" variable in Python
#> Processing:
#>  Switch backend from MsBackendMgf to MsBackendPy [Fri Apr 24 07:23:22 2026] 

## Without any data manipulation, both variables are valid
intensity(s)
#> NumericList of length 100
#> [[1]] 754969.625 1058878.75 20211204
#> [[2]] 15660 28164 8616 4516 855092 68248 6148 8996 10596 5176
#> [[3]] 20176 2073248 282228 21640
#> [[4]] 864 2548 10716 1756 72988 12436 1036
#> [[5]] 20176 2073248 282228 21640
#> [[6]] 864 2548 10716 1756 72988 12436 1036
#> [[7]] 20176 2073248 282228 21640
#> [[8]] 864 2548 10716 1756 72988 12436 1036
#> [[9]] 340 416 2580 412
#> [[10]] 360 588 15672 1420 384 3996 368 544 440 13060
#> ...
#> <90 more elements>
intensity(s_sub)
#> NumericList of length 10
#> [[1]] 754969.625 1058878.75 20211204
#> [[2]] 15660 28164 8616 4516 855092 68248 6148 8996 10596 5176
#> [[3]] 20176 2073248 282228 21640
#> [[4]] 864 2548 10716 1756 72988 12436 1036
#> [[5]] 20176 2073248 282228 21640
#> [[6]] 864 2548 10716 1756 72988 12436 1036
#> [[7]] 20176 2073248 282228 21640
#> [[8]] 864 2548 10716 1756 72988 12436 1036
#> [[9]] 340 416 2580 412
#> [[10]] 360 588 15672 1420 384 3996 368 544 440 13060

## However, any data manipulation on any of the two variables will affect
## the shared MS data in Python. Below we assign retention times to the
## data subset `s_sub`. This will cause the associated Python data in
## `"mgf_data"` to be subset to the first 10 spectra and the retention times
## be added.
s_sub$rtime <- 1:10 + 0.1
rtime(s_sub)
#>  [1]  1.1  2.1  3.1  4.1  5.1  6.1  7.1  8.1  9.1 10.1

## Calling `s` throws an error; the data in Python was restricted to the
## first 10 spectra, so, `s` is no longer *valid*.

## AVOID DATA COPRRUPTION WIHT COPY ON REPLACE

## If any MS data needs to be updated, replaced or added during an analysis,
## it is suggested to enable the *copy-on-replace* option to ensure that
## also the MS data in Python gets copied. We repeat the above analysis
## after enabling *copy-on-replace*:
s <- Spectra(fl, source = MsBackendMgf())
#> Start data import from 1 files ... 
#> done
s <- setBackend(s, MsBackendPy(), pythonVariableName = "mgf_data")

pyspec_copy_on_replace(TRUE)

## Subset `s` and add/replace retention times in the subset
s_sub <- s[1:10]
s_sub
#> MSn data (Spectra) with 10 spectra in a MsBackendPy backend:
#> Data stored in the "mgf_data" variable in Python
#> Processing:
#>  Switch backend from MsBackendMgf to MsBackendPy [Fri Apr 24 07:23:22 2026] 

s_sub$rtime <- 1:10 + 0.5

## `s` and `s_sub` are now pointing to two **different** variables in Python
s
#> MSn data (Spectra) with 100 spectra in a MsBackendPy backend:
#> Data stored in the "mgf_data" variable in Python
#> Processing:
#>  Switch backend from MsBackendMgf to MsBackendPy [Fri Apr 24 07:23:22 2026] 
s_sub
#> MSn data (Spectra) with 10 spectra in a MsBackendPy backend:
#> Data stored in the "mgf_data_1" variable in Python
#> Processing:
#>  Switch backend from MsBackendMgf to MsBackendPy [Fri Apr 24 07:23:22 2026] 

## Thus, the data from `s` was not changed
rtime(s)
#>   [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [26] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [51] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [76] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA

## While the one from `s_sub` was
rtime(s_sub)
#>  [1]  1.5  2.5  3.5  4.5  5.5  6.5  7.5  8.5  9.5 10.5

## Disable *copy-on-replace* after the data replacement
pyspec_copy_on_replace(FALSE)
```
