# Spectra similarity calculations using Python's matchms library

The `compareSpectriPy()` function allows to calculate spectral
similarity scores using the `calculate_scores()` function of the Python
[matchms.similarity](https://matchms.readthedocs.io/en/latest/api/matchms.similarity.html).
module.

Selection and configuration of the algorithm can be performed with one
of the *parameter* objects/functions:

- `CosineGreedy()`: calculate the *cosine similarity score* between
  spectra. The score is calculated by finding the best possible matches
  between peaks of two spectra. Two peaks are considered a potential
  match if their m/z ratios lie within the given `tolerance`. The
  underlying peak assignment problem is here solved in a *greedy* way.
  This can perform notably faster, but does occasionally deviate
  slightly from a fully correct solution (as with the `CosineHungarian`
  algorithm). In practice this will rarely affect similarity scores
  notably, in particular for smaller tolerances. The algorithm can be
  configured with parameters `tolerance`, `mz_power` and
  `intensity_power` (see parameter description for more details). See
  also [matchms
  CosineGreedy](https://matchms.readthedocs.io/en/latest/api/matchms.similarity.CosineGreedy.html)
  for more information.

- `CosineHungarian()`: calculate the *cosine similarity score* as with
  `CosineGreedy`, but using the Hungarian algorithm to find the best
  matching peaks between the compared spectra. The algorithm can be
  configured with parameters `tolerance`, `mz_power` and
  `intensity_power` (see parameter description for more details). See
  also [matchms
  CosingHungarian](https://matchms.readthedocs.io/en/latest/api/matchms.similarity.CosineHungarian.html)
  for more information.

- `ModifiedCosine()`: The modified cosine score aims at quantifying the
  similarity between two mass spectra. The score is calculated by
  finding the best possible matches between peaks of two spectra. Two
  peaks are considered a potential match if their m/z ratios lie within
  the given `tolerance`, or if their m/z ratios lie within the tolerance
  once a mass shift is applied. The mass shift is simply the difference
  in precursor-m/z between the two spectra. See also [matchms
  ModifiedCosine](https://matchms.readthedocs.io/en/latest/api/matchms.similarity.ModifiedCosine.html)
  for more information.

- `NeutralLossesCosine()`: The neutral losses cosine score aims at
  quantifying the similarity between two mass spectra. The score is
  calculated by finding the best possible matches between peaks of two
  spectra. Two peaks are considered a potential match if their m/z
  ratios lie within the given `tolerance` once a mass shift is applied.
  The mass shift is the difference in precursor-m/z between the two
  spectra. See also [matchms
  NeutralLossesCosine](https://matchms.readthedocs.io/en/latest/api/matchms.similarity.NeutralLossesCosine.html)
  for more information.

- `FingerprintSimilarity()`: Calculate similarity between molecules
  based on their fingerprints. For this similarity measure to work,
  fingerprints are expected to be derived by running
  *add_fingerprint()*. See also [matchms
  FingerprintSimilarity](https://matchms.readthedocs.io/en/latest/api/matchms.similarity.FingerprintSimilarity.html)
  for more information.

## Usage

``` r
CosineGreedy(tolerance = 0.1, mz_power = 0, intensity_power = 1)

CosineHungarian(tolerance = 0.1, mz_power = 0, intensity_power = 1)

ModifiedCosine(tolerance = 0.1, mz_power = 0, intensity_power = 1)

NeutralLossesCosine(
  tolerance = 0.1,
  mz_power = 0,
  intensity_power = 1,
  ignore_peaks_above_precursor = TRUE
)

# S4 method for class 'Spectra,Spectra,CosineGreedy'
compareSpectriPy(x, y, param, ...)

# S4 method for class 'Spectra,missing,CosineGreedy'
compareSpectriPy(x, y, param, ...)
```

## Arguments

- tolerance:

  `numeric(1)`: tolerated differences in the peaks' m/z. Peaks with m/z
  differences `<= tolerance` are considered matching.

- mz_power:

  `numeric(1)`: the power to raise m/z to in the cosine function. The
  default is 0, in which case the peak intensity products will not
  depend on the m/z ratios.

- intensity_power:

  `numeric(1)`: the power to raise intensity to in the cosine function.
  The default is 1.

- ignore_peaks_above_precursor:

  For `NeutralLossesCosine()`: `logical(1)`: if `TRUE` (the default),
  peaks with m/z values larger than the precursor m/z are ignored.

- x:

  A [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  object.

- y:

  A [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  object to compare against. If missing, spectra similarities are
  calculated between all spectra in `x`.

- param:

  One of the parameter classes listed above (such as `CosineGreedy`)
  defining the similarity scoring function in Python and its parameters.

- ...:

  ignored.

## Value

`compareSpectriPy()` Returns a `numeric` matrix with the scores, with
the number of rows equal to `length(x)` and the number of columns equal
to `length(y)`.

## Note

Parameters and algorithms are named as originally defined in the
*matchms* library (i.e. all parameters in *snake_case* while *CamelCase*
is used for the algorithms.

## See also

[`Spectra::compareSpectra()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
in the *Spectra* package for pure R implementations of spectra
similarity calculations.

## Author

Carolin Huber, Michael Witting, Johannes Rainer, Helge Hecht, Marilyn De
Graeve

## Examples

``` r

library(Spectra)
## Create some example Spectra.
DF <- DataFrame(
    msLevel = c(2L, 2L, 2L),
    name = c("Caffeine", "Caffeine", "1-Methylhistidine"),
    precursorMz = c(195.0877, 195.0877, 170.0924)
)
DF$intensity <- list(
    c(340.0, 416, 2580, 412),
    c(388.0, 3270, 85, 54, 10111),
    c(3.407, 47.494, 3.094, 100.0, 13.240)
)
DF$mz <- list(
    c(135.0432, 138.0632, 163.0375, 195.0880),
    c(110.0710, 138.0655, 138.1057, 138.1742, 195.0864),
    c(109.2, 124.2, 124.5, 170.16, 170.52)
)
sps <- Spectra(DF)

## Calculate pairwise similarity beween all spectra within sps with
## matchms' CosineGreedy algorithm
## Note: the first compareSpectriPy will take longer because the Python
## environment needs to be set up.
res <- compareSpectriPy(sps, param = CosineGreedy())
res
#>           [,1]      [,2] [,3]
#> [1,] 1.0000000 0.1948181    0
#> [2,] 0.1948181 1.0000000    0
#> [3,] 0.0000000 0.0000000    1

## Next we calculate similarities for all spectra against the first one
res <- compareSpectriPy(sps, sps[1], param = CosineGreedy())

## Calculate pairwise similarity of all spectra in sps with matchms'
## ModifiedCosine algorithm
res <- compareSpectriPy(sps, param = ModifiedCosine())
res
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.1948181 0.1384183
#> [2,] 0.1948181 1.0000000 0.8520549
#> [3,] 0.1384183 0.8520549 1.0000000

## Note that the ModifiedCosine method requires the precursor m/z to be
## known for all input spectra. Thus, it is advisable to remove spectra
## without precursor m/z before using this algorithm.
sps <- sps[!is.na(precursorMz(sps))]
compareSpectriPy(sps, param = ModifiedCosine())
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.1948181 0.1384183
#> [2,] 0.1948181 1.0000000 0.8520549
#> [3,] 0.1384183 0.8520549 1.0000000
```
