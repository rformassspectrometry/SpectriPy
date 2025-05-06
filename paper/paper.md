---
title: 'SpectriPy: Enhancing Cross-Language Mass Spectrometry Data Analysis with R and Python'
tags:
  - Mass Spectrometry
  - Cross-Language Integration
  - R
  - Python
  - Computational Efficiency
  - Reproducible Workflows
  - Proteomics
  - Metabolomics
authors:
  - name: Marilyn De Graeve
    orcid: 0000-0001-6916-401X
    affiliation: 1
  - name: Wout Bittremieux
    orcid: 0000-0002-3105-1359
    affiliation: 2
  - name: Thomas Naake
    orcid: 0000-0001-7917-5580
    affiliation: 3
  - name: Carolin Huber
    orcid: 0000-0002-9355-8948
    affiliation: 4
  - name: Matthias Anagho-Mattanovich
    orcid: 0000-0001-7561-7898
    affiliation: 5
  - name: Nils Hoffmann
    orcid: 0000-0002-6540-6875
    affiliation: 6
  - name: Pierre Marchal
    orcid: 0009-0006-6567-6257
    affiliation: 7
  - name: Victor Chrone
    orcid: 0009-0007-2121-4066
    affiliation: 8
  - name: Philippine Louail
    orcid: 0009-0007-5429-6846
    affiliation: 1
  - name: Helge Hecht
    orcid: 0000-0001-6744-996X
    affiliation: 9
  - name: Michael Witting
    orcid: 0000-0002-1462-4426
    affiliation: "10, 11"
  - name: Johannes Rainer
    orcid: 0000-0002-6977-7147
    corresponding: true
    affiliation: 1
affiliations:
  - name: Institute for Biomedicine, Eurac Research, Bolzano, Italy.
    index: 1
  - name: Department of Computer Science, University of Antwerp, Antwerpen, Belgium.
    index: 2
  - name: Genome Biology Unit, European Molecular Biology Laboratory (EMBL), Heidelberg, Germany.
    index: 3
  - name: Department of Exposure Science, Helmholtz Centre for Environmental Research - UFZ, Leipzig, Germany.
    index: 4
  - name: Novo Nordisk Foundation Center for Basic Metabolic Research, University of Copenhagen, Copenhagen, Denmark.
    index: 5
  - name: Institute for Bio- and Geosciences (IBG-5), Forschungszentrum Jülich GmbH, Jülich, Germany.
    index: 6
  - name: Department of Medical Oncology, University of Bern, Bern, Switzerland.
    index: 7
  - name: Alphalyse, Odense, Denmark.
    index: 8
  - name: RECETOX, Faculty of Science, Masaryk University, Brno, Czech Republic.
    index: 9
  - name: Metabolomics and Proteomics Core, Helmholtz Zentrum München, Munich, Germany.
    index: 10
  - name: Chair of Analytical Food Chemistry, TUM School of Life Sciences, Technical University of Munich, Freising-Weihenstephan, Germany.
    index: 11
date: 03 April 2025
bibliography: paper.bib
---

![SpectriPy package logo](../logo.png){height="150pt"}

# Summary

Mass spectrometry (MS) is a key technology used across multiple fields,
including biomedical research and life sciences. The data is often times large
and complex, and analyses must be tailored to the experimental and instrumental
setups. Excellent software libraries for such data analysis are available in
both R and Python, including R packages from the RforMassSpectrometry initiative
such as *Spectra*, *MsCoreUtils*, *MetaboAnnotation*, and *CompoundDb*
[@rainer_modular_2022], as well as Python libraries like *matchms*
[@huber_matchms_2020], *spectrum_utils* [@bittremieux_spectrum_utils_2020],
*Pyteomics* [@goloborodko_pyteomics_2013] and *pyOpenMS*
[@rost_pyopenms_2014]. The *reticulate* R package [@reticulate_2025] provides an
R interface to Python enabling interoperability between the two programming
languages. The open-source *SpectriPy* R package builds upon *reticulate* and
provides functionality to efficiently translate between R and Python MS data
structures. It can convert between R’s `Spectra::Spectra` and Python’s
`matchms.Spectrum` and `spectrum_utils.spectrum.MsmsSpectrum` objects and
includes functionality to directly apply spectral similarity, filtering,
normalization, etc. routines from the Python *matchms* library on MS data in
R. *SpectriPy* hence enables and simplifies the integration of R and Python for
MS data analysis, empowering data analysts to benefit from the full power of
algorithms in both programming languages. Furthermore, software developers can
reuse algorithms across languages rather than re-implementing them, enhancing
efficiency and collaboration.

# Statement of need

Over the past decade, tremendous efforts have been made to develop powerful
algorithms and excellent data analysis software for MS data analysis. Each of
these software covers different and in part complementary aspects in the
analysis of MS data, but their integration into a single workflow remains a
challenge, in particular across programming languages. To avoid the need for
repeated implementation of algorithms in different programming languages we
developed the *SpectriPy* package. By leveraging R's *reticulate* package, and
translating between R and Python MS data structures, this package enables
seamless cross-language integration of MS data analysis algorithms within
unified analysis workflows.

# Description

Reproducible examples and use case analyses on how to share and translate MS
data structures between R and Python, and combined Python and R-based analysis
workflows for LC-MS/MS data annotation enabled by *SpectriPy* can be found in
the package’s vignette and in one of the [example
workflows](https://rformassspectrometry.github.io/Metabonaut/articles/SpectriPy_tutorial_metabonaut.html)
of the Metabonaut resource [@louail_metabonaut_2025]. In this paper, we
primarily focus on the technical details and features of the package.

## Installation

During installation, *SpectriPy* automatically configures a Python environment
and installs all required libraries. The installation can be configured through
several environment variables that also allow users to disable the automatic
setup and instead use an available Python environment of the host
system. Detailed installation instructions can be found in the package’s GitHub
repository and in the package’s vignette.

## Translating MS data objects between R and Python

As its core functionality, *SpectriPy* allows translation between R and Python
MS data structures. In particular, *SpectriPy* provides the functions
`rspec_to_pyspec()` and `pyspec_to_rspec()` to convert between R’s
`Spectra::Spectra` and Python’s `matchms.Spectrum` and
`spectrum_utils.spectrum.MsmsSpectrum` objects. These functions also handle the
conversion, and any required renaming and reformatting of spectra metadata, such
as MS level, retention times, or any other arbitrary metadata available in the
MS data object. For more efficient integration of Python MS data objects into R,
*SpectriPy* implements a dedicated backend class for `Spectra::Spectra`
objects. Such backend classes handle the MS data for `Spectra::Spectra`
objects. Through different backend implementations, `Spectra::Spectra` can gain
support for additional data and file formats or memory-efficient on-disk or
remote data storage modes. *SpectriPy*’s `MsBackendPy` backend keeps only a
reference to the original data object in Python and retrieves and translates MS
data, or subsets thereof, only upon request from that object. This enables
seamless and memory-efficient integration of Python MS data objects into R for
more powerful cross-language analysis workflows. An example of such a combined
R-Python data analysis workflow, which can be realized e.g. using the Quarto
system, is provided in the following code snippets. In the Python code block
below, MS data are imported and processed.

```python
#' Python session:
#'  Import data and perform initial processing
import matchms
import matchms.filtering as mms_filt
from matchms.importing import load_from_mgf
mgf_py = list(load_from_mgf(<MGF file>))

#'  Scale intensities
for i in range(len(mgf_p)):
  mgf_py[i] = mms_filt.normalize_intensities(mgf_py[i])
```

To continue the analysis in R, a `Spectra::Spectra` object with a `MsBackendPy`
backend class is created, referring to the Python data object defined in the
associated Python session. All data from this data object is accessible in R,
with the entire or subsets of the data translated on-the-fly upon request. This
strategy ensures memory efficiency and minimizes the number of data copies.

```r
#' R session:
#'  Create an R data object for the MS data in the associated Python session
library(Spectra)
library(SpectriPy)
sps <- Spectra(“mgf_py”, source = MsBackendPy())

#'  Retrieve the MS peaks data for the 1st spectrum
peaksData(sps[1])
```

## Integrated functionality from the *matchms* Python library

*SpectriPy*’s `compareSpectriPy()` and `filterSpectriPy()` functions allow
spectra comparison, and filtering and processing routines, respectively, from
the *matchms* Python library to be called directly from R. These functions
internally translate the MS data from a `Spectra::Spectra` object to the
respective Python MS data structures, execute the Python functions, and collect
and convert the results to R data types, enabling the integration of
functionality from the *matchms* Python library directly into R-based analysis
workflows.

# Perspective

*SpectriPy* started as a collaboration of R and Python developers, with the
latest contributions added during the EuBIC-MS Developers Meeting
in 2025. Collaborative development will be further encouraged to extend
*SpectriPy* with additional functionality, support for additional libraries, and
data structures. New use cases will be integrated into larger interactive
tutorial frameworks such as the Metabonaut resource [@louail_metabonaut_2025],
enabling users to seamlessly integrate R and Python into their MS data analysis
pipelines.

Ultimately, the long-term goal is to promote cross-language compatibility and
reproducibility in computational mass spectrometry. By leveraging the strengths
of both R and Python, *SpectriPy* will help develop flexible and efficient MS
data analysis workflows, reduce redundancy and promote innovation in the field.

# Acknowledgements

The authors declare to not have any competing financial or personal interests
that could have influenced the work reported in this paper. Part of this work
was supported by the European Union HORIZON-MSCA-2021 project 101073062:
"HUMAN - Harmonizing and unifying blood metabolic analysis networks" to
Philippine Louail, Johannes Rainer and Micheal Witting. Helge Hecht thanks the
RECETOX Research Infrastructure (No LM2023069) financed by the Ministry of
Education, Youth and Sports for supportive background. Part of this project was
supported from the European Union’s Horizon 2020 research and innovation
programme under grant agreement 857560 (CETOCOEN Excellence) and from the
Horizon Europe programme under grant agreement 101079789 (EIRENE PPP). Views and
opinions expressed are however those of the author(s) only and do not
necessarily reflect those of the European Union or European Research Executive
Agency (REA). Neither the European Union nor the granting authority can be held
responsible for any use that may be made of the information it contains.


# References
