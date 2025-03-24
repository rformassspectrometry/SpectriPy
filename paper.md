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
    affiliation: 2
  - name: Thomas Naake
    affiliation: 3
  - name: Matthias Anagho-Mattanovich
    affiliation: 4
  - name: Nils Hoffmann
    affiliation: 5
  - name: Pierre Marchal
    affiliation: 6
  - name: Victor Chrone
    affiliation: 7
  - name: Philippine Louail
    orcid: 0009-0007-5429-6846    
    affiliation: 1
  - name: Carolin Huber
    affiliation: 8
  - name: Helge Hecht
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
 - name: Novo Nordisk Foundation Center for Basic Metabolic Research, University of Copenhagen, Copenhagen, Denmark.
   index: 4
 - name: Institute for Bio- and Geosciences (IBG-5), Forschungszentrum Jülich GmbH, Jülich, Germany.
   index: 5
 - name: Department of medical oncology, University of Bern, Bern, Switzerland.
   index: 6
 - name: Alphalyse, Odense, Denmark.
   index: 7
 - name: Department of Exposure Science, Helmholtz Centre for Environmental Research - UFZ, Leipzig, Germany.
   index: 8
 - name: Faculty of Science, Masaryk University, Brno, Czech Republic.
   index: 9
 - name: Metabolomics and Proteomics Core, Helmholtz Munich, Munich, Germany.
   index: 10
 - name: Chair of Analytical Food Chemistry, TUM School of Life Sciences, Technical University of Munich, Munich, Germany.
   index: 11
date: 06 March 2025
bibliography: paper.bib
---

![SpectriPy_logo](logo.png){height="150pt"}

- TODO-trackcahnges keep latest
- TODO max 80l
- TODO ref metabonaut *2

# Summary

Mass spectrometry (MS) is a key technology used across multiple fields, including biomedical research and life sciences. The data is generally large and complex and analyses need to be tailored to the experimental and instrumental setups. Excellent libraries for such data analysis are available in both R and Python, including R packages from the RforMassSpectrometry initiative such as Spectra [@rainer_modular_2022], MsCoreUtils [@rainer_modular_2022], MetaboAnnotation [@rainer_modular_2022] and CompoundDb [@rainer_modular_2022], as well as Python libraries like matchms [@huber_matchms_2020], spectrum_utils [@bittremieux_spectrum_utils_2020], Pyteomics [@goloborodko_pyteomics_2013] and pyOpenMS [@rost_pyopenms_2014]. The reticulate R package [@reticulate_2025] provides an R interface to Python and enables interoperability between the two programming languages. The open source SpectriPy R package builds upon reticulate and provides functionality to efficiently translate between R and Python MS data structures. For example, SpectriPy can leverage the spectral similarity, filtering, normalization etc. calculations from the Python matchms library and contains functions to convert between R's Spectra::Spectra objects and Python's matchms.Spectrum objects. SpectriPy hence enables and simplifies the integration of R and Python for MS data analysis, empowering data analysts to benefit from the full power of algorithms from both programming languages. Further, software developers can now reuse algorithms across languages rather than re-implementing them, enhancing efficiency and collaboration.

# Statement of need

Over the past decade, tremendous efforts have been made to develop powerful algorithms and excellent data analysis softwares for MS data analysis. Each of these softwares covers different and in part complementary aspects in the analysis of MS data, but their integration into a single workflow remains, in particular across programming languages, challenging. To avoid the need for repeated implementation of algorithms in different programming languages we developed the SpectriPy package. By leveraging R's `reticulate` package, and translating between Python and R MS data structures, this package enables a seamless cross-language integration of MS data analysis algorithms within unified analysis workflows.

# Description

Reproducible examples and use case analyses enabled by `SpectriPy` can be found in the package’s vignette and the Metabonaut resource [REF]. In this paper, we focus on technical details and features of the package.

## Installation

During installation, `SpectriPy` automatically configures a Python environment and installs all required libraries. The installation can be configured via several environment variables that also allow you to disable the automatic setup and instead use an available Python environment of the host system. Detailed installation instructions can be found in the package’s github repository and the package's vignette.

## Translating MS data objects between R and Python

As its core functionality, `SpectriPy` allows translation between R and Python data structures for MS data. In particular, `SpectriPy` provides the functions `rspec_to_pyspec()` and `pyspec_to_rspec()` to convert between R’s `Spectra::Spectra` and Python’s `matchms.Spectrum` objects. These functions also handle the conversion, and any required renaming and reformatting of spectra metadata, such as MS level, retention times or any other arbitrary metadata available in the MS data object. For more efficient integration of Python MS data objects into R, `SpectriPy` implements a dedicated backend class for `Spectra::Spectra` objects. Such backend classes handle the MS data for `Spectra::Spectra` objects. Through different backend implementations, `Spectra::Spectra` can gain support for additional data and file formats or memory-efficient on-disk- or remote data storage modes. `SpectriPy`’s `MsBackendPy` backend keeps only a reference to the original data object in Python and retrieves and translates MS data, or subsets thereof, only upon request from that object. This enables seamless and memory-efficient integration of Python MS data objects into R for more powerful cross-language analysis workflows. An example of such a combined R-Python data analysis workflow, for example using the Quarto system, is provided in the following code snippets. In the Python code block below, MS data are imported and processed.

```python
#' Python session:
#’ Import data and perform initial processing
import matchms
import matchms.filtering as mms_filt
from matchms.importing import load_from_mgf
mgf_py = list(load_from_mgf(<MGF file>))


#’ Scale intensities
for i in range(len(mgf_p)):
  mgf_py[i] = mms_filt.normalize_intensities(mgf_py[i])
```

Next, to continue the analysis in R, a `Spectra::Spectra` object with a `MsBackendPy` backend class is created, referring to the Python data object defined in the associated Python session above. All data from this data object is accessible in R, with the entire or subsets of the data translated on-the-fly upon request. This strategy ensures memory efficiency and minimizes the number of data copies.

```r
#' R session:
#’ Create an R data object for the MS data in the associated Python session
library(Spectra)
library(SpectriPy)
sps <- Spectra(“mgf_py”, source = MsBackendPy())


#’ Retrieve the MS peaks data for the 1st spectrum
peaksData(sps[1])
```

## Integrated functionality from the matchms Python library

`SpectriPy`’s `compareSpectriPy()` and `filterSpectriPy()` functions allow spectra comparison and filtering and processing routines, respectively, from the `matchms` Python library to be called directly from R. These functions internally translate the MS data from a `Spectra::Spectra` object to the respective Python MS data structures, execute the Python functions, and collect and convert the results to R data types, enabling the integration of functionality from the `matchms` Python library directly into R-based analysis workflows.

# Perspective
	
`SpectriPy` started as a collaboration of R and Python developers, with the latest contributions added during the EuBIC-MS hackathon in 2025. Collaborative development will be further encouraged to extend `SpectriPy` with additional functionality, support for additional libraries and data structures. New use cases will be integrated into larger interactive tutorial frameworks such as the Metabonaut resource [REF], enabling users to seamlessly integrate R and Python into their MS data analysis pipelines.
Ultimately, the long-term goal is to promote cross-language compatibility and reproducibility in computational mass spectrometry. By leveraging the strengths of both R and Python, `SpectriPy` will help develop flexible and efficient MS data analysis workflows, reduce redundancy and promote innovation in the field.

# Acknowledgements

The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

# References
