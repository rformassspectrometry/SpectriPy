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

# Summary

Mass spectrometry (MS) is a key technology used across multiple fields, including biomedical research and life sciences. The data is generally large and complex and analyses need to be tailored to the experimental and instrumental setups. Excellent libraries for such data analysis are available in both R and Python, including R packages from the RforMassSpectrometry initiative such as Spectra [@rainer_modular_2022], MsCoreUtils [@rainer_modular_2022], MetaboAnnotation [@rainer_modular_2022] and CompoundDb [@rainer_modular_2022], as well as Python libraries like matchms [@huber_matchms_2020], spectrum_utils [@bittremieux_spectrum_utils_2020], Pyteomics [@goloborodko_pyteomics_2013] and pyOpenMS [@rost_pyopenms_2014]. The reticulate R package [@reticulate_2025] provides an R interface to Python and enables interoperability between the two programming languages. The open source SpectriPy R package builds upon reticulate and provides functionality to efficiently translate between R and Python MS data structures. [some info on the bundled functionality from Python?] SpectriPy hence enables and simplifies the integration of R and Python for MS data analysis, empowering data analysts to benefit from the full power of algorithms from both programming languages. Further, software developers can now reuse algorithms across languages rather than re-implementing them, enhancing efficiency and collaboration.

# Statement of need

Over the past decade, tremendous efforts have been made to develop powerful algorithms and excellent data analysis softwares for MS data analysis. Each of these softwares covers different and in part complementary aspects in the analysis of MS data, but their integration into a single workflow remains, in particular across programming languages, challenging. To avoid the need for repeated implementation of algorithms in different programming languages we developed the SpectriPy package. By leveraging R's “reticulate” system, and translating between Python and R MS data structures, this package enables a seamless cross-language integration of MS data analysis algorithms within unified analysis workflows.

# Example workflow

The concepts and examples can be checked by performing some steps from the package’s vignette, from source https://rformassspectrometry.github.io/SpectriPy/articles/SpectriPy.html.

[Short example workflow. Maybe loading Python MGF file… example from the vignette? Something else? Ideally linking to actual reproducible results/vignettes]

# Perspective (outlook?)

[maybe some words on the hackathon? Future development?]

# Acknowledgements

The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

# References
