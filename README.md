
[![Build Status](https://travis-ci.com/abodein/netOmics.svg?branch=master)](https://travis-ci.com/abodein/netOmics)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# netOmics

With netOmics, we go beyond integration by introducing an interpretation tool.
netOmics is a package for the creation and exploration of multi-omics networks.

Depending on the provided dataset, it allows to create inference networks from expression data but also interaction networks from knowledge databases.
After merging the sub-networks to obtain a global multi-omics network, we propose network exploration methods using propoagation techniques to perform functional prediction or identification of molecular mechanisms.

Furthermore, the package has been developed for longitudinal multi-omics data and can be used in conjunction with our previously published package timeOmics.

for more examples, please visite https://github.com/abodein/netOmics-case-studies

## Installation

### Latest `GitHub` Version


Install the devtools package in R, then load it and install the latest stable version of `netOmics` from `GitHub`

```r
## install devtools if not installed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
## install netOmics
devtools::install_github("abodein/netOmics")
```

## Citing

*"Bodein, A., Scott-Boyer, M. P., Perin, O., Le Cao, K. A., & Droit, A. (2020). Interpretation of network-based integration from multi-omics longitudinal data. bioRxiv."*

## Maintainer
Antoine Bodein (<antoine.bodein.1@ulaval.ca>)

## Bugs/Feature requests

If you have any bugs or feature requests, [let us know](https://github.com/abodein/netOmics/issues). 
Thanks!

