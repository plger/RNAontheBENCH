# The RNAontheBENCH package

This package provides empirical and computational resources for benchmarking RNAseq analysis methods, harnessing a number of features:

- a RNAseq dataset from 12 human iPSC lines all including ExFold ERCC spike-ins (92 transcripts in known concentrations, some of which are differentially-expressed across mixes), as well as an additional validation dataset of 6 samples
- a panel of 150 genes (including multiple isoforms for some genes) measured, in the same samples, by the very precise Nanostring nCounter technology
- a set of internal, 'genetic controls' (genes with different copy-numbers and known to be expressed linearly with copy-numbers)
- accuracy metrics based on feature-wise z-scores across samples, in order to assess differential expression accuracy.

## Installation

To install the package, download it and install it using the following R command:
```
install.packages("path/to/RNAontheBENCH.tar.gz", repos=NULL)
```

Alternatively, if you have `devtools` installed you can install the package directly from the git repository using:
```
library(devtools)
install_git("https://github.com/plger/RNAontheBENCH")
```

## Documentation

Once the package is installed, you can load it and access the vignette for some examples:
```
library(RNAontheBENCH)
vignette("RNAontheBENCH")
```

Alternatively, you can view the vignette [here](inst/doc/RNAontheBENCH.html).
