[![Build Status](https://api.travis-ci.org/psobczyk/geneSLOPE.png)](https://travis-ci.org/psobczyk/geneSLOPE)

[<img src="http://www.ideal.rwth-aachen.de/wp-content/uploads/2013/08/banner1.png">](http://www.ideal.rwth-aachen.de/)

-------------

**geneSLOPE** -- Genome Wide Association Study with SLOPE
-------------------------

Package **geneSLOPE** can be used to perform Genome-wide association study with 
[SLOPE](http://statweb.stanford.edu/~candes/SortedL1/). 

Such an analysis  is split into three steps.

1. Data is read and immediately screened using marginal test for each SNP

2. SNPs are clumped based on correlations

3. SLOPE is performed on reduced data. Each clump is
represented by one SNP

-------------------------

### How do I get set up? ###

You can install the latest development version of the code using the `devtools` R package.

```R
# Install devtools, if you haven't already.
install.packages("devtools")

library(devtools)
install_github("psobczyk/geneSLOPE")
```

-------------------------

* You might need to install package dependencies:
    * SLOPE
    * ggplot2
    * bigmemory
    * grid
* Read vignette **"Tutorial for GWAS with SLOPE"** to get familiar with basic usage


#### Running GUI ####

GUI based on *shiny* R package is available

```R
library(geneSLOPE)
gui_geneSLOPE()
```

----------------------

### Who do I talk to? ###
* If help provided in the package documentation does not solve your problem
please contact Piotr.Sobczyk[at]pwr.edu.pl

-------------

Research reported in this software was supported by National Institutes of Health under award number R01 HG006695.

-------------

This project has received funding from the European Union’s
Seventh Framework Programme for research, technological
development and demonstration under grant agreement no 602552.
