# connector
It is R package that is able to fit and cluster growth data using four different fitting models (i.e. Functional Clustering Models, Malthus model, Gompertz model, and logistic model). 
It exploits an unsupervised clustering algorithm to cluster the fitted data, and the separation and tightness measures are provide to evaluete the quality of the derived  clusters.

### Required installed packages
The following R packages must be installed to use connector:
DEoptim, GenSA, cowplot, fda, funcy, ggforce, ggplot2, optimr, plyr, readxl.

### How to install connector
To install it you can use  **devtools**:

```
install.packages("devtools")
library(devtools)
install_github("qBioTurin/connector", ref="master",dependencies=TRUE)
```

### An example of connector analysis:
The **MainCommandsList.R** script contains the list of commands that can be used to reproduce the analysis described in the paper "CONNECTOR: fitting and clustering analysis of biological growth data".

To execute this script you can use ***demo()*
```
demo("MainCommandsList", package = "connector")
```


### Diclaimer:
connector developers have no liability for any use of docker4seq functions, including without limitation, any loss of data, incorrect results, or any costs, liabilities, or damages that result from use of connector. 
