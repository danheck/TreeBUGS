# TreeBUGS

TreeBUGS is an R package that facilitates using hierarchical multinomial processing tree (MPT) models that are often used in cognitive psychology (Erdfelder et al., 2009). Specifically, the Beta-MPT (Smith & Batchelder, 2010) and the latent-trait MPT model (Klauer, 2010) are implemented.

## General Procedure of Using TreeBUGS

In the most simple user scenario, the following steps are required:

1. Define path to existing MPT model file in .eqn format (cf. multiTree; Moshagen, 2010)
2. Define path to data set with individual frequencies (.csv file: comma separated, rows=persons, columns=labeled categories)
3. Call `betaMPT` or `traitMPT` (exact code in manual/vignette)
4. Check convergence of MCMC chains
5. Summarize and plot results using functions tailored to MPT models

These steps are explained in more detail in the package vignette, which can be opened in R by typing `vignette("TreeBUGS")`. Note that TreeBUGS requires a valid installation of the software JAGS (http://mcmc-jags.sourceforge.net/). 

To install TreeBUGS, run:
```
install.packages("devtools")
library(devtools)
install_github("denis-arnold/TreeBUGS")
```

### References

* Erdfelder, E., Auer, T.-S., Hilbig, B. E., Assfalg, A., Moshagen, M., & Nadarevic, L. (2009). Multinomial processing tree models: A review of the literature. Journal of Psychology, 217, 108–124. http://doi.org/10.1027/0044-3409.217.3.108

* Klauer, K. C. (2010). Hierarchical multinomial processing tree models: A latent-trait approach. Psychometrika, 75, 70–98. http://doi.org/10.1007/s11336-009-9141-0

* Matzke, D., Dolan, C. V., Batchelder, W. H., & Wagenmakers, E.-J. (2015). Bayesian estimation of multinomial processing tree models with heterogeneity in participants and items. Psychometrika, 80, 205–235. http://doi.org/10.1007/s11336-013-9374-9

* Moshagen, M. (2010). multiTree: A computer program for the analysis of multinomial processing tree models. Behavior Research Methods, 42, 42–54. http://doi.org/10.3758/BRM.42.1.42

* Smith, J. B., & Batchelder, W. H. (2010). Beta-MPT: Multinomial processing tree models for addressing individual differences. Journal of Mathematical Psychology, 54, 167–183. http://doi.org/10.1016/j.jmp.2009.06.007

