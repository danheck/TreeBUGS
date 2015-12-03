# TreeBUGS

TreeBUGS is an R package that facilitates using hierarchical multinomial processing tree (MPT) models that are often used in cognitive psychology (Erdfelder et al., 2009). Specifically, the Beta-MPT (Smith & Batchelder, 2010) and the latent-trait MPT model (Klauer, 2010) are implemented.

## General Procedure of Using TreeBUGS

In the most simple user scenario, the following steps are required:

1. Save MPT model file on disk in .eqn format
2. Save individual data to .csv file (comma separated, rows=persons, columns=labeled categories)
3. Call `mpt2BetaMPT` or `mpt2TraitMPT` (exact code below)
4. Check convergence of MCMC chains
5. Summarize and plot results

These steps are explained in more detail in the package vignette, which can be opened in R by typing `vignette("TreeBUGS")`. Note that TreeBUGS requires a valid installation of the software JAGS (http://mcmc-jags.sourceforge.net/).

### References

Erdfelder, E., Auer, T.-S., Hilbig, B. E., Assfalg, A., Moshagen, M., & Nadarevic, L. (2009). Multinomial processing tree models: A review of the literature. Journal of Psychology, 217, 108–124. http://doi.org/10.1027/0044-3409.217.3.108

Klauer, K. C. (2010). Hierarchical multinomial processing tree models: A latent-trait approach. Psychometrika, 75, 70–98. http://doi.org/10.1007/s11336-009-9141-0

Smith, J. B., & Batchelder, W. H. (2010). Beta-MPT: Multinomial processing tree models for addressing individual differences. Journal of Mathematical Psychology, 54, 167–183. http://doi.org/10.1016/j.jmp.2009.06.007

