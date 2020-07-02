# Seismic Inversion

Experimental codes for RF inversion to 1D earth model.  Work in progress techniques for inversion of receiver
functions and seismograms.

## Reverse Jump Monte Carlo Markov Chain

Experimental code in subfolder `mcmc` using Reverse Jump Monte Carlo Markov Chain methods.
Fortran solver application with thin Python scripts for collecting and plotting results. Code may benefit from
further validation study.

## Model inference from constraint minimization

Experimental method for automating of [Tao method][1] for solving 1D model within
 prescribed number of layers and parameter bounds. See subfolder `wavefield_decomp`.

[1]: https://doi.org/10.1093/gji/ggt515
