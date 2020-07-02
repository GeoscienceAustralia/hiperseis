# Wavefield decomposition methods for inversion

Direct seismic inversion or RF inversion techniques based on principles of wavefield decomposition (Kennett, 1978).
Generally lower cost than Bayesian techniques but also generates less detailed information about solution space ensemble.

The method implemented here wraps the method of [Tao][1] in a Monte Carlo Markov Chain (MCMC) solver that can return
multiple candidate solution points in N-dimensional space. The user specifies a number of layers L and material
properties of each layer, plus that of the underlying mantle, and the solver tries to solve for unknowns H (thickness)
and V<sub>S</sub> in each layer. Therefore when there are L layers given, the solver solves for N = 2L unknowns.

This method is experimental, but has been validated for 1 layer (over the mantle) and given good results in some
cases for 2 layers (in good agreement with H-k stacking). More constraints may be needed to get solutions in higher
dimensions (hence 'experimental' status).

## Using the solver

The main module for running cases using this technique is [`runners.py`][2]. The are two running modes:
- single station: See `python runners.py single-job --help` for details
- batch of stations: See `python runners.py batch-job --help` for details

See example JSON files in folder [`seismic/inversion/wavefield_decomp`][3] for details.

## References

1. Kai Tao, Tianze Liu, Jieyuan Ning, Fenglin Niu, "Estimating sedimentary and crustal structure using wavefield
   continuation: theory, techniques and applications", Geophysical Journal International, Volume 197, Issue 1,
   April, 2014, Pages 443-457, https://doi.org/10.1093/gji/ggt515
2. Charles A. Langston, "Wave-Field Continuation and Decomposition for Passive Seismic Imaging Under Deep
   Unconsolidated Sediments", Bulletin of the Seismological Society of America, Vol. 101, No. 5, pp. 2176-2190,
   October 2011, doi: 10.1785/0120100299
3. Kennett, B.L.N., Kerry, N.J. & Woodhouse, J.H., "Symmetries in the reflection and transmission of elastic waves",
   Geophys. J. R. astr. Soc., 52, 215-229, 1978. https://doi.org/10.1111/j.1365-246X.1978.tb04230.x

[1]: https://doi.org/10.1093/gji/ggt515

[2]: https://github.com/GeoscienceAustralia/hiperseis/blob/develop/seismic/inversion/wavefield_decomp/runners.py

[3]: https://github.com/GeoscienceAustralia/hiperseis/tree/develop/seismic/inversion/wavefield_decomp
