# Synthetic Seismogram Generation

This package provides modules for generating synthetic seismic waveforms for the purposes
of testing and algorithm validation. In particular, the [`synth`][3] module provides a standard
interface for producing synthetic seismograms from known earth model.

## Synth module

The `synth` module provides a single point of entry, function [`synth.synthesize_dataset`][3], for synthesizing
seismograms using provided underlying synthesis methods. This function takes the following inputs to generate
an output dataset in HDF5 format with event-based indexing:
- Synthesis method: 'propmatrix' or 'syngine'
- Name of output file
- Receiver station location (lat/lon)
- Source event location(s) (lat/lon)
- Sampling rate (Hz)
- Time window about onset (sec)

When method 'propmatrix' is selected, the propagation matrix method is used via the [`telewavesim`][1] library.
With this method, the user must also supply the 1D earth model to simulate. This is done by adding details to the
`kwargs` passed to `synth.synthesize_dataset` with the aid of [`class LayerProps`][6]. This is quite straightforward
and documented by example in file [`example_synth.py`][4]. This method is good for computing idealized seismograms
for an arbitrary 1D earth model provided by the user.

When method 'syngine' is used, the [IRIS Synthetics Engine (Syngine)][5] is used. This is a web service that uses
precalculated Green's Functions to compute a received seismogram according to the selected standard earth model
(e.g. IASP91). This method requires an internet connection. This method is good for computing more realistic
synthetic seismograms matching a standard 1D earth model.

## Custom synthesis methods

User-define synthesis methods can be added by creating a new class that inherits from `class Synthesizer`
and implements the abstract `synthesize` function signature. Any new synthesis methods can be made available
to all by adding support to the `synth.synthesize_dataset` function.

## References

Telewavesim [library][1] and [publication on JOSS][2].


[1]: (https://zenodo.org/badge/latestdoi/204565459)

[2]: (https://joss.theoj.org/papers/10.21105/joss.01818)

[3]: https://github.com/GeoscienceAustralia/hiperseis/blob/develop/seismic/synthetics/synth.py

[4]: https://github.com/GeoscienceAustralia/hiperseis/blob/develop/seismic/synthetics/example_synth.py

[5]: http://ds.iris.edu/ds/products/syngine/

[6]: https://github.com/GeoscienceAustralia/hiperseis/blob/develop/seismic/model_properties.py

[7]: https://github.com/GeoscienceAustralia/hiperseis/blob/develop/seismic/synthetics/backends/synthesizer_base.py

