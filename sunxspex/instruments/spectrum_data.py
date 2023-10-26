"""
Contains the spectrum data class.

This class is used to contain only the necessary information needed for 
spectral fitting.

Needed for fitting
	* counts, 1D, `counts`
	* error, 1D (same length as counts) `counts_error`
	* count bin edges, 1D `[keV]` or convert to` [keV]` via `<<`
	* photon bin edges, 1D `[keV]` or convert to `[keV]` via `<<`
	* SRM, 2D `[cm^2 count photon^-1]`
	* effective exposure `[s]`

The SRM will be provided by the instrument specific spectral response 
class that inherits from `base_spectral_response.SpectralResponseBase`.

All other information that is useful to have handy for spectral fitting, 
or general data-model comparison, will be calculated from the above list.
"""

class SpectrumData:
    pass
