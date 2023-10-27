"""
Contains the spectrum data class.

This class is used to contain only the necessary information needed for 
spectral fitting.

Needed for fitting
    * counts 
        * _1D, label:`counts`, unit:`counts`_
    * counts error 
        * _1D (same length as counts), label:`counts_error`,  unit:`counts`_
    * count energy bin edges 
        * _1D, label:`counts_energy edges`, unit:`keV`_
    * photon energy bin edges 
        * _1D, label:`photons_energy edges`, unit:`keV`_
    * spectral response matrix 
        * _2D, label:`spectral_response_matrix`, unit:`cm^2 count photon^-1`_
    * effective exposure 
        * _scalar or 1D (same length as counts), label:`effective_exposure`, unit:`s`_
    * distribution from which data is sampled
        * _string, label:`sample_distribution`_

The SRM will be provided by the instrument specific spectral response 
class that inherits from `base_spectral_response.BaseSpectralResponse`.

All other information that is useful to have handy for spectral fitting, 
or general data-model comparison, will be calculated from the above list.
"""

import astropy.units as u
import numpy as np

class SpectrumData:
    """ 
    A class to contain the information required for spectral data-model 
    comparison. 

    Instance Attributes
    -------------------
    counts : `numpy.ndarray`
        The 1D counts array.

    counts_error : `numpy.ndarray`
        The 1D array with the associated errors for `counts`.

    counts_energy_edges : `numpy.ndarray`
        The 1D array (values in `keV`) of the count channel bin edges.
        Should be size `counts.size+1`.

    photons_energy_edges : `numpy.ndarray`
        The 1D array (values in `keV`) of the photon channel bin edges.

    spectral_response_matrix : `numpy.ndarray`
        The 2D array (values in cm^2 counts photon^-1) of the photon to 
        counts conversion. The shape of the spectral response matrix (SRM) 
        should be `(photons_energy_edges.size-1,counts_energy_edges.size-1)` 
        where array rows represent photon axis and columns the count axis.

    effective_exposure : `int`, `float`, or `numpy.ndarray`
        The single value or 1D array of size `counts.size` describing the 
        effective exposure that resulted in the data (values in `s`).

    sample_distribution : `str`
        A string describing from where the data is likely sampled. Valid 
        options are given by `SAMPLE_OPTIONS` (see Class Attributes).

    counts_energy_widths : `numpy.ndarray`
        The 1D array (values in `keV`) of the count channel bin widths.
        Should be size `counts.size`.

    photons_energy_widths : `numpy.ndarray`
        The 1D array (values in `keV`) of the photon channel bin widths.
        Should be size `photons_energy_edges.size-1`.

    Class Attributes
    ----------------
    SAMPLE_OPTIONS : `list[str]`
        The categories used to describe from where the data is likely 
        sampled.

    Methods
    -------
    _check_valid_args()
        Checks all inputs are of a valid type, dimension, and size in 
        relation to each other.
    """
    
    SAMPLE_OPTIONS = ['Gaussian', 'Poissonian']

    @u.quantity_input
    def __init__(self, 
                 counts: u.ct, 
                 counts_error: u.ct, 
                 counts_energy_edges: u.keV, 
                 photons_energy_edges: u.keV, 
                 spectral_response_matrix: u.cm**2*u.ct/u.ph,
                 effective_exposure: u.s,
                 sample_distribution: str):
        """ 
        Inputs are checked to have compatible `astropy.units` with those 
        defined for each argument.

        Parameters
        ----------
        counts : `astropy.units.quantity.Quantity`
            The 1D counts array. 
            Units: `astropy.units.ct`
        counts_error : `astropy.units.quantity.Quantity`
            The 1D array with the associated errors for `counts`.
            Units: `astropy.units.ct`
        counts_energy_edges : `astropy.units.quantity.Quantity`
            The 1D array of the count channel bin edges.
            Should be size `counts.size+1`.
            Units: `astropy.units.keV`
        photons_energy_edges :`astropy.units.quantity.Quantity`
            The 1D array of the photon channel bin edges.
            Units: `astropy.units.keV`
        spectral_response_matrix : `astropy.units.quantity.Quantity`
            The 2D array of the photon to counts conversion. Shape of 
            the SRM should be `(photons_energy_edges.size-1,counts.size)` 
            where array rows represent photon axis and columns represent
            the count axis.
            Units: `astropy.units.cm**2 * astropy.units.ct/astropy.units.ph`
        effective_exposure : `astropy.units.quantity.Quantity`
            The single value or 1D array of size `counts.size` describing 
            the effective exposure that resulted in the data.
            Units: `astropy.units.s`
        sample_distribution : `str`
            A string describing from where the data is likely sampled. 
            Valid options are given by `SAMPLE_OPTIONS` (see Class
            Attributes).
        """
        
        # check units are convertable to what's needed, convert then strip 
        # them ready for the fitter
        self.counts = (counts << u.ct).value
        self.counts_error = (counts_error << u.ct).value
        self.counts_energy_edges = (counts_energy_edges << u.keV).value
        self.photons_energy_edges = (photons_energy_edges << u.keV).value
        self.spectral_response_matrix = (spectral_response_matrix << (u.cm**2 * u.ct/u.ph)).value
        self.effective_exposure = (effective_exposure << u.s).value

        # check `sample_distribution` can be a string
        self.sample_distribution = str(sample_distribution)

        # check dimensions are correct for inputs
        self._check_valid_args()

        # calculate the energy bin widths
        self.counts_energy_widths = np.diff(self.counts_energy_edges)
        self.photons_energy_widths = np.diff(self.photons_energy_edges)

    def _check_valid_args(self):
        """ Check all the arrays are valid and consistent with each other. """

        # check arrays which should be 1D
        # include effective exposure as it can be a single value OR a 1D array
        valid1d_dim = (self.effective_exposure.ndim<=self.counts.ndim==self.counts_error.ndim==self.counts_energy_edges.ndim==self.photons_energy_edges.ndim==1)
        if not valid1d_dim:
            raise ValueError("At least one expected 1D input has the wrong dimension.")
        
        # check arrays which should be 2D
        valid2d_dim = (self.spectral_response_matrix.ndim==2)
        if not valid2d_dim:
            raise ValueError("At least one expected 2D input has the wrong dimension.")
        
        # check `sample_distribution` value is allowed
        if self.sample_distribution not in self.SAMPLE_OPTIONS:
            raise ValueError(f"The `sample_distribution` should be from {self.SAMPLE_OPTIONS}.")
        
        # check counts-related arrays are consistent sizes
        valid_counts_size = (
            (self.effective_exposure.size==self.counts.size or self.effective_exposure.size == 1) and
            (self.counts.size==self.counts_error.size==(self.counts_energy_edges.size-1)))
        if not valid_counts_size:
            raise ValueError("At least one countsrelated array is inconsistent with the others.")
        
        # check spectral response matrix size/shape is consistent with count/photon bins
        # counts axis is across SRM columns, photon axis is across rows
        valid_srm_size = (self.spectral_response_matrix.shape==(self.photons_energy_edges.size-1,self.counts_energy_edges.size-1))
        if not valid_srm_size:
            raise ValueError("The spectral response matrix shape is not consistent with the count and photon bins.")