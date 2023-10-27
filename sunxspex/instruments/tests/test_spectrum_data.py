import astropy.units as u
import numpy as np
from numpy.testing import assert_array_equal

from sunxspex.instruments.spectrum_data import SpectrumData


def spectrum_data_compare(spec_data_obj, 
                          counts, 
                          counts_error, 
                          counts_energy_edges, 
                          photons_energy_edges, 
                          spectral_response_matrix,
                          effective_exposure,
                          sample_distribution,
                          counts_energy_widths,
                          photons_energy_widths):
    """ Easily compare all attributes to test data. """
    # given and converted to unitless 
    assert_array_equal(spec_data_obj.counts, counts.value)
    assert_array_equal(spec_data_obj.counts_error, counts_error.value)
    assert_array_equal(spec_data_obj.counts_energy_edges, counts_energy_edges.value)
    assert_array_equal(spec_data_obj.photons_energy_edges, photons_energy_edges.value)
    assert_array_equal(spec_data_obj.spectral_response_matrix, spectral_response_matrix.value)
    assert_array_equal(spec_data_obj.effective_exposure, effective_exposure.value)
    # check string 
    assert spec_data_obj.sample_distribution==sample_distribution, "The `sample_distribution` string mismatch"
    # check bin width calculation
    assert_array_equal(spec_data_obj.counts_energy_widths, counts_energy_widths)
    assert_array_equal(spec_data_obj.photons_energy_widths, photons_energy_widths)
    

def test_spectrum_data():
    """ Check the `SpectrumData` exhibits expected behaviour. """
    # define some count data
    cts = np.arange(1,11) << u.ct
    cts_err = np.sqrt(cts.value) << u.ct

    # TEST0:
    # define the count, photon, and other observed parameters
    cts_edges0 = np.arange(cts.size+1) << u.keV
    ph_edges0 = np.arange(cts.size+2) << u.keV
    srm0 = np.zeros((ph_edges0.size-1,cts_edges0.size-1)) << u.cm**2*u.ct/u.ph
    eff_exp0 = 0.5 << u.s
    samp_dist0 = "Poissonian"
    # pass the test data to the `SpectrumData` class
    data0 = SpectrumData(counts=cts, 
                         counts_error=cts_err, 
                         counts_energy_edges=cts_edges0, 
                         photons_energy_edges=ph_edges0, 
                         spectral_response_matrix=srm0,
                         effective_exposure=eff_exp0,
                         sample_distribution=samp_dist0)
    # manually say what the count/photon bin widths should be
    cts_e_widths0 = np.array([1]*cts.size)
    ph_e_widths0 = np.array([1]*(ph_edges0.size-1))
    # collect the known, test data and compare to what `SpectrumData` gets
    test0 = (cts, cts_err, cts_edges0, ph_edges0, srm0, eff_exp0, 
             samp_dist0, cts_e_widths0, ph_e_widths0)
    spectrum_data_compare(data0, *test0)

    # TEST1:
    # define the count, photon, and other observed parameters
    cts_edges1 = np.arange(0,2*cts.size+1,2) << u.keV
    ph_edges1 = np.arange(cts.size+2) << u.keV
    srm1 = np.zeros((ph_edges1.size-1,cts_edges1.size-1)) << u.cm**2*u.ct/u.ph
    eff_exp1 = np.arange(10,cts.size+10) << u.s
    samp_dist1 = "Gaussian"
    # pass the test data to the `SpectrumData` class
    data1 = SpectrumData(counts=cts, 
                         counts_error=cts_err, 
                         counts_energy_edges=cts_edges1, 
                         photons_energy_edges=ph_edges1, 
                         spectral_response_matrix=srm1,
                         effective_exposure=eff_exp1,
                         sample_distribution=samp_dist1)
    # manually say what the count/photon bin widths should be
    cts_e_widths1 = np.array([2]*cts.size)
    ph_e_widths1 = np.array([1]*(ph_edges1.size-1))
    # collect the known, test data and compare to what `SpectrumData` gets
    test1 = (cts, cts_err, cts_edges1, ph_edges1, srm1, eff_exp1, 
             samp_dist1, cts_e_widths1, ph_e_widths1)
    spectrum_data_compare(data1, *test1)