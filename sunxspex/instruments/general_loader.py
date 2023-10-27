"""
Not a permanent file.

This file only exists to provide a general layout of how instrument loaders 
are to be generated.

Contains three main parts:
    * Instrument loader class
    * Instrument spectral response class
    * Instrument file reading function
"""

from base_instrument import BaseInstrument
from base_spectral_response import BaseSpectralResponse
from spectrum_data import SpectrumData

class GeneralLoader(BaseInstrument):
    """
    A class to load in an instrument's data from a given format (file or 
    otherwise). Inspection and alteration of the loaded data should also 
    be possible such as time-range selection, rebinning, energy-range 
    selection, etc. 
    """

    def __init__(self, *args):
        """
        *args could be any pre-defined number of inputs. This will depend 
        on what files or inputs the General instrument needs.
        """
        # read in, or just use, the input files and/or arguments
        _ = self.read_general(*args)

        # some work, instrument specific most likely
        ...

        # instrument spectrum data for spectral fitting or data-model comparison
        spec = list()
        # make sure to get the instrument spectral response matrix
        srm = GeneralSpectralResponse(*args)

        # get the spectrum data to the SpectrumData class
        self.event_data = SpectrumData(spec, srm)
        # additionally, if there is background data you have from the files
        self.background_data = SpectrumData(...)

    def read_general(self, *args):
        """ 
        Use the stand-alone `io` function.
        
        Stand-alone is easier for testing and troubleshooting the reading 
        of instrument files.
        """
        _ = general_io(*args)
        ...
    
class GeneralSpectralResponse(BaseSpectralResponse):
    def __init__(self, *args):
        """
        *args could be any pre-defined number of inputs. This will depend 
        on what files or inputs the General instrument needs to work with
        and create the General instrument's spectral response matrix.

        Can contain methods specific to altering the specific instrument's 
        spectral response matrix such as time-range selection, rebinning, 
        energy-range selection, etc.
        """
        # some work, instrument specific most likely
        ...

        # the 2D spectral response matrix
        self.spectral_response_matrix = args

def general_io(file):
    """ Function to read in the General instrument's data files. """
    ...