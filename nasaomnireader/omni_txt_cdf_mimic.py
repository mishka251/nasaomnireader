import logging

import numpy as np
import pandas as pd

from geospacepy import special_datetime

from nasaomnireader import omnitxtcdf
from nasaomnireader.omni_txt_cdf_mimic_var import omni_txt_cdf_mimic_var

log = logging.getLogger(__name__)


class omni_txt_cdf_mimic(object):
    """
    A class to make reading from a text file emulate
    the inteface of pycdf.CDF instance, so I don't
    have to clutter up the rest of the code with
    alternate versions for txt or cdf
    """

    def __init__(self, omnitxt, cadence):
        self.txtfn = omnitxt
        self.cadence = cadence
        try:
            # self.data_old = np.genfromtxt(omnitxt)
            self.data = pd.read_csv(omnitxt, sep='\s+').to_numpy()
        except Exception as ex:
            log.error(str(ex))
            print(f"Reading from {omnitxt} error = {ex}")
            raise ex

        # Load the dictionaries that map CDF variable names in
        # the omni CDFs to columns in the text files
        cdfvars_meta = omnitxtcdf.metadata[cadence]['vars']
        self.vars = {varname: omni_txt_cdf_mimic_var(varname, cdfvars_meta[varname], self.data, cadence) for varname in
                     cdfvars_meta}
        self.attrs = omnitxtcdf.metadata[cadence]['attrs']
        # Compute the equivalent to the CDF variable'Epoch', i.e. the time
        # of each observation as an array of datetimes
        year, doy = self.vars['YR'][:], self.vars['Day'][:]
        if 'HR' in self.vars:
            doy += self.vars['HR'][:] / 24.
        if 'Minute' in self.vars:
            doy += self.vars['Minute'][:] / 24. / 60.
        epoch_vardict = {'column': -1, 'attrs': {'FILLVAL': np.nan}}
        epoch = special_datetime.doyarr2datetime(doy, year).flatten()
        self.vars['Epoch'] = omni_txt_cdf_mimic_var('Epoch', epoch_vardict, epoch, cadence, data_is_column=True)

    def __getitem__(self, var):
        try:
            data = self.vars[var]
        except KeyError:
            print(self.vars.keys())
            raise
        return data