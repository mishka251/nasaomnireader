import numpy as np


class omni_txt_cdf_mimic_var(object):
    """
    A class to mimic the interface to a CDF
    variable
    """

    def __init__(self, name, vardict, data, cadence, data_is_column=False):
        # Column of text data that
        # is the same as this variable
        self.name = name
        self.cadence = cadence
        self.column = vardict['column']

        if not data_is_column:
            self.data = data[:, int(vardict['column'])]
        else:
            self.data = data

        if 'attrs' in vardict:
            self.attrs = vardict['attrs']
        else:
            self.attrs = {'FILLVAL': np.nan}
            self.attrs['FILLVAL'] = self.identify_fill()

        # if 'FILLVAL' in self.attrs and np.count_nonzero(self.data==self.attrs['FILLVAL'])>2:
        #   pass
        # else:
        #   self.identify_fill()

    def identify_fill(self, debug=False):
        if debug:
            print("Defilling data from %s (column %d)..." % (self.name, self.column))
        # Convert fill to NaN by testing for presence of all possible fills
        possible_fills_no_decimal = ['999', '9999', '99999', '999999', '9999999', '99999999', \
                                     '999999999', '9999999999']
        possible_fills_no_decimal.reverse()  # To prevent false fill identification, we must go
        # from largest to smallest
        # last_time_column = 3 if self.cadence is 'hourly' else 4 #we don't want time columns to be defilled
        for fill_no_decimal in possible_fills_no_decimal:
            # Check just the integer
            this_fill = int(fill_no_decimal)
            isfill = self.data == this_fill
            nfill = np.count_nonzero(isfill)
            if nfill > 2.:
                if debug:
                    print("Found %d instances of integer fill value %f" % (nfill, this_fill))
                self.data[isfill] = np.nan
                continue
            # Check all possible decimal locations
            this_fill_chars = list(fill_no_decimal)  # split into list of characters
            for k in range(len(this_fill_chars)):
                # A possible fill must begin and end with 9,
                # but can have a decimal in any of the intermediate
                # values (e.g. 999 -> 9.9, .99, 99.)
                this_fill_chars_with_dec = [this_fill_chars[i] if i != k else '.' for i in range(len(this_fill_chars))]
                this_fill = float(''.join(this_fill_chars_with_dec))
                isfill = self.data == this_fill
                nfill = np.count_nonzero(isfill)
                if nfill > 2.:
                    if debug:
                        print("Found %d instances of float fill value %f for column %d" % (nfill, this_fill))
                    self.data[isfill] = np.nan
                    break
        expected_fill = '<missing>' if 'FILLVAL' not in self.attrs else self.attrs['FILLVAL']
        print("Fillval for %s (column %d) was identified as %f, tabluated as %s" % (
        self.name, self.column, this_fill, str(expected_fill)))
        return this_fill

    def _nan_fill_datapoints(self, vardata):
        fillval = self.attrs['FILLVAL']
        if np.isfinite(fillval):
            probably_fill = np.isclose(vardata, fillval, rtol=0., atol=1.)
            n_fill, n_total = np.count_nonzero(probably_fill), len(vardata)
            # print('{} NaNd {}/{} (fill was {})'.format(self.name,
            # 											n_fill,
            # 											n_total,
            # 											fillval))
            vardata[probably_fill] = np.nan
        return vardata

    def __getitem__(self, *args):
        vardata = self.data.__getitem__(*args)
        return self._nan_fill_datapoints(vardata)