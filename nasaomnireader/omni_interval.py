import datetime

import numpy as np

from nasaomnireader.omnireader import omni_downloader
from nasaomnireader.utils import juliandate, borovsky, newell, knippjh


class omni_interval(object):
    def __init__(self, startdt, enddt, cadence, yd_token, yd_dir, silent=False, cdf_or_txt='cdf', force_download=False, proxy_url=None,
                 proxy_key=None):
        # log.debug("omnireader.py:482")
        # Just handles the possiblilty of having a read running between two CDFs
        self.dwnldr = omni_downloader(yd_token, yd_dir, cdf_or_txt=cdf_or_txt, force_download=force_download)
        self.silent = silent  # No messages
        self.cadence = cadence
        self.startdt = startdt
        self.enddt = enddt
        # log.debug("omnireader.py:489")

        self.startdt, self.enddt = self.dwnldr.fix_interval_yadisk(self.startdt,self.enddt, cadence, proxy_url=proxy_url, proxy_key=proxy_key)

        self.cdfs = [self.dwnldr.get_cdf_from_ya_disk(self.startdt, cadence, proxy_url=proxy_url, proxy_key=proxy_key)]
        # log.debug("omnireader.py:491")
        self.attrs = self.cdfs[-1].attrs  # Mirror the global attributes for convenience
        self.transforms = dict()  # Functions which transform data automatically on __getitem__
        # Find the index corresponding to the first value larger than startdt
        self.si = np.searchsorted(self.cdfs[0]['Epoch'][:], self.startdt)
        # log.debug("omnireader.py:496")
        while self.cdfs[-1]['Epoch'][-1] < self.enddt:
            # Keep adding CDFs until we span the entire range
            # log.debug("omnireader.py:499")
            self.cdfs.append(self.dwnldr.get_cdf_from_ya_disk(self.cdfs[-1]['Epoch'][-1] + datetime.timedelta(days=1), cadence,
                                                 proxy_url=proxy_url, proxy_key=proxy_key))
        # Find the first index larger than the enddt in the last CDF
        # log.debug("omnireader.py:502")
        self.ei = np.searchsorted(self.cdfs[-1]['Epoch'][:], self.enddt)

        if not self.silent:
            print("Created interval between %s and %s, cadence %s, start index %d, end index %d" % (
            self.startdt.strftime('%Y-%m-%d'),
            self.enddt.strftime('%Y-%m-%d'), self.cadence, self.si, self.ei))
        self.add_transform('KP', ['hourly'], lambda x: x / 10., 'Hourly Kp*10 -> Kp')
        # Implement computed variables
        self.computed = dict()
        self.computed['juliandate'] = juliandate(self)
        self.computed['borovsky'] = borovsky(self)
        self.computed['newell'] = newell(self)
        self.computed['knippjh'] = knippjh(self)

        _vars = [
            'BX_GSE',
            'BY_GSM',
            'BZ_GSM',
        ]

        nan_var = None

        for var in _vars:
            values = self[var]
            if any([np.isnan(val) for val in  values]):
                nan_var = var
                break

        if nan_var is not None:
            if len(self.cdfs) == 1:
                fill = self.cdfs[-1][nan_var].attrs['FILLVAL']
                int_len = self.ei - self.si
                values = self.cdfs[0][nan_var][:]
                is_nan = np.isnan(values)

                ei = self.ei
                while is_nan[ei] and ei>0:
                    ei-=1

                if ei< int_len:
                    timespan = self.enddt - self.startdt

                    self.enddt = self.cdfs[0]['Epoch'][0] - datetime.timedelta(days=2)
                    self.startdt = self.enddt - timespan

                    dwnldr = omni_downloader(yd_token, yd_dir, cdf_or_txt=cdf_or_txt, force_download=True)

                    self.cdfs = [dwnldr.get_cdf_from_ya_disk(self.startdt, cadence, proxy_url=proxy_url, proxy_key=proxy_key)]
                    # log.debug("omnireader.py:491")
                    self.attrs = self.cdfs[-1].attrs  # Mirror the global attributes for convenience
                    self.transforms = dict()  # Functions which transform data automatically on __getitem__
                    # Find the index corresponding to the first value larger than startdt
                    self.si = np.searchsorted(self.cdfs[0]['Epoch'][:], self.startdt)
                    # log.debug("omnireader.py:496")
                    while self.cdfs[-1]['Epoch'][-1] < self.enddt:
                        # Keep adding CDFs until we span the entire range
                        # log.debug("omnireader.py:499")
                        self.cdfs.append(
                            self.dwnldr.get_cdf_from_ya_disk(self.cdfs[-1]['Epoch'][-1] + datetime.timedelta(days=1), cadence,
                                                proxy_url=proxy_url, proxy_key=proxy_key))
                    # Find the first index larger than the enddt in the last CDF
                    # log.debug("omnireader.py:502")
                    self.ei = np.searchsorted(self.cdfs[-1]['Epoch'][:], self.enddt)

                    if not self.silent:
                        print("Created interval between %s and %s, cadence %s, start index %d, end index %d" % (
                            self.startdt.strftime('%Y-%m-%d'),
                            self.enddt.strftime('%Y-%m-%d'), self.cadence, self.si, self.ei))
                    self.add_transform('KP', ['hourly'], lambda x: x / 10., 'Hourly Kp*10 -> Kp')
                    # Implement computed variables
                    self.computed = dict()
                    self.computed['juliandate'] = juliandate(self)
                    self.computed['borovsky'] = borovsky(self)
                    self.computed['newell'] = newell(self)
                    self.computed['knippjh'] = knippjh(self)

                else:
                    si = ei-int_len

                    self.si = si
                    self.ei = ei

                    self.startdt = self.cdfs[0]['Epoch'][self.si]
                    self.enddt = self.cdfs[0]['Epoch'][self.ei]

                    self.computed['juliandate'] = juliandate(self)
                    self.computed['borovsky'] = borovsky(self)
                    self.computed['newell'] = newell(self)
                    self.computed['knippjh'] = knippjh(self)

    def get_var_attr(self, var, att):
        """Get a variable attribute"""
        if var in self.computed:
            # print(f'1 {var}, {att}, {self.computed[var].attrs[att]}')
            return self.computed[var].attrs[att]
        elif att in self.cdfs[-1][var].attrs:
            return self.cdfs[-1][var].attrs[att]
        else:
            return None

    def __getitem__(self, cdfvar):
        # If it's a derived variable go get it
        # with it's own __call__ method
        if cdfvar in self.computed:
            # print('1')
            return self.computed[cdfvar]()
        # print(self.cdfs)
        # Attempt the getitem on all the cdfs in order
        if len(self.cdfs) > 1:
            ret = []
            for icdf, cdf in enumerate(self.cdfs):
                if icdf == 0:
                    ret.append(cdf[cdfvar][self.si:])
                elif icdf == len(self.cdfs) - 1:
                    ret.append(cdf[cdfvar][:self.ei])
                else:
                    ret.append(cdf[cdfvar][:])
            data = np.concatenate(ret)
        else:
            data = self.cdfs[-1][cdfvar][self.si:self.ei]
        # Fix the fill values
        try:
            if np.isfinite(self.cdfs[-1][cdfvar].attrs['FILLVAL']):
                filled = data == self.cdfs[-1][cdfvar].attrs['FILLVAL']
                if np.count_nonzero(filled) > 0:
                    data[filled] = np.nan
        except:
            print("Unhandled fill value %s for variable %s" % (self.cdfs[-1][cdfvar].attrs['FILLVAL'], cdfvar))
        # Check for transforms which need to be performed
        if cdfvar in self.transforms:
            transform = self.transforms[cdfvar]
            if self.cadence in transform['cadences']:
                if not self.silent:
                    print("Applying transform %s to omni %s variable %s" % (transform['desc'], self.cadence,
                                                                            cdfvar))
                # print "Data before", data
                data = transform['fcn'](data)
                # print "Data after", data
        # print('2')
        return data

    def add_transform(self, cdfvar, cadences, fcn, desc):
        """
            Call some function to manipulate the returned data
            whenever a particular CDF variable is '__getitem__'d
            Added to fix the obnoxious hourly 'KP' variable being Kp*10
            source of confusion.

            Arguments:
                cdfvar - str
                    name of cdf variable
                cadences - list
                    cadences which have this variable
                fcn - function
                    function to call to manipulate data
                desc - description of manipulation performed

        """
        self.transforms[cdfvar] = {'cadences': cadences, 'fcn': fcn, 'desc': desc}

    def __str__(self):
        return 'oi '+str(self.cdfs[0])