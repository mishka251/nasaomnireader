import numpy as np
from geospacepy import special_datetime

from nasaomnireader.omni_derived_var import omni_derived_var


class juliandate(omni_derived_var):
    """Julian date of timestamp"""

    def __init__(self, oi):
        omni_derived_var.__init__(self, oi)
        self.attrs['CATDESC'] = 'Julian date'
        self.attrs['UNITS'] = 'days'

    def __call__(self):
        if self.varvals is None:
            self.varvals = special_datetime.datetimearr2jd(self.oi['Epoch']).flatten()
        return self.varvals


class borovsky(omni_derived_var):
    """Borovsky solar wind coupling function"""

    def __init__(self, *args, **kwargs):
        omni_derived_var.__init__(self, *args, **kwargs)
        self.attrs['CATDESC'] = 'Borovsky Solar Wind Coupling Function'
        self.attrs['UNITS'] = 'nT km/s'

    def __call__(self):
        # Short circut if already computed
        if self.varvals is not None:
            return self.varvals

        oi = self.oi
        # Deal with names that differ between cadences
        (densvar, vswvar) = ('N', 'V') if oi.cadence == 'hourly' else ('proton_density', 'flow_speed')
        bx, by, bz = oi['BX_GSE'], oi['BY_GSM'], oi['BZ_GSM']
        n, vsw, mach = oi[densvar], oi[vswvar], oi['Mach_num']
        bt = np.sqrt(by ** 2 + bz ** 2)
        # Compute IMF clock angle
        ca = np.arctan2(by, bz)
        borovsky = 3.29e-2 * (np.sin(ca / 2) ** 2) * np.sqrt(n) * vsw ** 2 * mach ** (-.18) * np.exp(
            np.sqrt((mach / 3.42)))
        self.varvals = borovsky
        return borovsky


class newell(omni_derived_var):
    """Newell emperical solar wind coupling function"""

    def __init__(self, *args, **kwargs):
        omni_derived_var.__init__(self, *args, **kwargs)
        self.attrs['CATDESC'] = 'Newell Solar Wind Coupling Function'
        self.attrs['UNITS'] = 'm/s^(4/3) T^(2/3)'

    def __call__(self):
        # Short circut if already computed
        if self.varvals is not None:
            return self.varvals

        oi = self.oi
        # Deal with names that differ between cadences
        vswvar = 'V' if oi.cadence == 'hourly' else 'flow_speed'
        bx, by, bz = oi['BX_GSE'], oi['BY_GSM'], oi['BZ_GSM']
        vsw, mach = oi[vswvar], oi['Mach_num']
        bt = np.sqrt(by ** 2 + bz ** 2)
        # Compute IMF clock angle
        ca = np.arctan2(by, bz)
        neg_ca = bt * np.cos(ca) * bz < 0
        ca[neg_ca] = ca[net_ca] + np.pi
        sin_ca = np.abs(np.sin(ca / 2.))

        newell = (vsw * 1000.) ** (4. / 3) * (bt * 1.0e-9) ** (2. / 3) * (sin_ca) ** (8. / 3);
        self.varvals = newell
        return newell


class knippjh(omni_derived_var):
    """Knipp Joule Heating Index (Old Version)"""

    def __init__(self, *args, **kwargs):
        omni_derived_var.__init__(self, *args, **kwargs)
        self.attrs['UNITS'] = 'GW'
        self.attrs['CATDESC'] = 'Knipp Joule Heating Index'

    def __call__(self):
        # Short circut if already computed
        if self.varvals is not None:
            return self.varvals

        oi = self.oi
        # Computes the Joule heating index from Knipp, Tobiska, Emery,
        # Direct and Indirect Thermospheric Heating Sources For Solar Cycles 21-23,
        # Solar Physics
        dt = oi['Epoch'][:].flatten()
        doy = special_datetime.datetimearr2doy(dt).flatten()

        # Leap year modifier (adds one if this is a leap year)
        # ---NOTE: this implementation does not an edge cases:
        #   1. An omni interval runs between two years one of which is a leap
        # The workaround (computing a lymod value for each dt) would degrade performance
        # and the effect is small
        # so I chose not to address it. -LMK
        if np.mod(dt[0].year, 4) == 0:
            lymod = 1.
        else:
            lymod = 0.

        # Paper uses absolute value of pcn and dst
        PC = np.abs(self.oi['PC_N_INDEX']).flatten()
        Dst = np.abs(self.oi['DST' if oi.cadence == 'hourly' else 'SYM_H']).flatten()

        jhindex = np.zeros_like(PC)
        jhindex.fill(np.nan)

        # 4 Seasons
        annual = np.logical_or(doy > 335. + lymod, doy < 31.)
        jhindex[annual] = 24.89 * PC[annual] + 3.41 * PC[annual] ** 2 + .41 * Dst[annual] + .0015 * Dst[annual] ** 2

        # Winter is 21 October (294 for non-leap, 295 for leap) to 20
        # February (DOY 51 for both non-leap and leap)
        winter = np.logical_or(doy > 294. + lymod, doy < 51.)
        jhindex[winter] = 13.36 * PC[winter] + 5.08 * PC[winter] ** 2 + .47 * Dst[winter] + .0011 * Dst[winter] ** 2

        # Summer is 21 April (DOY 111, 112 leap) - 20 August (DOY 232, 233 leap)
        summer = np.logical_and(doy > 111. + lymod, doy < 232. + lymod)
        jhindex[summer] = 29.27 * PC[summer] + 8.18 * PC[summer] ** 2 - .04 * Dst[summer] + .0126 * Dst[summer] ** 2

        # Equinox is 21 Feb (DOY 51) - 20 Apr (DOY 110, 111 leap)
        # and 21 Aug (233, 234 leap) - 20 Oct (293 non-leap, 294 leap)
        equinox = np.logical_or(
            np.logical_and(doy > 51. + lymod, doy < 110. + lymod),
            np.logical_and(doy > 233. + lymod, doy < 293. + lymod)
        )

        jhindex[equinox] = 29.14 * PC[equinox] + 2.54 * PC[equinox] ** 2 + .21 * Dst[equinox] + .0023 * Dst[
            equinox] ** 2
        self.attrs['UNITS'] = 'GW'
        self.attrs['CATDESC'] = 'Knipp Joule Heating Index'
        self.varvals = jhindex
        return jhindex