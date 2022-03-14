import datetime

import numpy as np
from geospacepy import special_datetime

from nasaomnireader.omni_interval import omni_interval


class omni_interval_delay_smooth(object):
    def __init__(self, startdt, enddt, cadence, yd_token, yd_dir, delay_mins=10, avg_mins=45, proxy_url=None, proxy_key=None):
        """
        Create lagged and smoothed 1-minute omniweb solar wind data
        appropriately for driving CS10
        """
        self.startdt = startdt
        self.enddt = enddt
        self.delay_mins = delay_mins
        self.avg_mins = avg_mins
        total_lag = datetime.timedelta(minutes=delay_mins + avg_mins + 1)
        delayed_startdt = startdt - total_lag
        self.delayed_startdt = delayed_startdt
        self.oi = omni_interval(delayed_startdt, enddt, cadence, yd_token, yd_dir, proxy_url=proxy_url, proxy_key=proxy_key)
        self.dts = self.oi['Epoch']
        self.jds = special_datetime.datetimearr2jd(self.dts).flatten()

    def __getitem__(self, varname):
        if varname == 'Epoch':
            return self.dts
        else:
            return self._lagged_smoothed(self.oi[varname])

    def _mins2elements(self, jd, n_mins):
        """
        Calculate the number of elements
        of a time variable jd that has
        units of days that corresponds
        to n_mins minutes.
        """
        delta_t = np.nanmedian(np.diff(jd * 24. * 60.))
        if delta_t <= 0.:
            raise RuntimeError('Negative or zero delta-t'
                               + ' (%f minutes) for avg/lag' % (delta_t))
        n_elements = int(np.round(n_mins / delta_t))
        # print('%d minutes convertes to %d elements' % (n_mins,n_elements))
        return n_elements

    def _delay(self, jd, y, n_mins):
        """
        Lag an array n_mins minutes
        relative to a time variable in units of days (jd)
        """
        # Calculate number of elements to lag by
        n_elements = self._mins2elements(jd, n_mins)
        y_delay = np.roll(y, -1 * n_elements)
        y_delay[:n_elements] = np.nan
        return y_delay

    def _backward_smooth(self, jd, y, n_mins):
        """
        Backward smooth a variable y for n_mins minutes using
        a time variable jd in units of days
        """
        n_elements = self._mins2elements(jd, n_mins)
        y_smooth = np.zeros_like(y)
        lagged_ys = [np.roll(y, i) for i in range(n_elements)]
        lagged_ys_arr = np.column_stack(lagged_ys)
        y_smooth = np.nanmean(lagged_ys_arr, axis=1)
        y_smooth[:n_elements] = np.nan  # Clean up wraparound
        return y_smooth

    def _lagged_smoothed(self, y):
        """
        Use class variables to define
        lag and smoothing interval and
        jd
        """
        return self._backward_smooth(self.jds,
                                     self._delay(self.jds, y, self.delay_mins),
                                     self.avg_mins)