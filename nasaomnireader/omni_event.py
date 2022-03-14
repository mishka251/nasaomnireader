import numpy as np
from geospacepy import special_datetime
from scipy import interpolate as interpolate

from nasaomnireader.omni_interval import omni_interval


class omni_event(object):
    def __init__(self, startdt, enddt, yd_token, yd_dir, label=None, cadence='5min', cdf_or_txt='cdf', proxy_url=None, proxy_key=None):
        self.interval = omni_interval(startdt, enddt, cadence, yd_token, yd_dir, cdf_or_txt=cdf_or_txt, proxy_url=proxy_url,
                                      proxy_key=proxy_key)
        datetime2doy = lambda \
            dt: dt.timetuple().tm_yday + dt.hour / 24. + dt.minute / 24. / 60. + dt.second / 86400. + dt.microsecond / 86400. / 1e6
        self.doy = special_datetime.datetimearr2doy(self.interval['Epoch'])
        self.jd = special_datetime.datetimearr2jd(self.interval['Epoch'])
        self.label = '%s-%s' % (startdt.strftime('%m-%d-%Y'), enddt.strftime('%m-%d-%Y')) if label is None else label
        self.interpolants = dict()
        self.attrs = self.interval.attrs

    def __getitem__(self, *args):
        return self.interval.__getitem__(*args)

    def get_var_attr(self, var, att):
        """Get a variable attribute from the last CDF in the interval"""
        return self.interval.get_var_attr(var, att)

    def interpolate(self, var, jd, **kwargs):
        """Interpolate a variable to the julian dates in jd"""
        # Create an interpolant if we don't have one yet
        if var not in self.interpolants:
            print("No interpolant for variable %s, creating %d point interpolant" % (var, len(self.jd.flatten())))
            t, y = self.jd.flatten(), self.interval[var].flatten()
            g = np.isfinite(y)
            self.interpolants[var] = interpolate.PchipInterpolator(t[g], y[g])
        # Return the interpolated result
        return self.interpolants[var].__call__(jd, **kwargs)

    def close(self):
        """Close the CDFs"""
        for cdf in self.interval.cdfs:
            cdf.close()