import datetime
import os

import matplotlib as mpl
import numpy as np
from geospacepy import special_datetime
from matplotlib import pyplot as pp

from nasaomnireader.omni_event import omni_event


class omni_sea(object):
    def __init__(self, center_ymdhm_list, name=None, ndays=3, cadence='5min', cdf_or_txt='cdf'):
        """
        A class for superposed epoch analysis
        """
        self.ndays = ndays
        self.name = name  # Name of the list
        self.nevents = len(center_ymdhm_list)
        self.center_dts = [datetime.datetime(y, mo, d, h) + datetime.timedelta(minutes=m) for [y, mo, d, h, m] in
                           center_ymdhm_list]
        self.center_jds = special_datetime.datetimearr2jd(self.center_dts).flatten()
        self.cadence = cadence
        # Create an omnidata interval for each event
        self.events = [
            omni_event(center_dt - datetime.timedelta(days=ndays), center_dt + datetime.timedelta(days=ndays),
                       cadence=cadence, cdf_or_txt=cdf_or_txt) for center_dt in self.center_dts]
        # mirror the attributes of the first event's last CDF
        self.attrs = self.events[0].attrs

    def plot_individual(self, ax, var, show=False, cmap=None, text_kwargs=None, plot_kwargs=None):
        """ Plot all individual intervals labeled at their maxima"""
        if cmap == None:
            norm = mpl.colors.Normalize(vmin=0, vmax=len(self.center_jds))
            cmap = mpl.cm.get_cmap('hot')
        for i, (event, center_jd) in enumerate(zip(self.events, self.center_jds)):
            t, y = event.jd - center_jd, event[var]
            maxind = np.nanargmax(np.abs(y))
            # Allow for optional arguments to text and plot
            plot_kwargs = dict() if plot_kwargs is None else plot_kwargs
            ax.plot(t, y, '-', color=cmap(norm(i)), **plot_kwargs)
            text_kwargs = {'backgroundcolor': 'grey', 'alpha': .7,
                           'ha': 'center'} if text_kwargs is None else text_kwargs
            ax.text(t[maxind], y[maxind], event.label, **text_kwargs)
        if show:
            pp.show()
            pp.pause(10)

    def get_var_attr(self, var, att):
        """Get a variable attribute from the last CDF in the interval"""
        return self.events[0].get_var_attr(var, att)

    def plot_stats(self, ax, var, x=None, xstep=0.042, show=False, plot_events=False, **kwargs):
        """Default to having 1 hour steps for interpolation"""
        if x is None:
            x = np.arange(-1 * self.ndays, self.ndays + xstep, xstep)

        iy = np.zeros((len(self.center_jds), len(x)))
        for i, (event, center_jd) in enumerate(zip(self.events, self.center_jds)):
            t, y = event.jd - center_jd, event[var]
            iy[i, :] = event.interpolate(var, x + center_jd)

            # Plot points
            # Get len(x) random numbers from -.5 - .5
            jitter = np.random.rand(*x.flatten().shape) - .5
            # Scale
            jitter = jitter * xstep / 3
            # Plot jittered points
            if plot_events:
                ax.plot(x.flatten() + jitter.flatten(), iy[i, :].flatten(), '.',
                        color='b' if 'color' not in kwargs else kwargs['color'],
                        zorder=5., alpha=.1)
        y_med, y_lq, y_uq = np.nanmedian(iy, axis=0), np.nanpercentile(iy, 25, axis=0), np.nanpercentile(iy, 75, axis=0)
        lab = '' if self.name is None else '%s: ' % (self.name)
        lab += 'Median %s Response' % (var)
        ax.plot(x, y_med, label=lab, linestyle='-', zorder=10, **kwargs)
        ax.plot(x, y_lq, linestyle=':', zorder=10, **kwargs)
        ax.plot(x, y_uq, linestyle=':', zorder=10, **kwargs)
        # Put units on the y axis
        un = self.get_var_attr(var, 'UNITS')
        un = '' if un is None else '[%s]' % (un)
        ax.set_ylabel(var + un)
        ax.legend()
        ax.fill_between(x, y_lq, y_uq, color='b' if 'color' not in kwargs else kwargs['color'], alpha=.1, zorder=7)
        if show:
            pp.show()
            pp.pause(10)

    def dump_stats(self, var, csvdir, csvfn=None, x=None, xstep=0.042):
        """
        Writes superposed epoch analysis results to a CSV file
        """
        if x is None:
            x = np.arange(-1 * self.ndays, self.ndays + xstep, xstep)

        iy = np.zeros((len(self.center_jds), len(x)))
        for i, (event, center_jd) in enumerate(zip(self.events, self.center_jds)):
            t, y = event.jd - center_jd, event[var]
            iy[i, :] = event.interpolate(var, x + center_jd)

        y_med, y_lq, y_uq = np.nanmedian(iy, axis=0), np.nanpercentile(iy, 25, axis=0), np.nanpercentile(iy, 75, axis=0)
        header = '' if self.name is None else '# %s: \n' % (self.name)
        header += '# Omni Cadence: %s\n' % (self.cadence)
        header += '# First Event: %s\n' % (self.events[0].label)
        header += '# Last Event: %s\n' % (self.events[-1].label)
        header += '# Generated: %s\n' % (datetime.datetime.now().strftime('%c'))
        header += '# Column 1: Time since center time / zero epoch hour [days] \n'
        header += '# Column 2: 25th Percentile / 1st Quartile of %s \n' % (var)
        header += '# Column 3: 50th Percentile / Median of %s \n' % (var)
        header += '# Column 4: 75th Percentile / 3rd Quartile of %s \n' % (var)

        # Characters to remove from filename
        whitepunc = [' ', ',', '/', ':', ';']
        if csvfn is None:
            csvfn = '%s_%s_stats.csv' % (self.name, var)

        # Remove non-filesystem characters
        for ch in whitepunc:
            csvfn = csvfn.replace(ch, '_')

        with open(os.path.join(csvdir, csvfn), 'w') as f:
            f.write(header)
            # print(header)
            for i in range(len(t.flatten())):
                ln = '%.5f,%e,%e,%e\n' % (t[i], y_lq[i], y_med[i], y_uq[i])
                f.write(ln)
                # print(ln)