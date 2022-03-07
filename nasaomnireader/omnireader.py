# (C) 2020 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam M. Kilcommons
import datetime
import logging
import os
import requests
import textwrap

import yadisk as yadisk
from requests import ReadTimeout
import re
import calendar

from nasaomnireader.omni_txt_cdf_mimic import omni_txt_cdf_mimic

log = logging.getLogger(__name__)

# Attempt to pull in spacepy to get pycdf interface
# to use faster CDFs
try:
    from spacepy import pycdf

    spacepy_is_available = True
except ImportError:
    # print(traceback.format_exc())
    print(textwrap.dedent("""
        ------------IMPORTANT----------------------------
        Unable to import spacepy. Will fall back to
        using Omni text files, which may have slightly
        different data and incomplete metadata
        -------------------------------------------------
        """))
    spacepy_is_available = False

# Variables in 1 Hour CDFS
# ABS_B: CDF_REAL4 [4344]
# AE: CDF_INT4 [4344]
# AL_INDEX: CDF_INT4 [4344]
# AP_INDEX: CDF_INT4 [4344]
# AU_INDEX: CDF_INT4 [4344]
# BX_GSE: CDF_REAL4 [4344]
# BY_GSE: CDF_REAL4 [4344]
# BY_GSM: CDF_REAL4 [4344]
# BZ_GSE: CDF_REAL4 [4344]
# BZ_GSM: CDF_REAL4 [4344]
# Beta: CDF_REAL4 [4344]
# DST: CDF_INT4 [4344]
# Day: CDF_INT4 [4344]
# E: CDF_REAL4 [4344]
# Epoch: CDF_EPOCH [4344]
# F: CDF_REAL4 [4344]
# F10_INDEX: CDF_REAL4 [4344]
# HR: CDF_INT4 [4344]
# IMF: CDF_INT4 [4344]
# IMF_PTS: CDF_INT4 [4344]
# KP: CDF_INT4 [4344]
# MFLX: CDF_INT4 [4344]
# Mach_num: CDF_REAL4 [4344]
# Mgs_mach_num: CDF_REAL4 [4344]
# N: CDF_REAL4 [4344]
# PC_N_INDEX: CDF_REAL4 [4344]
# PHI-V: CDF_REAL4 [4344]
# PHI_AV: CDF_REAL4 [4344]
# PLS: CDF_INT4 [4344]
# PLS_PTS: CDF_INT4 [4344]
# PR-FLX_1: CDF_REAL8 [4344]
# PR-FLX_10: CDF_REAL4 [4344]
# PR-FLX_2: CDF_REAL4 [4344]
# PR-FLX_30: CDF_REAL4 [4344]
# PR-FLX_4: CDF_REAL4 [4344]
# PR-FLX_60: CDF_REAL4 [4344]
# Pressure: CDF_REAL4 [4344]
# R: CDF_INT4 [4344]
# Ratio: CDF_REAL4 [4344]
# Rot#: CDF_INT4 [4344]
# SIGMA-ABS_B: CDF_REAL4 [4344]
# SIGMA-B: CDF_REAL4 [4344]
# SIGMA-Bx: CDF_REAL4 [4344]
# SIGMA-By: CDF_REAL4 [4344]
# SIGMA-Bz: CDF_REAL4 [4344]
# SIGMA-N: CDF_REAL4 [4344]
# SIGMA-PHI-V: CDF_REAL4 [4344]
# SIGMA-T: CDF_REAL4 [4344]
# SIGMA-THETA-V: CDF_REAL4 [4344]
# SIGMA-V: CDF_REAL4 [4344]
# SIGMA-ratio: CDF_REAL4 [4344]
# T: CDF_REAL4 [4344]
# THETA-V: CDF_REAL4 [4344]
# THETA_AV: CDF_REAL4 [4344]
# V: CDF_REAL4 [4344]
# YR: CDF_INT4 [4344]

# Variables in 5 Minute CDFs
# <CDF:
# AE_INDEX: CDF_INT4 [8928]
# AL_INDEX: CDF_INT4 [8928]
# ASY_D: CDF_INT4 [8928]
# ASY_H: CDF_INT4 [8928]
# AU_INDEX: CDF_INT4 [8928]
# BSN_x: CDF_REAL4 [8928]
# BSN_y: CDF_REAL4 [8928]
# BSN_z: CDF_REAL4 [8928]
# BX_GSE: CDF_REAL4 [8928]
# BY_GSE: CDF_REAL4 [8928]
# BY_GSM: CDF_REAL4 [8928]
# BZ_GSE: CDF_REAL4 [8928]
# BZ_GSM: CDF_REAL4 [8928]
# Beta: CDF_REAL4 [8928]
# Day: CDF_INT4 [8928]
# E: CDF_REAL4 [8928]
# Epoch: CDF_EPOCH [8928]
# F: CDF_REAL4 [8928]
# HR: CDF_INT4 [8928]
# IMF: CDF_INT4 [8928]
# IMF_PTS: CDF_INT4 [8928]
# Mach_num: CDF_REAL4 [8928]
# Mgs_mach_num: CDF_REAL4 [8928]
# Minute: CDF_INT4 [8928]
# PC_N_INDEX: CDF_REAL4 [8928]
# PLS: CDF_INT4 [8928]
# PLS_PTS: CDF_INT4 [8928]
# PR-FLX_10: CDF_REAL4 [8928]
# PR-FLX_30: CDF_REAL4 [8928]
# PR-FLX_60: CDF_REAL4 [8928]
# Pressure: CDF_REAL4 [8928]
# RMS_SD_B: CDF_REAL4 [8928]
# RMS_SD_fld_vec: CDF_REAL4 [8928]
# RMS_Timeshift: CDF_INT4 [8928]
# SYM_D: CDF_INT4 [8928]
# SYM_H: CDF_INT4 [8928]
# T: CDF_REAL4 [8928]
# Time_btwn_obs: CDF_INT4 [8928]
# Timeshift: CDF_INT4 [8928]
# Vx: CDF_REAL4 [8928]
# Vy: CDF_REAL4 [8928]
# Vz: CDF_REAL4 [8928]
# YR: CDF_INT4 [8928]
# flow_speed: CDF_REAL4 [8928]
# percent_interp: CDF_INT4 [8928]
# proton_density: CDF_REAL4 [8928]
# x: CDF_REAL4 [8928]
# y: CDF_REAL4 [8928]
# z: CDF_REAL4 [8928]
# >

from nasaomnireader import config

localdir = config['omnireader']['local_cdf_dir']



class omni_downloader(object):
    def __init__(self, yd_token, yd_dir, cdf_or_txt='cdf', force_download=False):
        self.localdir = localdir
        self.cdf_or_txt = cdf_or_txt if spacepy_is_available else 'txt'  # is set at top of file in imports
        self.force_download = force_download
        self.ftpserv = 'spdf.gsfc.nasa.gov'
        self.ftpdir = '/pub/data/omni'

        self.yd_token = yd_token
        self.yd_dir = yd_dir


        # Hourly CDF are every six months, 5 minute are every month as are 1 min
        if self.cdf_or_txt == 'cdf':
            self.cadence_subdir = {
                'hourly': 'omni_cdaweb/hourly',
                '5min': 'omni_cdaweb/hro_5min',
                '1min': 'omni_cdaweb/hro_1min'
            }
            self.filename_gen = {
                'hourly': lambda dt: '%d/omni2_h0_mrg1hr_%d%.2d01_v01.cdf' % (
                    dt.year, dt.year, 1 if dt.month < 7 else 7),
                '5min': lambda dt: '%d/omni_hro_5min_%d%.2d01_v01.cdf' % (dt.year, dt.year, dt.month),
                '1min': lambda dt: '%d/omni_hro_1min_%d%.2d01_v01.cdf' % (dt.year, dt.year, dt.month)
            }
            self.filepatterns = {
                'hourly': '\d{4}',
                '5min': '\d{4}',
                '1min': '\d{4}'
            }
            self.fileformats = {
                'hourly': '%Y',
                '5min': '%Y',
                '1min': '%Y'
            }
        elif self.cdf_or_txt == 'txt':
            self.cadence_subdir = {
                'hourly': 'low_res_omni',
                '5min': 'high_res_omni',
                '1min': 'high_res_omni/monthly_1min'
            }
            self.filename_gen = {
                'hourly': lambda dt: 'omni2_%d.dat' % (dt.year),
                '5min': lambda dt: 'omni_5min%d.asc' % (dt.year),
                '1min': lambda dt: 'omni_min%d%.2d.asc' % (dt.year, dt.month)
            }

            self.filepatterns = {
                'hourly': 'omni_m\d*.dat',
                '5min': 'omni_min\d*.asc',
                '1min': 'omni_min\d*.asc'
            }
            self.fileformats = {
                'hourly': 'omni_m%Y.dat',
                '5min': 'omni_min%Y%m.asc',
                '1min': 'omni_min%Y%m.asc'
            }

            self.filepatterns_yd = {
                'hourly': 'omni2_\d*.dat',
                '5min': 'omni_5min\d*.asc',
                '1min': 'omni_min\d*.asc'
            }

            self.fileformats_yd = {
                'hourly': 'omni2_%Y.dat',
                '5min': 'omni_5min%Y.asc',
                '1min': 'omni_min%Y%m.asc'
            }

        else:
            raise ValueError('Invalid value of cdf_or_txt argument. Valid values are "txt" and "cdf"')

    def get_response(self, url, proxy_url, proxy_key):
        response = None
        if proxy_url is not None and proxy_key is not None:
            api_key = proxy_key
            headers = {
                "apikey": api_key
            }
            params = [("url", url)]

            # log.debug("try get response content with proxy")
            get_timeout = 62.75
            try:
                response = requests.get(proxy_url, headers=headers, params=params, timeout=get_timeout)
                # print(response.status_code)

                if response.status_code >= 400:
                    raise RuntimeError(f"Ошибка запроса - ответ пришел с кодом {response.status_code}")
            except ReadTimeout as e:
                msg = f"TimeOut {str(e)} then try to get data from NASA server in {get_timeout} seconds"
                log.error(msg)
                raise RuntimeError(msg)
                # log.debug("get response content with proxy")

        else:
            # log.debug("try get response content without proxy")
            get_timeout = 12.75
            try:
                response = requests.get(url, timeout=get_timeout)
                # print(response.status_code)
            except ReadTimeout as e:
                msg = f"TimeOut {str(e)} then try to get data from NASA server in {get_timeout} seconds"
                log.error(msg)
                raise RuntimeError(msg)
        return response

    def fix_interval(self, start_dt, end_dt, cadence, proxy_url=None, proxy_key=None):
        remotefn = self.ftpdir + '/' + self.cadence_subdir[cadence] + '/'
        url = 'https://' + self.ftpserv + remotefn
        response = self.get_response(url, proxy_url, proxy_key)
        tmp = set(re.findall(self.filepatterns[cadence], response.text))
        tmp2 = [datetime.datetime.strptime(d, self.fileformats[cadence]) for d in tmp]

        min_date = min(tmp2)
        max_date = max(tmp2)

        max_date += datetime.timedelta(days=calendar.monthrange(max_date.year, max_date.month)[1])

        if start_dt < min_date:
            delta = min_date - start_dt
            return (start_dt + delta, end_dt + delta)
        elif end_dt > max_date:
            delta = end_dt - max_date + datetime.timedelta(days=1)
            return (start_dt - delta, end_dt - delta)
        else:
            return start_dt, end_dt

    def fix_interval_yadisk(self, start_dt, end_dt, cadence, **kwargs):
        y = yadisk.YaDisk(token=self.yd_token)
        yadisk_base_dir = self.yd_dir
        files = [file_inf for file_inf in y.listdir(yadisk_base_dir, fields=['name'])]
        file_names = [file_inf.name for file_inf in files]
        filtered_files = [file for file in file_names if re.match(self.filepatterns_yd[cadence], file)]
        dates = [datetime.datetime.strptime(d, self.fileformats_yd[cadence]) for d in filtered_files]

        min_date = min(dates)
        max_date = max(dates)

        max_date += datetime.timedelta(days=calendar.monthrange(max_date.year, max_date.month)[1])

        if start_dt < min_date:
            delta = min_date - start_dt
            return (start_dt + delta, end_dt + delta)
        elif end_dt > max_date:
            delta = end_dt - max_date + datetime.timedelta(days=1)
            return (start_dt - delta, end_dt - delta)
        else:
            return start_dt, end_dt

    def get_cdf(self, dt, cadence, proxy_url=None, proxy_key=None):
        # print(f"{cadence=}, {dt=}")
        remotefn = self.ftpdir + '/' + self.cadence_subdir[cadence] + '/' + self.filename_gen[cadence](dt)
        remote_path, fn = '/'.join(remotefn.split('/')[:-1]), remotefn.split('/')[-1]
        localfn = os.path.join(self.localdir, fn)
        # log.debug(f"omnireader.py:292, localfn={localfn}, remote={remote_path}")
        if not os.path.exists(localfn) or self.force_download:
            url = 'https://' + self.ftpserv + remotefn
            # log.debug(url)
            response = self.get_response(url, proxy_url, proxy_key)

            with open(localfn, 'wb') as f:
                f.write(response.content)

        if self.cdf_or_txt == 'txt':
            return omni_txt_cdf_mimic(localfn, cadence)
        elif self.cdf_or_txt == 'cdf':
            return pycdf.CDF(localfn)

    def get_cdf_from_ya_disk(self, dt, cadence, **kwargs):
        y = yadisk.YaDisk(token=self.yd_token)
        yadisk_base_dir = self.yd_dir
        fn = self.filename_gen[cadence](dt)
        remotefn = yadisk_base_dir + '/' + fn
        localfn = os.path.join(self.localdir, fn)
        # log.debug(f"omnireader.py:292, localfn={localfn}, remote={remote_path}")
        if not os.path.exists(localfn) or self.force_download:
            y.download(remotefn, localfn)

        if self.cdf_or_txt == 'txt':
            return omni_txt_cdf_mimic(localfn, cadence)
        elif self.cdf_or_txt == 'cdf':
            return pycdf.CDF(localfn)


if __name__ == '__main__':
    pass
