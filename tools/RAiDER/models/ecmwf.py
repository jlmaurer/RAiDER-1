import datetime

import numpy as np
from pyproj import CRS

from RAiDER.logger import *
from RAiDER import utilFcns as util
from RAiDER.models.weatherModel import WeatherModel


class ECMWF(WeatherModel):
    '''
    Implement ECMWF models
    '''

    def __init__(self):
        # initialize a weather model
        WeatherModel.__init__(self)

        # model constants
        self._k1 = 0.776   # [K/Pa]
        self._k2 = 0.233   # [K/Pa]
        self._k3 = 3.75e3  # [K^2/Pa]

        self._lon_res = 0.2
        self._lat_res = 0.2

        self._proj = CRS.from_epsg(4326)

        self._expver = '1'
        self._classname = 'od'

    def load_weather(self, *args, **kwargs):
        '''
        Consistent class method to be implemented across all weather model types.
        As a result of calling this method, all of the variables (x, y, z, p, q,
        t, wet_refractivity, hydrostatic refractivity, e) should be fully
        populated.
        '''
        self._load_model_level(*self.files)


    def _load_model_level(self, fname):

        # read data from file
        lats, lons, xs, ys, t, q, p, h = self._makeDataCubes(
            fname,
            verbose=False
        )

        self._p = p
        self._q = q
        self._t = t
        self._lats = lats
        self._lons = lons
        self._xs = xs.copy()
        self._ys = ys.copy()
        self._zs = h


    def _fetch(self, lats, lons, time, out, Nextra=2):
        '''
        Fetch a weather model from ECMWF
        '''
        # bounding box plus a buffer
        lat_min, lat_max, lon_min, lon_max = self._get_ll_bounds(lats, lons, Nextra)

        # execute the search at ECMWF
        try:
            self._get_from_ecmwf(
                lat_min,
                lat_max,
                self._lat_res,
                lon_min,
                lon_max,
                self._lon_res,
                time,
                out
            )
        except Exception as e:
            logger.warning('Query point bounds are {}/{}/{}/{}'.format(lat_min, lat_max, lon_min, lon_max))
            logger.warning('Query time: {}'.format(time))
            logger.exception(e)

    def _get_from_ecmwf(self, lat_min, lat_max, lat_step, lon_min, lon_max,
                        lon_step, time, out):
        import ecmwfapi

        server = ecmwfapi.ECMWFDataServer()

        corrected_date = util.round_date(time, datetime.timedelta(hours=6))

        server.retrieve({
            "class": self._classname,  # ERA-Interim
            'dataset': self._dataset,
            "expver": "{}".format(self._expver),
            # They warn me against all, but it works well
            "levelist": 'all',
            "levtype": "ml",  # Model levels
            "param": "lnsp/q/z/t",  # Necessary variables
            "stream": "oper",
            # date: Specify a single date as "2015-08-01" or a period as
            # "2015-08-01/to/2015-08-31".
            "date": datetime.datetime.strftime(corrected_date, "%Y-%m-%d"),
            "type": "an", # should be reanalysis ("an")
            # time: With type=an, time can be any of
            # "00:00:00/06:00:00/12:00:00/18:00:00".  With type=fc, time can
            # be any of "00:00:00/12:00:00",
            "time": datetime.time.strftime(corrected_date.time(), "%H:%M:%S"),
            # step: With type=an, step is always "0". With type=fc, step can
            # be any of "3/6/9/12".
            "step": "0",
            # grid: Only regular lat/lon grids are supported.
            "grid": '{}/{}'.format(lat_step, lon_step),
            "area": '{}/{}/{}/{}'.format(lat_max, lon_min, lat_min, lon_max),  # area: N/W/S/E
            "format": "netcdf",
            "resol": "av",
            "target": out,    # target: the name of the output file.
        })

    def _get_from_cds(self, lat_min, lat_max, lat_step, lon_min, lon_max,
                      lon_step, acqTime, outname):
        import cdsapi

        pls = ['1', '2', '3', '5', '7', '10', '20', '30', '50', '70', '100', '125', '150', '175', '200', '225', '250', '300', '350', '400', '450', '500', '550', '600', '650', '700', '750', '775', '800', '825', '850', '875', '900', '925', '950', '975', '1000']
        mls = np.arange(137) + 1

        c = cdsapi.Client(verify=0)
        # corrected_date = util.round_date(time, datetime.timedelta(hours=6))
        if self._model_level_type == 'pl':
            var = ['geopotential', 'relative_humidity', 'specific_humidity', 'temperature']
            levels = 'all'
            levType = 'pressure_level'
        else:
            var = ['lnsp', 'q', 'z', 't']
            levels = mls
            levType = 'model_level'

        bbox = [lat_max, lon_min, lat_min, lon_max]

        dataDict = {
            "product_type": "reanalysis",
            "{}".format(levType): levels,
            "levtype": "{}".format(self._model_level_type),  # 'ml' for model levels or 'pl' for pressure levels
            'variable': var,
            "stream": "oper",
            "type": "an",
            "year": "{}".format(acqTime.year),
            "month": "{}".format(acqTime.month),
            "day": "{}".format(acqTime.day),
            "time": "{}".format(datetime.time.strftime(acqTime.time(), '%H:%M')),
            # step: With type=an, step is always "0". With type=fc, step can
            # be any of "3/6/9/12".
            "step": "0",
            "area": bbox,
            "format": "netcdf"}

        try:
            c.retrieve('reanalysis-era5-pressure-levels', dataDict, outname)
        except Exception as e:
            logger.warning('Query point bounds are {}/{} latitude and {}/{} longitude'.format(lat_min, lat_max, lon_min, lon_max))
            logger.warning('Query time: {}'.format(acqTime))
            logger.exception(e)
            raise Exception


    def _makeDataCubes(self, fname, verbose=False):
        '''
        Create a cube of data representing temperature and relative humidity
        at specified pressure levels
        '''
        # get ll_bounds
        S, N, W, E = self._ll_bounds

        with xr.open_dataset(fname) as ds:
            # Fix longitudes to be -180 - 180. 
            # Note that if lons are already in this format it will not change them, except 180 -> -180.
            ds = ds.assign_coords(longitude=(((ds.longitude + 180) % 360) - 180))
    
            # mask based on query bounds
            m1 = (S <= ds.latitude) & (N >= ds.latitude)
            m2 = (W <= ds.longitude) & (E >= ds.longitude)
            block = ds.where(m1 & m2, drop=True)
    
            # Pull the data
            z = np.squeeze(block['z'].values)[0, ...]
            t = np.squeeze(block['t'].values)
            q = np.squeeze(block['q'].values)
            lnsp = np.squeeze(block['lnsp'].values)[0, ...]
            lats = np.squeeze(block.latitude.values)
            lons = np.squeeze(block.longitude.values)
    
            xs = lons.copy()
            ys = lats.copy()
    
        if z.size == 0:
            raise RuntimeError('There is no data in z, '
                               'you may have a problem with your mask')
    
            # ECMWF appears to give me this backwards
            if lats[0] > lats[1]:
                z = z[::-1]
                lnsp = lnsp[::-1]
                t = t[:, ::-1]
                q = q[:, ::-1]
                lats = lats[::-1]
            # Lons is usually ok, but we'll throw in a check to be safe
            if lons[0] > lons[1]:
                z = z[..., ::-1]
                lnsp = lnsp[..., ::-1]
                t = t[..., ::-1]
                q = q[..., ::-1]
                lons = lons[::-1]
            # pyproj gets fussy if the latitude is wrong, plus our
            # interpolator isn't clever enough to pick up on the fact that
            # they are the same
            lons[lons > 180] -= 360
    
            geo_hgt, pres, hgt = self._calculategeoh(z, lnsp)

            # re-assign lons, lats to match heights
            _lons = np.broadcast_to(lons[np.newaxis, np.newaxis, :], hgt.shape)
            _lats = np.broadcast_to(lats[np.newaxis, :, np.newaxis], hgt.shape)

            # ys is latitude
            h = self._get_heights(_lats, hgt)
    
            # We want to support both pressure levels and true pressure grids.
            # If the shape has one dimension, we'll scale it up to act as a
            # grid, otherwise we'll leave it alone.
            if len(pres.shape) == 1:
                p = np.broadcast_to(pres[:, np.newaxis, np.newaxis], _zs.shape)
            else:
                p = pres
    
            # Re-structure everything from (heights, lats, lons) to (lons, lats, heights)
            p = np.transpose(p).swapaxes(0, 1)
            t = np.transpose(t).swapaxes(0, 1)
            q = np.transpose(q).swapaxes(0, 1)
            h = np.transpose(h).swapaxes(0, 1)
            _lats = np.transpose(_lats).swapaxes(0, 1)
            _lons = np.transpose(_lons).swapaxes(0, 1)
    
            # Flip all the axis so that zs are in order from bottom to top
            p = np.flip(p, axis=2)
            t = np.flip(t, axis=2)
            q = np.flip(q, axis=2)
            h = np.flip(h, axis=2)
            _lats = np.flip(_lats, axis=2)
            _lons = np.flip(_lons, axis=2)
    
        return lats, lons, xs, ys, t, q, p, h
