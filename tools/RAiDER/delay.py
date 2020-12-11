#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Authors: Jeremy Maurer, Raymond Hogenson, David Bekaert & Yang Lei
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import h5py
import os

import numpy as np

import RAiDER.delayFcns

from RAiDER.constants import _STEP, _ZREF, Zenith
from RAiDER.interpolator import interp_along_axis
from RAiDER.dem import getHeights
from RAiDER.logger import *
from RAiDER.losreader import getLookVectors
from RAiDER.processWM import prepareWeatherModel
from RAiDER.ioFcns import (
    make_weather_model_filename, 
    writeDelays, 
    writePnts2HDF5
)


def computeDelay(
        weather_model_file_name, 
        pnts_file_name, 
        useWeatherNodes=False,
        zlevels=None,
        zref=_ZREF, 
        step=_STEP,
        out=None, 
    ):
    """
    Calculate troposphere delay using a weather model file and query 
    points file. 
    """
    logger.debug('Beginning delay calculation')
    logger.debug('Max integration height is {:1.1f} m'.format(zref))
    logger.debug('Reference integration step is {:1.1f} m'.format(step))

    # If weather model nodes only are desired, the calculation is very quick
    if useWeatherNodes:
        # Get the weather model data
        with h5py.File(weather_model_file_name, 'r') as f:
            zs_wm = f['z'][()].copy()
            total_wet = f['wet_ztd'][()].copy()
            total_hydro = f['hydro_ztd'][()].copy()
        if zlevels is None:
            return total_wet, total_hydro
        else:
            wet_delays = interp_along_axis(zs_wm, zlevels, total_wet, axis=-1)
            hydro_delays = interp_along_axis(zs_wm, zlevels, total_hydro, axis=-1)
            return wet_delays, hydro_delays

    else:
        RAiDER.delayFcns.calculate_rays(
            pnts_file_name, 
            step
        )

        wet, hydro = RAiDER.delayFcns.get_delays(
            step,
            pnts_file_name, 
            weather_model_file_name,
        )

        logger.debug('Finished delay calculation')

        return wet, hydro


def tropo_delay(args):
    """
    raiderDelay main function.
    """

    # unpacking the dictionairy
    los = args['los']
    lats = args['lats']
    lons = args['lons']
    ll_bounds = args['ll_bounds']
    heights = args['heights']
    flag = args['flag']
    weather_model = args['weather_model']
    wmLoc = args['wmLoc']
    zref = args['zref']
    outformat = args['outformat']
    time = args['times']
    out = args['out']
    download_only = args['download_only']
    wetFilename = args['wetFilenames']
    hydroFilename = args['hydroFilenames']
    
    # logging
    logger.debug('Starting to run the weather model calculation')
    logger.debug('Time type: {}'.format(type(time)))
    logger.debug('Time: {}'.format(time.strftime('%Y%m%d')))
    logger.debug('Flag type is {}'.format(flag))
    logger.debug('DEM/height type is "{}"'.format(heights[0]))
    
    # Flags
    useWeatherNodes = flag == 'bounding_box'
    delayType = ["Zenith" if los is Zenith else "LOS"]

    # location of the weather model files
    logger.debug('Beginning weather model pre-processing')
    logger.debug('Download-only is {}'.format(download_only))

    if wmLoc is None:
        wmLoc = os.path.join(out, 'weather_files')
        
    # weather model calculation    
    wm_filename = make_weather_model_filename(weather_model['name'], time, ll_bounds)   
    weather_model_file = os.path.join(wmLoc, wm_filename)

    if not os.path.exists(weather_model_file):
        weather_model, lats, lons = prepareWeatherModel(
            weather_model, wmLoc, lats=lats, lons=lons, los=los, zref=zref,
            time=time, download_only=download_only, makePlots=True
        )
        
        if download_only:
            return None, None

        try:
            weather_model.write2HDF5(weather_model_file)
        except Exception:
            logger.exception("Unable to save weathermodel to file")

        del weather_model
    else:
        logger.warning(
            'Weather model already exists, please remove it ("%s") if you want '
            'to create a new one.', weather_model_file
        )



    if download_only:
        return None, None



    # Pull the DEM.
    logger.debug('Beginning DEM calculation')
    in_shape = lats.shape
    
    lats, lons, hgts = getHeights(lats, lons, heights, useWeatherNodes)

    pnts_dir = os.path.join(out, 'geom')
    os.makedirs(os.path.abspath(pnts_dir, exists_ok = True)
    pnts_file = os.path.join(pnts_dir, 'query_points.h5')
    if not useWeatherNodes:
        zlevels = None
        if not os.path.exists(pnts_file):
            # Convert the line-of-sight inputs to look vectors
            logger.debug('Lats shape is {}'.format(lats.shape))
            logger.debug(
                'lat/lon box is %f/%f/%f/%f (SNWE)',
                np.nanmin(lats), np.nanmax(lats), np.nanmin(lons), np.nanmax(lons)
            )
            logger.debug(
                'DEM height range is %.2f-%.2f m',
                np.nanmin(hgts), np.nanmax(hgts)
            )
            logger.debug('Beginning line-of-sight calculation')
            los = getLookVectors(los, lats, lons, hgts, zref)

            # write to an HDF5 file
            query_shape = lats.shape
            checkLOS(los, np.prod(query_shape))
            chunkSize = getChunkSize(query_shape)
            NDV = 0.

            write(
                {
                    'lat': {
                        'data': lats,
                        'attrs': {
                            'Shape': query_shape,
                            'standard_name': np.string_('latitude'),
                            'units': np.string_('degrees_east'),
                            'ChunkSize': chunkSize,
                        },
                    },
                    'lon': {
                        'data': lons,
                        'attrs': {
                            'Shape': query_shape,
                            'standard_name': np.string_('longitude'),
                            'units': np.string_('degrees_north'),
                            'ChunkSize': chunkSize,
                        },
                    },
                    'hgt': {
                        'data': hgts, 
                        'attrs': {
                            'Shape': query_shape,
                            'standard_name': np.string_('height'),
                            'units': np.string_('m'),
                            'ChunkSize': chunkSize,
                        },
                    },
                    'LOS': {
                        'data': los, 
                        'attrs': {
                            'Shape': query_shape + (3,),
                            'ChunkSize': chunkSize + (3,),
                        },
                    },
                    'Rays_SP': {
                        'data': None, 
                        'attrs': {
                            'Shape': query_shape + (3,),
                            'ChunkSize': chunkSize + (3,),
                            'dtype': '<f8',
                        },
                    },
                    'Rays_SLV': {
                        'data': None, 
                        'attrs': {
                            'Shape': query_shape + (3,),
                            'ChunkSize': chunkSize + (3,),
                            'dtype': '<f8',
                        },
                    },
                    'Rays_len': {
                        'data': None,
                        'attrs': {
                            'Shape': query_shape + (3,),
                            'ChunkSize': chunkSize + (3,),
                            'dtype': '<f8',
                        },
                    },
                },
                attrs = {
                    'ChunkSize': chunkSize,
                    'NoDataValue': NDV,
                    'NumRays': np.prod(query_shape),
                },
                outName = pnts_file,
                chunkSize = chunkSize,
                NoDataValue = NDV
            )

    elif heights[0] == 'lvs':
        zlevels = hgts

    else:
        zlevels = None

    wetDelay, hydroDelay = computeDelay(
        weather_model_file, 
        pnts_file, 
        useWeatherNodes, 
        zlevels,
        zref, 
        out = out,
    )

    if heights[0] == 'lvs':
        outName = wetFilename[0].replace('wet', 'delays')
        writeDelays(flag, wetDelay, hydroDelay, lats, lons,
                    outName, zlevels=hgts, outformat=outformat, delayType=delayType)
        logger.info('Finished writing data to %s', outName)
    elif useWeatherNodes:
        logger.info(
            'Delays have been written to the weather model file; see %s',
            weather_model_file
        )
    else:
        writeDelays(flag, wetDelay, hydroDelay, lats, lons,
                    wetFilename, hydroFilename, outformat=outformat,
                    proj=None, gt=None, ndv=0.)
        logger.info('Finished writing data to %s', wetFilename)

    return wetDelay, hydroDelay


def weather_model_debug(los, lats, lons, ll_bounds, weather_model, wmLoc, zref,
                        time, out, download_only):
    """
    raiderWeatherModelDebug main function.
    """

    log.debug('Starting to run the weather model calculation with debugging plots')
    log.debug('Time type: %s', type(time))
    log.debug('Time: %s', time.strftime('%Y%m%d'))

    # location of the weather model files
    logger.debug('Beginning weather model pre-processing')
    logger.debug('Download-only is %s', download_only)
    if wmLoc is None:
        wmLoc = os.path.join(out, 'weather_files')

    # weather model calculation
    wm_filename = make_weather_model_filename(weather_model['name'], time, ll_bounds)
    weather_model_file = os.path.join(wmLoc, wm_filename)

    if not os.path.exists(weather_model_file):
        weather_model, lats, lons = prepareWeatherModel(
            weather_model, wmLoc, lats=lats, lons=lons, los=los, zref=zref,
            time=time, download_only=download_only, makePlots=True)
        try:
            weather_model.write2HDF5(weather_model_file)
        except Exception:
            logger.exception("Unable to save weathermodel to file")

        del weather_model
    else:
        logger.warning(
            'Weather model already exists, please remove it ("%s") if you want '
            'to create a new one.', weather_model_file
        )

    return 1

 
def getChunkSize(in_shape):
    ''' 
    Try to estimate a reasonable chunk size for writing data 

    Parameters
    ----------
    in_shape    - shape of the array to be chunked

    Returns
    -------
    Tuple of chunk sizes for each dimension in in_shape
    '''
    minChunkSize = 100
    maxChunkSize = 10000
    cpu_count = mp.cpu_count()
    chunkSize = tuple(
            max(
                min(maxChunkSize, s // cpu_count),
                min(s, minChunkSize)
            ) for s in in_shape
        )
    return chunkSize
