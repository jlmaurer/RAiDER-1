#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  Author: Jeremy Maurer, Raymond Hogenson & David Bekaert
#  Copyright 2019, by the California Institute of Technology. ALL RIGHTS
#  RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import contextlib
import os
import sys
from datetime import datetime


def getWMFilename(weather_model_name, time, outLoc, verbose=False):
    '''
    Check whether the output weather model exists, and
    if not, download it.
    '''
    with contextlib.suppress(FileExistsError):
        os.mkdir('weather_files')
    f = os.path.join(
        outLoc,
        '{}_{}.nc'.format(
            weather_model_name,
            datetime.strftime(time, '%Y_%m_%d_T%H_%M_%S')
        )
    )

    if verbose:
        print('Storing weather model at: {}'.format(f))

    download_flag = True
    if os.path.exists(f):
        print('WARNING: Weather model already exists, skipping download')
        download_flag = False

    return download_flag, f


def prepareWeatherModel(weatherDict, wmFileLoc, out, lats=None, lons=None,
                        los=None, zref=None, time=None, verbose=False,
                        download_only=False, makePlots=False):
    '''
    Parse inputs to download and prepare a weather model grid for interpolation
    '''
    import numpy as np

    from RAiDER.models.allowed import checkIfImplemented
    from RAiDER.utilFcns import getTimeFromFile

    # Make weather
    weather_model, weather_files, weather_model_name = \
    weatherDict['type'], weatherDict['files'], weatherDict['name']
    checkIfImplemented(weather_model_name.upper().replace('-', ''))

    # check whether weather model files are supplied
    if weather_files is None:
        download_flag, f = getWMFilename(weather_model.Model(), time, wmFileLoc, verbose)
    else:
        download_flag = False
        time = getTimeFromFile(weather_files[0])

    # if no weather model files supplied, check the standard location
    if download_flag:
        try:
            weather_model.fetch(lats, lons, time, f)
        except Exception as e:
            print('ERROR: Unable to download weather data')
            print('Exception encountered: {}'.format(e))
            sys.exit(0)

        # exit on download if download_only requested
        if download_only:
            print('WARNING: download_only flag selected. I will only '
                  'download the weather model, '
                  ' without doing any further processing.')
            return None, None, None

    # Load the weather model data
    if weather_files is not None:
        weather_model.load(*weather_files, outLats=lats, outLons=lons, los=los, zref=zref)
        download_flag = False
    else:
        weather_model.load(f, outLats=lats, outLons=lons, los=los, zref=zref)

    # weather model name
    if verbose:
        print('Number of weather model nodes: {}'.format(np.prod(weather_model.getWetRefractivity().shape)))
        print('Shape of weather model: {}'.format(weather_model.getWetRefractivity().shape))
        print('Bounds of the weather model: {}/{}/{}/{} (SNWE)'
              .format(np.nanmin(weather_model._ys), np.nanmax(weather_model._ys),
                      np.nanmin(weather_model._xs), np.nanmax(weather_model._xs)))
#        print('Using weather nodes only? (true/false): {}'.format(uwn))
        print('Weather model: {}'.format(weather_model.Model()))
        print('Mean value of the wet refractivity: {}'
              .format(np.nanmean(weather_model.getWetRefractivity())))
        print('Mean value of the hydrostatic refractivity: {}'
              .format(np.nanmean(weather_model.getHydroRefractivity())))
        # If the verbose option is called, write out the weather model to a pickle file
        print('Saving weather model object to pickle file')
        import pickle
        pickleFilename = os.path.join(out, 'pickledWeatherModel.pik')
        with open(pickleFilename, 'wb') as f:
            pickle.dump(weather_model, f)
        print('Weather Model Name: {}'.format(weather_model.Model()))
        print(weather_model)

    if makePlots:
        p = weather_model.plot('wh', True)
        p = weather_model.plot('pqt', True)

    return weather_model, lats, lons
