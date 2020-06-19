#!/usr/bin/env python3
#
import copy
import datetime as dt
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from urllib.request import urlopen 


def checkInt(int2test):
    '''
    checks for an integer
    '''
    if not isinstance(int2test, int):
        try:
            return int(int2test)
        except:
            raise RuntimesError('not an integer')
    else: 
        return int2test
    

def doy2date(doy, year):
    '''
    Convert a year + day-of-year combo to date.
    doy - integer day-of-year
    year- year of interest
    '''
    doy = checkInt(doy)
    year = checkInt(year)

    return dt.datetime(year, 1, 1) + dt.timedelta(doy - 1)


def getDate(stationFile):
    '''
    extract the date from a station delay file
    '''
    import re

    # find the date info
    p = re.compile(r'\d{4}\/\d{3}')
    try:
        year, doy = [int(t) for t in p.search(stationFile).group().split('/')]
    except:
        return None

    date = doy2date(doy, year)
    return date, year, doy


def getDelays(stationFile, returnTime = None):
    '''
    Parses and returns a dictionary containing either (1) all
    the GPS delays, if returnTime is None, or (2) only the delay
    at the closest times to to returnTime. 
    Inputs: 
         stationFile - a .gz station delay file 
         returnTime  - a .gz station delay file 
    Outputs:
         a dict containing the times and delay information
         (delay in mm, delay uncertianty,  delay gradients)
    '''
    import gzip

    # get the date of the file
    time, yearFromFile, doyFromFile = getDate(stationFile)

    d, ngrad, egrad,timesList,Sig = [], [], [],[],[]
    flag = False
    with gzip.open(stationFile, 'rb') as f:
        for line in f.readlines():
            line = line.decode('utf-8')
            if flag:
                if 'SITE' in line:
                   continue
                try:
                    de, deSD, ng, ngSD,  eg, egSD = [float(t) for t in line.split()[2:]]
                except:
                    continue
                site = line.split()[0]
                year, doy, seconds = [int(n) for n in line.split()[1].split(':')]
                if doy != doyFromFile:
                   continue

                d.append(de)
                ngrad.append(ng)
                egrad.append(eg)
                timesList.append(seconds)
                Sig.append(deSD)

            if 'TROP/SOLUTION' in line:
                flag = True

    # check for missing times
    true_times = list(range(0,86400, 300))
    if len(timesList)!=len(true_times):
       missing = [True if t not in timesList else False for t in true_times]
       mask = np.array(missing)
       delay, sig, east_grad, north_grad = [np.full((288,), np.nan)]*4 
       delay[~mask] = d
       sig[~mask] = Sig
       east_grad[~mask] = egrad
       north_grad[~mask] = ngrad
       times = true_times.copy()
    else:
       delay     = np.array(d)
       times     = np.array(timesList)
       sig       = np.array(Sig)
       east_grad = np.array(egrad)
       north_grad= np.array(ngrad)

 
    if returnTime == None:
       return {'StatName': site, 'Date': time, 'ZTD': delay, 'north_grad': north_grad, 'east_grad': east_grad, 'Datetimes': times, 'sigZTD': sig}
    else:
       index = np.argmin(np.abs(np.array(timesList) - returnTime))
       return {'StatName': site, 'Date': times, 'ZTD': delay[index], 'north_grad': ngrad[index], 'east_grad': egrad[index], 'Datetimes': timesList[index], 'sigZTD': sig[index]}


def convertMM2Rad(v, lam = 0.056):
    return -(v/1000)*4*np.pi/lam

def cosd(x):
    """Return the cosine of x when x is in degrees."""
    return np.cos(np.radians(x))


def initializeNetcdf(filename, stationName, metadata = None, outformat = 'NETCDF4'):
    '''
    Open a NetCDF file to prepare for writing GNSS delays
    '''
    import netCDF4

    with netCDF4.Dataset(filename, mode='w',format=outformat) as ncfile:
       times_dim = ncfile.createDimension('times', None) # unlimited axis 
       seconds_dim = ncfile.createDimension('seconds', 288)
       ncfile.station_name = stationName
       if metadata is not None:
          for key, item in metadata.items():
              ncfile.setncattr(key,item)

       crs = ncfile.createVariable('spatial_ref', 'i4')
       crs.spatial_ref='GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'

       # create variables to write
       times = ncfile.createVariable('times', 'f8', ('times',))
       times.units = "hours since 0001-01-01 00:00:00.0"
       times.long_name = 'times_in_days'
       times.calendar = "gregorian"

       # 
       seconds = ncfile.createVariable('seconds', 'int', ('seconds',))
       seconds.units = "Seconds since midnight"
       seconds.long_name = 'seconds_since_midnight'
       seconds[:] = list(range(0,86400, 300))

       # note: unlimited dimension is leftmost
       ztd = ncfile.createVariable('ztd','f8',('times','seconds')) 
       ztd.units = 'mm' 
       ztd.long_name = 'zenith_delay' 

       sig_ztd = ncfile.createVariable('sig_ztd','f8',('times','seconds')) 
       sig_ztd.units = 'mm' 
       sig_ztd.long_name = 'zenith_delay_1sigma_uncertainty' 

       east_grad = ncfile.createVariable('east_grad','f8',('times','seconds')) 
       east_grad.units = '' 
       east_grad.long_name = 'zenith_east_gradient'

       north_grad = ncfile.createVariable('north_grad','f8',('times','seconds')) 
       north_grad.units = '' 
       north_grad.long_name = 'zenith_north_gradient'


def writeStationDataToNetcdf(filename, data):
    '''
    Open a NetCDF file to prepare for writing GNSS delays
    '''
    import netCDF4
    from netCDF4 import num2date, date2num
    import datetime

    with netCDF4.Dataset(filename, mode='a') as ncfile:
       times      = ncfile.variables["times"]
       ztd        = ncfile.variables["ztd"]
       sig_ztd    = ncfile.variables["sig_ztd"]
       north_grad = ncfile.variables["north_grad"]
       east_grad  = ncfile.variables["east_grad"]

       date = datetime.datetime.combine(data['Date'].date(), datetime.time())
       newDate = date2num(date, units=times.units, calendar=times.calendar)
       nt = len(times[:])
       times[nt] = newDate
       ztd[nt,:] = data['ZTD'][np.newaxis,...]
       sig_ztd[nt,:] = data['sigZTD'][np.newaxis,...]
       east_grad[nt,:] = data['east_grad'][np.newaxis,...]
       north_grad[nt,:] = data['north_grad'][np.newaxis,...]


def getStationData(statName, dirLoc = None, times = None,outDir = None, 
                   inc = None, metadata = None, outformat = 'NETCDF4'):
    '''
    Pull tropospheric delay data for a given station name
    inc is inclination in degrees. 
    '''

    if dirLoc is None:
        dirLoc = os.getcwd()
    if outDir is None:
        outDir = os.getcwd()

    stationFiles = glob.glob(dirLoc+ os.sep + '*' + 
                             os.sep + '*' + os.sep + 
                             statName + '*.gz')

    stationFiles.sort()

    print('Found {} delay files for station {}'.format(len(stationFiles), statName))
    if len(stationFiles)>0:
       name = os.path.join(outDir, statName + '_ztd.nc')
       initializeNetcdf(name, statName, metadata, outformat)
       statList = []
       for sf in stationFiles:
           result = getDelays(sf, returnTime = times)
           writeStationDataToNetcdf(name, result)


def extractDelayFromNetCDF(filename, timeindex, delayThresh=0.8): 
    '''
    Extract delay information from a ZTD delay netcdf file at a given timeindex
    '''
    import netCDF4
    from netCDF4 import date2index, num2date
    import numpy as np

    with netCDF4.Dataset(filename) as f: 
        times = f.variables['times']  
        delays = f.variables['ztd']  
        delaySig = f.variables['sig_ztd']  

        outDates, outDelay, outSig = [],[],[]
        for k,t in enumerate(times):
           dest = delays[k].copy()[timeindex]/1000
           if dest < delayThresh:
              dest = np.nan
           outDelay.append(dest)
           outSig.append(delaySig[k].copy()[timeindex]/1000)
           outDates.append(num2date(t, times.units))
    return outDates, outDelay, outSig


def getZTDAllStations(stationFiles, timeindex, delayThresh=0.8, baseDate=None, lastDate=None):
    '''
    Finds all netcdf files in dirLoc and pulls the date and time
    from each one.
    timeindex    - number of five-minute intervals since midnight
    '''
    import datetime
    import glob
    import re
    import pandas as pd

    # set up an initial dataframe to populate
    if baseDate is None:
        baseDate = datetime.datetime(2010, 1, 1)
    if lastDate is None:
        lastDate = datetime.datetime.today()
    masterDates = [baseDate + k*datetime.timedelta(days=1) for k in range((lastDate - baseDate).days)]
    df = pd.DataFrame(index=masterDates)
    dfSig = pd.DataFrame(index=masterDates)


    # Pull all the delays at the given time
    p = re.compile(r'[A-Z0-9]{4}')
    delays, sigs = [],[]
    for sf in stationFiles:
        ID = p.search(sf).group()
        outDates, outDelay, outSig = extractDelayFromNetCDF(sf, timeindex)
        df.loc[outDates, ID] = outDelay
        dfSig.loc[outDates, ID] = outSig

    # now tidy the result
    df['Datetime']    = df.index
    dfSig['Datetime'] = dfSig.index
    df    = df.reset_index().drop('index', axis=1)
    dfSig = dfSig.reset_index().drop('index', axis=1)
    outdf    = pd.melt(frame=df,    id_vars=['Datetime'], var_name = 'ID', value_name = 'ZTD_m')
    outdfSig = pd.melt(frame=dfSig, id_vars=['Datetime'], var_name = 'ID', value_name = 'ZTD_Sig_m')
    outdf = outdf.merge(outdfSig)

    return outdf


