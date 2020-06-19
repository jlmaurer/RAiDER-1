#!/usr/bin/env python3
import os
import sys
import requests
import argparse


UNR_URL = "http://geodesy.unr.edu/"


def download_UNR(statID, year, check = True, writeDir = '.', baseURL = UNR_URL, verbose = False):
    '''
    Download a zip file containing tropospheric delays for a given station and year
    The URL format is http://geodesy.unr.edu/gps_timeseries/trop/<ssss>/<ssss>.<yyyy>.trop.zip
    Inputs:
        statID   - 4-character station identifier
        year     - 4-numeral year
    '''
    URL = "{baseURL}gps_timeseries/trop/{statID}/{statID}.{year}.trop.zip".format(baseURL, statID.lower(), year)
    if check:
        flag = check_url(URL, verbose=verbose)
    else:
        saveLoc = os.path.join(writeDir, '{statID}.{year}.trop.zip'.format(statID.lower(), year))
        flag = download_url(URL, saveLoc, verbose=verbose)
    return flag


def getStatsByllh(llhBox = None, baseURL = UNR_URL):
    '''
    Function to pull lat, lon, height, beginning date, end date, and number of solutions for stations inside the bounding box llhBox. 
    llhBox should be a tuple with format (lat1, lat2, lon1, lon2), where lat1, lon1 define the lower left-hand corner and lat2, lon2 
    define the upper right corner. 
    '''
    import pandas as pd
    from urllib.request import urlopen 

    if llhBox is None:
        llhBox = [-90, 90, 0, 360]

    stationHoldings = '{baseURL}NGLStationPages/DataHoldings.txt'.format(baseURL)
    data = urlopen(stationHoldings) # it's a file like object and works just like a file
    stations = []
    for ind, line in enumerate(data): # files are iterable
        if ind == 0:
            continue
        statID, lat, lon = getID(line)
        if inBox(lat, lon, llhBox):
            stations.append({'ID': statID, 'Lat': lat, 'Lon': lon})
    
    print('{} stations were found'.format(len(stations)))
        
    return(pd.DataFrame(stations))


def check_url(url, verbose = False):
    '''
    Check whether a file exists at a URL. Modified from 
    https://stackoverflow.com/questions/9419162/download-returned-zip-file-from-url
    '''
    r = requests.get(url, stream=True)
    if r.status_code==404:
        return False
    else:
        return True


def download_url(url, save_path, verbose = False, chunk_size=128):
    '''
    Download a file from a URL. Modified from 
    https://stackoverflow.com/questions/9419162/download-returned-zip-file-from-url
    '''
    r = requests.get(url, stream=True)
    if r.status_code==404:
        return False
    else:
        if verbose:
            print('Beginning download of {} to {}'.format(url, save_path))
        with open(save_path, 'wb') as fd:
            for chunk in r.iter_content(chunk_size=chunk_size):
                fd.write(chunk)
        if verbose:
            print('Completed download of {} to {}'.format(url, save_path))
        return True


def readStationFile(filename):
    '''
    Read a list of GNSS station names from a plain text file
    '''
    statList = []
    with open(filename, 'r') as f:
        for line in f:
            statList.append(line.strip())
    return statList 


def writeStationList(statList, filename = 'gnssStationList.txt'):
    '''
    Write a python list of GNSS station names to a simple text file
    '''
    with open(filename, 'w') as f:
        for stat in statList:
            f.write('{}\n'.format(stat))


def getStationList(bbox = None, writeLoc = None):
    '''
    Creates a list of stations inside a lat/lon bounding box from a source
    Inputs: 
        bbox    - length-4 list of floats that describes a bounding box. Format is 
                  S N W E
    '''
    if bbox is not None:
        statList = getStatsByllh(bbox)
    else:
        statList = getStatsByllh()
    stations = list(statList.keys())

    if writeLoc is not None:
        writeStationList(stations, writeLoc)

    return stations


def inBox(lat, lon, llhbox):
    '''
    Checks whether the given lat, lon pair are inside the bounding box llhbox
    '''
    lattest = lat < llhbox[1] and lat > llhbox[0]    
    lontest = lon < llhbox[3] and lon > llhbox[2]
    if lattest and lontest:
        return True
    else:
        return False


def getID(line):
    '''
    Pulls the station ID, lat, and lon for a given entry in the UNR text file
    '''
    statID = str(line.split()[0].decode('utf-8'))
    lat = float(line.split()[1].decode('utf-8'))
    lon = float(line.split()[2].decode('utf-8'))
    return statID, lat, lon


def downloadTropoDelays(stats, years, check = True, writeDir = '.', verbose = False):
    '''
    Check for and download GNSS tropospheric delays from an archive. If check is True then 
    I will only return which stations have data in which years. 
    '''

    # argument checking
    if not isinstance(stats, list):
        if not isinstance(stats, str):
            raise TypeError('stats should be a string or a list of strings')
        stats = [stats]
    if not isinstance(years, list):
        if not isinstance(years, int):
            try:
                years = list(int(years))
            except TypeError:
                raise TypeError('years should be an int or a list of ints')
        years = [years]

    # Iterate over stations and years and check or download data
    results = []
    stat_year_tup = itertools.product(stats, years)
    for s, y in stat_year_tup:
        flag = download_UNR(s, y, check = check, writeDir = writeDir, verbose = verbose)
        results.append({'ID': s, 'year': y, 'flag': flag})
    return pd.DataFrame(results).set_index('ID')
            

def parse_years(input):
    '''
    Takes string input and returns a list of years as integers
    '''
    test = len(input.split(',')) > 1
    if test:
        years = [int(y) for y in input.split(',')]
        if years[1] > years[0] + 1:
            years = [years[0] + k for k in range(years[1] - years[0])]
    else:
        years = [int(input)]
    return years
    

def parse_args():
    """Parse command line arguments using argparse."""
    p = argparse.ArgumentParser(
          formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Check for and download tropospheric zenith delays for a set of GNSS stations from UNR
Usage examples: 
downloadGNSSdelay.py -y 2010 --check
downloadGNSSdelay.py -y 2010 -b 40 -79 39 -78 -v
downloadGNSSdelay.py -y 2010 -f station_list.txt --out products
""")
 
    p.add_argument(
        '--years', '-y', dest='years', nargs='+',
        help="""Year to check or download delays (format YYYY).
Can be a single value or a comma-separated list. If two years non-consequtive years are given, I will download each year in between as well. 
""", type=parse_years, required=True)

    # Stations to check/download
    area = p.add_argument_group('Stations to check/download. Can be a lat/lon bounding box or file, or will run the whole world if not specified')
    area.add_argument(
        '--station_file', '-f', default = None, dest='station_file',
        help=('Text file containing a list of 4-char station IDs separated by newlines'))
    area.add_argument(
        '--BBOX', '-b', nargs=4, dest='bounding_box', default=None,
        help="""Lat/Lon Bounding box in NWSE format""",
        metavar=('S', 'N', 'W', 'E'))

    misc = p.add_argument_group("Run parameters")
    misc.add_argument(
        '--outformat', 
        help='GDAL-compatible file format if surface delays are requested.',
        default=None)

    misc.add_argument(
        '--out', dest='out',
        help='Directory to download products', 
        default='.')

    misc.add_argument(
        '--check',
        help='Only check if data exists, do not download',
        action='store_true',dest='check', default = False)

    misc.add_argument(
        '--verbose', '-v',
        help='Run in verbose (debug) mode? Default False',
        action='store_true',dest='verbose', default = False)

    return p.parse_args(), p


def parseCMD():
    """
    Parse command-line arguments and pass to tropo_delay
    We'll parse arguments and call delay.py.
    """
    args, p = parse_args()

    # Handle different station requests 
    if args.station_file:
        stats = readStatNames(args.station_file)
    elif args.bbox:
        bbox = [float(d) for d in args.bbox]
        bbox[2] +=360
        bbox[3] +=360
        stats = getStationList()
    else:
        stats = getStationList()

    # iterate over years
    for yr in args.years:
        statDF = downloadTropoDelays(stats, yr, check = args.check, verbose = args.verbose)
        statDF.to_csv(os.path.join(args.out, 'stationList.csv'))
    if verbose:
        print('Completed processing')

