#!/usr/bin/env python3
import os
import sys
import requests


UNR_URL = "http://geodesy.unr.edu/"

def download_UNR(statID, year, writeDir = '.', baseURL = UNR_URL):
    '''
    Download a zip file containing tropospheric delays for a given station and year
    The URL format is http://geodesy.unr.edu/gps_timeseries/trop/<ssss>/<ssss>.<yyyy>.trop.zip
    Inputs:
        statID   - 4-character station identifier
        year     - 4-numeral year
    '''
    URL = "{baseURL}gps_timeseries/trop/{statID}/{statID}.{year}.trop.zip".format(baseURL, statID.lower(), year)
    saveLoc = os.path.join(writeDir, '{statID}.{year}.trop.zip'.format(statID.lower(), year))
    download_url(URL, saveLoc)
    

def getStatsByllh(llhBox, baseURL = UNR_URL):
    '''
    Function to pull lat, lon, height, beginning date, end date, and number of solutions for stations inside the bounding box llhBox. 
    llhBox should be a tuple with format (lat1, lat2, lon1, lon2), where lat1, lon1 define the lower left-hand corner and lat2, lon2 
    define the upper right corner. 
    '''
    import pandas as pd
    from urllib.request import urlopen 

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


def download_url(url, save_path, chunk_size=128):
    '''
    Download a file from a URL. Modified from 
    https://stackoverflow.com/questions/9419162/download-returned-zip-file-from-url
    '''
    r = requests.get(url, stream=True)
    with open(save_path, 'wb') as fd:
        for chunk in r.iter_content(chunk_size=chunk_size):
            fd.write(chunk)


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


def getStationList(bbox, stationFileName = 'gnssStationList.txt'):
    '''
    Creates a list of stations inside a lat/lon bounding box from a source
    Inputs: 
        bbox    - length-4 list that describes a bounding box. Format is 
                  S N W E
    '''
    bbox = [float(p) for p in bbox]
    statList = getStatsByllh(bbox)
    stations = list(statList.keys())
    writeStationList(stations, stationFileName)
    return stations


def inBox(self, llhbox):
    '''
    Checks whether the given lat, lon pair are inside the bounding box llhbox
    '''
    lattest = self._Lat < llhbox[1] and self._Lat > llhbox[0]    
    lontest = self._Lon < llhbox[3] and self._Lon > llhbox[2]
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


def downloadTropoDelays(stats, years, writeDir = '.'):
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

    for stat in stats:
        for year in years:
            download_UNR(stat, year, writeDir = writeDir)
            

    
if __name__=='__main__':

    if len(sys.argv)==6:
       year = int(sys.argv[1])
       stats = getStationList([float(d) for d in sys.argv[2:]])
    elif len(sys.argv)==3:
       year = int(sys.argv[1])
       stats = readStatNames(sys.argv[2])
    else:
       print('USAGE: ')
       print('downloadGNSSDelays.py <year> [<bounding box> OR <stationFileName>]')
       print(' ')
       print('Example: ')
       print('downloadGNSSDelays.py 2016 15.5 21.5 258 266')
       print('downloadGNSSDelays.py 2016 gnssStationList.txt')
       sys.exit(0)

    downloadTropoDelays(stats, year)


