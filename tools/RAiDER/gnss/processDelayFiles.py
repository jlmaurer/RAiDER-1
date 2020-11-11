import datetime
import glob
import os

import pandas as pd


def combineDelayFiles(outName, loc = os.getcwd(), ext='.csv'):
    files = glob.glob(loc + '*' + ext)
    addDateTimeToFiles(files)
    df = concatDelayFiles(files)
    df.sort_values(by=['Datetime', 'ID'], inplace=True).to_csv(os.path.join(loc, outName), index=False)


def addDateTimeToFiles(fileList, force=False):
    ''' Run through a list of files and add the datetime of each file as a column '''
    for f in fileList:
        data = pd.read_csv(f)

        if 'Datetime' in data.columns and not force:
            print('Files already have been processed, pass "force = True" if you want to continue')
            return
        dt = getDateTime(f)
        data['Datetime'] = dt
        data.to_csv(f, index=False)

def getDateTime(filename):
    ''' Parse a datetime from a RAiDER delay filename '''
    parts = filename.split('_')
    dt = parts[2]
    return datetime.datetime.strptime(dt, '%Y%m%dT%H%M%S')


def concatDelayFiles(fileList):
    dfList = []
    for f in fileList:
        dfList.append(pd.read_csv(f))
    return pd.concat(dfList)



