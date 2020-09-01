import argparse
from textwrap import dedent

from RAiDER.checkArgs import checkArgs
from RAiDER.cli.parser import add_bbox, add_out, add_verbose
from RAiDER.cli.validators import DateListAction, LOSFileTypeAction, date_type, time_type
from RAiDER.constants import _ZREF
from RAiDER.delay import tropo_delay, weather_model_debug
from RAiDER.logger import *
from RAiDER.models.allowed import ALLOWED_MODELS



def create_parser():
    """Parse command line arguments using argparse."""
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=dedent("""\
            Calculate tropospheric delay from a weather model.
            Usage examples:
            raiderDelay.py --date 20200103 --time 23:00:00 -b 39 40 -79 -78 --model ERA5 --zref 15000 -v
            raiderDelay.py --date 20200103 --time 23:00:00 -b 39 40 -79 -78 --model ERA5 --zref 15000 --heightlvs 0 100 200 -v
            raiderDelay.py --date 20200103 --time 23:00:00 --latlon test/scenario_1/geom/lat.dat test/scenario_1/geom/lon.dat --model ERA5 --zref 20000 -v --out test/scenario_1/
            """)
    )

    datetime = p.add_argument_group('Datetime')
    datetime.add_argument(
        '--date', dest='dateList',
        help=dedent("""\
            Date to calculate delay.
            Can be a single date, a list of two dates (earlier, later) with 1-day interval, or a list of two dates and interval in days (earlier, later, interval).
            Example accepted formats:
               YYYYMMDD or
               YYYYMMDD YYYYMMDD
               YYYYMMDD YYYYMMDD N
            """),
        nargs="+",
        action=DateListAction,
        type=date_type,
        required=True
    )

    datetime.add_argument(
        '--time', dest='time',
        help=dedent('''\
        Calculate delay at this time.
        Example formats:
           THHMMSS,
           HHMMSS, or
           HH:MM:SS'''),
        type=time_type, required=True)

    # Area
    area = p.add_argument_group('Area of Interest (Supply one)').add_mutually_exclusive_group(required=True)
    area.add_argument(
        '--latlon',
        '-ll',
        nargs=2,
        dest='query_area',
        help='GDAL-readable latitude and longitude raster files (2 single-band files)',
        metavar=('LAT', 'LONG')
    )
    add_bbox(area)
    area.add_argument(
        '--station_file',
        default=None,
        type=str,
        dest='query_area',
        help=('CSV file with a list of stations, containing at least '
              'the columns "Lat" and "Lon"')
    )

    # Line of sight
    los = p.add_argument_group(
        'Specify a Line-of-sight or state vector file. If neither is supplied, the Zenith delay will be returned'
    )
    los.add_argument(
        '--LOS_file_option', '-l',
        help=dedent('''\
        Can be:
            A GDAL-readable two-band file (B1: inclination, B2: heading), 
            An ESA orbit file (time must correspond to the input query time), or
            A 7-column text file containing state vectors
            '''),
        metavar='LOS',
        action=LOSFileTypeAction,
        dest='lineofsight',
        default=None,
    )
    los.add_argument(
        '--zref', '-z',
        help='Reference vertical integration height limit (meters) (default: {} m)'.format(_ZREF),
        type=float,
        default=_ZREF
    )

    # heights
    heights = p.add_argument_group('Height data. Default is ground surface for specified lat/lons, height levels otherwise')
    heights.add_argument(
        '--dem', '-d', default=None,
        help="""Specify a DEM to use with lat/lon inputs."""
    )
    heights.add_argument(
        '--heightlvs',
        help=("""A space-deliminited list of heights"""),
        default=None, nargs='+', type=float
    )

    # Weather model
    weather = p.add_argument_group("Weather model. See documentation for details")
    weather.add_argument(
        '--model',
        help="Weather model option to use.",
        type=lambda s: s.upper().replace("-", ""),
        choices=ALLOWED_MODELS,
        default='ERA5T')
    weather.add_argument(
        '--files',
        help="""OUT/PLEV or HDF5 file(s) """,
        default=None, nargs='+', type=str, metavar="FILES")
    weather.add_argument(
        '--weatherFiles', '-w',
        help='Directory location of/to write weather model files',
        default=None, dest='wmLoc')

    misc = p.add_argument_group("Run parameters")
    misc.add_argument(
        '--zref', '-z',
        help='Height limit when integrating (meters) (default: {} m)'.format(_ZREF),
        type=float,
        default=_ZREF)
    misc.add_argument(
        '--outformat',
        help='GDAL-compatible file format if surface delays are requested.',
        default=None)

    add_out(misc)

    misc.add_argument(
        '--download_only',
        help='Download weather model only without processing? Default False',
        action='store_true', dest='download_only', default=False)

    add_verbose(misc)

    return p


def parseCMD():
    """
    Parse command-line arguments and pass to tropo_delay
    We'll parse arguments and call delay.py.
    """

    p = create_parser()
    args = p.parse_args()

    # Argument checking
    los, lats, lons, ll_bounds, heights, flag, weather_model, wmLoc, zref, outformat, \
        times, out, download_only, verbose, \
        wetNames, hydroNames = checkArgs(args, p)

    if verbose:
        logger.setLevel(logging.DEBUG)

    # Loop over each datetime and compute the delay
    for t, wfn, hfn in zip(times, wetNames, hydroNames):
        try:
            (_, _) = tropo_delay(los, lats, lons, ll_bounds, heights, flag, weather_model, wmLoc, zref,
                                 outformat, t, out, download_only, wfn, hfn)

        except RuntimeError:
            logger.exception("Date %s failed", t)
            continue


def parseCMD_weather_model_debug():
    """
    Parse command-line arguments and pass to prepareWeatherModel
    We'll parse arguments and call delay.py.
    """

    p = create_parser()
    args = p.parse_args()

    # Argument checking
    los, lats, lons, ll_bounds, heights, flag, weather_model, wmLoc, zref, outformat, \
        times, out, download_only, verbose, \
        wetNames, hydroNames = checkArgs(args, p)

    if verbose:
        logger.setLevel(logging.DEBUG)

    # Loop over each datetime
    for t in times:
        try:
            weather_model_debug(los, lats, lons, ll_bounds, weather_model, wmLoc, zref, t, out, download_only)

        except RuntimeError:
            logger.exception("Date %s failed", t)
            continue
