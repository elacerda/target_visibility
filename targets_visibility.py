import io
import sys
import warnings
import argparse as ap
import astropy.units as u
from os.path import basename
from astropy.io import ascii
from datetime import datetime
from astropy.time import Time
from contextlib import redirect_stdout
from astroplan.utils import time_grid_from_range
from astropy.coordinates import EarthLocation, Angle, SkyCoord
from astroplan import (Observer, FixedTarget, observability_table,
                       AltitudeConstraint, AirmassConstraint, AtNightConstraint, 
                       MoonSeparationConstraint, is_event_observable)

# Ignore warnings
warnings.filterwarnings('ignore')

__script_name__ = basename(sys.argv[0])
__script_desc__ = """Check targets visibility at the sky at some location of the Earth for a time interval.

It uses three constraints: maximal airmass, minimal moon separation and minimal altitude.
The night is defined by the interval between evening and morning twilights (either civil,
nautical or astronomical). 

The scripts produces a CSV file containing the following information at each line:
    DATE,RA,DEC,ATNIGHT,AIRMASS,ALTITUDE,MOONSEP,TOTAL

    DATE: The night defined by DATE 00:00:00
    RA,DEC: The coordinates of the target using the input units.
    ATNIGHT: 0/1 - observability of the target at night
    AIRMASS: 0/1 - constraint of maximal airmass achieved
    ALTITUDE: 0/1 - constraint of minimal altitude achieved
    MOONSEP: 0/1 - constraint of the maximal moon separation achieved
    TOTAL: 0/1 - The AND between all constraint and ATNIGHT.
"""

args_help = {
    'filename': """CSV file containing, at least, two columns with ICRS RA and DEC coordinates. 
Other column will be ignored. The unit of the RA could be entered in degrees or hourangle 
(set by --unit_ra argument). The DEC coordinate should be in degrees. Each target must be an 
entry line: 

    RA,DEC[,NAME]
""",
    'twilight': """It defines the observable night: 
    (evening twilight)<>(obs. night)<>(morning twilight).
Values are: 
    "civil": sun crossing the -6 deg altitude horizon;
    "nautical": sun crossing the -12 deg altitude horizon;
    "astronomical": sun crossing the -18 deg altitude horizon. 
Defaults to astronomical.""",
    'min_moonsep': 'Minimal separation angle to the moon, in degrees. Defaults to 40 deg.',
    'max_airmass': 'Maximal airmass. Defaults to 2',
    'min_alt': 'Minimal altitude. Defaults to 40 deg.',
    'unit_ra': 'Helps astropy.coordinates.SkyCoord to parse the RA from FILENAME. Defaults to deg.',
    'tel_lat': 'Telescope geodetic latitude in degrees. Could be in decimal or dms. Defaults to T80-South latitude (-30d10m04.31s).',
    'tel_lon': 'Telescope geodetic longitude in degrees. Could be in decimal or dms. Defaults to T80-South longitude (-70d48m20.48s).',
    'tel_hei': 'Telescope geodetic height in meters. Defaults to T80-South height (2178 m).',
    'tel_name': 'Telescope name. Defaults to T80-South.',
    'tel_tz': """Telescope timezone. Defaults to America/Santiago. 
See: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones""",
    'start_date': """Initial date to the timeline to be averiguated. If --end_date is not setted, the program will
run for this night only.""",
    'end_date': 'Final date of the timeline to be averiguated. Defaults to None.',
    'output': 'Outputs to file. Defaults to None.',
}

def parse_arguments():
    t_ch = ['astronomical', 'nautical', 'civil']

    parser = ap.ArgumentParser(
        prog=__script_name__, 
        description=__script_desc__, 
        formatter_class=ap.RawTextHelpFormatter
    )
    
    # POSITIONAL ARGUMENTS
    parser.add_argument('filename', metavar='FILENAME', help=args_help['filename'])
    parser.add_argument('start_date', metavar='YYYY-MM-DD', help=args_help['start_date'])

    # OPTIONAL ARGUMENTS
    parser.add_argument('--end_date', metavar='YYYY-MM-DD', default=None, help=args_help['end_date'])
    parser.add_argument('--output', '-O', metavar='FILENAME', default=None, help=args_help['output'])
    parser.add_argument('--twilight', choices=t_ch, default=t_ch[0], help=args_help['twilight'])
    parser.add_argument('--min_moonsep', default=40, type=float, help=args_help['min_moonsep'])
    parser.add_argument('--max_airmass', default=2, type=float, help=args_help['max_airmass'])
    parser.add_argument('--min_alt', default=40, type=float, help=args_help['min_alt'])
    parser.add_argument('--unit_ra', choices=['deg', 'hourangle'], default='deg', help=args_help['unit_ra'])
    parser.add_argument('--tel_lat', default='-30d10m04.31s', help=args_help['tel_lat'])
    parser.add_argument('--tel_lon', default='-70d48m20.48s', help=args_help['tel_lon'])
    parser.add_argument('--tel_hei', default=2178, type=float, help=args_help['tel_hei'])
    parser.add_argument('--tel_name', default='T80-South', help=args_help['tel_name'])
    parser.add_argument('--tel_tz', default='America/Santiago', help=args_help['tel_tz'])
    args = parser.parse_args(args=sys.argv[1:])
    _atnight = {
        'astronomical': AtNightConstraint.twilight_astronomical(),
        'nautical': AtNightConstraint.twilight_nautical(),
        'civil': AtNightConstraint.twilight_civil(),
    }
    args.constraints = {
        'airmass': AirmassConstraint(max=args.max_airmass),
        'atnight': _atnight[args.twilight],
        'altitud': AltitudeConstraint(min=args.min_alt*u.deg),
        'moonsep': MoonSeparationConstraint(min=args.min_moonsep*u.deg),
    }
    args.observer = Observer(
        location=EarthLocation.from_geodetic(
            lon=Angle(args.tel_lon, 'deg'), 
            lat=Angle(args.tel_lat, 'deg'), 
            height=args.tel_hei*u.m,
        ), 
        timezone=args.tel_tz, 
        name=args.tel_name
    )
    if args.output is None:
        args.output = sys.stdout

    args.targets = []
    # READ TARGETS
    T = ascii.read(args.filename)
    n_cols = len(T.columns)
    name_prefix = 'TARGET'
    targets = []
    for i, t in enumerate(T):
        ra = t[0]
        dec = t[1]
        units = (args.unit_ra, 'deg')
        t_coord = SkyCoord(ra=ra, dec=dec, unit=units)
        name = '{}_{:03d}'.format(name_prefix, i) if n_cols == 2 else t[2]
        targets.append(FixedTarget(coord=t_coord, name=name))
    args.targets = targets
    return args

def ieo(constraints, observer, target, jd_night, time_resolution=10):
    str_local_time = f'{jd_night.strftime("%Y-%m-%d")} 23:59:00'
    local_dt = datetime.fromisoformat(str_local_time)
    utc_Time = observer.datetime_to_astropy_time(local_dt)
    start_time = observer.twilight_evening_astronomical(utc_Time, which='nearest')
    end_time = observer.twilight_morning_astronomical(utc_Time, which='nearest')
    time_grid = time_grid_from_range(
        [start_time, end_time],
        time_resolution=time_resolution*u.min
    )
    oc = {}
    for cname, c in constraints.items():
        ieobs = is_event_observable(c, observer, target.coord, times=time_grid)
        oc[cname] = ieobs.any()
    oc['ieo'] = all(x for x in oc.values())
    return oc

if __name__ == '__main__':
    args = parse_arguments()

    start_date = Time(args.start_date)
    end_date = Time(args.end_date) if args.end_date is not None else start_date + 1*u.day
    
    ot = observability_table(
        constraints=[v for v in args.constraints.values()],
        observer=args.observer,
        targets=args.targets,
        time_range=[start_date, end_date]
    )
    n_possible_nights = sum([int(x['ever observable']) for x in ot])

    time_grid = time_grid_from_range([start_date, end_date], time_resolution=1*u.day)
    if n_possible_nights > 0:
        close_f = True
        if isinstance(args.output, io.TextIOWrapper):
            f = args.output
            close_f = False
        else:
            f = open(args.output, 'w')
        with redirect_stdout(f):
            print(f'# TELESCOPE: {args.observer.name}')
            formatted_loc = ', '.join(["{:.4f} {}".format(i.value, i.unit) for i in args.observer.location.to_geodetic()])
            print(f'# TELESCOPE LOCATION (lon, lat, hei): ({formatted_loc})')
            print(f'# TELESCOPE TIMEZONE: {args.observer.timezone}')
            if args.end_date is not None:
                print(f'# START DATE: {args.start_date}')
                print(f'# END DATE: {args.end_date}')
            else:
                print(f'# DATE: {args.start_date}')
            print('DATE,RA,DEC,ATNIGHT,AIRMASS,ALTITUDE,MOONSEP,TOTAL')
            for i_t, t in enumerate(args.targets):
                if ot[i_t]:
                    ra = t.coord.ra.to(args.unit_ra).value
                    dec = t.coord.dec.value
                    # CHECK TARGET NIGHT OBSERVABILITY
                    for i_day, _jd in enumerate(time_grid):
                        date_str = datetime.strftime(_jd.datetime, '%Y-%m-%d')
                        oc = ieo(args.constraints, args.observer, t, _jd)
                        avail_airmass = int(oc['airmass'])
                        avail_atnight = int(oc['atnight'])
                        avail_altitud = int(oc['altitud'])
                        avail_moonsep = int(oc['moonsep'])
                        avail = oc['ieo']
                        print(
                            '{},{},{},{},{},{},{},{}'.format(
                                date_str, ra, dec, avail_airmass, 
                                avail_atnight, avail_altitud, avail_moonsep, 
                                avail,
                            )
                        )
        if close_f:
            f.close()
    else:
        sys.exit(f'{__script_name__}: Zero targets are observable during the requested timeline within the requested constraints.')