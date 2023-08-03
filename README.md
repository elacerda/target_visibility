Targets visibility
==================

Check targets visibility at the sky at some location of the Earth for a time interval.

It uses three constraints: maximal airmass, minimal moon separation and minimal altitude.
The night is defined by the interval between evening and morning twilights (either civil,
nautical or astronomical). 

The scripts produces a CSV file containing the following information at each line:
```
    DATE,RA,DEC,ATNIGHT,AIRMASS,ALTITUDE,MOONSEP,TOTAL

    DATE: The night defined by DATE 00:00:00
    RA,DEC: The coordinates of the target using the input units.
    ATNIGHT: 0/1 - observability of the target at night
    AIRMASS: 0/1 - constraint of maximal airmass achieved
    ALTITUDE: 0/1 - constraint of minimal altitude achieved
    MOONSEP: 0/1 - constraint of the maximal moon separation achieved
    TOTAL: 0/1 - The AND between all constraint and ATNIGHT.
```

Usage
-----

**targets_visibility.py** usage:
```
    usage: targets_visibility.py [-h] [--end_date YYYY-MM-DD] [--output FILENAME] [--twilight {astronomical,nautical,civil}]
                                 [--min_moonsep MIN_MOONSEP] [--max_airmass MAX_AIRMASS] [--min_alt MIN_ALT] [--unit_ra {deg,hourangle}]
                                 [--tel_lat TEL_LAT] [--tel_lon TEL_LON] [--tel_hei TEL_HEI] [--tel_name TEL_NAME] [--tel_tz TEL_TZ]
                                 FILENAME YYYY-MM-DD
    Check targets visibility at the sky of some location of the Earth.

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

    positional arguments:
      FILENAME              CSV file containing, at least, two columns with ICRS RA and DEC coordinates. 
                            Other column will be ignored. The unit of the RA could be entered in degrees or hourangle 
                            (set by --unit_ra argument). The DEC coordinate should be in degrees. Each target must be an 
                            entry line: 
                        
                                RA,DEC[,NAME]
      YYYY-MM-DD            Initial date to the timeline to be averiguated. If --end_date is not setted, the program will
                            run for this night only.

    options:
      -h, --help            show this help message and exit
      --end_date YYYY-MM-DD
                            Final date of the timeline to be averiguated. Defaults to None.
      --output FILENAME, -O FILENAME
                            Outputs to file. Defaults to None.
      --twilight {astronomical,nautical,civil}
                            It defines the observable night: 
                                (evening twilight)<>(obs. night)<>(morning twilight).
                            Values are: 
                                "civil": sun crossing the -6 deg altitude horizon;
                                "nautical": sun crossing the -12 deg altitude horizon;
                                "astronomical": sun crossing the -18 deg altitude horizon.
                            Defaults to astronomical.
      --min_moonsep MIN_MOONSEP
                            Minimal separation angle to the moon, in degrees. Defaults to 40 deg.
      --max_airmass MAX_AIRMASS
                            Maximal airmass. Defaults to 2.
      --min_alt MIN_ALT     Minimal altitude. Defaults to 40 deg.
      --unit_ra {deg,hourangle}
                            Helps astropy.coordinates.SkyCoord to parse the RA from FILENAME. Defaults to deg.
      --tel_lat TEL_LAT     Telescope geodetic latitude in degrees. Could be in decimal or dms. Defaults to T80-South latitude (-30d10m04.31s).
      --tel_lon TEL_LON     Telescope geodetic longitude in degrees. Could be in decimal or dms. Defaults to T80-South longitude (-70d48m20.48s).
      --tel_hei TEL_HEI     Telescope geodetic height in meters. Defaults to T80-South height (2178 m).
      --tel_name TEL_NAME   Telescope name. Defaults to T80-South.
      --tel_tz TEL_TZ       Telescope timezone. Defaults to America/Santiago. 
                            See: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
                            
```

Example
-------
```bash
    $ python3 targets_visibility.py targets.csv 2023-01-01 --unit_ra hourangle
    # TELESCOPE: T80-South
    # TELESCOPE LOCATION (lon, lat, hei): (-70.8057 deg, -30.1679 deg, 2178.0000 m)
    # TELESCOPE TIMEZONE: America/Santiago
    # DATE: 2023-01-01
    DATE,RA,DEC,ATNIGHT,AIRMASS,ALTITUDE,MOONSEP,TOTAL
    2023-01-01,2.105175095121998,-0.3815939553815652,1,0,1,0,0
    2023-01-01,17.921878822876767,-70.95404058902848,0,0,0,1,0
    2023-01-01,6.201913813027653,-70.92900805509113,1,0,0,1,0

    $ python3 targets_visibility.py targets.csv 2023-01-01 --min_alt 80
    targets_visibility.py: Zero targets are observable during the requested timeline within the requested constraints.
```

Contact
-------
	
Contact us: [dhubax@gmail.com](mailto:dhubax@gmail.com).
