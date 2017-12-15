import sys
import argparse


def get_gnrl_path(season):
    '''GET_GNRL_PATH'''
    '''season: season of treatment'''
    '''return the general path for data directory on RAIJIN'''

    custom_path = '/g/data2/rr5/vhl548/season_' + season + '/'
    custom_path_bis_1213 = '/g/data2/rr5/vhl548/season_1213_cfradial/'
    uf_path = '/g/data2/rr5/arm/data/cpol_UF/' + season + '/ppi/'
    lassen_path = '/g/data2/rr5/arm/data/cpol_lassen/cpol_' + season + '/PPI/'

    the_path = {"9899": {'path': lassen_path, 'subdir': False},
                "9900": {'path': lassen_path, 'subdir': False},
                "0102": {'path': lassen_path, 'subdir': False},
                "0203": {'path': lassen_path, 'subdir': False},
                "0304": {'path': uf_path, 'subdir': False},
                "0405": {'path': uf_path, 'subdir': False},
                "0506": {'path': uf_path, 'subdir': True},
                "0607": {'path': uf_path, 'subdir': True},
                "0910": {'path': custom_path, 'subdir': True},
                "1011": {'path': custom_path, 'subdir': True},
                "1112": {'path': custom_path, 'subdir': True},
                "1213": {'path': custom_path_bis_1213, 'subdir': True},
                "1314": {'path': custom_path, 'subdir': True},
                "1415": {'path': custom_path, 'subdir': True},
                "1516": {'path': custom_path, 'subdir': False}
                }

    return the_path[season]  # type: Dict[str, bool]


def get_date(season_str):
    '''
    GET_YEAR
    Does not work for season BEFORE 1998 or AFTER 2020 (just add dates)

    Parameter
    ==========
        season_str: str
            string of season, e.g. "1213" for 2012/2013

    Returns
    =======
        season_start: str
            Beginning of the season
        season_end: str
            End of the season
    '''

    the_years = {"9899": ['19981206', '19990507'],
                 "9900": ['19991104', '20000403'],
                 "0102": ['20011026', '20020407'],
                 "0203": ['20021029', '20030714'],
                 "0304": ['20031020', '20040504'],
                 "0405": ['20041201', '20050228'],
                 "0506": ['20051110', '20060430'],
                 "0607": ['20061012', '20070419'],
                 "0910": ['20091124', '20100430'],
                 "1011": ['20101105', '20110329'],
                 "1112": ['20111120', '20120420'],
                 "1213": ['20121109', '20130506'],
                 "1314": ['20131011', '20140504'],
                 "1415": ['20141203', '20150602'],
                 "1516": ['20151006', '20160215']
                 }

    return the_years[season_str]  # type: List[int, int]


def write_example():

    ex = """[general]
# Number of CPU for multiprocessing
# start and end date in YYYYMMDD format
ncpu = 8
start_date = 20121109
end_date = 20130506

[path]
ground_radar = /path/to/ground/radar/data/
satellite = /path/to/satellite/radar/data/
output = /path/to/ouptut/directory/

[radar]
# rmin is the minimum radar range (in meter) we start looking for data.
# must be at least 15000m. rmax is the radar maximum range (in m).
# offset is the offset, in dB, by which we want to correct the reflectivity.
# Units in meters and degrees
radar_name = MYSUPERRADAR
radar_id = MYSUPERRADAR_ID1
rmin = 15000
rmax = 150000
longitude = 12.8
latitude = 48.0
altitude = 42
beamwidth = 1.0
offset = 0
sat_offset = 0

[thresholds]
# Threshold on satellite reflectivity
# Minimum number of pair
# Minimum number of satellite profiles
# Maximum time diffenrece between radar and satellite, in seconds
# Threshold on ground radar reflectivity
min_sat_reflec = 17
min_pair = 10
min_profiles = 10
max_time_delta = 300
min_gr_reflec = 17

[switch]
# Case insenstive, can be yes/no, y/n, true/false, 1/0
# Using dBZ or natural units for the statistical calculations
# Satellite is GPM (false for TRMM)
# Ground radar is C-Band (false for S-Band)
# Writing results in output directory
# Correct ground radar attenuation using pyart
dbz = False
gpm = False
cband = True
write = True
correct_gr_attenuation = False
"""

    with open('matchvol_configuration.ini', 'w') as fid:
        fid.write(ex)

    print("Example configuration file written in: matchvol_configuration.ini")
    print("Once you have modified the configuration file, type: ")
    print("matchvol.py -c matchvol_configuration.ini")
    return None


def main():

    st_date, end_date = get_date(season)
    gr_path = get_gnrl_path(season)['path']

    pfix = '''[general]
# Number of CPU for multiprocessing
# start and end date in YYYYMMDD format
ncpu = %i
start_date = %s
end_date = %s

[path]
ground_radar = %s
satellite = /g/data2/rr5/vhl548/TRMM
output = /home/548/vhl548/data/

''' % (ncpu, st_date, end_date, gr_path)

    sfix = '''[radar]
# rmin is the minimum radar range we start looking for data
# Units in meters and degrees
radar_name = CPOL
radar_id = IDR59
rmin = 15000
rmax = 150000
longitude = 131.04530334
latitude = -12.24880028
altitude = 42
beamwidth = 1.0
offset = 0

[satellite]
# Satellite reflectivity offset.
sat_offset = 0

[thresholds]
# Threshold on satellite reflectivity
# Minimum number of pair
# Minimum number of satellite profiles
# Maximum time diffenrece between radar and satellite, in seconds
# Threshold on ground radar reflectivity
min_sat_reflec = 17
min_pair = 10
min_profiles = 10
max_time_delta = 300
min_gr_reflec = 17

[switch]
# Case insenstive, can be yes/no, y/n, true/false, 1/0
# Using dBZ or natural units for the statistical calculations
# Satellite is GPM (false for TRMM)
# Ground radar is C-Band (false for S-Band)
# Writing results in output directory
# Correct ground radar attenuation using pyart
dbz = False
gpm = False
cband = True
write = True
correct_gr_attenuation = True
intermediary = False
'''

    with open(output_fname, 'w') as fid:
        fid.write(pfix)
        fid.write(sfix)

    return None


if __name__ == '__main__':
    welcome_msg = "Automatic creation of the configuration file for MSGR for CPOL data on NCI RAIJIN (only)."

    parser = argparse.ArgumentParser(description=welcome_msg)
    parser.add_argument('-j', '--cpu', dest='ncpu', default=16, type=int, help='Number of process')
    parser.add_argument('-o', '--output', dest='output', default=None, type=str, help='Name of the output file.')
    parser.add_argument('-s', '--season', dest='season', default=None, type=str, help='Season to parse.')

    feature_parser = parser.add_mutually_exclusive_group(required=False)
    feature_parser.add_argument('--example', dest='example', action='store_true', help='Create example file.')
    feature_parser.add_argument('--no-example', dest='example', action='store_false', help="Don't create example file.")
    parser.set_defaults(example=False)

    args = parser.parse_args()
    ncpu = args.ncpu
    season = args.season
    example = args.example
    output_fname = args.output

    if output_fname is None:
        output_fname = "config.ini"

    if ".ini" not in output_fname:
        output_fname += ".ini"

    # No season provided and this is not an example.
    if season is None:
        if not example:
            print("You need to provide an argument, type `generate_config_matchvol -h` for help.")
            sys.exit()

    if example:
        write_example()
    else:
        main()

    # End of __main__
