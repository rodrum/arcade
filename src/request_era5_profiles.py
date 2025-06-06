#!/usr/bin/env python
import cdsapi
from os import mkdir
from os.path import join, isfile
import toml
import datetime 
import subprocess
import numpy as np

def format_csv(in_file, out_file):
    """Function to format CSV"""
    with open(in_file, 'r') as f_in:
        with open(out_file, 'w') as f_out:
            f_out.write("Latitude,Longitude,Value\n")
            for line in f_in.readlines():
                if line[0:8] != "Latitude":
                    f_out.write(','.join(line.split()) + '\n')
    print("Formatted {0} --> {1}".format(in_file, out_file))

def check_downloads(output_name):
    """Check if profiles have already been downloaded"""
    is_downloaded = []
    temp_name = f"{output_name}_temperature.grib"
    is_downloaded.append(isfile(temp_name))
    zonal_name = f"{output_name}_zonalWinds.grib"
    is_downloaded.append(isfile(zonal_name))
    merid_name = f"{output_name}_meridionalWinds.grib"
    is_downloaded.append(isfile(merid_name))
    return is_downloaded

def main():
    c = cdsapi.Client()

    output_dir = "../output/ecmwf"
    try:
        mkdir(output_dir)
        print(f"Folder {output_dir} created")
    except FileExistsError:
        print(f"Folder {output_dir} already there")

    config = toml.load("../input/config.toml")
    year = config['discretization']['year']
    doys = config['discretization']['doys']
    secs = config['discretization']['sec']
    for sec in secs:
        for doy in doys:
            date_time = datetime.datetime.strptime(f"{year:04d} {doy:03d}", "%Y %j")
            date_str = date_time.strftime('%Y-%m-%d')  # "2016-06-10"
            hours = int(sec/3600)
            mins = int((sec-hours*3600)/60)
            seconds = int(sec-hours*3600-mins*60)
            time_str = f"{hours:02d}:{mins:02d}:{seconds:02d}"  # "14:00:00"

            dlon = config['discretization']['ecmwf']['dlon']  # grid step in degrees
            dlat = config['discretization']['ecmwf']['dlat']
            min_lon, max_lon, min_lat, max_lat = 0, 0, 0, 0
            if config['discretization']['ecmwf']['auto_area'] is False:
                min_lon = config['discretization']['ecmwf']['min_lon']
                max_lon = config['discretization']['ecmwf']['max_lon']
                min_lat = config['discretization']['ecmwf']['min_lat']
                max_lat = config['discretization']['ecmwf']['max_lat']
            else:
                sou_pos = config['discretization']['sources']['pos_latlon']
                sta_pos = config['discretization']['stations']['pos_latlon']
                min_lon = int(np.min([s[1] for s in (sou_pos+sta_pos)])) - 2*dlon
                max_lon = int(np.max([s[1] for s in (sou_pos+sta_pos)])) + 2*dlon
                min_lat = int(np.min([s[0] for s in (sou_pos+sta_pos)])) - 2*dlat
                max_lat = int(np.max([s[0] for s in (sou_pos+sta_pos)])) + 2*dlat

            output_name = join(
                output_dir,
                f"{date_str}_{time_str}"
                )

            print("")
            print("Requesting ECMWF ERA 5 data")
            print("===========================")
            print(f"Date        : {date_str}")
            print(f"Time        : {time_str}")
            print(f"NW corner   : ({max_lat}, {min_lon})")
            print(f"SE corner   : ({min_lat}, {max_lon})")
            print(f"dlon, dlat  : {dlon}, {dlat}")
            print(f"Output file : {output_name}")
            print("")

            # check if profiles have already been downloaded
            print("-> Checking for downloaded profiles...")
            temp_name = f"{output_name}_temperature.grib"
            zonal_name = f"{output_name}_zonalWinds.grib"
            merid_name = f"{output_name}_meridionalWinds.grib"
            is_downloaded = check_downloads(output_name)
            if is_downloaded[0] is True:
                print(f"--> Skipping download of {output_name}_temperature.grib")
            else:
                # temperature [K]
                print("-> Requesting temperatures...")
                c.retrieve('reanalysis-era5-complete', {
                    'class': 'ea',
                    'date': date_str,
                    'expver': '0001',
                    'levelist': '1/to/137',
                    'levtype': 'ml',
                    'param': '130',
                    'stream': 'oper',
                    'time': time_str,
                    'type': 'an',
                    'area': f"{max_lat}/{min_lon}/{min_lat}/{max_lon}", # N/W/S/E
                    'grid': f"{dlat}/{dlon}",
                    'format': 'grib'
                }, temp_name)

            if is_downloaded[1] is True:
                print(f"--> Skipping download of {output_name}_zonalWinds.grib")
            else:
                # u-wind (zonal) [m/s]
                print("-> Requesting zonal winds...")
                c.retrieve('reanalysis-era5-complete', {
                    'class': 'ea',
                    'date': date_str,
                    'expver': '0001',
                    'levelist': '1/to/137',
                    'levtype': 'ml',
                    'param': '131',
                    'stream': 'oper',
                    'time': time_str,
                    'type': 'an',
                    'area': f"{max_lat}/{min_lon}/{min_lat}/{max_lon}",
                    'grid': f"{dlat}/{dlon}",
                    'format': 'grib'
                }, zonal_name)

            if is_downloaded[2] is True:
                print(f"--> Skipping download of {output_name}_meridionalWinds.grib")
            else:
                # v-wind (meridonial) [m/s]
                print("-> Requesting meridional winds...")
                c.retrieve('reanalysis-era5-complete', {
                    'class': 'ea',
                    'date': date_str,
                    'expver': '0001',
                    'levelist': '1/to/137',
                    'levtype': 'ml',
                    'param': '132',
                    'stream': 'oper',
                    'time': time_str,
                    'type': 'an',
                    'area': f"{max_lat}/{min_lon}/{min_lat}/{max_lon}",
                    'grid': f"{dlat}/{dlon}",
                    'format': 'grib'
                }, merid_name)

            # do this only if end output are not already present
            new_files = isfile(f"{output_name}_temperature_new.csv") \
                        and \
                        isfile(f"{output_name}_zonalWinds_new.csv") \
                        and \
                        isfile(f"{output_name}_meridionalWinds_new.csv")
            if not new_files:
                #=== convert from grib to pre-csv
                f = open(f"{output_name}_temperature.csv", "w")
                subprocess.call(['grib_get_data', temp_name], stdout=f)
                f.close()

                f = open(f"{output_name}_zonalWinds.csv", "w")
                subprocess.call(['grib_get_data', zonal_name], stdout=f)
                f.close()

                f = open(f"{output_name}_meridionalWinds.csv", "w")
                subprocess.call(['grib_get_data', merid_name], stdout=f)
                f.close()

                #=== convert to proper CSV
                format_csv(
                    f"{output_name}_temperature.csv",
                    f"{output_name}_temperature_new.csv"
                    )

                format_csv(
                    f"{output_name}_zonalWinds.csv",
                    f"{output_name}_zonalWinds_new.csv"
                    )

                format_csv(
                    f"{output_name}_meridionalWinds.csv",
                    f"{output_name}_meridionalWinds_new.csv"
                    )
            else:
                print("-> CSV files present, skipping...")

if __name__ == '__main__':
    main()
