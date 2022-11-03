"""
Based on 'interpolate_ecmwf_notes.py' from Puyehue&Calbuco

2020-02-09
RSDNL @ West Campus
"""

import pandas
from os.path import join
import numpy as np
from obspy.geodetics.base import gps2dist_azimuth, kilometer2degrees
from scipy import interpolate 
import toml
import datetime

def main_rng_ind():
    #=== load levels 1/to/137
    levels_file = join("../input", 'ECMWF - L137 model level definitions.csv')
    levels = pandas.read_csv(levels_file, header=1)
    levels.info()

    #=== get altitudes
    # note that the levels are from top to bottom, and it's in [m]
    geom_alt = levels['Geometric Altitude [m]'].to_numpy() 
    geom_alt_flip = np.flip(geom_alt)/1000.  # to km

    #=== get temperature and winds
    #>>> from 'request_era5_profiles.py'
    output_dir = "../output/ecmwf"
    disc_param = toml.load("../input/discretize_parameters.toml")
    year = disc_param['year']
    #doy =  disc_param['doys'][0] # NOTE: only one day implemented for now
    doys =  disc_param['doys']
    sec =  disc_param['sec']
    for doy in doys:
        date_time = datetime.datetime.strptime(f"{year:04d} {doy:03d}", "%Y %j")
        date_str = date_time.strftime('%Y-%m-%d')  # "2016-06-10"
        hours = int(sec/3600)
        mins = int((sec-hours*3600)/60)*0      # NOTE: multiplied by zero because it seems the analysis are each hour 
        secs = int(sec-hours*3600-mins*60)*0
        time_str = f"{hours:02d}:{mins:02d}:{secs:02d}"  # "14:00:00"
        base_dir = join(
            output_dir,
            f"{date_str}_{time_str}"
            )

        lvls_temp = pandas.read_csv(base_dir+"_temperature_new.csv")
        lvls_zWin = pandas.read_csv(base_dir+"_zonalWinds_new.csv")
        lvls_mWin = pandas.read_csv(base_dir+"_meridionalWinds_new.csv")

        #=== discretization points from 'request_era5_profiles.py'
        dlon = disc_param['ecmwf']['dlon']  # grid step in degrees
        dlat = disc_param['ecmwf']['dlat']
        #lon0 = int(np.min([l[1] for l in (disc_param['sou_pos']+disc_param['sta_pos'])])) - 2*dlon
        #lon1 = int(np.max([l[1] for l in (disc_param['sou_pos']+disc_param['sta_pos'])])) + 2*dlon
        #lat0 = int(np.min([l[0] for l in (disc_param['sou_pos']+disc_param['sta_pos'])])) - 2*dlat
        #lat1 = int(np.max([l[0] for l in (disc_param['sou_pos']+disc_param['sta_pos'])])) + 2*dlat
        lon0, lon1, lat0, lat1 = 0, 0, 0, 0
        if disc_param['ecmwf']['auto_area'] == False:
            lon0 = disc_param['ecmwf']['min_lon']
            lon1 = disc_param['ecmwf']['max_lon']
            lat0 = disc_param['ecmwf']['min_lat']
            lat1 = disc_param['ecmwf']['max_lat']
        else:
            lon0 = int(np.min([l[1] for l in (disc_param['sou_pos']+disc_param['sta_pos'])])) - 2*dlon
            lon1 = int(np.max([l[1] for l in (disc_param['sou_pos']+disc_param['sta_pos'])])) + 2*dlon
            lat0 = int(np.min([l[0] for l in (disc_param['sou_pos']+disc_param['sta_pos'])])) - 2*dlat
            lat1 = int(np.max([l[0] for l in (disc_param['sou_pos']+disc_param['sta_pos'])])) + 2*dlat

        nlon, nlat = int(((lon1-lon0)/dlon+1)), int(((lat1-lat0)/dlat+1))
        npts = nlon * nlat

        #=== source, receiver, discretization along... 
        stations = disc_param['sta_pos']
        sources = disc_param['sou_pos']

        #=== calculate profile
        dl = kilometer2degrees(disc_param['ds']) # deg, discretization along prof.
        lons = np.arange(lon0-360., lon1-360.+dlon, dlon)
        lats = np.arange(lat1, lat0-dlat, -dlat)  # decreasing in latitude following data

        for nsou, (soulat, soulon) in enumerate(sources):
            print("\nSource #: {0} @ ({1:4.7f}, {2:4.7f})".format(
                nsou, soulat, soulon
                ))
            for nsta, (stalat, stalon) in enumerate(stations):
                print("\n -> Station #: {0} @ ({1:4.7f}, {2:4.7f})".format(
                    nsta, stalat, stalon
                    ))
                dist_m, az12_deg, _ = gps2dist_azimuth(soulat, soulon, stalat, stalon)
                dist_deg = kilometer2degrees(dist_m/1000.) # deg
                az12_rad = np.deg2rad(az12_deg)
                dist_vec = np.linspace(0, dist_deg, int(dist_deg/dl)+1)
                lon_new = soulon + dist_vec*np.sin(az12_rad)
                lat_new = soulat + dist_vec*np.cos(az12_rad)

                # here is where the average values will be saved
                av_temp = []
                av_zonal = []
                av_mer = []

                # Load each layer by height 
                start_ind = 0
                end_ind = npts
                for lvl in range(geom_alt.shape[0]):
                    temp_val = lvls_temp['Value'][start_ind:end_ind].to_numpy()
                    zWin_val = lvls_zWin['Value'][start_ind:end_ind].to_numpy()
                    mWin_val = lvls_mWin['Value'][start_ind:end_ind].to_numpy()

                    # Temperature 
                    z_temp = np.reshape(temp_val, (nlat, nlon))
                    f = interpolate.interp2d(lons, lats, z_temp)
                    znew_temp = f(lon_new, lat_new)

                    # Zonal winds 
                    z_zonal = np.reshape(zWin_val, (nlat, nlon))
                    f = interpolate.interp2d(lons, lats, z_zonal)
                    znew_zonal = f(lon_new, lat_new)

                    # Meridonial winds 
                    z_mer = np.reshape(mWin_val, (nlat, nlon))
                    f = interpolate.interp2d(lons, lats, z_mer)
                    znew_mer = f(lon_new, lat_new)

                    # Average on level 
                    # note: diagonal is along arc (xi, yi) i = 0,...,n
                    znew_temp_av = np.mean(np.diagonal(znew_temp))
                    znew_zonal_av = np.mean(np.diagonal(znew_zonal))
                    znew_mer_av = np.mean(np.diagonal(znew_mer))
                    av_temp.append(znew_temp_av)
                    av_zonal.append(znew_zonal_av)
                    av_mer.append(znew_mer_av)

                    start_ind = (lvl+1)*npts
                    end_ind = start_ind+npts


                # flip levels (low to high)
                av_temp = np.flip(av_temp)
                av_zonal = np.flip(av_zonal)
                av_mer = np.flip(av_mer)

                # fill data table
                table_name = f"ecmwf_{doy:03d}_{nsou+1:05d}_{nsta+1:04d}.csv"
                headerstr = ['Height [km]', 'Temperature [K]', 'Zonal wind [m/s]', 'Meridonial wind [m/s]']
                data_mat = np.zeros((geom_alt.shape[0], 4))
                data_mat[:, 0] = geom_alt_flip
                data_mat[:, 1] = av_temp
                data_mat[:, 2] = av_zonal
                data_mat[:, 3] = av_mer
                df = pandas.DataFrame(data=data_mat, columns=headerstr)
                df.to_csv(join("../output/ecmwf", table_name), index=False)

                # save in strict formatting too for reading with fortran
                table_name_txt = f"ecmwf_{doy:03d}_{nsou+1:05d}_{nsta+1:04d}.dat"
                with open(join("../output/ecmwf", table_name_txt), 'w') as f:
                    for lvl in range(geom_alt.shape[0]):
                        f.write("{0:8.3f} {1:8.3f} {2:8.3f} {3:8.3f}\n".format(
                            geom_alt_flip[lvl],
                            av_temp[lvl],
                            av_zonal[lvl],
                            av_mer[lvl]
                            ))

def main_rng_dep():
    #=== load levels 1/to/137
    levels_file = join("../input", 'ECMWF - L137 model level definitions.csv')
    levels = pandas.read_csv(levels_file, header=1)
    levels.info()

    #=== get altitudes
    # note that the levels are from top to bottom, and it's in [m]
    geom_alt = levels['Geometric Altitude [m]'].to_numpy() 
    geom_alt_flip = np.flip(geom_alt)/1000.  # to km

    #=== get temperature and winds
    #>>> from 'request_era5_profiles.py'
    output_dir = "../output/ecmwf"
    disc_param = toml.load("../input/discretize_parameters.toml")
    year = disc_param['year']
    #doy =  disc_param['doys'][0] # NOTE: only one day implemented for now
    doys =  disc_param['doys']
    sec =  disc_param['sec']
    for doy in doys:
        date_time = datetime.datetime.strptime(f"{year:04d} {doy:03d}", "%Y %j")
        date_str = date_time.strftime('%Y-%m-%d')  # "2016-06-10"
        hours = int(sec/3600)
        mins = int((sec-hours*3600)/60)*0      # NOTE: multiplied by zero because it seems the analysis are each hour 
        secs = int(sec-hours*3600-mins*60)*0
        time_str = f"{hours:02d}:{mins:02d}:{secs:02d}"  # "14:00:00"
        base_dir = join(
            output_dir,
            f"{date_str}_{time_str}"
            )

        lvls_temp = pandas.read_csv(base_dir+"_temperature_new.csv")
        lvls_zWin = pandas.read_csv(base_dir+"_zonalWinds_new.csv")
        lvls_mWin = pandas.read_csv(base_dir+"_meridionalWinds_new.csv")

        #=== discretization points from 'request_era5_profiles.py'
        lon0, lon1, lat0, lat1 = 0, 0, 0, 0
        dlon = disc_param['ecmwf']['dlon']  # grid step in degrees
        dlat = disc_param['ecmwf']['dlat']
        if disc_param['ecmwf']['auto_area'] == False:
            lon0 = disc_param['ecmwf']['min_lon']
            lon1 = disc_param['ecmwf']['max_lon']
            lat0 = disc_param['ecmwf']['min_lat']
            lat1 = disc_param['ecmwf']['max_lat']
        else:
            lon0 = int(np.min([l[1] for l in (disc_param['sou_pos']+disc_param['sta_pos'])])) - 2*dlon
            lon1 = int(np.max([l[1] for l in (disc_param['sou_pos']+disc_param['sta_pos'])])) + 2*dlon
            lat0 = int(np.min([l[0] for l in (disc_param['sou_pos']+disc_param['sta_pos'])])) - 2*dlat
            lat1 = int(np.max([l[0] for l in (disc_param['sou_pos']+disc_param['sta_pos'])])) + 2*dlat
        lons = np.arange(lon0-360., lon1-360.+dlon, dlon)
        lats = np.arange(lat1, lat0-dlat, -dlat)  # decreasing in latitude following data

        nlon, nlat = int(((lon1-lon0)/dlon+1)), int(((lat1-lat0)/dlat+1))
        npts = nlon * nlat

        #=== source, receiver, discretization along... 
        stations = disc_param['sta_pos']
        sources = disc_param['sou_pos']

        #=== calculate columns on specified discretized points
        dlat_new = disc_param['range_dependent']['dlat'] # NOTE: different from above
        dlon_new = disc_param['range_dependent']['dlon'] #   in lines 191 and 192
        all_lats = [l[0] for l in (stations+sources)]
        all_lons = [l[1] for l in (stations+sources)]
        min_lat = np.min(all_lats) - 2*dlat_new
        min_lon = np.min(all_lons) - 2*dlon_new
        max_lat = np.max(all_lats) + 2*dlat_new
        max_lon = np.max(all_lons) + 2*dlon_new
        lats_new = np.arange(min_lat, max_lat+dlat_new, dlat_new)
        lons_new = np.arange(min_lon, max_lon+dlon_new, dlon_new)

        for nsou, (soulat, soulon) in enumerate(sources):
            print("\nSource #: {0} @ ({1:4.7f}, {2:4.7f})".format(
                nsou, soulat, soulon
                ))
            for nsta, (stalat, stalon) in enumerate(stations):
                print("\n -> Station #: {0} @ ({1:4.7f}, {2:4.7f})".format(
                    nsta, stalat, stalon
                    ))

                for ilon, lon in enumerate(lons_new):
                    for ilat, lat in enumerate(lats_new):
                        file_out = f"vals_ecmwf_{doy:03d}_{nsou+1:05d}_{nsta+1:04d}"\
                                    +f"_{ilat+1:04d}_{ilon+1:04d}"
                        print(f"    --> ilat, ilon={ilat+1:04d}, {ilon+1:04d}")

                        # here is where the values will be saved for all heights
                        temp = []
                        zonw = []
                        merw = []

                        # Load each layer by height 
                        start_ind = 0
                        end_ind = npts
                        for lvl in range(geom_alt.shape[0]):
                            temp_val = lvls_temp['Value'][start_ind:end_ind].to_numpy()
                            zWin_val = lvls_zWin['Value'][start_ind:end_ind].to_numpy()
                            mWin_val = lvls_mWin['Value'][start_ind:end_ind].to_numpy()

                            # Temperature 
                            z_temp = np.reshape(temp_val, (nlat, nlon))
                            f = interpolate.interp2d(lons, lats, z_temp)
                            znew_tem = f(lon, lat)[0]
                            
                            # Zonal winds 
                            z_zonal = np.reshape(zWin_val, (nlat, nlon))
                            f = interpolate.interp2d(lons, lats, z_zonal)
                            znew_zon = f(lon, lat)[0]

                            # Meridonial winds 
                            z_mer = np.reshape(mWin_val, (nlat, nlon))
                            f = interpolate.interp2d(lons, lats, z_mer)
                            znew_mer = f(lon, lat)[0]
                            
                            temp.append(znew_tem)
                            zonw.append(znew_zon)
                            merw.append(znew_mer)

                            start_ind = (lvl+1)*npts
                            end_ind = start_ind+npts

                        # flip levels (low to high)
                        temp = np.flip(temp)
                        zonw = np.flip(zonw)
                        merw = np.flip(merw)

                        # fill data table
                        table_name = f"{file_out}.csv"
                        headerstr = [
                            'Height [km]', 'Temperature [K]', 
                            'Zonal wind [m/s]', 'Meridonial wind [m/s]'
                            ]
                        data_mat = np.zeros((geom_alt.shape[0], 4))
                        data_mat[:, 0] = geom_alt_flip
                        data_mat[:, 1] = temp
                        data_mat[:, 2] = zonw
                        data_mat[:, 3] = merw
                        df = pandas.DataFrame(data=data_mat, columns=headerstr)
                        df.to_csv(join("../output/ecmwf", table_name), index=False)
                        
                        # save in strict formatting too for reading with fortran
                        table_name_txt = f"{file_out}.dat"
                        with open(join("../output/ecmwf", table_name_txt), 'w') as f:
                            for lvl in range(geom_alt.shape[0]):
                                f.write(
                                    f"{geom_alt_flip[lvl]:8.3f} "
                                    f"{temp[lvl]:8.3f} "
                                    f"{zonw[lvl]:8.3f} "
                                    f"{merw[lvl]:8.3f}\n"
                                    )

if __name__ == '__main__':
    par_disc = toml.load("../input/discretize_parameters.toml")
    if par_disc['range_dependent']['use_rng_dep'] == False:
        main_rng_ind()
    elif par_disc['range_dependent']['use_rng_dep'] == True:
        main_rng_dep()
