"""
Group of functions to discretize different cases to calculate profiles

Input GTD8D (NRLMSIS2.0)
========================
iyd     integer, year and day YYDDD
sec     integer, universal time in seconds
alt     float, geodetic altitude [km]
glat    float, geodetic latitude [deg]
glon    float, geodetic longitude [deg]
stl     float, local solar time (IGNORED)
f107a   float, 81 day average, centered on input time, of f10.7 solar activity
               index
f107    float, daily f10.7 for previous day
apd     float, daily geomagnetic activity
mass    float, mass number (IGNORED)


Input HWM14
===========
iyd     integer, year and day as yyddd
sec     integer, ut(sec)
alt     float, altitude(km)
glat    float, geodetic latitude(deg)
glon    float, geodetic longitude(deg)
stl     Not used
f107a   Not used
f107    Not used
aph     current 3hr AP index
"""

import shutil
from os import mkdir, chdir, listdir
from os.path import join, isdir
import toml
from geographiclib.geodesic import Geodesic
import numpy as np
from subprocess import Popen, PIPE
import pandas
from itertools import product
from datetime import datetime, timedelta

geod = Geodesic.WGS84

def pres(dens, temp, R=287.058):
    """
    R    : specific gas constant for dry air, in J/kg*K
    dens : density in g/cm^3
    temp : temperature in K

    Returns
    -------
    pressure in mbar (10^-2 Pa)
    """
    return dens*1.0e1*R*temp

def rng_ind_clim(ds, alts, clim_params, all_comb, sou_pos, sta_pos, out_path):
    stl     = clim_params['stl']
    f107a   = clim_params['f107a']
    f107    = clim_params['f107']
    apd     = clim_params['apd']
    aph     = clim_params['aph']

    sec_doy_sou_sta = []
    for sec, doy, isou, ista in all_comb:
        iyd = (year%1000)*1000+doy
        # print(f"-> iyd={iyd}")
        sou_lat, sou_lon = sou_pos[isou][0], sou_pos[isou][1]
        sta_lat, sta_lon = sta_pos[ista][0], sta_pos[ista][1]
        file_out = f"nodes_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.txt"
        # print(f"--> ({sou_lat:.2f}, {sou_lon:.2f}) to ({sta_lat:.2f}, {sta_lon:.2f})")
        print(f"-> {sec} {doy} {isou} {ista}")
        with open(join(out_path, file_out), 'w') as f:
            # calculate nodes along arc
            line = geod.InverseLine(sou_lat, sou_lon, sta_lat, sta_lon)
            # along arc discretization
            nl = int(np.ceil(line.s13/ds/1000))  # number of points
            # to save for later
            sec_doy_sou_sta.append([sec, doy, isou, ista, nl])
            arc = np.linspace(0, line.s13, nl)
            for alt in alts:
                for s in arc:
                    g = line.Position(s, Geodesic.STANDARD | Geodesic.LONG_UNROLL)
                    f.write(
                        f"{int(iyd):>5d} "
                        f"{int(sec):>5d} "
                        f"{alt:>6.1f} "
                        f"{g['lat2']:>6.1f} "
                        f"{g['lon2']:>6.1f} "
                        f"{stl:>6.2f} "
                        f"{f107a:>6.1f} "
                        f"{f107:>6.1f} "
                        f"{apd:>6.1f} "
                        f"{aph:6.1f}\n"
                        )
            # print(f"--> saved {join(out_path, file_out)}")

    # =========================================================
    # Run calculate_sph_nodes
    # =========================================================
    print("\nRunning: calculate_sph_modes")
    print("============================")
    outs, errs = Popen(['bash', '-c', "./calculate_sph_nodes"],
        stdout=PIPE, stderr=PIPE).communicate()
    outs = outs.decode('UTF-8')
    print(outs)

    for sec, doy, isou, ista, nl in sec_doy_sou_sta:
        file_in = f"descr_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.txt"
        file_path_in = join(out_path, file_in)
        prof = np.loadtxt(file_path_in)
        print(f"Loaded\n\t{file_path_in}")

        # =======================================================
        # Two launch angles will require two intependent profiles
        # =======================================================
        file_out_1 = f"1_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.met"
        file_out_2 = f"2_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.met"
        file_path_out_1 = join(out_path_prof, file_out_1)
        file_path_out_2 = join(out_path_prof, file_out_2)

        for file_path_out in [file_path_out_1, file_path_out_2]:
            lay_ind = 0
            print(f"Writing to\n\t{file_path_out}")
            with open(file_path_out, 'w') as f:
                for i in range(alts.shape[0]):
                    # Average density [g/m^3]
                    ave_dens = np.average(prof[lay_ind:lay_ind+nl, 5])
                    # Average tempereture [K]
                    ave_temp = np.average(prof[lay_ind:lay_ind+nl, 6])
                    # Meridional winds [m/s]
                    ave_merw = np.average(prof[lay_ind:lay_ind+nl, 7])
                    # Zonal winds [m/s]
                    ave_zonw = np.average(prof[lay_ind:lay_ind+nl, 8])
                    # Average pressure [Pa]
                    ave_pres = pres(ave_dens, ave_temp)
                    # Write in file
                    f.write(f"{alts[i]:>5.1f} {ave_temp:>7.2f} {ave_zonw:>7.2f} "
                            f"{ave_merw:>7.2f} {ave_dens:>10.2e} {ave_pres:>10.2e}\n")
                    lay_ind += nl
    print("Done.")

def rng_ind_ecmwf(year, ds, h1, h2, alts, clim_params, all_comb, sou_pos, sta_pos, out_path, skip_download=False):
    stl     = clim_params['stl']
    f107a   = clim_params['f107a']
    f107    = clim_params['f107']
    apd     = clim_params['apd']
    aph     = clim_params['aph']

    print("-> Hybrid profile with ECMWF ERA 5 lower ~80 km values")
    if skip_download is False:
        import request_era5_profiles
        request_era5_profiles.main()
    import interpolate_ecmwf
    interpolate_ecmwf.main_rng_ind()

    #=== get fixed altitudes from levels file (taken from 'interpolate_ecmwf.py')
    # load levels 1/to/137
    levels_file = join("../input", 'ECMWF - L137 model level definitions.csv')
    levels = pandas.read_csv(levels_file, header=1)
    # get altitudes (levels are from top to bottom, and it's in [m])
    geom_alt = levels['Geometric Altitude [m]'].to_numpy()
    geom_alt_flip = np.flip(geom_alt)/1000.  # to km

    sec_doy_sou_sta = []

    for sec, doy, isou, ista in all_comb:
        iyd = (year%1000)*1000+doy
        # print(f"-> idy={iyd}")
        sou_lat, sou_lon = sou_pos[isou][0], sou_pos[isou][1]
        sta_lat, sta_lon = sta_pos[ista][0], sta_pos[ista][1]
        climt_out = f"nodes_climt_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.txt"
        ecmwf_out = f"nodes_ecmwf_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.txt"
        # print(f"--> ({sou_lat:.2f}, {sou_lon:.2f}) to ({sta_lat:.2f}, {sta_lon:.2f})")
        print(f"-> {sec} {doy} {isou} {ista}")
        for file_out, alts_new in [(ecmwf_out, geom_alt_flip), (climt_out, alts)]:
            with open(join(out_path, file_out), 'w') as f:
                # calculate nodes along arc
                line = geod.InverseLine(sou_lat, sou_lon, sta_lat, sta_lon)
                # along arc discretization
                nl = int(np.ceil(line.s13/ds/1000))  # number of points
                # to save for later
                sec_doy_sou_sta.append([sec, doy, isou, ista, nl])
                arc = np.linspace(0, line.s13, nl)
                for alt in alts_new:
                    for s in arc:
                        g = line.Position(s, Geodesic.STANDARD | Geodesic.LONG_UNROLL)
                        f.write(
                            f"{int(iyd):>5d} "
                            f"{int(sec):>5d} "
                            f"{alt:>9.4f} "
                            f"{g['lat2']:>6.1f} "
                            f"{g['lon2']:>6.1f} "
                            f"{stl:>6.2f} "
                            f"{f107a:>6.1f} "
                            f"{f107:>6.1f} "
                            f"{apd:>6.1f} "
                            f"{aph:6.1f}\n"
                            )
                # print(f"--> saved {join(out_path, file_out)}")

    # =========================================================
    # Run complete_ecmwf
    # NOTE: range dependent (sph)
    # =========================================================
    print("")
    print("Running: complete_ecmwf")
    print("=======================")
    outs, errs = Popen(['bash', '-c', "./complete_ecmwf"],
        stdout=PIPE, stderr=PIPE).communicate()
    outs = outs.decode('UTF-8')
    print(outs)

    for sec, doy, isou, ista, nl in sec_doy_sou_sta:
        ecmwf_in = f"descr_ecmwf_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.txt"
        climt_in = f"descr_climt_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.txt"
        file_ecmwf_in = join(out_path, ecmwf_in)
        file_climt_in = join(out_path, climt_in)

        for fnum, file_path_in in enumerate([file_ecmwf_in, file_climt_in]):
            prof = np.loadtxt(file_path_in)
            print(f"Loaded\n\t{file_path_in}")

            ptype = 'ecmwf' if fnum == 0 else 'climt'
            alts_new = geom_alt_flip if fnum == 0 else alts
            # Idea of creating this profiles for two beams before running arcade
            file_out = f"prof_{ptype}_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.txt"
            file_path_out = join(out_path_prof, file_out)

            lay_ind = 0
            print(f"Writing to\n\t{file_path_out}")
            with open(file_path_out, 'w') as f:
                for i in range(alts_new.shape[0]):
                    ave_dens = np.average(prof[lay_ind:lay_ind+nl, 5])
                    ave_temp = np.average(prof[lay_ind:lay_ind+nl, 6])
                    # note: for ecmwf the average is done for same values,
                    #       giving same result. This could be avoided, although
                    #       it works as it is.
                    # meridional winds [m/s]
                    ave_merw = np.average(prof[lay_ind:lay_ind+nl, 7])
                    # zonal winds [m/s]
                    ave_zonw = np.average(prof[lay_ind:lay_ind+nl, 8])
                    ave_pres = pres(ave_dens, ave_temp)
                    # Write in file
                    f.write(f"{alts_new[i]:>9.4f} {ave_temp:>7.2f} {ave_zonw:>7.2f} "
                            f"{ave_merw:>7.2f} {ave_dens:>10.2e} {ave_pres:>10.2e}\n")
                    lay_ind += nl
            print("Done.")

        #=== load both output profiles (climat and ecmwf) and merge
        ecmwf_file = join(
            out_path_prof,
            f"prof_ecmwf_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.txt"
            )
        ecmwf_data = np.loadtxt(ecmwf_file)
        climt_file = join(
            out_path_prof,
            f"prof_climt_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.txt"
            )
        climt_data = np.loadtxt(climt_file)
        mixed_data = []
        for i in range(ecmwf_data.shape[0]):
            height = ecmwf_data[i,0]
            if height <= h1:
                mixed_data.append(ecmwf_data[i])
            else:
                mixed_data.append(
                    (h2-height)/(h2-h1)*ecmwf_data[i]+ \
                    (height-h1)/(h2-h1)*climt_data[int(height/0.5)]
                    )
        for i in range(int(h2/0.5), climt_data.shape[0]):
            mixed_data.append(climt_data[i])
        mixed_data = np.asarray(mixed_data)
        for i in range(1,3):
            mixed_file = join(
                "../output/profiles",
                f"{i}_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}_mix.met"
                )
            np.savetxt(mixed_file, mixed_data, fmt='%9.4E')
        print("Saved {0}".format(mixed_file))

def rng_dep_clim(year, dlat, dlon, alts, clim_params, all_comb, sou_pos, sta_pos, out_path):
    stl     = clim_params['stl']
    f107a   = clim_params['f107a']
    f107    = clim_params['f107']
    apd     = clim_params['apd']
    aph     = clim_params['aph']

    all_lats = [d[0] for d in (sou_pos+sta_pos)]
    all_lons = [d[1] for d in (sou_pos+sta_pos)]
    min_lat = np.min(all_lats) - 2*dlat
    min_lon = np.min(all_lons) - 2*dlon
    max_lat = np.max(all_lats) + 2*dlat
    max_lon = np.max(all_lons) + 2*dlon
    lats = np.arange(min_lat, max_lat+dlat, dlat)
    lons = np.arange(min_lon, max_lon+dlon, dlon)
    # Save nodes to use by range dependent calculation
    np.savetxt("../output/profiles/nodes-lon.loc", lons)
    print("../output/profiles/nodes-lon.loc saved.")
    np.savetxt("../output/profiles/nodes-lat.loc", lats)
    print("../output/profiles/nodes-lat.loc saved.")

    sec_doy_sou_sta = []
    for sec, doy, isou, ista in all_comb:
        iyd = (year%1000)*1000+doy
        # print(f"-> idy={iyd}")
        # # sou_lat, sou_lon = sou_pos[isou][0], sou_pos[isou][1]
        # sta_lat, sta_lon = sta_pos[ista][0], sta_pos[ista][1]
        # print(f"--> ({sou_lat:.2f}, {sou_lon:.2f}) to ({sta_lat:.2f}, {sta_lon:.2f})")
        sec_doy_sou_sta.append([sec, doy, isou, ista])
        print(f"-> {sec} {doy} {isou} {ista}")
        for ilon, ilat in product(range(len(lons)), range(len(lats))):
            lon = lons[ilon]
            lat = lats[ilat]
            file_out = f"nodes_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}"\
                +f"_{ilat+1:04d}_{ilon+1:04d}.txt"
            # print(f"    --> ilat, ilon={ilat+1:04d}, {ilon+1:04d}")
            with open(join(out_path, file_out), 'w') as f:
                for alt in alts:
                    f.write(
                        f"{int(iyd):>5d} "
                        f"{int(sec):>5d} "
                        f"{alt:>9.4f} "
                        f"{lat:>6.1f} "
                        f"{lon:>6.1f} "
                        f"{stl:>6.2f} "
                        f"{f107a:>6.1f} "
                        f"{f107:>6.1f} "
                        f"{apd:>6.1f} "
                        f"{aph:6.1f}\n"
                        )
                # print(f"--> saved {join(out_path, file_out)}")

    # =========================================================
    # Run calculate_rngdep_nodes
    # =========================================================
    print("\nRunning: calculate_rngdep_modes")
    print("===============================")
    outs, errs = Popen(['bash', '-c', "./calculate_rngdep_nodes"],
        stdout=PIPE, stderr=PIPE).communicate()
    outs = outs.decode('UTF-8')
    print(outs)

    out_path_prof = "../output/profiles/"
    if not isdir(out_path_prof):
        mkdir(out_path_prof)
    for sec, doy, isou, ista in sec_doy_sou_sta:
        prof_num = 0
        for ilon, ilat in product(range(len(lons)), range(len(lats))):
            lon = lons[ilon]
            lat = lats[ilat]
            file_in = f"descr_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}"\
                    + f"_{ilat+1:04d}_{ilon+1:04d}.txt"
            file_path_in = join(out_path, file_in)
            prof = np.loadtxt(file_path_in)
            print(f"Loaded\n\t{file_path_in}")

            file_out_1 = f"1_prof_{sec:05d}_{doy:03d}_{isou+1:05d}"\
                    + f"_{ista+1:04d}"\
                    + f"_{prof_num:d}.met"
            file_out_2 = f"2_prof_{sec:05d}_{doy:03d}_{isou+1:05d}"\
                    + f"_{ista+1:04d}"\
                    + f"_{prof_num:d}.met"
            file_path_out_1 = join(out_path_prof, file_out_1)
            file_path_out_2 = join(out_path_prof, file_out_2)

            for file_path_out in [file_path_out_1, file_path_out_2]:
                print(f"Writing to\n\t{file_path_out}")
                with open(file_path_out, 'w') as f:
                    for i in range(alts.shape[0]):
                        ave_dens = prof[i, 5]
                        ave_temp = prof[i, 6]
                        ave_merw = prof[i, 7] # meridional winds [m/s]
                        ave_zonw = prof[i, 8] # zonal winds [m/s]
                        ave_pres = pres(ave_dens, ave_temp)
                        # Write in file
                        f.write(
                            f"{alts[i]:>5.1f} "
                            f"{ave_temp:>7.2f} "
                            f"{ave_zonw:>7.2f} "
                            f"{ave_merw:>7.2f} "
                            f"{ave_dens:>10.2e} "
                            f"{ave_pres:>10.2e}\n")
            prof_num += 1
            print("Done.")

def rng_dep_ecmwf(year, dlat, dlon, dh, h1, h2, alts, clim_params, all_comb, sou_pos, sta_pos, out_path, skip_download=False):
    stl     = clim_params['stl']
    f107a   = clim_params['f107a']
    f107    = clim_params['f107']
    apd     = clim_params['apd']
    aph     = clim_params['aph']

    if skip_download is False:
        import request_era5_profiles
        request_era5_profiles.main()
    import interpolate_ecmwf
    interpolate_ecmwf.main_rng_dep()

    all_lats = [l[0] for l in (sou_pos+sta_pos)]
    all_lons = [l[1] for l in (sou_pos+sta_pos)]
    min_lat = np.min(all_lats) - 2*dlat
    min_lon = np.min(all_lons) - 2*dlon
    max_lat = np.max(all_lats) + 2*dlat
    max_lon = np.max(all_lons) + 2*dlon
    lats = np.arange(min_lat, max_lat+dlat, dlat)
    lons = np.arange(min_lon, max_lon+dlon, dlon)
    # Save nodes to use by range dependent calculation
    np.savetxt("../output/profiles/nodes-lon.loc", lons)
    print("../output/profiles/nodes-lon.loc saved.")
    np.savetxt("../output/profiles/nodes-lat.loc", lats)
    print("../output/profiles/nodes-lat.loc saved.")

    sec_doy_sou_sta = []
    #=== get fixed altitudes from levels file (taken from 'interpolate_ecmwf.py')
    # load levels 1/to/137
    levels_file = join("../input", 'ECMWF - L137 model level definitions.csv')
    levels = pandas.read_csv(levels_file, header=1)
    # get altitudes (levels are from top to bottom, and it's in [m])
    geom_alt = levels['Geometric Altitude [m]'].to_numpy()
    geom_alt_flip = np.flip(geom_alt)/1000.  # to km

    for sec, doy, isou, ista in all_comb:
        iyd = (year%1000)*1000+doy
        print(f"-> idy={iyd}")
        sou_lat, sou_lon = sou_pos[isou][0], sou_pos[isou][1]
        sta_lat, sta_lon = sta_pos[ista][0], sta_pos[ista][1]
        print(f"--> ({sou_lat:.2f}, {sou_lon:.2f}) to ({sta_lat:.2f}, {sta_lon:.2f})")
        sec_doy_sou_sta.append([sec, doy, isou, ista])
        for ilon, ilat in product(range(len(lons)), range(len(lats))):
            lon = lons[ilon]
            lat = lats[ilat]
            climt_out = f"nodes_climt_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}"\
                +f"_{ilat+1:04d}_{ilon+1:04d}.txt"
            ecmwf_out = f"nodes_ecmwf_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}"\
                +f"_{ilat+1:04d}_{ilon+1:04d}.txt"
            print(f"    --> ilat, ilon={ilat+1:04d}, {ilon+1:04d}")
            for file_out, alts_new in [(ecmwf_out, geom_alt_flip), (climt_out, alts)]:
                with open(join(out_path, file_out), 'w') as f:
                    for alt in alts_new:
                        f.write(
                            f"{int(iyd):>5d} "
                            f"{int(sec):>5d} "
                            f"{alt:>9.4f} "
                            f"{lat:>6.1f} "
                            f"{lon:>6.1f} "
                            f"{stl:>6.2f} "
                            f"{f107a:>6.1f} "
                            f"{f107:>6.1f} "
                            f"{apd:>6.1f} "
                            f"{aph:6.1f}\n"
                            )
                    print(f"--> saved {join(out_path, file_out)}")


    # =========================================================
    # Run calculate_rngdep_nodes_ecmwf
    # =========================================================
    print("")
    print("Running: calculate_rngdep_modes_ecmwf")
    print("=====================================")
    print("")
    outs, errs = Popen(['bash', '-c', "./calculate_rngdep_nodes_ecmwf"],
        stdout=PIPE, stderr=PIPE).communicate()
    outs = outs.decode('UTF-8')
    print(outs)

    """
    1. Load 'descr_climt_{doy}_{isou}_{ista}_{ilat}_{ilon}.txt and
       'descr_ecmwf_{doy}_{isou}_{ista}_{ilat}_{ilon}.txt.
    2. Load 'descr_climt_...' and 'descr_ecmwf_...' and stich together,
       saving
        "{1,2}_prof_{doy:03d}_{isou+1:05d}_{ista+1:04d}_{prof_num}.met"
    """
    out_path_prof = "../output/profiles/"
    if not isdir(out_path_prof):
        mkdir(out_path_prof)

    for sec, doy, isou, ista in sec_doy_sou_sta:
        prof_num = 0
        for ilon, ilat in product(range(len(lons)), range(len(lats))):
            lon = lons[ilon]
            lat = lats[ilat]
            climt_in = f"descr_climt_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}"\
                    + f"_{ilat+1:04d}_{ilon+1:04d}.txt"
            ecmwf_in = f"descr_ecmwf_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}"\
                    + f"_{ilat+1:04d}_{ilon+1:04d}.txt"
            climt_path_in = join(out_path, climt_in)
            ecmwf_path_in = join(out_path, ecmwf_in)
            climt_data = np.loadtxt(climt_path_in)
            ecmwf_data = np.loadtxt(ecmwf_path_in)
            print(f"Loaded\n\t{climt_path_in}")
            print(f"Loaded\n\t{ecmwf_path_in}")

            #=== stitch them
            mixed_data = []
            for i in range(ecmwf_data.shape[0]):
                ecmwf_data_row = ecmwf_data[i].copy()
                height = ecmwf_data_row[0]
                climt_data_row = climt_data[int(height*2)].copy()
                if height <= h1:
                    mixed_data.append([
                        height,
                        ecmwf_data_row[2],
                        ecmwf_data_row[4],
                        ecmwf_data_row[3],
                        ecmwf_data_row[1],
                        pres(ecmwf_data_row[1], ecmwf_data_row[2])
                    ])
                else:
                    new_vals = (h2-height)/(h2-h1)*ecmwf_data_row+ \
                        (height-h1)/(h2-h1)*climt_data_row
                    mixed_data.append([
                        height,
                        new_vals[2],
                        new_vals[4],
                        new_vals[3],
                        new_vals[1],
                        pres(new_vals[1], new_vals[2])
                    ])
            for i in range(int(h2/dh), climt_data.shape[0]):
                climt_data_row = climt_data[i].copy()
                mixed_data.append([
                    climt_data_row[0],
                    climt_data_row[2],
                    climt_data_row[4],
                    climt_data_row[3],
                    climt_data_row[1],
                    pres(climt_data_row[1],climt_data_row[2])
                ])

            mixed_data = np.asarray(mixed_data)
            #=== save both copies
            for i in range(1,3):
                mixed_file = join(
                    "../output/profiles",
                    f"{i}_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}"\
                    f"_{prof_num:d}.met"
                    )
                np.savetxt(mixed_file, mixed_data, fmt='%9.4E')
            print("Saved {0}".format(mixed_file))
            prof_num += 1
            print("Done.")

def rng_ind_ncpag2s(year, ds, all_comb, sou_pos, sta_pos, path_ncpag2s, out_path):
    for sec, doy, isou, ista in all_comb:
        # print(f"-> iyd={iyd}")
        sou_lat, sou_lon = sou_pos[isou][0], sou_pos[isou][1]
        sta_lat, sta_lon = sta_pos[ista][0], sta_pos[ista][1]

        # Transform year-doy-sec to YYYY-MM-DD HOUR
        this_date = datetime(year, 1, 1) + timedelta(doy - 1)
        this_date_str = f"{this_date.year:d}-{this_date.month:02d}-{this_date.day:02d}"
        this_hour_str = f"{int(sec/3600):d}"

        # Number of points based on ds
        # calculate nodes along arc
        line = geod.InverseLine(sou_lat, sou_lon, sta_lat, sta_lon)
        # along arc discretization
        nl = int(np.ceil(line.s13/ds/1000))  # number of points


        # =========================================================
        # Run ncpag2s.py
        # =========================================================
        print("")
        print("Running: ncpag2s")
        print("================")
        print("")
        print([f"python {join('..',path_ncpag2s,'ncpag2s.py')} "
                            f"line --date {this_date_str} --hour {this_hour_str} "
                            f"--startlat {np.min([sou_lat, sta_lat])} "
                            f"--startlon {np.min([sou_lon, sta_lon])} "
                            f"--endlat {np.max([sou_lat, sta_lat])} "
                            f"--endlon {np.max([sou_lon, sta_lon])} "
                            f"--points {nl} "
                            f"--output {out_path} "
                            f"--outputformat ncpaprop"])
        outs, errs = Popen([f"python {join('..',path_ncpag2s,'ncpag2s.py')} "
                            f"line --date {this_date_str} --hour {this_hour_str} "
                            f"--startlat {np.min([sou_lat, sta_lat])} "
                            f"--startlon {np.min([sou_lon, sta_lon])} "
                            f"--endlat {np.max([sou_lat, sta_lat])} "
                            f"--endlon {np.max([sou_lon, sta_lon])} "
                            f"--points {nl} "
                            f"--output {out_path} "
                            f"--outputformat ncpaprop"],
            shell=True, stdout=PIPE, stderr=PIPE).communicate()
        outs = outs.decode('UTF-8')
        print(outs)

        # Read summary file to know the order of the names
        summ_files = []
        with open(join(out_path, 'summary.dat'), 'r') as f:
            for line in f.readlines():
                summ_files.append(line.split())
        # Get altitudes from one of the profile files
        aux_prof = np.loadtxt(join(out_path, summ_files[0][1]))
        alts = aux_prof[:,0]

        # =======================================================
        # Two launch angles will require two independent profiles
        # =======================================================
        # NOTE: both profiles are the same. It's a good approximation
        #       but it would be better to redefine them depending on the
        #       launch angle...
        file_out_1 = f"1_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.met"
        file_out_2 = f"2_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.met"
        file_path_out_1 = join('../output/profiles', file_out_1)
        file_path_out_2 = join('../output/profiles', file_out_2)

        for file_path_out in [file_path_out_1, file_path_out_2]:
            print(f"Writing to\n\t{file_path_out}")
            with open(file_path_out, 'w') as f:
                for i in range(alts.shape[0]):
                    ave_dens = 0
                    ave_temp = 0
                    ave_merw = 0
                    ave_zonw = 0
                    ave_pres = 0
                    for dist_in, file_in in summ_files:
                        prof_in = np.loadtxt(join(out_path, file_in))
                        # density [g/m^3]
                        ave_dens += prof_in[i, 4]
                        # tempereture [K]
                        ave_temp += prof_in[i, 1]
                        # Meridional winds [m/s]
                        ave_merw += prof_in[i, 3]
                        # Zonal winds [m/s]
                        ave_zonw += prof_in[i, 2]
                        # pressure [Pa]
                        ave_pres += prof_in[i, 5]
                    ave_dens = ave_dens/nl
                    ave_temp = ave_temp/nl
                    ave_merw = ave_merw/nl
                    ave_zonw = ave_zonw/nl
                    ave_pres = ave_pres/nl
                    f.write(f"{alts[i]:>5.1f} {ave_temp:>7.2f} {ave_zonw:>7.2f} "
                            f"{ave_merw:>7.2f} {ave_dens:>10.2e} {ave_pres:>10.2e}\n")

        # Move output of NCPAG2S away to avoid FileExistsError
        this_date_str = f"{this_date.year:d}{this_date.month:02d}{this_date.day:02d}_{sec:05d}"
        shutil.move(join(out_path, 'summary.dat'), join(out_path, f"{this_date_str}_summary.dat"))
    print("Done.")


def rng_dep_ncpag2s(year, dlat, dlon, all_comb, sou_pos, sta_pos, path_ncpag2s, out_path):
    # create the grid
    all_lats = [d[0] for d in (sou_pos+sta_pos)]
    all_lons = [d[1] for d in (sou_pos+sta_pos)]
    min_lat = np.min(all_lats) - 2*dlat
    min_lon = np.min(all_lons) - 2*dlon
    max_lat = np.max(all_lats) + 2*dlat
    max_lon = np.max(all_lons) + 2*dlon
    lats = np.arange(min_lat, max_lat+dlat, dlat)
    lons = np.arange(min_lon, max_lon+dlon, dlon)

    for sec, doy, isou, ista in all_comb:
        # Transform year-doy-sec to YYYY-MM-DD HOUR
        this_date = datetime(year, 1, 1) + timedelta(doy - 1)
        this_date_str = f"{this_date.year:d}-{this_date.month:02d}-{this_date.day:02d}"
        this_hour_str = f"{int(sec/3600):d}"

        # =========================================================
        # Run calculate_rngdep_nodes_ecmwf
        # =========================================================
        print("")
        print("Running: ncpag2s rng-dep")
        print("========================")
        print("")
        command =[f"python {join('..', path_ncpag2s, 'ncpag2s.py')} "
                f"grid --date {this_date_str} --hour {this_hour_str} "
                f"--startlat {min_lat} "
                f"--startlon {min_lon} "
                f"--endlat {max_lat} "
                f"--endlon {max_lon} "
                f"--latpoints {int(len(lats))} "
                f"--lonpoints {int(len(lons))} "
                f"--output {out_path} "
                "--outputformat infraga"]
        print(command)
        outs, errs = Popen([f"python {join('..', path_ncpag2s, 'ncpag2s.py')} "
                            f"grid --date {this_date_str} --hour {this_hour_str} "
                            f"--startlat {min_lat} "
                            f"--startlon {min_lon} "
                            f"--endlat {max_lat} "
                            f"--endlon {max_lon} "
                            f"--latpoints {int(len(lats))} "
                            f"--lonpoints {int(len(lons))} "
                            f"--output {out_path} "
                            "--outputformat infraga"],
            shell=True, stdout=PIPE, stderr=PIPE).communicate()
        outs = outs.decode('UTF-8')
        print(outs)


        # saving
        #    "{1,2}_prof_{doy:03d}_{isou+1:05d}_{ista+1:04d}_{prof_num}.met"

        # use lats.dat and lons.dat to rename the files
        lats_dat = np.loadtxt(join(out_path, 'lats.dat'))
        lons_dat = np.loadtxt(join(out_path, 'lons.dat'))

        out_files = [f for f in listdir(out_path) if f.split('.')[-1]=='met']
        prof_num = 0
        print("Saving profiles...")
        for ilon, ilat in product(range(len(lons_dat)), range(len(lats_dat))):
            # match file
            for ifile in out_files:
                ifile_yy = int(float(ifile[4:6])) == this_date.year%1000
                ifile_mm = int(float(ifile[6:8])) == this_date.month
                ifile_dd = int(float(ifile[8:10])) == this_date.day
                ifile_hr = int(float(ifile[10:12])) == int(sec/3600)
                ifile_num = int(float(ifile.split('_')[-1].split('.')[0])) == prof_num
                if ifile_yy and ifile_mm and ifile_dd and ifile_hr and ifile_num:
                    shutil.copy(join(out_path, ifile),
                        join("../output/profiles",
                        f"1_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}"\
                        f"_{prof_num:d}.met"))
                    shutil.copy(join(out_path, ifile),
                        join("../output/profiles",
                        f"2_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}"\
                        f"_{prof_num:d}.met"))
            prof_num += 1
        print("Done.")

if __name__ == '__main__':

    # setting working directory from ../bin
    chdir("./bin")

    params = toml.load("../input/config.toml")

    year = params['discretization']['year'][0]  # NOTE: one year
    doys = params['discretization']['doys']
    doy_step = params['discretization']['doy_step']
    # Option of writing a [start, stop] and setting a doy step to create a 
    # list as in [1, 2,..., 365] for whole year or long time intervals
    if doy_step > 0:
        if len(doys) == 2:
            doys = np.arange(int(float(doys[0])), 
                             int(float(doys[1]))+doy_step, 
                             doy_step)
    secs = params['discretization']['sec']
    hmin = params['discretization']['hmin']
    hmax = params['discretization']['hmax']
    dh   = params['discretization']['dh']
    ds   = params['discretization']['ds']

    # Climatological model parameters (solad and geomagnetic activity)
    f107a   = params['atmospheric_model']['clim']['f107a']
    f107    = params['atmospheric_model']['clim']['f107']
    apd     = params['atmospheric_model']['clim']['apd']
    aph     = params['atmospheric_model']['clim']['aph']
    ## Nominal values for unused variables
    stl = 0.0
    mass = 0
    clim_params = {
        'f107a' : f107a,
        'f107'  : f107,
        'apd'   : apd,
        'aph'   : aph,
        'stl'   : stl
        }

    # Load sources
    sou_nam = params['discretization']['sou_name']
    sou_alt = params['discretization']['sou_alt']
    sou_pos = params['discretization']['sou_pos']

    # Load stations
    sta_nam = params['discretization']['sta_name']
    sta_pos = params['discretization']['sta_pos']

    print("Times of the day (hrs.):")
    for sec in secs:
        print(f"-> {sec/3600}")
    print(f"Number of days: {len(doys)}")
    print(f"-> {doys}")
    print(f"Number of sources: {len(sou_pos)}")
    for lat, lon in sou_pos:
        print(f"-> ({lat:.2f}, {lon:.2f})")
    print(f"Number of stations: {len(sta_pos)}")
    for lat, lon in sta_pos:
        print(f"-> ({lat:.2f}, {lon:.2f})")


    sec_num = len(secs)
    doy_num = len(doys)
    sou_num = len(sou_pos)
    sta_num = len(sta_pos)
    alts = np.arange(hmin, hmax+dh, dh)

    header_st = \
        "{0:>5} {1:>5} {2:>6} {3:>6} {4:>6} {5:>6} {6:>6} {7:>6} {8:>6} {9:>6}"\
        .format(
            "IYD", "SEC", "ALT", "GLAT", "GLON", "STL", "F107A", "F107", "APD", "APH"
            )
    print(f"\nSaving profiles with\n{header_st}")

    # Check that output folder and its subfolders exists or create them
    if not isdir("../output"):
        mkdir("../output")
    out_path = "../output/nodes"
    if not isdir(out_path):
        mkdir(out_path)
    out_path_prof = "../output/profiles/"
    if not isdir(out_path_prof):
        mkdir(out_path_prof)

    # ==============================
    # To facilitate read of files in
    #    calculate_sph_nodes.f90
    # ==============================
    all_comb = [combi for combi in product(secs,
                                           doys,
                                           range(len(sou_pos)),
                                           range(len(sta_pos)))]
    with open("../input/secs_doys_sources_stations.txt", "w") as f:
        for nsec, ndoy, isou, ista in all_comb:
            f.write(f"{nsec:5d} {ndoy:3d} {isou+1:5d} {ista+1:4d}\n")


    use_rngdep      = params['atmospheric_model']['prop_model'] == 'range_dep'
    recycle         = params['atmospheric_model']['ecmwf']['recycle']
    use_ecmwf       = params['atmospheric_model']['type'] == 'hybrid'
    use_ncpag2s     = params['atmospheric_model']['type'] == 'ncpag2s'

    if use_rngdep is False:
        print("\n-> Range Independent (infraga-sph) 3D ray-tracing")
        if use_ncpag2s is True:
            print("-> NCPA-G2S atmospheric descriptions")
            if not isdir(join(out_path_prof, 'ncpag2s')):
                mkdir(join(out_path_prof, 'ncpag2s'))
            path_ncpag2s = params['atmospheric_model']['ncpag2s']['path']
            rng_ind_ncpag2s(year, ds, all_comb, sou_pos, sta_pos, path_ncpag2s,
                join(out_path_prof, 'ncpag2s'))
        elif use_ecmwf is True:
            print("-> Hybrid ECMWF ERA 5 reanalysis (~0-80 km)"
                  "+ HWM14/MSIS2.0 (>~80 km) atmospheric descriptions")
            h1 = params['atmospheric_model']['ecmwf']['h1']
            h2 = params['atmospheric_model']['ecmwf']['h2']
            if recycle is True:
                print("--> WARNING: Assuming the ECMWF ERA 5 descriptions "
                      "for this area have already been downloaded.")
                print("    These profiles should be manually put in output/ecmwf.")
            rng_ind_ecmwf(year, ds, h1, h2, alts, clim_params, all_comb,
                          sou_pos, sta_pos, out_path, recycle == True)
        else:
            print("-> HWM14/MSIS2.0 atmospheric descriptions")
            rng_ind_clim(ds, alts, clim_params, all_comb, sou_pos, sta_pos, out_path)
    else:
        print("\n-> Range Dependent (infraga-sph-rngdep) 3D ray tracing")
        if use_ncpag2s is True:
            dlat = params['discretization']['range_dep']['dlat']
            dlon = params['discretization']['range_dep']['dlon']
            path_ncpag2s = params['atmospheric_model']['ncpag2s']['path']
            if not isdir(join(out_path_prof, 'ncpag2s')):
                mkdir(join(out_path_prof, 'ncpag2s'))
            rng_dep_ncpag2s(year, dlat, dlon, all_comb, sou_pos, sta_pos,
                            path_ncpag2s, join(out_path_prof, 'ncpag2s'))
        elif use_ecmwf is True:
            print("-> Hybrid ECMWF ERA 5 reanalysis (~0-80 km)"
                  " + HWM14/MSIS2.0 (>~80 km) atmospheric descriptions")
            h1 = params['atmospheric_model']['ecmwf']['h1']
            h2 = params['atmospheric_model']['ecmwf']['h2']
            dlat = params['discretization']['range_dep']['dlat']
            dlon = params['discretization']['range_dep']['dlon']
            if recycle is True:
                print("--> WARNING: Assuming the ECMWF ERA 5 descriptions "
                      "for this area have already been downloaded.")
                print("    These profiles should be manually put in output/ecmwf.")
            rng_dep_ecmwf(year, dlat, dlon, dh, h1, h2, alts, clim_params,
                          all_comb, sou_pos, sta_pos, out_path, recycle == True)
        else:
            print("-> HSM14+MSIS2.0 atmospheric descriptions")
            dlat = params['discretization']['range_dep']['dlat']
            dlon = params['discretization']['range_dep']['dlon']
            rng_dep_clim(year, dlat, dlon, alts, clim_params, all_comb, sou_pos,
                         sta_pos, out_path)
