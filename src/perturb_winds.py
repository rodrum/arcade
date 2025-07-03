"""
Script to add a Gaussian perturbation centered at 40 km by default that 
increases the "along" winds in 30% of the maximum current value.
"""

import numpy as np
from os.path import join
from obspy.geodetics.base import gps2dist_azimuth
import toml
from itertools import product

config = toml.load("./input/config.toml")

# Parameters that we could define from outside later
max_pert = 0.3
pert_hgt = 40 # km


# Get DOY, number of sources, and number of stations
#   (this means we already run 'discretize.py') 
secs = config['discretization']['sec']
doys = config['discretization']['doys']
doy_step = config['discretization']['doy_step']
# Option of writing a [start, stop] and setting a doy step to create a
# list as in [1, 2,..., 365] for whole year or long time intervals
if doy_step > 0:
    if len(doys) == 2:
        doys = np.arange(int(float(doys[0])),
                         int(float(doys[1]))+doy_step,
                         doy_step)
sources  = config['discretization']['sou_pos']
stations = config['discretization']['sta_pos']
all_comb = product(secs, doys, range(len(sources)), range(len(stations)))

# Create the profile names following 'discretize.py'
if config['atmospheric_model']['prop_model'] == 'range_ind':
    end_str = '.met' # default
    if config['atmospheric_model']['type'] == 'hybrid':
        end_str = '_mix.met'
    prof_names = []

    for sec, doy, isou, ista in all_comb:
        sou_lat, sou_lon = sources[isou][0], sources[isou][1]
        sta_lat, sta_lon = stations[ista][0], stations[ista][1]
        file_in = join(
            "./output/profiles", 
            f"1_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}{end_str}"
        )
        print(f"Reading {file_in}")
        # Azimuth from source to station
        _, a12, _ = gps2dist_azimuth(sou_lat, sou_lon, sta_lat, sta_lon)
        print(f"a12={a12:.2f} deg.")
        a12 = np.deg2rad(a12) # pass to radians
        # Load profiles, modify them, and save the perturbed version
        prof = np.loadtxt(file_in)
        h = prof[:, 0] # altitude (km)
        u = prof[:, 2] # Zonal winds, East (+)
        v = prof[:, 3] # Meridional winds, North (+)
        # Along winds
        kx = np.sin(a12)
        ky = np.cos(a12)
        print(f"kx={kx:.2f}")
        print(f"ky={ky:.2f}")
        c_along = u*kx + v*ky
        bool1 = h<pert_hgt+10 
        bool2 = h>pert_hgt-10 
        bool3 = [b1 and b2 for (b1, b2) in zip(bool1, bool2)]
        c_along_max = np.max(c_along[bool3])  # maximum that will serve to perturb
                                                # the right amount (30%)
        if c_along_max < 0:
            print("Note: along winds in perturbation heights are negative.")
            print("      perturbed profiles will be empty.")
            # Save winds
            file_out_1 = join(
                "./output/profiles", 
                f"1_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}_pert.met"
            )
            file_out_2 = join(
                "./output/profiles", 
                f"2_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}_pert.met"
            )
            for file_out in [file_out_1, file_out_2]:
                with open(file_out, 'w') as f:
                    f.write("")
        else:
            print(f"c_along_max={c_along_max:.2f}")

            # Perturbed winds
            c0 = max_pert*c_along_max
            wz = np.sqrt(2.0)*10.0 # this defines the widthas +-10 km around pert_hgt
            cz = c0*np.exp(-(h-pert_hgt)**2/wz**2)
            # So u_pert*kx + v_pert*kx = u*kx + v*ky + cz;
            # while keeping meridional winds intact
            u_pert = u+cz/kx/1.0
            v_pert = v+cz/ky/1.0
            # Save winds
            file_out_1 = join(
                "./output/profiles", 
                f"1_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}_pert.met"
            )
            file_out_2 = join(
                "./output/profiles", 
                f"2_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}_pert.met"
            )
            data_out = np.zeros(shape=prof.shape, dtype=float)
            data_out[:, 0] = prof[:, 0]
            data_out[:, 1] = prof[:, 1]
            data_out[:, 2] = u_pert
            data_out[:, 3] = v_pert
            data_out[:, 4] = prof[:, 4]
            data_out[:, 5] = prof[:, 5]
            #fmt = ["%5.1f", "%7.2f", "%7.2f", "%7.2f", "%10.2e", "%10.2e"]
            fmt = '%9.4E'
            for file_out in [file_out_1, file_out_2]:
                np.savetxt(file_out, data_out, fmt=fmt)
                print(f"{file_out} saved.")
else:
    print("Note: perturbing profiles has not been implemented")
    print("      in range dependent calculations yet.")
