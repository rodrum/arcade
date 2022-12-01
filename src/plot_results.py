"""
Code to plot resulting arrivals and raypaths after running ARCADE with infraGA.
"""

import numpy as np 
import matplotlib.pyplot as plt
from os.path import join
import toml
from itertools import product
    
def main():
    arrv_path = "../output/profiles"
    rays_path = "../output/profiles"
    fig_path = "../output/figures/"

    plot_arrivals = toml.load("../input/arcade_config.toml")['plot_arrivals']
    params = toml.load("../input/discretize_parameters.toml")

    if plot_arrivals == True:
        secs = params['sec']
        doys = params['doys']
        sources = params['sou_pos']
        stations = params['sta_pos']
        all_comb = [combi for combi in product(secs, doys, range(len(sources)), range(len(stations)))]
        run_type = toml.load("../input/arcade_config.toml")['run_type']
        rng_dep = params['range_dependent']['use_rng_dep']
        use_ecmwf = params['ecmwf']['use_ecmwf']


        end_str = ".arrivals.dat"
        if rng_dep == False:
            if run_type == "pert":
                end_str = "_pert.arrivals.dat"
            elif use_ecmwf == True:
                end_str = "_mix.arrivals.dat"   
        else:
            end_str = "_.arrivals.dat"


        # Arrivals (Ground intercepts)

        prof_num = 0
        for sec, doy, isou, ista in all_comb:
            sou_lat, sou_lon = sources[isou][0], sources[isou][1]
            sta_lat, sta_lon = stations[ista][0], stations[ista][1]
            arrv1_name = join(
                arrv_path,
                f"1_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}{end_str}"
                )
            arrv2_name = join(
                arrv_path,
                f"2_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}{end_str}"
                )
            arrv_1, arrv_2 = [], []
            try:
                arrv_1 = np.loadtxt(arrv1_name)
                print(f"File {arrv1_name} loaded")
            except FileNotFoundError:
                print(f"ERROR: File {arrv1_name} not found!")

            try:
                arrv_2 = np.loadtxt(arrv2_name)
                print(f"File {arrv2_name} loaded")
            except FileNotFoundError:
                print(f"ERROR: File {arrv2_name} not found!")

            """
            Fields are:
            0: incl [deg]
            1: az [deg]
            2: n_b
            3: lat_0 [deg] 
            4: lon_0 [deg] 
            5: time [s]    
            6: cel [km/s]  
            7: turning ht [km] 
            8: inclination [deg]   
            9: back azimuth [deg]  
            10: trans. coeff. [dB]  
            11: absorption [dB]
            """

            if len(arrv_1) != 0 and len(arrv_2):
                # Plot spatially 
                # lat_0 = arrv_1[:, 3]
                # lon_0 = arrv_1[:, 4]
                lat_0 = np.append(arrv_1[:, 3], arrv_2[:, 3])
                lon_0 = np.append(arrv_1[:, 4], arrv_2[:, 4])
                # wrap negative longitude values to +180 degrees
                lon_0 = [ilon if ilon>0 else ilon+360 for ilon in lon_0]
                sou_lon = sou_lon if sou_lon>0 else sou_lon+360
                sta_lon = sta_lon if sta_lon>0 else sta_lon+360

                # trn_c = 10**(-arrv_1[:, 10]/20)
                # absor = 10**(-arrv_1[:, 11]/20)
                trn_c = 10**(-np.append(arrv_1[:, 10], arrv_2[:, 10])/20)
                absor = 10**(-np.append(arrv_1[:, 11], arrv_2[:, 11])/20)
                t_loss = +20*np.log10(trn_c+absor)
                # time = arrv_1[:, 5]/3600.0 # s->h
                # cel = arrv_1[:, 6]
                # turn_h = arrv_1[:, 7]
                time = np.append(arrv_1[:, 5], arrv_2[:, 5])/3600.0 # s->h
                cel = np.append(arrv_1[:, 6], arrv_2[:, 6])
                turn_h = np.append(arrv_1[:, 7], arrv_2[:, 7])

                fig, axs = plt.subplots(
                    num=prof_num,
                    nrows=2,
                    ncols=2,
                    figsize=(8,5),
                    # gridspec_kw = {"width_ratios" : [1, 0.05]},
                    )

                ax0 = axs[0, 0]
                ax1 = axs[0, 1]
                ax2 = axs[1, 0]
                ax3 = axs[1, 1]

                ax0.tick_params(direction='in')
                ax0.set_facecolor('xkcd:black')
                ax1.tick_params(direction='in')
                ax1.set_facecolor('xkcd:black')
                ax2.tick_params(direction='in')
                ax2.set_facecolor('xkcd:black')
                ax3.tick_params(direction='in')
                ax3.set_facecolor('xkcd:black')

                # take y-axis out
                ax1.set_yticklabels([])
                ax3.set_yticklabels([])


                sc0 = ax0.scatter(
                    x=lon_0,
                    y=lat_0,
                    c=t_loss,
                    marker='.',
                    cmap='viridis',
                    )
                sc1 = ax1.scatter(
                    x=lon_0,
                    y=lat_0,
                    c=time,
                    marker='.',
                    cmap='plasma',
                    )
                sc2 = ax2.scatter(
                    x=lon_0,
                    y=lat_0,
                    c=cel,
                    marker='.',
                    cmap='cividis',
                    )
                sc3 = ax3.scatter(
                    x=lon_0,
                    y=lat_0,
                    c=turn_h,
                    vmin=0,
                    vmax=170,
                    marker='.',
                    cmap='tab20c',
                    )

                # Add source and station
                ax0.plot(sou_lon, sou_lat, '^r')
                ax0.plot(sta_lon, sta_lat, 'vb')
                ax1.plot(sou_lon, sou_lat, '^r')
                ax1.plot(sta_lon, sta_lat, 'vb')
                ax2.plot(sou_lon, sou_lat, '^r')
                ax2.plot(sta_lon, sta_lat, 'vb')
                ax3.plot(sou_lon, sou_lat, '^r')
                ax3.plot(sta_lon, sta_lat, 'vb')

                ax0.set_ylabel('Lat. ($^\circ$)')
                ax2.set_ylabel('Lat. ($^\circ$)')

                ax2.set_xlabel('Lon. ($^\circ$)')
                ax3.set_xlabel('Lon. ($^\circ$)')

                ax0.grid() 
                ax1.grid() 
                ax2.grid() 
                ax3.grid() 
                # ax0.axis('equal')
                # ax1.axis('equal')
                # ax2.axis('equal')
                
                plt.colorbar(sc0, ax=ax0, label='TL [dB]')
                plt.colorbar(sc1, ax=ax1, label='Time [hr]')
                plt.colorbar(sc2, ax=ax2, label='Cel. [km/s]')
                plt.colorbar(sc3, ax=ax3, label='Turn. h. [km]')

                # title
                title = f"Arrivals: {sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}{end_str}"
                plt.suptitle(title, fontsize=14)

                plt.tight_layout(rect=[0, 0.03, 1, 0.95])

                name_fig = join(
                    fig_path,
                    f"arrv_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}.png"
                    )

                plt.savefig(name_fig, dpi=300, bbox_inches='tight')
                prof_num += 1
                print(f"   > Figure {name_fig} saved")

            
            # Rays (Ground intercepts) NOTE: incomplete
        prof_num = 0
        for sec, doy, isou, ista in all_comb:
            sou_lat, sou_lon = sources[isou][0], sources[isou][1]
            sta_lat, sta_lon = stations[ista][0], stations[ista][1]
            rays1_name = join(
                rays_path,
                f"1_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}{end_str}"
                )
            rays2_name = join(
                rays_path,
                f"2_prof_{sec:05d}_{doy:03d}_{isou+1:05d}_{ista+1:04d}{end_str}"
                )
            rays_1, rays_2 = [], []
            try:
                rays_1 = np.loadtxt(rays1_name)
                print(f"File {rays1_name} loaded")
            except FileNotFoundError:
                print(f"ERROR: File {rays1_name} not found!")

            try:
                rays_2 = np.loadtxt(rays2_name)
                print(f"File {rays2_name} loaded")
            except FileNotFoundError:
                print(f"ERROR: File {rays2_name} not found!")

            """
            Fields are:
            0: lat [deg] 
            1: lon [deg]   
            3: z [km]  
            4: trans. coeff. [dB]  
            5: absorption [dB]
            6: time [s]
            """
    else:
        print(f"Note: skipping arrival plots (results)...")


