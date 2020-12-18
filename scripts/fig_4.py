#!/usr/bin/env python3
from modules_plot import *

########## Read in literature data

def literature_comparison():

    Goldblatt13_Ts  = []
    Goldblatt13_OLR = []
    with open("../data/fig4_radiation_limits/Goldblatt13_data.txt", 'r') as data_file:
        for line in data_file:
            if not line.startswith('#'):
                line = line.rstrip('\n')
                line = line.split(",")
                Goldblatt13_Ts.append(float(line[0]))
                Goldblatt13_OLR.append(float(line[1]))
    Kopparapu13_Ts  = []
    Kopparapu13_OLR = []
    with open("../data/fig4_radiation_limits/Kopparapu13_data.txt", 'r') as data_file:
        for line in data_file:
            if not line.startswith('#'):
                line = line.rstrip('\n')
                line = line.split(",")
                Kopparapu13_Ts.append(float(line[0]))
                Kopparapu13_OLR.append(float(line[1]))
    Hamano15_Ts  = []
    Hamano15_OLR = []
    with open("../data/fig4_radiation_limits/Hamano15_data.txt", 'r') as data_file:
        for line in data_file:
            if not line.startswith('#'):
                line = line.rstrip('\n')
                line = line.split(",")
                Hamano15_Ts.append(float(line[0]))
                Hamano15_OLR.append(float(line[1]))

    ### Plot and annotate literature comparison
    ax1.plot(Goldblatt13_Ts, Goldblatt13_OLR, color=vol_colors["qgray"], ls=":", lw=1.0, zorder=0.1)
    ax1.text(1900, 320, "Goldblatt+ 13", va="bottom", ha="right", fontsize=7, color=vol_colors["qgray"], bbox=dict(fc='white', ec="white", alpha=0.5, pad=0.05, boxstyle='round'))
    # ax1.plot(Kopparapu13_Ts, Kopparapu13_OLR, color=vol_colors["qgray"], ls="-.", lw=1.0, zorder=0.1)
    ax1.plot(Hamano15_Ts, Hamano15_OLR, color=vol_colors["qgray"], ls="-.", lw=1.0, zorder=0.1)
    # ax1.text(2180, 330, "Hamano+ 15", va="top", ha="left", fontsize=7, color=vol_colors["qgray"], bbox=dict(fc='white', ec="white", alpha=0.5, pad=0.05, boxstyle='round'))
    ax1.text(2180, 350, "Kopparapu et al.(2013)\nHamano et al. (2015)", va="top", ha="left", fontsize=7, color=vol_colors["qgray"], bbox=dict(fc='white', ec="white", alpha=0.5, pad=0.05, boxstyle='round'))
    # ax1.text(1500, 330, "Kopparapu+ 13", va="top", ha="left", fontsize=7, color=vol_colors["qgray"], bbox=dict(fc='white', ec="white", alpha=0.5, pad=0.05, boxstyle='round'))

def define_mixing_ratios(vol, vol_list):

    # Set current volatile to 1, others to zero
    for vol1 in vol_list.keys():

        ### Pure cases
        if vol1 == vol:
            vol_list[vol1] = 1.0
        else:
            vol_list[vol1] = 0.0

        ### Mixed  cases
        if vol == "H2O-CO2":
            if vol1 == "H2O" or vol1 == "CO2":
                vol_list[vol1] = 0.5
            else:
                vol_list[vol1] = 0.0
        if vol == "H2-CO":
            if vol1 == "H2" or vol1 == "CO":
                vol_list[vol1] = 0.5
            else:
                vol_list[vol1] = 0.0
        if vol == "H2-CH4":
            if vol1 == "H2" or vol1 == "CH4":
                vol_list[vol1] = 0.5
            else:
                vol_list[vol1] = 0.0
        if vol == "H2O-H2":
            if vol1 == "H2O" or vol1 == "H2":
                vol_list[vol1] = 0.5
            else:
                vol_list[vol1] = 0.0
        if vol == "H2-N2":
            if vol1 == "H2" or vol1 == "N2":
                vol_list[vol1] = 0.5
            else:
                vol_list[vol1] = 0.0
        if vol == "CO2-N2":
            if vol1 == "CO2" or vol1 == "N2":
                vol_list[vol1] = 0.5
            else:
                vol_list[vol1] = 0.0

    # Color ranges for pure species
    if vol in vol_list.keys(): 
        vol_color = vol_colors[vol][5]
    # Specific ones for mixtures
    else:
        vol_color = vol_colors[vol]

    return vol_list, vol_color

### Initial conditions

# Check if data dirs exists, otherwise create
if not os.path.exists(dirs["data_dir"]):
    os.makedirs(dirs["data_dir"])
    print("--> Create data directory:", dirs["data_dir"])

# Planet age and orbit
time = { "planet": 0., "star": 0.100e+9 } # yr,

# Star age range, yr
star_age_range = [ 0.100e+9 ]

# Star mass range, M_sun
Mstar_range = [ 1.0 ]

# Planet-star distance range, au
distance_range = [ 1.0, 0.4 ]

# Surface pressure range (Pa) for plot A
prs_rangeA    = [ 260e+5, 1e+5 ]

# Surface pressure range (Pa) for plot B
prs_rangeB    = [ 10e+5 ]

# Surface temperature range (K)
tmp_range   = np.arange(200, 3001, 50)
tmp_range   = [ int(round(Ts)) for Ts in tmp_range ]

# Volatile molar concentrations: ! must sum to one !
vol_list    = { 
              "H2O" : .0, 
              "CO2" : .0,
              "H2"  : .0, 
              "N2"  : .0,  
              "CH4" : .0, 
              "O2"  : .0, 
              "CO"  : .0 
            }

ls_list = [ "-", "--", ":", "-." ]
lw      = 1.5
col_idx = 6


# Define volatile combinations plotted, options: 
#   Single species: "H2O", "CO2", "H2", "CH4", "N2", "CO", "O2"
#   Mixtures: "H2O-CO2", "H2-CO", "H2-CH4", "H2O-H2", "H2-N2", "CO2-N2"
# vol_array = [ "H2O", "CO2", "H2", "CH4", "N2", "CO", "O2" ]
vol_array = [ "H2O", "CO2", "H2", "CH4", "N2", "CO", "O2" ]

##### PLOT A
for setting in [ "trpp" ]:

    print("----------->>>>> Setting: ", setting)

    # Purge previous plot settings
    legendA1_handles = []
    legendA2_handles = []
    legendB1_handles = []
    legendB2_handles = []
    a_ymax = 0
    a_ymin = 0
    b_ymax = 0
    b_ymin = 0
    plt.close("all")

    # Set up new plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,6))

    # Implement setting
    if setting == "trpp":
        trpp    = True
        cp_dry  = False
    if setting == "moist":
        trpp    = False
        cp_dry  = False
    if setting == "tstep":
        trpp    = False
        cp_dry  = True

    # Loop through volatiles
    for vol_idx, vol in enumerate(vol_array): 

        # Define mixing ratios and color based on function
        vol_list, vol_color = define_mixing_ratios(vol, vol_list)

        # Loop through surface pressures
        for prs_idx, P_surf in enumerate(prs_rangeA+prs_rangeB):

            # Loop through distance range
            for dist_idx, dist in enumerate(distance_range): 

                # Loop through star masses
                for Mstar_idx, Mstar in enumerate(Mstar_range): 

                    # Loop through star ages
                    for tstar_idx, tstar in enumerate(star_age_range):

                        data_file_name = '../data/fig4_radiation_limits/TF_'+str(vol)+"_p"+str(round(P_surf/1e+5))+"_a"+str(dist)+".txt"

                        time["star"] = tstar

                        OLR_array    = []
                        NET_array    = []
                        tmp_range    = []

                        print("Read:", data_file_name)

                        with open(data_file_name, "r") as filestream:
                            next(filestream) # skip header
                            filestream = [line.rstrip('\n') for line in filestream]
                            for line in filestream:

                                data_values = line.split(" ")
                                tmp_range.append(float(data_values[0]))
                                OLR_array.append(float(data_values[1]))
                                NET_array.append(float(data_values[2]))
                                # print(data_values[0], data_values[1], data_values[2])

                        ##### Plot A: OLR
                        if P_surf in prs_rangeA and dist_idx == 0:

                            l1, = ax1.plot(tmp_range, OLR_array, color=vol_color, ls=ls_list[prs_idx], lw=lw, label=vol_latex[vol])
                            
                            # Fill color and P_surf legends each only once
                            if prs_idx == 0: 
                                legendA1_handles.append(l1)
                            if vol_idx == 0: 
                                l2, = ax1.plot([0],[0], color="gray", ls=ls_list[prs_idx], lw=lw, label=r"$a$ = "+str(dist)+" au, $P_\mathrm{surf}$ = "+str(round(P_surf/1e+5))+" bar")
                                legendA2_handles.append(l2)

                            # Literature comparison for corresponding correct settings
                            if P_surf == 260e+5 and vol == "H2O":
                                literature_comparison()

                            # Set ylim range for subplot A
                            a_ymin = np.min([ a_ymin, np.min(OLR_array) ])
                            a_ymax = np.max([ a_ymax, np.max(OLR_array) ])

                        ##### Plot B: NET flux
                        if P_surf in prs_rangeB:

                            l3, = ax2.plot(tmp_range, NET_array, color=vol_color, ls=ls_list[dist_idx], lw=lw, label=vol_latex[vol])
                            
                            # Fill color and P_surf legends each only once
                            if dist_idx == 0: 
                                legendB1_handles.append(l3)
                            if Mstar_idx == 0 and vol_idx == 0:
                                l4, = ax2.plot([0],[0], color=vol_colors["qgray"], ls=ls_list[dist_idx], lw=lw, label=r"$a$ = "+str(dist)+" au, $P_\mathrm{surf}$ = "+str(round(P_surf/1e+5))+" bar")
                                legendB2_handles.append(l4)

                            # Set ylim range for subplot B
                            b_ymin = np.min([ b_ymin, np.min(NET_array) ])
                            b_ymax = np.max([ b_ymax, np.max(NET_array) ])

    ########## GENERAL PLOT SETTINGS

    label_fs    = 12
    legend_fs   = 12
    ticks_fs    = 12
    annotate_fs = 12

    ##### PLOT A settings
    # Legend for the main volatiles
    legendA1 = ax1.legend(handles=legendA1_handles, loc=2, ncol=2, fontsize=legend_fs, title="Volatiles")
    ax1.add_artist(legendA1)
    # Legend for the line styles
    legendA2 = ax1.legend(handles=legendA2_handles, loc=4, ncol=1, fontsize=legend_fs, title="Fixed orbit")

    ax1.set_xlabel(r'Surface temperature, $T_\mathrm{surf}$ (K)', fontsize=label_fs)
    ax1.set_ylabel(r'Outgoing longwave radiation, $F^{\uparrow}_\mathrm{LW}$ (W m$^{-2}$)', fontsize=label_fs)
    ax1.set_yscale("log")
    ax1.set_xlim(left=np.min(tmp_range), right=np.max(tmp_range))
    ax1.set_ylim(top=a_ymax*10)
    ax1.set_xticks([np.min(tmp_range), 500, 1000, 1500, 2000, 2500, np.max(tmp_range)])
    ax1.tick_params(axis='both', which='major', labelsize=ticks_fs)
    ax1.tick_params(axis='both', which='minor', labelsize=ticks_fs)

    ##### PLOT B settings
    legendB2 = ax2.legend(handles=legendB2_handles, loc=4, ncol=1, fontsize=legend_fs, title="Fixed surface pressure")

    ax2.set_xlabel(r'Surface temperature, $T_\mathrm{s}$ (K)', fontsize=label_fs)
    ax2.set_ylabel(r'Outgoing net radiation, $F^{\uparrow}_\mathrm{atm}$ (W m$^{-2}$)', fontsize=label_fs)
    ax2.set_yscale("symlog")
    ax2.set_xlim(left=np.min(tmp_range), right=np.max(tmp_range))
    ax2.set_ylim(bottom=b_ymin*2., top=b_ymax*10)
    ax2.set_xticks([np.min(tmp_range), 500, 1000, 1500, 2000, 2500, np.max(tmp_range)])
    ax2.tick_params(axis='both', which='major', labelsize=ticks_fs)
    ax2.tick_params(axis='both', which='minor', labelsize=ticks_fs)

    # Annotate subplot numbering
    ax1.text(0.98, 0.985, 'A', color="k", rotation=0, ha="right", va="top", fontsize=20, transform=ax1.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))
    ax2.text(0.98, 0.985, 'B', color="k", rotation=0, ha="right", va="top", fontsize=20, transform=ax2.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))

    # Indicate cooling/heating regimes
    ax2.fill_between(tmp_range, 0, +1e+10, alpha=0.05, color="blue")
    ax2.text(np.max(tmp_range)*0.99, 0.3, "Planet cools down", va="bottom", ha="right", fontsize=legend_fs, color=vol_colors["qblue_dark"])
    ax2.text(np.max(tmp_range)*0.99, -0.3, "Planet heats up", va="top", ha="right", fontsize=legend_fs, color=vol_colors["qred_dark"])
    ax2.fill_between(tmp_range, 0, -1e+10, alpha=0.05, color="red")

    plt.savefig( "../figures/fig_4.pdf", bbox_inches="tight")
    plt.close(fig)
