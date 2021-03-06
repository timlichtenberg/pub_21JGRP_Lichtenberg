#!/usr/bin/env python3
from modules_plot import *

### Initial conditions

# Planet age and orbit
time = { "planet": 0., "star": 100e+6 } # yr,

# Star age, yr
Tstar_range = [ 0.100e+9 ]          # yr , 4.567e+9

# Star mass, M_sun
Mstar       = 1.0 

# Planet-star distance, au
distance    = 1.0

# Surface pressure range (Pa)
P_surf    = 260e+5

# Surface temperature range (K)
T_surf    = 1000

# Volatiles considered
vol_dict    = { 
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
col_idx = 5

# Font sizes 
fs_l = 16
fs_m = 14
fs_s = 12

legend1_handles = []
legend2_handles = []

fig, ax1 = plt.subplots(1, 1, figsize=(9,6))
sns.set_style("ticks")
sns.despine()

dirs =  {
           "output":   os.getcwd()+"/../figures/", 
           "dat_dir":  os.getcwd()+"/../data/int_atm/", 
           "rad_conv": "/Users/tim/bitbucket/pcd_couple-interior-atmosphere/atm_rad_conv"
        }

# Loop through volatiles, options: "H2O", "CO2", "H2", "N2", "CH4", "CO", "O2"
vol_list = [ "H2O", "CO2", "CH4", "O2", "CO", "N2", "H2" ]
for vol_idx, vol in enumerate(reversed(vol_list)): 

    atm_file = dirs["dat_dir"]+"atm_snapshots/"+vol+"_atm.pkl"
    print(atm_file)

    # Read pickle file
    if os.path.isfile(atm_file):
      atm_file_stream = open(atm_file,'rb')
      atm = pkl.load(atm_file_stream)
      atm_file_stream.close()

    # When soc-rad-conv environment is set
    else:

        # Set current volatile to 1, others to zero
        for vol1 in vol_dict.keys():
            if vol1 == vol:
                vol_dict[vol1] = 1.0
            else:
                vol_dict[vol1] = 0.0

        # Create atmosphere object
        atm                = atmos(T_surf, P_surf, vol_dict)
        atm.toa_heating    = SocRadConv.InterpolateStellarLuminosity(Mstar, time, distance, atm.albedo_pl)
        atm_dry, atm = SocRadConv.RadConvEqm(dirs, time, atm, [], [], standalone=False, cp_dry=False, trpp=True) 
        
        # Save pkl file
        with open(atm_file, "wb") as atm_stream: pkl.dump(atm, atm_stream)

    # Temperature vs. pressure
    l1, = ax1.semilogy(atm.tmp,(atm.p)/1e+5, color=vol_colors[vol][col_idx], ls="-", lw=lw, label=vol_latex[vol])

    # Saturation vapor pressure for given temperature
    Psat_array = [ p_sat(vol, T) for T in atm.tmp  ]
    ax1.semilogy( atm.tmp, np.array(Psat_array)/1e+5, lw=lw*0.8, ls="--", color=vol_colors[vol][col_idx])

    # Find intersection
    idx = np.argwhere(np.diff(np.sign(Psat_array-atm.p))).flatten()
    if len(idx) >= 1:
      idx = idx[-1]+1
      ax1.axhline( atm.p[idx]/1e+5, xmin=atm.tmp[idx]/np.max(atm.tmp), xmax=1, lw=lw*0.5, ls=":", color=vol_colors[vol][col_idx], alpha=0.5)

    # Add to legend 1
    legend1_handles.append(l1)

# Annotate
ax1.text(0.995, 0.24, 'Saturated: moist convection', color=vol_colors["H2O"][col_idx], rotation=0, ha="right", va="bottom", fontsize=fs_s, transform=ax1.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))
ax1.text(0.995, 0.225, 'Unsatured: dry convection', color=vol_colors["H2O"][col_idx], rotation=0, ha="right", va="top", fontsize=fs_s, transform=ax1.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))

# Legends 1
legend1 = ax1.legend(handles=legend1_handles, loc=1, ncol=3, fontsize=fs_m, framealpha=0.99, title="Volatile cases")
ax1.add_artist(legend1)
title = legend1.get_title()
title.set_fontsize(fs_m)

# Legend 2
l2a, = ax1.plot([0], [0], lw=lw, ls="-", color=vol_colors["qgray"], label="Adiabat")
l2b, = ax1.plot([0], [0], lw=lw, ls="--", color=vol_colors["qgray_dark"], label=r"$p^{i}_\mathrm{sat}$")
legend2_handles.append(l2a)
legend2_handles.append(l2b)
legend2 = ax1.legend(handles=legend2_handles, loc=5, ncol=1, fontsize=fs_m, framealpha=0.0)
ax1.add_artist(legend2)

ax1.invert_yaxis()
ax1.set_xlabel(r'Temperature, $T$ (K)', fontsize=fs_m)
ax1.set_ylabel(r'Pressure, $P$ (bar)', fontsize=fs_m)
ax1.set_xlim(left=0, right=T_surf)
ax1.set_ylim(top=1e-5, bottom=atm.ps*1.01/1e+5)
ax1.set_xticks([10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
ax1.set_yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 260])
ax1.set_xticklabels([10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000], fontsize=fs_s)
ax1.set_yticklabels([r"$10^{-5}$", r"$10^{-4}$", r"$10^{-3}$", r"$10^{-2}$", r"$10^{-1}$", r"$10^{0}$", r"$10^{1}$", r"$10^{2}$", "260"], fontsize=fs_s)

plt.savefig(dirs["output"]+"fig_3.pdf", bbox_inches="tight")
plt.close(fig)
