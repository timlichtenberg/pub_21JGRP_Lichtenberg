#!/usr/bin/env python3
from modules_plot import *

#====================================================================
def plot_global( host_dir, sub_dirs ):

    # Plot settings
    fscale     = 1.1
    fsize      = 18
    fs_title   = 18
    fs_legend  = 11
    width      = 12.00 #* 3.0/2.0
    height     = 10.0 #/ 2.0
    # Subplot titles
    title_fs   = 12
    title_xy   = (0.07, 0.02)
    title_x    = title_xy[0]
    title_y    = title_xy[1]
    title_xycoords = 'axes fraction'
    title_ha   = "left"
    title_va   = "bottom"
    title_font = 'Arial'
    txt_alpha  = 0.5
    txt_pad    = 0.1
    label_fs   = 11
    zorder_txt = 100

    fig_o = FigureData( 3, 2, width, height, '../figures/fig_5', units='yr' )
    fig_o.fig.subplots_adjust(wspace=0.05,hspace=0.1)

    ax0 = fig_o.ax[0][0]
    ax1 = fig_o.ax[1][0]
    ax2 = fig_o.ax[2][0]
    ax3 = fig_o.ax[0][1]
    ax4 = fig_o.ax[1][1]
    ax5 = fig_o.ax[2][1]

    title_xNumbering = 0.02
    title_yNumbering = 0.05
    fsplus = 5
    ax0.text(title_xNumbering, title_yNumbering, 'A', color="k", rotation=0, ha="left", va="bottom", fontsize=label_fs+fsplus, transform=ax0.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))
    ax1.text(title_xNumbering, title_yNumbering, 'B', color="k", rotation=0, ha="left", va="bottom", fontsize=label_fs+fsplus, transform=ax1.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))
    ax2.text(title_xNumbering, title_yNumbering, 'C', color="k", rotation=0, ha="left", va="bottom", fontsize=label_fs+fsplus, transform=ax2.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))
    ax3.text(title_xNumbering, title_yNumbering, 'D', color="k", rotation=0, ha="left", va="bottom", fontsize=label_fs+fsplus, transform=ax3.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))
    ax4.text(title_xNumbering, title_yNumbering, 'E', color="k", rotation=0, ha="left", va="bottom", fontsize=label_fs+fsplus, transform=ax4.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))
    ax5.text(title_xNumbering, title_yNumbering, 'F', color="k", rotation=0, ha="left", va="bottom", fontsize=label_fs+fsplus, transform=ax5.transAxes, bbox=dict(fc='white', ec="white", alpha=0.01, pad=0.1, boxstyle='round'))

    # Loop over all subdirectories
    for sub_dir in sub_dirs:

        if sub_dir == "N2_reduced":
            color_strength = 3
            lw             = 1.5
            ls             = "--"
        else:
            color_strength = 6
            lw             = 1.5
            ls             = "-"

        output_dir = host_dir+sub_dir
        print(output_dir)

        fig_o.time = get_all_output_times(output_dir)

        print("---------------------------------------------------------------")
        print(sub_dir, "times:", len(fig_o.time), ",", np.min(fig_o.time)/1e+6, "–", np.max(fig_o.time)/1e+6, "Myr")
        print("---------------------------------------------------------------")

        # Read in runtime helpfile and separate in atmosphere and interior params
        df = pd.read_csv(output_dir+"/runtime_helpfile.csv", sep=" ")
        df_int = df.loc[df['Input']=='Interior']
        df_atm = df.loc[df['Input']=='Atmosphere']
        
        # Remove duplicate atm entries for one timestep
        for idx, row in df_atm.iterrows():
            if len(df_atm.loc[df_atm["Time"] == int(row["Time"])]) > 1:
                df_atm = df_atm.drop(idx)

        ########## Global properties
        keys_t = ( ('atmosphere','mass_liquid'),
                   ('atmosphere','mass_solid'),
                   ('atmosphere','mass_mantle'),
                   ('atmosphere','mass_core'),
                   ('atmosphere','temperature_surface'),
                   ('atmosphere','emissivity'),
                   ('rheological_front_phi','phi_global'),
                   ('atmosphere','Fatm'),
                   ('atmosphere','pressure_surface'),
                   ('rheological_front_dynamic','depth'),
                   ('rheological_front_dynamic','mesh_index')
                   )
        data_a = get_dict_surface_values_for_times( keys_t, fig_o.time, output_dir )
        mass_liquid         = data_a[0,:]
        mass_solid          = data_a[1,:]
        mass_mantle         = data_a[2,:]
        mass_core           = data_a[3,:]
        T_surf              = data_a[4,:]
        emissivity          = data_a[5,:]
        phi_global          = data_a[6,:]
        Fatm                = data_a[7,:]
        P_surf              = data_a[8,:]
        rheol_front         = data_a[9,:]
        rheol_front_idx     = data_a[10,:]


        xlabel = r'Time, $t$ (yr)'
        xlim = (5e1,1e7)

        red = (0.5,0.1,0.1)
        blue = (0.1,0.1,0.5)
        black = 'black'

        xcoord_l = -0.10
        ycoord_l = 0.5
        xcoord_r = 1.11
        ycoord_r = 0.5

        nsteps       = 10

        # Replace NaNs
        for idx, val in enumerate(T_surf):
            # print(idx, val)
            if np.isnan(val):
                json_file_time = MyJSON( output_dir+'/{}.json'.format(fig_o.time[idx]) )
                int_tmp   = json_file_time.get_dict_values(['data','temp_b'])
                # print("T_surf:", idx, val, "-->", round(int_tmp[0],3), "K")
                T_surf[idx] = int_tmp[0]

        ##########
        # figure a
        ##########
        title = r'Heat flux to space'  
        if nsteps > 1:

            Fatm_atm_rolling = np.convolve(df_atm["F_atm"], np.ones((nsteps,))/nsteps, mode='valid')
            Time_atm_rolling = np.convolve(df_atm["Time"], np.ones((nsteps,))/nsteps, mode='valid')
            ax0.loglog( Time_atm_rolling, Fatm_atm_rolling, color=vol_colors[sub_dir][color_strength], lw=lw, ls=ls, label=vol_latex[sub_dir], zorder=vol_zorder[sub_dir] )
        else:
            ax0.loglog( df_atm["Time"], df_atm["F_atm"], color=vol_colors[sub_dir][color_strength], lw=lw, ls=ls, alpha=1.0, label=vol_latex[sub_dir], zorder=vol_zorder[sub_dir] )
            
        # fig_o.set_myaxes(ax0)
        ax0.set_ylabel(r'$F_\mathrm{atm}^{\uparrow}$ (W m$^{-2}$)', fontsize=label_fs)
        ax0.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20) )
        ax0.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8), numticks=20))
        ax0.set_xlim( *xlim )
        ax0.set_xticklabels([])
        ax0.yaxis.set_label_coords(xcoord_l,ycoord_l)
        ax0.set_title(title, fontname=title_font, fontsize=title_fs, x=title_x, y=title_y, ha=title_ha, va=title_va, bbox=dict(fc='white', ec="white", alpha=txt_alpha, pad=txt_pad), zorder=zorder_txt)
        ax0.set_ylim(8e-1, 1e+7)
        yticks = [ 1e+0, 1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7 ]
        ax0.set_yticks( yticks )

        # # # SHOW LEGEND
        ax0.legend(ncol=2, loc=1, frameon=1, fancybox=True, framealpha=0.9, fontsize=fs_legend)

        ##########
        # figure b
        ##########
        # T_surf = T_surf[np.logical_not(np.isnan(T_surf))]
        title = r'Surface temperature'
        if nsteps > 1:
            Time_int_rolling = np.convolve(df_int["Time"], np.ones((nsteps,))/nsteps, mode='valid')
            Ts_int_rolling = np.convolve(df_int["T_surf"], np.ones((nsteps,))/nsteps, mode='valid')
            Time_atm_rolling = np.convolve(df_atm["Time"], np.ones((nsteps,))/nsteps, mode='valid')
            Ts_atm_rolling = np.convolve(df_atm["T_surf"], np.ones((nsteps,))/nsteps, mode='valid')

            h1, = ax1.semilogx(Time_atm_rolling, Ts_atm_rolling, ls=ls, lw=lw, color=vol_colors[sub_dir][color_strength], label=vol_latex[sub_dir], zorder=vol_zorder[sub_dir])
        else:
            h1, = ax1.semilogx(df_atm["Time"], df_atm["T_surf"], ls=ls, lw=lw, color=vol_colors[sub_dir][color_strength], label=vol_latex[sub_dir], zorder=vol_zorder[sub_dir])

        ax1.set_ylabel(r'$T_\mathrm{surf}$ (K)', fontsize=label_fs)
        ax1.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20) )
        ax1.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8), numticks=20))
        ax1.xaxis.set_minor_formatter(ticker.NullFormatter())
        ax1.set_xlim( *xlim )
        ax1.set_xticklabels([])
        ax1.yaxis.set_label_coords(xcoord_l,ycoord_l)
        ax1.set_ylim(200, 3000)
        yticks = [ 200, 500, 1000, 1500, 2000, 2500, 3000]
        ax1.set_yticks( yticks )
        ax1.set_yticklabels( [ str(int(i)) for i in yticks ] )
        ax1.set_title(title, fontname=title_font, fontsize=title_fs, x=title_x, y=title_y, ha=title_ha, va=title_va, bbox=dict(fc='white', ec="white", alpha=txt_alpha, pad=txt_pad), zorder=zorder_txt)

        ##########
        # figure c
        ##########

        # Rheological front
        RF_depth_crit = 0.01 # normalized
        RF_depth_crit_num, RF_depth_crit_idx = find_nearest(rheol_front/np.max(rheol_front), RF_depth_crit)
        RF_depth_crit_time   = fig_o.time[RF_depth_crit_idx]
        Phi_global_intersect = phi_global[RF_depth_crit_idx]

        ax2.arrow(RF_depth_crit_time, 0, 0, Phi_global_intersect, head_width=0, head_length=0, fc=vol_colors[sub_dir][color_strength], ec=vol_colors[sub_dir][color_strength], lw=1.2, alpha=0.5, ls=":")

        if sub_dir == sub_dirs[0]:
            legend_handles_ax2 = []
            rf_label = "Rheological front\nreaches surface"
            l2, = ax2.plot( 0, 0, lw=1.2, color=qgray, label=rf_label, ls=":")
            legend_handles_ax2.append(l2)
            legend_ax2 = ax2.legend(handles=legend_handles_ax2, loc=1, ncol=1, fontsize=fs_legend, framealpha=0.3)
            ax2.add_artist(legend_ax2)

        # Mante melt + solid fraction
        ax2.semilogx( fig_o.time, phi_global, color=vol_colors[sub_dir][color_strength], ls=ls, lw=lw, label=r'Melt, $\phi_{\mathrm{mantle}}$', zorder=vol_zorder[sub_dir])

        ax2.set_xlim( *xlim )
        ax2.set_ylim( bottom=0, top=1.01 )
        ax2.set_xlabel(xlabel, fontsize=label_fs)
        ax2.set_ylabel(r'$\phi_{\mathrm{mantle}}$ (wt)', fontsize=label_fs)
        ax2.yaxis.set_label_coords(xcoord_l,ycoord_l)
        handles, labels = ax2.get_legend_handles_labels()
        ax2.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20) )
        ax2.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8), numticks=20))

        title = r'Mantle melt fraction'
        ax2.set_title(title, fontname=title_font, fontsize=title_fs, x=title_x, y=title_y, ha=title_ha, va=title_va, bbox=dict(fc='white', ec="white", alpha=txt_alpha, pad=txt_pad), zorder=zorder_txt)

        ### Plot axes setup
        title_ax3 = r'Surface volatile pressure equivalent'
        title_ax4 = r'Atmosphere volatile mass fraction'
        title_ax5 = r'Interior volatile mass fraction'

        ########## Volatile species-specific plots
        
        # Check for times when volatile is present in data dumps
        for vol in volatile_species:

            vol_times   = []

            print(vol, end = " ")

            # For all times
            for sim_time in fig_o.time:

                # Define file name
                json_file = output_dir+"/"+str(int(sim_time))+".json"

                # For string check
                vol_str = '"'+vol+'"'

                # # Find the times with the volatile in the JSON file
                with open(json_file, 'rb', 0) as file, \
                    mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ) as s:
                    if s.find(vol_str.encode()) != -1:
                        keys_t = ( ('atmosphere', vol, 'liquid_kg'),
                                   ('atmosphere', vol, 'solid_kg'),
                                   ('atmosphere', vol, 'initial_kg'),
                                   ('atmosphere', vol, 'atmosphere_kg'),
                                   ('atmosphere', vol, 'atmosphere_bar')
                                   )
                        vol_times.append(sim_time)

            # Only for volatiles that are present at some point 
            if vol_times:

                # Get the data for these files
                data_vol = get_dict_surface_values_for_times( keys_t, vol_times, output_dir )
                vol_liquid_kg       = data_vol[0,:]
                vol_solid_kg        = data_vol[1,:]
                vol_initial_kg      = data_vol[2,:]
                vol_atm_kg          = data_vol[3,:]
                vol_atm_pressure    = data_vol[4,:]
                vol_interior_kg     = vol_liquid_kg   + vol_solid_kg
                vol_total_kg        = vol_interior_kg + vol_atm_kg

                ##########
                # figure d
                ##########
                ax3.semilogx( vol_times, vol_atm_pressure, color=vol_colors[vol][color_strength], ls=ls, lw=lw, label=vol_latex[sub_dir], zorder=vol_zorder[sub_dir])
                ##########
                # figure e
                ##########
                ax4.semilogx( vol_times, vol_atm_kg/vol_total_kg, lw=lw, color=vol_colors[vol][color_strength], ls=ls, label=vol_latex[vol], zorder=vol_zorder[sub_dir])
                ##########
                # figure f
                ##########
                # ax5.semilogx( fig_o.time, vol_mass_interior/vol_mass_total, lw=lw, color="gray", linestyle='-', label=r'Total')
                ax5.semilogx( vol_times, vol_interior_kg/vol_total_kg, lw=lw, color=vol_colors[vol][color_strength], ls=ls, label=vol_latex[vol], zorder=vol_zorder[sub_dir] )

        ##########
        # figure d
        ##########
        fig_o.set_myaxes( ax3)
        ax3.set_ylabel('$p^{\mathrm{i}}_{\mathrm{surf}}$ (bar)', fontsize=label_fs)
        ax3.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20) )
        ax3.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8), numticks=20))
        ax3.xaxis.set_minor_formatter(ticker.NullFormatter())
        ax3.set_xlim( *xlim )
        ax3.set_xticklabels([])
        ax3.yaxis.tick_right()
        ax3.yaxis.set_label_position("right")
        ax3.yaxis.set_label_coords(xcoord_r,ycoord_r)
        handles, labels = ax3.get_legend_handles_labels()
        ax3.set_title(title_ax3, fontname=title_font, fontsize=title_fs, x=title_x, y=title_y, ha=title_ha, va=title_va, bbox=dict(fc='white', ec="white", alpha=txt_alpha, pad=txt_pad), zorder=zorder_txt)

        fig_o.set_myaxes( ax4)
        ax4.set_ylabel(r'$X_{\mathrm{atm}}^{\mathrm{i}}/X_{\mathrm{tot}}^{\mathrm{i}}$', fontsize=label_fs)
        ax4.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20) )
        ax4.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8), numticks=20))
        ax4.xaxis.set_minor_formatter(ticker.NullFormatter())
        ax4.yaxis.tick_right()
        ax4.set_xticklabels([])
        ax4.yaxis.set_label_position("right")
        ax4.yaxis.set_label_coords(xcoord_r,ycoord_r)
        ax4.set_xlim( *xlim )
        handles, labels = ax4.get_legend_handles_labels()
        ax4.set_title(title_ax4, fontname=title_font, fontsize=title_fs, x=title_x, y=title_y, ha=title_ha, va=title_va, bbox=dict(fc='white', ec="white", alpha=txt_alpha, pad=txt_pad), zorder=zorder_txt)


        ##########
        # figure f
        ##########
        fig_o.set_myaxes( ax5, title=title_ax5)
        ax5.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20) )
        ax5.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8), numticks=20))
        ax5.xaxis.set_minor_formatter(ticker.NullFormatter())
        ax5.set_xlim( *xlim )
        ax5.yaxis.tick_right()
        ax5.yaxis.set_label_coords(xcoord_r,ycoord_r)
        ax5.yaxis.set_label_position("right")
        handles, labels = ax5.get_legend_handles_labels()
        ax5.set_xlabel(xlabel, fontsize=label_fs)
        ax5.set_ylabel(r'$X_{\mathrm{mantle}}^{\mathrm{i}}/X_{\mathrm{tot}}^{\mathrm{i}}$', fontsize=label_fs)
        ax5.set_title(title_ax5, fontname=title_font, fontsize=title_fs, x=title_x, y=title_y, ha=title_ha, va=title_va, bbox=dict(fc='white', ec="white", alpha=txt_alpha, pad=txt_pad), zorder=zorder_txt)


    fig_o.savefig(6)
    plt.close()

#====================================================================
def main():

    # Define specific one
    output_dir  = "../data/int_atm/coupler_runs/"
    
    # Define volatiles considered
    sub_dirs    = [ "H2", "H2O", "CO2", "CH4", "CO", "O2", "N2", "N2_reduced" ]

    print("Host directory:", output_dir)

    plot_global(host_dir=output_dir, sub_dirs=sub_dirs)

#====================================================================

if __name__ == "__main__":

    # Import utils- and plot-specific modules
    from modules_plot import *
    main()
