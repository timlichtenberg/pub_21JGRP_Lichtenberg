import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import seaborn as sns
from matplotlib import cm

# plt.figure()
plt.figure(figsize=(10, 8)) 
plt.clf()
plt.subplot(111)
ax = plt.gca()

# Seaborn axis style
sns.set_style("ticks")
sns.despine()

# Color definitions: 
# https://www.codecademy.com/articles/seaborn-design-ii
# https://python-graph-gallery.com/python-colors/
# https://matplotlib.org/tutorials/colors/colormaps.html
no_colors   = 7
vol_colors = {
    "H2O"            : cm.get_cmap('PuBu', no_colors)(range(no_colors)),
    "CO2"            : cm.get_cmap("Reds", no_colors)(range(no_colors)),
    "H2"             : cm.get_cmap("Greens", no_colors)(range(no_colors)),
    "N2"             : cm.get_cmap("Purples", no_colors)(range(no_colors)),
    "O2"             : cm.get_cmap("Wistia", no_colors+2)(range(no_colors+2)),
    "CH4"            : cm.get_cmap("RdPu", no_colors)(range(no_colors)),
    "CO"             : cm.get_cmap("pink_r", no_colors)(range(no_colors)),
    "S"              : cm.get_cmap("YlOrBr", no_colors)(range(no_colors)),
    "He"             : cm.get_cmap("Greys", no_colors)(range(no_colors)),
    "NH3"            : cm.get_cmap("cool", no_colors)(range(no_colors)),
    "black_1"        : "#000000",
    "black_2"        : "#323232",
    "black_3"        : "#7f7f7f",
    "qgray"          : "#768E95",
    "qgray2"         : "#888888",
    "qblue"          : "#4283A9", # http://www.color-hex.com/color/4283a9
    "qgreen"         : "#62B4A9", # http://www.color-hex.com/color/62b4a9
    "qred"           : "#E6767A",
    "qturq"          : "#2EC0D1",
    "qorange"        : "#ff7f0e",
    "qmagenta"       : "#9A607F",
    "qyellow"        : "#EBB434",
    "qgray_dark"     : "#465559",
    "qblue_dark"     : "#274e65",
    "qgreen_dark"    : "#3a6c65",
    "qred_dark"      : "#b85e61",
    "qturq_dark"     : "#2499a7",
    "qmagenta_dark"  : "#4d303f",
    "qyellow_dark"   : "#a47d24",
    "qgray_light"    : "#acbbbf",
    "qblue_light"    : "#8db4cb",
    "qgreen_light"   : "#a0d2cb",
    "qred_light"     : "#eb9194",
    "qturq_light"    : "#57ccda",
    "qmagenta_light" : "#c29fb2",
    "qyellow_light"  : "#f1ca70",
}

# Volatile Latex names
vol_latex = {
    "H2O"   : r"H$_2$O",
    "CO2"   : r"CO$_2$",
    "H2"    : r"H$_2$" ,
    "CH4"   : r"CH$_4$",
    "CO"    : r"CO" ,
    "N2"    : r"N$_2$" ,
    "S"     : r"S"  ,
    "O2"    : r"O$_2$" ,
    "He"    : r"He" ,
    "NH3"   : r"NH$_3$"
}

data_dir = "../data/fig2_solubilities/"

# Define fitting functions
def henrys_law_linear(p_x, alpha):
    return p_x * alpha
def henrys_law(p_x, alpha, beta):
    X_v = (p_x**(1./beta)) * alpha
    return X_v

######################################################################################
######################################################################################
################################ Plotting settings ###################################
######################################################################################
######################################################################################

# Show literature fits that we use?
show_literature_fits = True

# Show fit coefficients in legend?
show_fit_coefficients = False

# Show literature comparison fits for used volatiles?
show_literature_comparison = False

# Show S?
show_S = False

# High-pressure H2O Silver+1990 data is ignored in fit, worsens fit at low pressures
# For comparison set ignore_silver90_data = False
ignore_silver90_data = True

# Same for some CO2 data
ignore_blank93_data  = True

# Min/max X/Y range
prs_min_Pa  = 1e3               # Pa
prs_max_Pa  = 1e12              # Pa
prs_min_bar = prs_min_Pa*1e-5   # bar
prs_max_bar = prs_max_Pa*1e-5   # bar
prs_min     = prs_min_Pa
prs_max     = prs_max_Pa
X_min_ppm   = 1e-6              # ppmw
X_max_ppm   = 1e+6              # ppmw

### Plot handles for legend
data_handles_title, = plt.plot([0], [0], color='k', ls='-', alpha=0.0, label=r'$\mathbf{Data}$')
# fits_handles_title, = plt.plot([0], [0], color='k', ls='-', alpha=0.0, label=r'$\mathbf{Applied \; fits}$ ($\alpha$, $\beta$)')
fits_handles_title, = plt.plot([0], [0], color='k', ls='-', alpha=0.0, label=r'$\mathbf{Volatile \; species}$')
lit_handles_title,  = plt.plot([0], [0], color='k', ls='-', alpha=0.0, label=r'$\mathbf{Literature \; comparison}$  ($\alpha$, $\beta$)')
handles_data = [data_handles_title] # data
handles_fits = []#[fits_handles_title] # fits
handles_lit  = [lit_handles_title]  # literature fits

# Plotting styles
ms      = 4.5      # marker size
lw      = 2.5      # fit line width
lw_lit  = lw-1.0   # line width literature estimates
col_idx = 4        # color index from color map to be used

######################################################################################
######################################################################################
######################################## H2O #########################################
######################################################################################
######################################################################################

### Import data 

# Moore et al. (1998) data
h2o_moore98_data             = np.loadtxt(data_dir+"H2O_Moore1998_data.csv", delimiter=",")
h2o_moore98_fH2O             = [] 
h2o_moore98_solubility       = [] 

# Liu et al. (2005) data
h2o_liu2005_data             = np.loadtxt(data_dir+"H2O_Liu2005_data.txt", delimiter=",")
h2o_liu2005_fH2O             = [] 
h2o_liu2005_solubility       = [] 

# Yamashita1999 
h2o_yamashita1999_data       = np.loadtxt(data_dir+"H2O_Yamashita_1999.txt", delimiter=",")
h2o_yamashita1999_fH2O       = [] 
h2o_yamashita1999_solubility = [] 

# Silver1990 
h2o_silver1990_data          = np.loadtxt(data_dir+"H2O_Silver_1990.txt", delimiter=",")
h2o_silver1990_fH2O          = [] 
h2o_silver1990_solubility    = [] 

# Holtz1995
h2o_holtz1995_data           = np.loadtxt(data_dir+"H2O_Holtz_1995.txt", delimiter=",")
h2o_holtz1995_fH2O           = [] 
h2o_holtz1995_solubility     = [] 

# Gardner 1999
h2o_gardner1999_data         = np.loadtxt(data_dir+"H2O_Gardner_1999.txt", delimiter=",")
h2o_gardner1999_fH2O         = [] 
h2o_gardner1999_solubility   = [] 

### Adjust units
for i in range(1, len(h2o_moore98_data)):
    h2o_moore98_fH2O.append(h2o_moore98_data[i][4]*1e5)                     # bar -> Pa
    h2o_moore98_solubility.append(h2o_moore98_data[i][2]*1e4)               # wt % -> ppm wt
for i in range(1, len(h2o_liu2005_data)):
    h2o_liu2005_fH2O.append(h2o_liu2005_data[i][0]*1e6)                     # MPa --> Pa
    h2o_liu2005_solubility.append(h2o_liu2005_data[i][1]*1e4)               # wt % --> ppm
for i in range(1, len(h2o_yamashita1999_data)):
    h2o_yamashita1999_fH2O.append(h2o_yamashita1999_data[i][0]*1e5)         # bar -> Pa
    h2o_yamashita1999_solubility.append(h2o_yamashita1999_data[i][1]*1e4)   # wt% -> ppm
for i in range(1, len(h2o_silver1990_data)):
    h2o_silver1990_fH2O.append(h2o_silver1990_data[i][0]*1e3*1e5)           # kbar -> Pa
    h2o_silver1990_solubility.append(h2o_silver1990_data[i][1]*1e4)         # wt% -> ppm
for i in range(1, len(h2o_holtz1995_data)):
    h2o_holtz1995_fH2O.append(h2o_holtz1995_data[i][0]*1e3*1e5)             # kbar -> Pa
    h2o_holtz1995_solubility.append(h2o_holtz1995_data[i][1]*1e4)           # wt% -> ppm
for i in range(1, len(h2o_gardner1999_data)):
    h2o_gardner1999_fH2O.append(h2o_gardner1999_data[i][0]*1e6)             # MPa -> Pa
    h2o_gardner1999_solubility.append(h2o_gardner1999_data[i][1]*1e4)       # wt% -> ppm

### Plot H2O data individually
data_h2o, = plt.plot(h2o_moore98_fH2O, h2o_moore98_solubility, '^', markersize=ms, color=vol_colors["H2O"][col_idx], label=vol_latex["H2O"])
plt.plot(h2o_gardner1999_fH2O, h2o_gardner1999_solubility, '^', markersize=ms, color=vol_colors["H2O"][col_idx])
plt.plot(h2o_holtz1995_fH2O, h2o_holtz1995_solubility, '^', markersize=ms, color=vol_colors["H2O"][col_idx])
plt.plot(h2o_liu2005_fH2O, h2o_liu2005_solubility, '^', markersize=ms, color=vol_colors["H2O"][col_idx])
plt.plot(h2o_yamashita1999_fH2O, h2o_yamashita1999_solubility, '^', markersize=ms, color=vol_colors["H2O"][col_idx])

### Combine data sets for fitting
h2o_pressure        = h2o_moore98_fH2O + h2o_gardner1999_fH2O + h2o_holtz1995_fH2O + h2o_yamashita1999_fH2O + h2o_liu2005_fH2O
h2o_solubility      = h2o_moore98_solubility + h2o_gardner1999_solubility + h2o_holtz1995_solubility + h2o_yamashita1999_solubility + h2o_liu2005_solubility

# Silver90 switch
if ignore_silver90_data == False:
    plt.plot(h2o_silver1990_fH2O, h2o_silver1990_solubility, '^', markersize=ms+4, color=vol_colors["H2O"][col_idx+2]) # highlight this data
    h2o_pressure    = h2o_pressure + h2o_silver1990_fH2O
    h2o_solubility  = h2o_solubility + h2o_silver1990_solubility

### Plot H2O fit
popt_h2o, pcov_h2o  = curve_fit(henrys_law, h2o_pressure, h2o_solubility)
pressure_h2o        = np.linspace(prs_min, prs_max, 100)
fit_h2o             = np.zeros(len(pressure_h2o))
for i in range(0, 100):
    fit_h2o[i]      = henrys_law(pressure_h2o[i], *popt_h2o)

if show_fit_coefficients == True:
    fit_h2o, = plt.plot(pressure_h2o, fit_h2o, color=vol_colors["H2O"][col_idx], ls='-', lw=lw, label=vol_latex["H2O"]+r' (%1.3e, %1.3e)' % tuple(popt_h2o))
else:
    fit_h2o, = plt.plot(pressure_h2o, fit_h2o, color=vol_colors["H2O"][col_idx], ls='-', lw=lw, label=vol_latex["H2O"])

### Add handles for legend
handles_data.append(data_h2o)
handles_fits.append(fit_h2o)

######################################################################################
######################################################################################
######################################## H2 ##########################################
######################################################################################
######################################################################################

### Import data 

# Hirschmann et al. (2012), Tab. 2; andesite+ O-H; fH2 (=p_H2; bar) and ln(K) (=ln(X_H2/alpha); ln(mole fraction/bar)), P (GPa), T (C), H2-wt%
h2_hirschmann_andesite_data = np.loadtxt(data_dir+"H2_Hirschmann12_andesite.txt", delimiter=", ")
h2_andesite_tab_fH2         = []
h2_andesite_tab_solubility  = []

# Hirschmann et al. (2012), Tab. 2; basalt+C-O-H: fH2 (=p_H2; bar) and ln(K) (=ln(X_H2/alpha); ln(mole fraction/bar)), P (GPa), T (°C), H2-wt%
h2_hirschmann_basalt_data   = np.loadtxt(data_dir+"H2_Hirschmann12_basalt.txt", delimiter=", ")
h2_basalt_tab_fH2           = []
h2_basalt_tab_solubility    = []

# Gaillard et al. (2003), Tab. 4
# No., fH2 (bar), SH2 (g/m^3), SH2(ppm-wt), log K (measured), log K (calculated), alpha, log K (calculated)
h2_gaillard_tab4            = np.loadtxt(data_dir+"H2_Gaillard_2003_tab4_data.txt", delimiter=", ")
h2_gaillard03_fh2           = []
h2_gaillard03_solubility    = []

### Adjust units
for i in range(1, len(h2_hirschmann_basalt_data)):
    h2_andesite_tab_fH2.append(h2_hirschmann_andesite_data[i][2]*1e9)           # GPa -> Pa
    h2_andesite_tab_solubility.append(h2_hirschmann_andesite_data[i][4]*1e4)    # wt% -> ppm wt
for i in range(1, len(h2_hirschmann_basalt_data)):
    h2_basalt_tab_fH2.append(h2_hirschmann_basalt_data[i][2]*1e9)               # GPa -> Pa
    h2_basalt_tab_solubility.append(h2_hirschmann_basalt_data[i][4]*1e4)        # wt% -> ppm 
for i in range(2, len(h2_gaillard_tab4)):
    h2_gaillard03_fh2.append(h2_gaillard_tab4[i][1]*1e5)                        # bar -> Pa
    h2_gaillard03_solubility.append(h2_gaillard_tab4[i][3])                     # ppm wt

### Plot data

data_h2, = plt.plot(h2_gaillard03_fh2, h2_gaillard03_solubility, '^', markersize=ms, color=vol_colors["H2"][col_idx], label=vol_latex["H2"])
plt.plot(h2_andesite_tab_fH2, h2_andesite_tab_solubility, '^', markersize=ms, color=vol_colors["H2"][col_idx])
plt.plot(h2_basalt_tab_fH2, h2_basalt_tab_solubility, '^', markersize=ms, color=vol_colors["H2"][col_idx])

### Combine data sets for fitting

h2_pressure             = h2_andesite_tab_fH2 + h2_basalt_tab_fH2 + h2_gaillard03_fh2
h2_solubility           = h2_andesite_tab_solubility + h2_basalt_tab_solubility + h2_gaillard03_solubility

### Plot fit

popt_h2, pcov_h2        = curve_fit(henrys_law_linear, h2_pressure, h2_solubility)
pressure_h2             = np.linspace(prs_min, prs_max, 100)
fit_h2                  = np.zeros(len(pressure_h2))
for i in range(0, 100):
    fit_h2[i] = henrys_law_linear(pressure_h2[i], *popt_h2)

if show_fit_coefficients == True:
    fit_h2, = plt.plot(pressure_h2, fit_h2, color=vol_colors["H2"][col_idx], ls='-', lw=lw, label=vol_latex["H2"]+r' (%1.3e, 1.0)' % popt_h2)
else:
    fit_h2, = plt.plot(pressure_h2, fit_h2, color=vol_colors["H2"][col_idx], ls='-', lw=lw, label=vol_latex["H2"])

### Add handles for legend
handles_data.append(data_h2)
handles_fits.append(fit_h2)

######################################################################################
######################################################################################
####################################### CO2 ##########################################
######################################################################################
######################################################################################

### Import data

# Pan et al. (1991)
co2_pan91_data_ftir = np.loadtxt(data_dir+"CO2_Pan_1991_tab2_data_FTIR.txt", delimiter=", ")
co2_pan91_fCO2_ftir_lowP        = [] 
co2_pan91_solubility_ftir_lowP  = [] 
co2_pan91_fCO2_ftir             = [] 
co2_pan91_solubility_ftir       = [] 
co2_pan91_data_sims = np.loadtxt(data_dir+"CO2_Pan_1991_tab2_data_SIMS.txt", delimiter=", ")
co2_pan91_fCO2_sims             = [] 
co2_pan91_solubility_sims       = []
co2_pan91_data_bulk = np.loadtxt(data_dir+"CO2_Pan_1991_tab2_data_bulk.txt", delimiter=", ")
co2_pan91_fCO2_bulk             = [] 
co2_pan91_solubility_bulk       = [] 

# Stolper & Holloway (1988)
co2_stolper88_data       = np.loadtxt(data_dir+"CO2_Stolper_1988_tab1.txt", delimiter=", ")
co2_stolper88_fCO2       = [] 
co2_stolper88_solubility = [] 

# Blank 1993
co2_blank93_data         = np.loadtxt(data_dir+"CO2_Blank_1993.txt", delimiter=", ")
co2_blank93_fCO2         = [] 
co2_blank93_solubility   = [] 

# Dixon 1995
co2_dixon95_data         = np.loadtxt(data_dir+"CO2_Dixon_1995.txt", delimiter=", ")
co2_dixon95_fCO2         = [] 
co2_dixon95_solubility   = [] 

# Mysen et al. (1975) 
co2_mysen75_data         = np.loadtxt(data_dir+"CO2_Mysen_1975_tab3.txt", delimiter=", ")
co2_mysen75_fCO2         = [] 
co2_mysen75_solubility   = [] 


### Adjust units
for i in range(0, 2):
    co2_pan91_fCO2_ftir_lowP.append(co2_pan91_data_ftir[i][0]*1e3*1e5)      # kbar -> Pa
    co2_pan91_solubility_ftir_lowP.append(co2_pan91_data_ftir[i][1]*1e4)    # wt % -> ppm wt
for i in range(2, len(co2_pan91_data_ftir)):
    co2_pan91_fCO2_ftir.append(co2_pan91_data_ftir[i][0]*1e3*1e5)           # kbar -> Pa
    co2_pan91_solubility_ftir.append(co2_pan91_data_ftir[i][1]*1e4)         # wt % -> ppm wt
for i in range(0, len(co2_pan91_data_sims)):    
    co2_pan91_fCO2_sims.append(co2_pan91_data_sims[i][0]*1e3*1e5)           # kbar -> Pa
    co2_pan91_solubility_sims.append(co2_pan91_data_sims[i][1]*1e4)         # wt % -> ppm wt
for i in range(0, len(co2_pan91_data_bulk)):    
    co2_pan91_fCO2_bulk.append(co2_pan91_data_bulk[i][0]*1e3*1e5)           # kbar -> Pa
    co2_pan91_solubility_bulk.append(co2_pan91_data_bulk[i][1]*1e4)         # wt % -> ppm wt
co2_pan91_fCO2_highP         = co2_pan91_fCO2_ftir + co2_pan91_fCO2_sims + co2_pan91_fCO2_bulk
co2_pan91_solubility_highP   = co2_pan91_solubility_ftir + co2_pan91_solubility_sims + co2_pan91_solubility_bulk

for i in range(0, len(co2_stolper88_data)):
    co2_stolper88_fCO2.append(co2_stolper88_data[i][0]*1e5)             # bar --> Pa
    co2_stolper88_solubility.append(co2_stolper88_data[i][1])           # ppm wt

for i in range(0, len(co2_blank93_data)):
    co2_blank93_fCO2.append(co2_blank93_data[i][0]*1e5)                       # bar -> Pa
    co2_blank93_solubility.append(co2_blank93_data[i][1]*1e-2*1e+6)           # wt% -> ppm

for i in range(0, len(co2_dixon95_data)):
    co2_dixon95_fCO2.append(co2_dixon95_data[i][0]*1e5)                 # bar -> Pa
    co2_dixon95_solubility.append(co2_dixon95_data[i][1]*1e-2*1e+6)     # wt% -> ppm

for i in range(0, len(co2_mysen75_data)):
    co2_mysen75_fCO2.append(co2_mysen75_data[i][0]*1e3*1e5)             # kbar -> Pa
    co2_mysen75_solubility.append(co2_mysen75_data[i][1]*1e-2*1e+6)     # wt.% -> ppm

### Plot data
data_co2, = plt.plot(co2_pan91_fCO2_ftir_lowP, co2_pan91_solubility_ftir_lowP, '^', markersize=ms, color=vol_colors["CO2"][col_idx], label=vol_latex["CO2"])
plt.plot(co2_pan91_fCO2_highP, co2_pan91_solubility_highP, '^', markersize=ms, color=vol_colors["CO2"][col_idx])
plt.plot(co2_stolper88_fCO2, co2_stolper88_solubility, '^', markersize=ms, color=vol_colors["CO2"][col_idx])
plt.plot(co2_dixon95_fCO2, co2_dixon95_solubility, '^', markersize=ms, color=vol_colors["CO2"][col_idx])
plt.plot(co2_mysen75_fCO2, co2_mysen75_solubility, '^', markersize=ms, color=vol_colors["CO2"][col_idx])

### Combine data sets for fit
co2_pressure        = co2_pan91_fCO2_ftir_lowP + co2_pan91_fCO2_highP + co2_stolper88_fCO2 + co2_dixon95_fCO2 + co2_mysen75_fCO2
co2_solubility      = co2_pan91_solubility_ftir_lowP + co2_pan91_solubility_highP + co2_stolper88_solubility + co2_dixon95_solubility + co2_mysen75_solubility

# Blank93 switch
if ignore_blank93_data == False:
    # Highlight this data
    plt.plot(co2_blank93_fCO2, co2_blank93_solubility, '^', markersize=8, color=vol_colors["CO2"][col_idx])
    co2_pressure    = co2_pressure + h2o_silver1990_fH2O
    co2_solubility  = co2_solubility + h2o_silver1990_solubility

### Plot fit
popt_co2, pcov_co2  = curve_fit(henrys_law, co2_pressure, co2_solubility)
pressure_co2        = np.linspace(prs_min, prs_max, 100)
fit_co2             = np.zeros(len(pressure_co2))
for i in range(0, 100):
    fit_co2[i]      = henrys_law(pressure_co2[i], *popt_co2)

if show_fit_coefficients == True:
    fit_co2, = plt.plot(pressure_co2, fit_co2, color=vol_colors["CO2"][col_idx], ls='-', lw=lw, label=vol_latex["CO2"]+r' (%1.3e, %1.3e)' % tuple(popt_co2))
else:
    fit_co2, = plt.plot(pressure_co2, fit_co2, color=vol_colors["CO2"][col_idx], ls='-', lw=lw, label=vol_latex["CO2"])

### Add handles for legend
handles_data.append(data_co2)
handles_fits.append(fit_co2)

######################################################################################
######################################################################################
####################################### CH4 ##########################################
######################################################################################
######################################################################################

### Import data

# Ardia et al. (2013), combination of Tabs 5, 6 (Method 6), and read-out from Fig. 8
# Experiment, P (GPa), fCH4 (GPa), a_C, C (ppm wt), S_CH4 (ppm wt)
ch4_ardia13_data        = np.loadtxt(data_dir+"CH4_Ardia_2013_data.txt", delimiter=", ")
ch4_ardia13_fCH4        = []
ch4_ardia13_solubility  = []

## Adjust units
for i in range(0, len(ch4_ardia13_data)):
    ch4_ardia13_fCH4.append(ch4_ardia13_data[i][1])                 # GPa 
    ch4_ardia13_solubility.append(ch4_ardia13_data[i][4])           # ppm wt
ch4_ardia13_fCH4        = [ i*1e4*1e5 for i in ch4_ardia13_fCH4 ]   # GPa --> Pa
ch4_ardia13_solubility  = [ i for i in ch4_ardia13_solubility ]     # ppm wt

### Plot data
data_ch4, = plt.plot(ch4_ardia13_fCH4, ch4_ardia13_solubility, '^', markersize=ms, color=vol_colors["CH4"][col_idx], label=vol_latex["CH4"])

### PLot fit
popt_ch4, pcov_ch4  = curve_fit(henrys_law_linear, ch4_ardia13_fCH4, ch4_ardia13_solubility)
pressure_ch4        = np.linspace(prs_min, prs_max, 100)
fit_ch4             = np.zeros(len(pressure_ch4))
for i in range(0, 100):
    fit_ch4[i] = henrys_law_linear(pressure_ch4[i], *popt_ch4)

if show_fit_coefficients == True:
    fit_ch4, = plt.plot(pressure_ch4, fit_ch4, color=vol_colors["CH4"][col_idx], ls='-', lw=lw, label=vol_latex["CH4"]+r' (%1.3e, 1.0)' % popt_ch4)
else:
    fit_ch4, = plt.plot(pressure_ch4, fit_ch4, color=vol_colors["CH4"][col_idx], ls='-', lw=lw, label=vol_latex["CH4"])

### Add handles for legend
handles_data.append(data_ch4)
handles_fits.append(fit_ch4)

######################################################################################
######################################################################################
######################################## N2 ##########################################
######################################################################################
######################################################################################

### Import data

# Libourel et al. (2003), Table 2
# fN2 (atm), N2 content (10e-9 mol g^-1), N content (ppm wt atm^-1), log(FO2) (atm)
n2_libourel03_data               = np.loadtxt(data_dir+"N2_Libourel_2003_tab2_data.txt", delimiter=", ")
n2_libourel03_fN2_lowfO2         = []
n2_libourel03_solubility_lowfO2  = []
n2_libourel03_fN2_highfO2        = []
n2_libourel03_solubility_highfO2 = []

# Li et al. (2013), Table 1
# p (kbar), N (ppm wt)
n2_li13_data                     = np.loadtxt(data_dir+"N2_Li_2013_tab1.txt", delimiter=", ")
n2_li13_fN2                      = []
n2_li13_solubility               = []

n2_mosenfelder19_data            = np.loadtxt(data_dir+"N2_Mosenfelder_2019.txt", delimiter=", ")
n2_mosenfelder19_fN2             = []
n2_mosenfelder19_solubility      = []

n2_dalou17_data                  = np.loadtxt(data_dir+"NC_Dalou_2017.txt", delimiter=", ")
n2_dalou17_fN2                   = []
n2_dalou17_solubility            = []

### Adjust units

for i in range(0, len(n2_libourel03_data)):
    # Separate by fO2 – ! ATTENTION !
    # Transition (from low to high fO2) must be masked, otherwise slope overshoot
    if n2_libourel03_data[i][3] < -9.5:
        n2_libourel03_fN2_lowfO2.append(n2_libourel03_data[i][0]*101325)    # atm -> Pa
        n2_libourel03_solubility_lowfO2.append(n2_libourel03_data[i][2])    # ppm wt
    else:
        n2_libourel03_fN2_highfO2.append(n2_libourel03_data[i][0]*101325)   # atm -> Pa
        n2_libourel03_solubility_highfO2.append(n2_libourel03_data[i][2])   # ppm wt
for i in range(0, len(n2_li13_data)):
    n2_li13_fN2.append(n2_li13_data[i][0]*1e3*1e5)                          # kbar -> Pa
    n2_li13_solubility.append(n2_li13_data[i][1])                           # ppm wt
for i in range(0, len(n2_mosenfelder19_data)):
    n2_mosenfelder19_fN2.append(n2_mosenfelder19_data[i][0]*1e6)            # MPa -> Pa
    n2_mosenfelder19_solubility.append(n2_mosenfelder19_data[i][1])         # ppm wt
for i in range(0, len(n2_dalou17_data)):
    n2_dalou17_fN2.append(n2_dalou17_data[i][0]*1e6)                        # MPa -> Pa
    n2_dalou17_solubility.append(n2_dalou17_data[i][1])                     # ppm wt

# Low fO2 data
n2_fN2_lowfO2         = n2_libourel03_fN2_lowfO2
n2_solubility_lowfO2  = n2_libourel03_solubility_lowfO2

# Very high fO2 data 
# (cf. buffer data in Li+13 Tab. 1 vs. Libourel+03)
n2_fN2_highfO2        = n2_libourel03_fN2_highfO2 + n2_li13_fN2
n2_solubility_highfO2 = n2_libourel03_solubility_highfO2 + n2_li13_solubility

# Combined data from several sources, ~mid range of fO2
n2_fN2_medfO2         = n2_mosenfelder19_fN2 + n2_dalou17_fN2 + n2_libourel03_fN2_lowfO2
n2_solubility_medfO2  = n2_mosenfelder19_solubility + n2_dalou17_solubility + n2_libourel03_solubility_lowfO2

# Fit all fO2 ranges for comparison
popt_n2_lowfO2, pcov_n2_lowfO2   = curve_fit(henrys_law, n2_fN2_lowfO2, n2_solubility_lowfO2)
popt_n2_medfO2, pcov_n2_medfO2   = curve_fit(henrys_law, n2_fN2_medfO2, n2_solubility_medfO2)
popt_n2_highfO2, pcov_n2_highfO2 = curve_fit(henrys_law, n2_fN2_highfO2, n2_solubility_highfO2)

# Create and fill approximated arrays
pressure_n2           = np.linspace(prs_min, prs_max, 100)
fit_n2_lowfO2         = np.zeros(len(pressure_n2))
fit_n2_medfO2         = np.zeros(len(pressure_n2))
fit_n2_highfO2        = np.zeros(len(pressure_n2))
for i in range(0, 100):
    fit_n2_lowfO2[i]  = henrys_law(pressure_n2[i], *popt_n2_lowfO2)
    fit_n2_medfO2[i]  = henrys_law(pressure_n2[i], *popt_n2_medfO2)
    fit_n2_highfO2[i] = henrys_law(pressure_n2[i], *popt_n2_highfO2)

# Combined fO2 data + fit
data_Nlow, = plt.plot(n2_fN2_medfO2, n2_solubility_medfO2, 'v', markersize=ms, color=vol_colors["N2"][col_idx-2], label=r"N, $f$O$_2 \lesssim$IW")

if show_fit_coefficients == True:
    fit_Nlow, = plt.plot(pressure_n2, fit_n2_medfO2, color=vol_colors["N2"][col_idx-2], ls='-', lw=lw, label=r'$N$, $f$O$_2 \lesssim$ IW (%1.3e, %1.3e)' % tuple(popt_n2_medfO2))
else:
    fit_Nlow, = plt.plot(pressure_n2, fit_n2_medfO2, color=vol_colors["N2"][col_idx-2], ls='-', lw=lw, label=r'$N$, $f$O$_2 \lesssim$ IW')

# High fO2 data for very oxidized systems
data_Nhigh, = plt.plot(n2_fN2_highfO2, n2_solubility_highfO2, '^', markersize=ms, color=vol_colors["N2"][col_idx+0], label=r"N, $f$O$_2 \gtrsim$IW" )

### ! ATTENTION ! – Hand fit for high fO2 data

# # ! This is high-P heavy, don't use !
# fit_Nhigh2, = plt.plot(pressure_n2, fit_n2_highfO2, color=vol_colors["N2"][col_idx+2], ls='--', label=r'$N$ fit, $fO_2 \gtrsim$ IW: $\alpha$=%1.3e, $\beta$=%1.3e' % tuple(popt_n2_highfO2))
# handles_fits.append(fit_Nhigh2)

# This looks better but is no automatic fit
pressure        = np.linspace(prs_min, prs_max, 100)
solubility      = np.zeros(len(pressure))
alpha           = 0.00007      # ppm/Pa
beta            = 1.8
for i in range(0, 100):
    solubility[i] = henrys_law(pressure[i], alpha, beta)


if show_fit_coefficients == True:
    fit_Nhigh, = plt.plot(pressure, solubility, color=vol_colors["N2"][col_idx+0], ls='-', lw=lw, label=r'$N$, $f$O$_2 \gtrsim$ IW (%1.3e, %1.2f)' % tuple([alpha, beta]))
else:
    fit_Nhigh, = plt.plot(pressure, solubility, color=vol_colors["N2"][col_idx+0], ls='-', lw=lw, label=r'$N$, $f$O$_2 \gtrsim$ IW')

### Add handles for legend
handles_data.append(data_Nlow)
handles_fits.append(fit_Nlow)
handles_data.append(data_Nhigh)
handles_fits.append(fit_Nhigh)

######################################################################################
######################################################################################
######################################## He ##########################################
######################################################################################
######################################################################################
### ! ATTENTION ! – He data needs to be double checked
### Data is inconsistent, comment following block in to visualize

'''###### Comment block start

### Import data
he_paonita2000_data= np.loadtxt(data_dir+"He_Paonita_2000.txt", delimiter=", ")
he_paonita2000_pressure   = []
he_paonita2000_solubility = []

### Adjust units
for i in range(0, len(he_paonita2000_data)):
    he_paonita2000_pressure.append(he_paonita2000_data[i][0])              # MPa 
    he_paonita2000_solubility.append(he_paonita2000_data[i][1])            # mol/mol
he_paonita2000_pressure    = [ i*1e6 for i in he_paonita2000_pressure ]    # MPa --> Pa
he_paonita2000_solubility  = [ i*1e6 for i in he_paonita2000_solubility ]  # mol/mol--> ppm wt 

### Plot data
data_he, = plt.plot(he_paonita2000_pressure, he_paonita2000_solubility, '^', markersize=ms, color="darkturquoise", label="$He$ data")

### Plot fit
popt_he, pcov_he  = curve_fit(henrys_law_linear, he_paonita2000_pressure, he_paonita2000_solubility)
pressure_he        = np.linspace(prs_min, prs_max, 100)
fit_he             = np.zeros(len(pressure_he))
for i in range(0, 100):
    fit_he[i] = henrys_law_linear(pressure_he[i], *popt_he)
fit_he, = plt.plot(pressure_he, fit_he, color="darkturquoise", ls='-', lw=lw, label=r'$He$ fit: $\alpha$=%1.3e, $\beta$=1.0' % popt_he)

### Add handles for legend
handles_data.append(data_he)
handles_fits.append(fit_he)

'''###### / Comment block end

#####################################################################################
############################# LITERATURE FITS #######################################
#####################################################################################

if show_literature_comparison == True:

    #####################
    ######## H2O ########
    #####################

    ######### Lebrun+ 2013 H2O fit
    pressure_h2o        = np.linspace(prs_min, prs_max, 100)
    fit_h2o             = np.zeros(len(pressure_h2o))
    henry_alpha_h2o     = 6.8e-2  # ppm / Pa
    henry_beta_h2o      = 1./0.7
    for i in range(0, 100):
        fit_h2o[i]      = henrys_law(pressure_h2o[i], henry_alpha_h2o, henry_beta_h2o)
    if show_fit_coefficients == True:
        lit_h2o_lebrun13,   = plt.plot(pressure_h2o, fit_h2o, color=vol_colors["H2O"][col_idx], ls='--', lw=lw_lit, label=vol_latex["H2O"]+r', Lebrun+ 2013 (%1.3e, %1.3e)' % tuple([henry_alpha_h2o, henry_beta_h2o]))
    else:
        lit_h2o_lebrun13,   = plt.plot(pressure_h2o, fit_h2o, color=vol_colors["H2O"][col_idx], ls='--', lw=lw_lit, label=vol_latex["H2O"])

    handles_lit.append(lit_h2o_lebrun13)

    #### Elkins-Tanton 2008 H2O fit
    pressure_h2o        = np.linspace(prs_min, prs_max, 100) # Pa
    fit_h2o             = np.zeros(len(pressure_h2o))
    alpha_et            = 2.08e-4               # wt%/Pa
    alpha_et            = alpha_et*1e-2*1e+6    # ppm/Pa
    beta_et             = 1./0.52
    for i in range(0, 100):
        fit_h2o[i]      = ((((pressure_h2o[i])**(1./beta_et)) * (alpha_et)) + 0.30) # Pa
    lit_h2o_et08, = plt.plot(pressure_h2o, fit_h2o, color=vol_colors["H2O"][col_idx], ls=':', lw=lw_lit, label=vol_latex["H2O"]+r', Elkins-Tanton 2008 (non-Henrian)')
    handles_lit.append(lit_h2o_et08)

    #####################
    ######### H #########
    #####################

    #### Hirschmann16 fit, H, very reduced
    pressure        = np.linspace(prs_min, prs_max, 100)
    solubility      = np.zeros(len(pressure))
    alpha   = 5               # ppm wt / MPa
    beta    = 1.
    alpha   = alpha * 1e-6    # --> ppm wt / Pa
    beta    = beta
    for i in range(0, 100):
        solubility[i] = henrys_law(pressure[i], alpha, beta)
    if show_fit_coefficients == True:
        lit_h_hirschmann16, = plt.plot(pressure, solubility, color=vol_colors["H2"][col_idx], ls='--', lw=lw_lit, label=r'H (IW-3.5), Hirschmann 2016 (%1.3e, 1.0)' % tuple([alpha]))
    else:
        lit_h_hirschmann16, = plt.plot(pressure, solubility, color=vol_colors["H2"][col_idx], ls='--', lw=lw_lit, label=r'H (IW-3.5)')
    handles_lit.append(lit_h_hirschmann16)

    ########### Gaillard 2003 H2 fit
    # C, oxidized
    pressure        = np.linspace(prs_min, prs_max, 100)
    solubility      = np.zeros(len(pressure))
    alpha           = 6.29e-8             # ppm wt / Pa
    beta            = 0.781
    for i in range(0, 100):
        solubility[i] = henrys_law(pressure[i], alpha, beta)
    lit_c1_gaillard03, = plt.plot(pressure, solubility, color=vol_colors["qgreen_light"], ls=':', lw=2, label=vol_latex["H2"]+r', Gaillard+2003 (%1.3e, %1.1f)' % tuple([alpha, beta]))
    handles_lit.append(lit_c1_gaillard03)


    #####################
    ######## CO2 ########
    #####################

    ##### Lebrun+13 CO2
    pressure_co2        = np.linspace(prs_min, prs_max, 100)
    fit_co2             = np.zeros(len(pressure_co2))
    henry_alpha_co2     = 4.4e-7                    #  wt / bar
    henry_alpha_co2     = henry_alpha_co2*1e+6*1e-5 #  -> ppm / Pa
    henry_beta_co2      = 1.0
    for i in range(0, 100):
        fit_co2[i] = henrys_law(pressure_co2[i], henry_alpha_co2, henry_beta_co2)
    
    if show_fit_coefficients == True:
        lit_co2_lebrun13, = plt.plot(pressure_co2, fit_co2, color=vol_colors["CO2"][col_idx], ls='--', lw=lw_lit, label=vol_latex["CO2"]+r', Lebrun+13 (%1.3e, %1.1f)' % tuple([henry_alpha_co2, henry_beta_co2]))
    else:
        lit_co2_lebrun13, = plt.plot(pressure_co2, fit_co2, color=vol_colors["CO2"][col_idx], ls='--', lw=lw_lit, label=vol_latex["CO2"])
    handles_lit.append(lit_co2_lebrun13)

    ##### Ni & Keppler 2013 – CO2"
    pressure        = np.linspace(prs_min, prs_max, 100) # Pa
    solubility      = np.zeros(len(pressure))
    alpha   = 0.155         # ppm wt / bar
    alpha   = alpha*1e-5    # -> ppm / Pa
    beta    = 1.
    for i in range(0, 100):
        solubility[i] = henrys_law(pressure[i], alpha, beta)
    if show_fit_coefficients == True:
        lit_co2_ni13, = plt.plot(pressure, solubility, color=vol_colors["CO2"][col_idx], ls=':', lw=lw_lit, label=vol_latex["CO2"]+r', Ni & Keppler 2013 (%1.3e, %1.1f)' % tuple([alpha, beta]))
    else:
        lit_co2_ni13, = plt.plot(pressure, solubility, color=vol_colors["CO2"][col_idx], ls=':', lw=lw_lit, label=vol_latex["CO2"])
    handles_lit.append(lit_co2_ni13)

    #####################
    ######## CH4 ########
    #####################

    ######### CH4 from Keppler & Golabek 2019 fit; data from Ardia et al. 2013
    pressure            = np.linspace(prs_min, prs_max, 1000)
    solubility          = np.zeros(len(pressure))
    alpha               = 0.0138       # ppm wt / bar
    alpha               = 0.0138*1e-5  # -> ppm wt / Pa
    beta                = 1.
    for i in range(0, 1000):
        solubility[i] = henrys_law(pressure[i], alpha, beta)
    if show_fit_coefficients == True:
        lit_ch4_keppler19, = plt.plot(pressure, solubility, color=vol_colors["CH4"][col_idx], ls='--', lw=lw_lit, label=vol_latex["CH4"]+r', Keppler & Golabek 2019 (%1.3e, %1.1f)' % tuple([alpha, beta]))
    else:
        lit_ch4_keppler19, = plt.plot(pressure, solubility, color=vol_colors["CH4"][col_idx], ls='--', lw=lw_lit, label=vol_latex["CH4"])
    handles_lit.append(lit_ch4_keppler19)

    #####################
    ######### C #########
    #####################

    # ########### Hirschmann 2016, American Mineralogist
    # # C, oxidized
    # pressure        = np.linspace(prs_min, prs_max, 100)
    # solubility      = np.zeros(len(pressure))
    # alpha           = 1.6                 # ppm wt / MPa
    # beta            = 1.
    # alpha           = alpha * 1e-6        # --> ppm wt / Pa
    # beta            = beta
    # for i in range(0, 100):
    #     solubility[i] = henrys_law(pressure[i], alpha, beta)
    # lit_c1_hirschmann16, = plt.plot(pressure, solubility, color=vol_colors["qgreen_light"], ls=':', lw=2, label=r'$C$, oxidized: Hirschmann 2016: $\alpha$=%1.3e' % tuple([alpha]))
    # handles_lit.append(lit_c1_hirschmann16)
    # # C, reduced
    # pressure        = np.linspace(prs_min, prs_max, 100)
    # solubility      = np.zeros(len(pressure))
    # alpha           = 0.55                # ppm wt / MPa
    # beta            = 1.
    # alpha           = alpha * 1e-6        # --> ppm wt / Pa
    # beta            = beta
    # for i in range(0, 100):
    #     solubility[i] = henrys_law(pressure[i], alpha, beta)
    # lit_c2_hirschmann16, = plt.plot(pressure, solubility, color=vol_colors["qgreen_light"], ls=':', lw=2, label=r'$C$, reduced: Hirschmann 2016: $\alpha$=%1.3e' % tuple([alpha]))
    # handles_lit.append(lit_c2_hirschmann16)
    # # C, very reduced
    # pressure        = np.linspace(prs_min, prs_max, 100)
    # solubility      = np.zeros(len(pressure))
    # alpha           = 0.22                # ppm wt / MPa
    # beta            = 1.
    # alpha           = alpha * 1e-6        # --> ppm wt / Pa
    # beta            = beta
    # for i in range(0, 100):
    #     solubility[i] = henrys_law(pressure[i], alpha, beta)
    # lit_c3_hirschmann16, = plt.plot(pressure, solubility, color=vol_colors["qgreen_light"], ls=':', lw=2, label=r'$C$, very reduced: Hirschmann 2016: $\alpha$=%1.3e' % tuple([alpha]))
    # handles_lit.append(lit_c3_hirschmann16)

if show_literature_fits == True:

    #####################
    ######### CO ########
    #####################

    ######### Yoshioka et al. 2019 – CO", from Keppler & Golabek 2019
    pressure        = np.linspace(prs_min, prs_max, 100)
    solubility      = np.zeros(len(pressure))
    alpha   = 0.016          # ppm wt / bar
    alpha   = alpha * 1e-5   # -> ppm / Pa
    beta    = 1.
    for i in range(0, 100):
        solubility[i] = henrys_law(pressure[i], alpha, beta)
    
    if show_fit_coefficients == True:
        lit_co_yoshioka19, = plt.plot(pressure, solubility, color=vol_colors["CO"][col_idx], ls='--', lw=lw, label=vol_latex["CO"]+r', Yoshioka+ 2019 (%1.3e, %1.1f)' % tuple([alpha, beta]))
    else:
        lit_co_yoshioka19, = plt.plot(pressure, solubility, color=vol_colors["CO"][col_idx], ls='--', lw=lw, label=vol_latex["CO"])
        
    handles_fits.append(lit_co_yoshioka19)

    if show_S == True:

        #####################
        ######### S #########
        #####################

        # Hirschmann 16
        # S, all fugacities
        pressure    = np.linspace(prs_min, prs_max, 100)
        solubility  = np.zeros(len(pressure))
        alpha       = 5000                          # ppm wt / MPa
        alpha       = alpha * 1e-6                  # -> ppm / Pa
        beta        = 1.
        for i in range(0, 100):
            solubility[i] = henrys_law(pressure[i], alpha, beta)
        if show_fit_coefficients == True:
            lit_s_hirschmann16, = plt.plot(pressure, solubility, color=vol_colors["S"][col_idx], ls='--', lw=lw, label=vol_latex["S"]+r', Hirschmann 2016 (%1.3e, 1.0)' % tuple([alpha]))
        else:
            lit_s_hirschmann16, = plt.plot(pressure, solubility, color=vol_colors["S"][col_idx], ls='--', lw=lw, label=vol_latex["S"])
        handles_fits.append(lit_s_hirschmann16)

######################## PLOT SETTINGS ########################

### Axes settings

plt.xscale("log")
plt.yscale("log")

plt.xlim((prs_min, prs_max))
plt.ylim((X_min_ppm, X_max_ppm))

ax.tick_params(axis='both', which='major', labelsize=16)
ax.tick_params(axis='both', which='minor', labelsize=16)

# Set xtick labels to bar
xticks      = [1e3, 1e5, 1e7, 1e9, 1e11]
ax.set_xticks(xticks)
xticklabels = [ "0.01", "1", "100", r"$10^{4}$", r"$10^{6}$" ]
ax.set_xticklabels(xticklabels)

plt.xlabel('Gas phase pressure, $p^{\mathrm{i}}_{\mathrm{vapor}}$ (bar)', fontsize=19)
plt.ylabel('Volatile solubility, $X_{\mathrm{magma}}^{\mathrm{i}}$ (ppmw)', fontsize=19)

## Legends

# ## Data legend individual
# legend_data = plt.legend(handles=handles_data, fontsize=13, ncol=2, loc=2)
# ax = plt.gca().add_artist(legend_data) # loc=[0.70, 0.01]

## Data legend generic
data_generic_handles_title, = plt.plot([0], [0], '^', markersize=ms, color='k', alpha=0.99, label=r'Experimental data')
handles_generic_data = [data_generic_handles_title] 
legend_generic_data = plt.legend(handles=handles_generic_data, fontsize=17, ncol=2, loc=2)
ax = plt.gca().add_artist(legend_generic_data) # loc=[0.70, 0.01]


## Fits legend individual
if show_literature_fits == True: fsl=15
else: fsl=17
legend_fits = plt.legend(handles=handles_fits, ncol=2, loc=4, fontsize=fsl, title='Solubility fits') 
plt.setp(legend_fits.get_title(),fontsize=17)
ax = plt.gca().add_artist(legend_fits)

if show_literature_comparison == True:
    legend_literature = plt.legend(handles=handles_lit, fontsize=11, ncol=2, loc=[0.01, 1.01])
    for text in legend_literature.get_texts():
        text.set_color(vol_colors["qgray_dark"])
    ax = plt.gca().add_artist(legend_literature)

### Output figure

plt.savefig("../figures/fig_2.pdf", bbox_inches="tight")
