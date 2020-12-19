# Plotting scripts + simulation data (Lichtenberg et al., 2021, JGRP)

This repository is part of the combined data + script container at https://osf.io/m4jh7. 
For questions, help and notifications surrounding this repository please contact [Tim Lichtenberg](https://timlichtenberg.net/).

# References

When you make use of material from this repository, please acknowledge:

[![DOI](https://zenodo.org/badge/321733308.svg)](https://zenodo.org/badge/latestdoi/321733308)

> Lichtenberg, T., Bower, D. J., Hammond, M., Boukrouche, R., Sanan, P., Tsai, S.-M., Pierrehumbert R. T. (2021) Vertically resolved magma ocean–protoatmosphere evolution: H2, H2O, CO2, CH4, CO, O2, and N2 as primary absorbers. Journal of Geophysical Research: Planets.

# Repository structure

    .
    ├── /data                           # Simulation data files, top-level directory
    │   ├── /fig2_solubilities          # Data for Fig. 2
    │   ├── /fig4_radiation_limits      # Data for Fig. 4
    │   ├── /fig8_cff                   # Data for Fig. 8
    │   └── /int_atm                    # Data for all other figs, **see note 'Simulation data' below**
    ├── /figures                        # Figures
    ├── /scripts                        # Plotting scripts, help files, and illustration source files
    ├── LICENSE                         # License agreement
    └── README.md                       # This README file

# Simulation data

In order to execute the plotting scripts, please download *int_atm.zip* from https://doi.org/10.17605/OSF.IO/M4JH7 and unzip it in */data/*. Plotting scripts for Figs. 3, 5, 6, 7, and 10 only work when this simulation data set is present in the data folder.

# Software requirements

Python 3.x, see the individual plotting scripts and */scripts/modules_plot.py* for the required packages.
