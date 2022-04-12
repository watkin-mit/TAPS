Tool for Air Pollution Scenarios (TAPS)
Version 1.0
Last updated: 08 April 2022

Contact: Will Atkinson (watkin@mit.edu)
Alternative Contacts: Noelle E. Selin (selin@mit.edu), Sebastian Eastham (seastham@mit.edu)
MIT License: https://github.com/watkin-mit/TAPS/blob/main/LICENSE

Welcome! This is the repository for TAPS v1.0, as submitted to Geoscientific Model Development (GMD). You can report bugs, suggest features or view the current version on GitHub: https://github.com/watkin-mit/TAPS/

Below is a summary of the main files. See the user manual at https://github.com/watkin-mit/TAPS/wiki for a full description of dependencies, methods, and input/output files. 

~~~~~~~~~~~~~~~~~~~~~~

~~~ Main directory: scale_emissions.py

Combines emissions inventories, activity scenarios, and intensity scenarios (in input_files) generate emissions scaling scenarios (stored in scaling_output).

~~~ Main directory: analyze_emissions.py

Uses the scale_emissions.py exports to calculate emissions scenarios and create the figures submitted to GMD.

~~~ Main directory: output_for_CTM.py

Uses the scale_emissions.py exports to create gridded netCDF scaling files (stored in netcdf_output) for global chemical transport models (CTMs). 

~~~ Main directory: TAPS_install_packages.py

Installs (or gives instructions for) all packages used, via dependencies.txt.

~~~ Subdirectory: input_files

Inputs required for running the main scripts. The necessary inventory and GAINS files are available at separate links (as described in https://github.com/watkin-mit/TAPS/wiki#2-external-data-sources). A subfolder (regional_mapping) includes files to edit grid-to-region mappings if needed for future applications.

~~~ Subdirectory: scaling_output

Outputs mainly from scale_emissions.py, including inventory emissions (em_CEDS.csv, em_GFED.csv), scaling trends (CEDS_scaling.csv, GFED_scaling.csv), and other exports (including data for Table A1).

~~~ Subdirectory: netcdf_output

Example folder to hold the files created in output_for_CTM.py. 

~~~ Subdirectory: Figures

Figure outputs from analyze_emissions.py, as submitted to GMD.







