# spectra-config: Generate configuration files for MARBL-SPECTRA

MARBL-SPECTRA is an extension of the base marine ecosystem in the Marine Biogeochemistry Library (MARBL) to include 15 plankton: 9 phytoplankton and 6 zooplankton with broad size and trait diversity. These codes generate the runtime-configurable namelist settings for the SPECTRA model, but are written so that one could generate a starting set of parameters for any number of plankton greater than 4 phytoplankton and 2 zooplankton. The codes here also generate initial conditions for the model.

### To generate the initial conditions:
```
cd initial_conditions
./gen_input_dataset.py
```
This generates a netcdf file that will need to be called in `user_nl_pop`:
```
init_ecosys_init_file = 'YOUR_PATH/ecosys_jan_IC_gx1v6_9p6z_20180124_3impCalc.nc'
```

### To generate the namelist settings:
```
./gen_input_data.R
./gen_marbl_input_data.R
```
This then generates a namelist file called `user_nl_marbl`, which needs to be placed in your `CASEROOT` in order for these settings to take effect.

You can modify parameters in `input_parms.R`, and generate plots that visualizes the parameters by running `analyze_model_setup.R`.

The `size_structured_ecosystem.R` script is a wrapper script that takes one input, a casename, and creates a new copy of the four main R scripts in a new case folder and runs the scripts. This can assist in setting up new parameter cases.

### Reference:
Gabriela Negrete-García, Jessica Y. Luo,  Matthew C. Long,  Keith Lindsay,  Michael Levy, Andrew D. Barton. Plankton energy flows using a global size-structured and trait-based model. Preprint. doi: https://doi.org/10.1101/2022.02.01.478546

### Contact:
Jessica Luo: jessica.luo [at] noaa.gov<br>
Gaby Negrete-García: g1negret [at] ucsd.edu

