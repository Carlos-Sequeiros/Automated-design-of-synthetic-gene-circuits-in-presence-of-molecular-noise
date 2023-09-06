# Automated-design-of-synthetic-gene-circuits-in-presence-of-molecular-noise
Example code for the article _Automated design of synthetic gene circuits in presence of molecular noise_
by Carlos Sequeiros, Carlos Vázquez, Julio R. Banga and Irene Otero-Muras
## Requirements
- The code has been tested on Windows 10 64 bits
- A CUDA capable GPU
- MATLAB optimization toolbox
- MATLAB parallel computing toolbox
- [MEIGO toolbox](http://gingproc.iim.csic.es/meigo.html)
## Instructions
- Install MEIGO in the current MATLAB session by running `install_MEIGO.m` in MEIGO folder
- Go inside the desired example's folder where `.m` files are located
- Load the MEIGO settings by opening `MEIGO_configuration.mat`
- Run the line `Results=MEIGO(problem,opts,'ESS');`. In the examples `Target_desing` or `bimodal_probabilities_stationary` load respectively a `Target.mat` file or `domains_sym_split.mat` and include the variable inside as last argument 
