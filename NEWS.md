# dev

# v.1.5.11.9000
  - restore functionality to restart a simulation from a saved state
  - if a config object is used to run a simulation, it is now required
    to define the gen3sis$general$config_name variable with the name
    of the subdirectory where the outputs of an individual simulation
    will be saved (if the config object is created from a config file,
    this variable will be automatically set from the file name)
  - fix bug in verify_config() where missing settings were not
    properly identified

# v.1.5.11 release 11.2023
  - fix comb phylogeny
  - color deficient, blind and B&W safe colours
  - speed-up of loop_ecology function
  - fix nexus phylo file format
  - fix accounting of extinctions times 
  - fix phylo checks if simulation did not end at t=0
  - fix package build notes citation and class comparisons

# v.1.4 release 10.2021
  - fix bracket compatibility with new R version
  - added new abundance plotting function 

# v.1.3 release 07.2021

# v.1.2 release	12.2020 

# v.1.1 release	08.2020
  
# v.1.0 release 06.2020
