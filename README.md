# Mass_Separation_icecubegen2_scripts
## This repository contains all the scripts related to the analysis: "Estimating the sensitivity of IceCube-Gen2 to Cosmic Ray mass separation"
Datafiles (basically CORSIKA showers) used are located at: /data/sim/IceCubeUpgrade/CosmicRay/Radio/coreas/data/continuous/star-pattern

In this analysis, the combined response of scintillators array and optical array is studied towards cosmic ray mass separation. 
Surface array simulations are done using surface array meta project.
Variables studied corresponding to scintillators array: beta (slope of LDF fitting well the scintillator signal) and number of CORSIKA muons at surface level. Here, beta is a reconstructed variable, reconstructed using "Rockbottom" project.
Reconstructed beta is read using the script: beta_reader_scints_New.py
CORSIKA Muons at surface are read using the script: 
