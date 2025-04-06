## script to calculate angular resolution

import matplotlib.pyplot as plt
import numpy as np
import math
from math import erf
from scipy import special
from scipy.special import erf

from icecube import icetray, dataclasses, clsim, photonics_service, simclasses
from icecube.icetray import I3Units
from icecube import dataio
from glob import glob

##
##np.loadtxt("proton_hits_beta.txt", unpack = True)
SIGMA = special.erf(1/np.sqrt(2)) # this is just the exact value of the 0.68....
angular_resolution = np.percentile((np.loadtxt("Iron_opening_angle.txt", unpack = True)), 100*SIGMA)
	
print(angular_resolution, "Iron_resolution_contained_tracks")


