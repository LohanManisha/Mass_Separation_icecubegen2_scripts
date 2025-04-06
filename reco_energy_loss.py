## To read total reconstructed energy loss from reconstructed I3 file

import matplotlib.pyplot as plt
import numpy as np
import math
from icecube import icetray, dataclasses, clsim, photonics_service, simclasses
from icecube import dataio
from glob import glob

files = glob("/home/user/inice/src/millipede/resources/examples/MPEFit_testing/*.i3.gz")
for filename in files:
        file = dataio.I3File(filename)
        print(filename)
        while file.more():
                frame = file.pop_frame()
                if 'MillipedeHighEnergy' in frame:
                        MillipedeHighEnergy = frame['MillipedeHighEnergy']
                        dEdX = sum([x.energy for x in frame['MillipedeHighEnergy']])
                        print(dEdX)
			#print(dEdx/10)
