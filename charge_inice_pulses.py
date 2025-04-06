## Script to read the charge of in-ice pulses.

import matplotlib.pyplot as plt
import numpy as np
import math
from icecube import icetray, dataclasses, clsim, photonics_service, simclasses
from icecube import dataio
from glob import glob

files = glob("/home/user/inice/src/Gen2-Scripts/python/segments/Helium_trigger_PDOM_1/helium*_trigger.i3.gz")
for filename in files:
        file = dataio.I3File(filename)
        print(filename)
        while file.more():
                frame = file.pop_frame()
                if 'I3RecoPulseSeriesMapGen2' in frame:
                        I3RecoPulseSeriesMapGen2 = frame['I3RecoPulseSeriesMapGen2']
                        chargesum = sum([sum([y.charge for y in x]) for x in I3RecoPulseSeriesMapGen2.values()])
                        print(chargesum)
