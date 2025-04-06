## Script to read the number of DOM hits.

import matplotlib.pyplot as plt
import numpy as np
import math
from icecube import icetray, dataclasses, clsim, photonics_service, simclasses
from icecube import dataio
from glob import glob

files = glob("/home/user/inice/src/Gen2-Scripts/python/segments/Detector_SIM_Iron_Resample_4_times_lgE_16.9_sin2_0.0/Iron_Trigger_*.i3.gz")
for filename in files:
	file = dataio.I3File(filename)
	print(filename)	
	while file.more():
		frame = file.pop_frame()
		if 'I3RecoPulseSeriesMapGen2' in frame:
			I3RecoPulseSeriesMapGen2 = frame["I3RecoPulseSeriesMapGen2"]
			modules = []
			for OMKey, pulses in I3RecoPulseSeriesMapGen2:
				if any([p.flags & dataclasses.I3RecoPulse.PulseFlags.LC for p in pulses]):					modules.append((OMKey.string, OMKey.om))
				n_lc_doms = len(set(modules))

			print(n_lc_doms)



