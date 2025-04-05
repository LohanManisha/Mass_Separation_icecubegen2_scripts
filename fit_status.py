import matplotlib.pyplot as plt
import numpy as np
import math
from icecube import icetray, dataclasses, clsim, photonics_service, simclasses
from icecube.icetray import I3Units
from icecube import dataio
from glob import glob

files = glob("/home/manisha/t3store3/inice/src/millipede/resources/examples/Millipede_Iron_SplineMPE/*.i3.gz")
for filename in files:
	file = dataio.I3File(filename)
	print(filename)
	while file.more():
		frame = file.pop_frame()
		if 'SplineMPEOut' in frame:
			SplineMPEOut = frame['SplineMPEOut']
			status = SplineMPEOut.fit_status
			print(status==0)

