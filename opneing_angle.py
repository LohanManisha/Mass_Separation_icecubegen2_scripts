import matplotlib.pyplot as plt
import numpy as np
import math
from icecube import icetray, dataclasses, clsim, photonics_service, simclasses
from icecube.icetray import I3Units
from icecube import dataio
from glob import glob

files = glob("/home/manisha/t3store3/inice/src/millipede/resources/examples/hit_cleaning_COG_hits_millipede_Proton_SplineMPE_400m_R_2000ns_T/Proton_millipede_*_SplineMPE.i3.gz")
for filename in files:
	file = dataio.I3File(filename)
	print(filename)
	while file.more():
		frame = file.pop_frame()
		if 'I3MCTree' in frame:
			I3MCTree = frame['I3MCTree']
			truth = dataclasses.get_most_energetic_primary(frame['I3MCTree'])
			#truth1 = dataclasses.get_most_energetic_muon(frame['I3MCTree'])
			#particle = frame['Fe56Nucleus'] 
			direction1 = (truth.dir)
			#direction3 = (truth1.dir)
#			print(direction1)
			#print(direction3)
		if 'SplineMPEOut' in frame:
			SplineMPEOut = frame['SplineMPEOut']
			direction2 = (SplineMPEOut.dir)
#	print(direction2)
#	print(math.acos(direction1*direction2))
	print(np.degrees(math.acos(direction1*direction2)))

