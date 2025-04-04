## To read zenith angle in radian and number of muons from I3MCTree


import matplotlib.pyplot as plt
import numpy as np
import math
from icecube import icetray, dataclasses, clsim, photonics_service, simclasses
from icecube import dataio
from glob import glob
ZenithAngle = []
NumberOFMuonsPLus = []
files = glob("/user/iron/lgE_18.5/sin2_0.0/000*.i3.gz")
for filename in files:
	file = dataio.I3File(filename)
	print(filename)
	Tray_Info = file.pop_frame()
	Shower = file.pop_frame()
	I3MCTreeIT = Shower["I3MCTreeIT"]
	MCPrimary = Shower["MCPrimary"]
	ZenithAngle.append(MCPrimary.dir.zenith*180/math.pi)
	totalMu = 0
	for particle in I3MCTreeIT:
		if "Mu" in particle.type_string: 
			totalMu+=1
	NumberOFMuonsPLus.append(totalMu)
	

print(NumberOFMuonsPLus, ZenithAngle)
lgMu = np.log10(NumberOFMuonsPLus)
plt.scatter(ZenithAngle, lgMu)
plt.xlabel("ZenithAngle")
plt.ylabel("Number of Muons")
plt.savefig("MuonsZenith.png")

