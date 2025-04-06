## script to read the number of muons and electrons from I3MCTree

import matplotlib.pyplot as plt
import numpy as np
from icecube import icetray, dataclasses, clsim, photonics_service, simclasses
from icecube import dataio
from glob import glob
NumberofElectrons = []
NumberOFMuons = []
files = glob("/user/iron/lgE_18.5/sin2_0.0/000*.i3.gz")
for filename in files:
	file = dataio.I3File(filename)
	print(filename)
	Tray_Info = file.pop_frame()
	Shower = file.pop_frame()
	I3MCTreeIT = Shower["I3MCTreeIT"]
	I3MCTreeIT = Shower["I3MCTreeIT"]
	totalMu = 0
	for particle in I3MCTreeIT:
		if "Mu" in particle.type_string: 
			totalMu+=1
	NumberOFMuons.append(totalMu)
	
	totalE = 0
	for particle in I3MCTreeIT:
		if "E" in particle.type_string:
			totalE+=1
	NumberofElectrons.append(totalE)


print(NumberOFMuons)
lgE = np.log10(NumberofElectrons)
lgMu = np.log10(NumberOFMuons)

plt.scatter(lgMu, lgE)
plt.xlabel("NumberOFMuons")
plt.ylabel("NumberofElectrons")
plt.savefig("MuonsElectrons.png")

