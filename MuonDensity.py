import matplotlib.pyplot as plt
import numpy as np
from icecube import icetray, dataclasses, clsim, photonics_service, simclasses
from icecube import dataio
from glob import glob
energy = []
#NumberOFMuonsPLus = []
#radius = []
#BinEdges = []
Density = [] 
#files = glob("./*.i3.gz")
files = glob("/home/manisha/t3store3/surface16August/proton/lgE_16.*/sin2_0.0/000*.i3.gz")
for filename in files:
	file = dataio.I3File(filename)
	radius = []
	print(filename)
	Tray_Info = file.pop_frame()
	Shower = file.pop_frame()
	#MCPrimaryInfo = Shower["MCPrimaryInfo"]
	I3MCTreeIT = Shower["I3MCTreeIT"]
	MCPrimary = Shower["MCPrimary"]
	#energy.append(MCPrimary.energy)
#	NumberOFMuonsPLus.append(I3MCTreeIT.MuPlus)
	energy.append(MCPrimary.energy)
	for p in I3MCTreeIT:
#		print(p.type)
		if "Mu" in p.type_string:
		#if (p.type==5) or (p.type==6):
#	for i in radius:	
#		print(i)
			x = p.pos.x
			y = p.pos.y	
#		z = p.pos.z
			r = np.sqrt(x*x + y*y)
			radius.append(r)
#	np.histogram(radius, bins=np.arange(300, 1000, 20))
	counts, bin_edges = np.histogram(radius, bins=np.arange(300, 1000, 20))
#	print(radius)
#	BinEdges.append(bin_edges)
	area = np.pi*(bin_edges[1:]**2)
	density = np.cumsum(counts)/area
	Density.append(density)
	
Density = np.array(Density)
print(Density)
Number_of_Rows = Density.shape[1]
for row in range(Number_of_Rows):
#	plt.scatter(energy, Density[:,row])
	plt.plot(energy, Density[:,row], linestyle="", marker="")
#print(
	print(Density[:,row])
#lgE = mp.log10(energy)
#lgD = 
plt.xscale("log")
plt.yscale("log")
#plt.scatter(energy, Density)
#energy = np.log10(energy)
#plt.scatter(lgMu, lgE)
plt.xlabel("Energy")
plt.ylabel("Muon Density")
plt.savefig("MuonDensity24.png")

