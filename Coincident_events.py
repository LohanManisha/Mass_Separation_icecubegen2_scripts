import matplotlib.pyplot as plt
import numpy as np
from icecube import icetray, dataclasses, clsim, photonics_service, simclasses, dataio
from icecube import MuonGun
from glob import glob

from icecube.phys_services import I3ScaleCalculator
I3ScaleCalculator.IT_CUSTOM
#icecube._phys_services.IceTopConfig.IT_CUSTOM
I3ScaleCalculator.IC_CUSTOM
#icecube._phys_services.IceCubeConfig.IC_CUSTOM
#print(dir(I3ScaleCalculator))
# Grab the geometry

gcdfile = "/home/manisha/t3store3/inice/src/Gen2-Scripts/python/segments/IceCubeHEX_Sunflower_240m_v4.0beta_ExtendedDepthRange.GCDScintAnt-shift.i3.gz"
gcd = dataio.I3File(gcdfile)
frame = gcd.pop_frame()
if "I3Geometry" in frame.keys(): frame = gcd.pop_frame()

geo = frame["I3Geometry"]

#icgen2 = list(range(1000, 1120))
#itgen2 = list(range(1000, 1120))
icgen2 = [1117, 1014, 1027, 1040, 1052, 1064, 1075, 1085, 1095, 1104, 1113, 1099, 1108, 1115, 1102, 1111,1097, 1106, 1091, 1100, 1109, 1116, 1103, 1112, 1098, 1107, 1114, 1101, 1110, 1096, 1105, 1090, 1073, 1055, 1035, 1119, 1005, 1120]
#itgen2 = [1117, 14, 27, 40, 52, 64, 75, 85, 95, 104, 113, 99, 108, 115, 102, 111, 97, 106, 91, 100, 109, 116, 103, 112, 98, 107, 114, 101, 110, 96, 105, 90, 1073, 55, 35, 119, 5, 120]

#icgen2 = [14, 27, 40, 52, 64, 75, 85, 95, 104, 113, 99, 108, 115, 102, 111, 97, 106, 91, 100, 109, 116, 103, 112, 98, 107, 114, 101, 110, 96, 105, 90, 55, 35, 119, 5, 120]
#itgen2 = [14, 27, 40, 52, 64, 75, 85, 95, 104, 113, 99, 108, 115, 102, 111, 97, 106, 91, 100, 109, 116, 103, 112, 98, 107, 114, 101, 110, 96, 105, 90, 55, 35, 119, 5, 120]


#icgen2 = ["1000-1120"]
#itgen2 = ["1000-1120"]

calc = I3ScaleCalculator(geo, I3ScaleCalculator.IC_CUSTOM, I3ScaleCalculator.IT_CUSTOM, icgen2)

files = glob("/home/manisha/t3store3/inice/src/Gen2-Scripts/python/segments/CLSIM_Iron_Resample_4_times_lgE_17.2_sin2_0.0/Iron_*.i3.gz")
for filename in files:
	file = dataio.I3File(filename)
	print(filename)
	while file.more():
		frame = file.pop_frame()
		if 'I3MCTree' in frame:
#			MMCTrackList = frame["MMCTrackList"]	
#			track = frame["MMCTrackList"]
			track = dataclasses.get_most_energetic_muon(frame["I3MCTree"])	
			volume_containment1 = calc.scale_inice(track)
			#area_containment = calc.scale_icetop(track)
			print("I3MCTree", volume_containment1)
			#print(area_containment)

		if 'I3MCTreeIT' in frame:
#                       MMCTrackList = frame["MMCTrackList"]    
#                       track = frame["MMCTrackList"]
                        track = dataclasses.get_most_energetic_muon(frame["I3MCTreeIT"])
                        volume_containment2 = calc.scale_inice(track)
                        #area_containment = calc.scale_icetop(track)
                        print("I3MCTreeIT", volume_containment2)
                        #print(area_containment)








