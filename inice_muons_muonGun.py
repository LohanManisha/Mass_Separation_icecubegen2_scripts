## Script to estimate in-ice muons using MuonGun

import matplotlib.pyplot as plt
import numpy as np
from icecube import icetray, dataclasses, clsim, photonics_service, simclasses, dataio
from icecube import MuonGun
from glob import glob
gcdfile = "/home/user/inice/src/Gen2-Scripts/python/segments/MuonGun/IceCubeHEX_Sunflower_240m_v4.0beta_ExtendedDepthRange.GCDScintAnt-shift.i3.gz"
gcd = dataio.I3File(gcdfile)
frame = gcd.pop_frame()
if "I3Geometry" in frame.keys(): frame = gcd.pop_frame()
omgeomap = frame['I3Geometry'].omgeo

positions = [omgeo.position for omgeo in omgeomap.values()] 
gen2_inice_surface = MuonGun.ExtrudedPolygon(positions, 300)

files = glob("/home/user/inice/src/Gen2-Scripts/python/segments/helium_clsim/helium_clsim_9L_h5.i3.gz")
for filename in files:
	file = dataio.I3File(filename)
	while file.more():
		frame = file.pop_frame()
#	print(frame)
#	Shower = file.pop_frame()
		if 'I3MCTree' in frame:
			I3MCTree = frame["I3MCTree"]
			MMCTrackList = frame["MMCTrackList"]	
			nMuonsEntering, nMuonsExiting = 0, 0
			for track in MuonGun.Track.harvest(frame['I3MCTree'], frame['MMCTrackList']):
    # Find distance to entrance and exit from sampling volume
				intersections = gen2_inice_surface.intersection(track.pos, track.dir)

    # Get the corresponding energies at entrance/exit in volume
				e0, e1 = track.get_energy(intersections.first), track.get_energy(intersections.second)

				if e0 > 0: nMuonsEntering += 1
				if e1 > 0: nMuonsExiting += 1
				#print(e0)
			print("Found {} muons entering the detector and {} exiting".format(nMuonsEntering, nMuonsExiting))
#				#print(e0)


