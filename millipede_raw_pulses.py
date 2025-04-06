#!/usr/bin/env python

## To reconstruct energy loss of in-ice muons using raw pulses

from I3Tray import *
import sys
import numpy as np
from icecube import icetray, dataio, dataclasses, photonics_service, finiteReco, lilliput, gulliver, rpdf
import icecube.lilliput.segments
from icecube.dataclasses import I3RecoPulseSeriesMap, I3TimeWindow
load('millipede')
from icecube import linefit
#from icecube import MPEFit
from icecube import phys_services
#from filters_InIceSplit_2015 import MPEFit
##if len(sys.argv) < 3:
##	print('Usage: %s output.i3 input1.i3 [input2.i3] ...' % sys.argv[0])
##	sys.exit(1)
##files = sys.argv[2:]


import argparse
parser = argparse.ArgumentParser()
GCD = '/home/user/IceCubeHEX_Sunflower_240m_v4.0beta_ExtendedDepthRange.GCDScintAnt-shift.i3.gz'
parser.add_argument('--gcd', type=str, metavar='FILE', default=GCD, help='GCD file')
#parser.add_argument('--input', type=str, metavar='FILE', default='/home/user/inice/src/Gen2-Scripts/python/segments/Proton_trigger_PDOM_1/*.i3.gz', help='Input')
parser.add_argument('--output', type=str, metavar='FILE', default='/home/user/inice/src/millipede/resources/examples/Millipede_HLC_pulses/Proton_HLC_millipede_28_linefit.i3.gz', help='Output')


parser.add_argument('result', type=str, nargs='+', metavar='FILE', help='Simulation Files')
args = parser.parse_args()

#table_base = os.path.expandvars('/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/cascade_single_spice_bfr-v2_flat_z20_a5.%s.fits')
#muon_service = photonics_service.I3PhotoSplineService(table_base % 'abs', table_base % 'prob', timingSigma=0)
table_base = os.path.expandvars('/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/cascade_single_spice_bfr-v2_flat_z20_a5.%s.fits')
cascade_service = photonics_service.I3PhotoSplineService(table_base % 'abs', table_base % 'prob', timingSigma=0)

tray = I3Tray()
tray.AddModule('I3Reader', 'reader', FilenameList = [args.gcd] + args.result)
tray.Add("I3NullSplitter")

def hlc_cleaning(frame):
    frame['Gen2InIcePulsesHLC'] = dataclasses.I3RecoPulseSeriesMapMask(frame, 'I3RecoPulseSeriesMapGen2',lambda omkey, index, pulse: pulse.flags & pulse.PulseFlags.LC)

tray.Add(hlc_cleaning, Streams=[icetray.I3Frame.Physics,])

tray.Add(linefit.simple, "linefit",
         inputResponse = "Gen2InIcePulsesHLC", # or whatever pulses you have
         fitName = 'linefit')

#tray.Add(lilliput.segments.I3SinglePandelFitter, "SPEFitSingle",
#	pulses = "I3RecoPulseSeriesMapGen2",
#	seeds = ['linefit',],
#	fitname = "SPEFitPulses")

#tray.Add(lilliput.segments.I3IterativePandelFitter, "SPEFit4",
#	n_iterations=4,
#        pulses = "I3RecoPulseSeriesMapGen2",
#        seeds = ['SPEFitPulses',],
#        fitname = "SPEFit4Out")

#tray.Add(lilliput.segments.I3SinglePandelFitter, "MPEFit",
#        pulses = "I3RecoPulseSeriesMapGen2",
#        seeds = ['SPEFit4Out',],
#	fitname = "MPEFitPulses")


def readout(frame, name_of_pulsemap = 'Gen2InIcePulsesHLC'):
    timing = []
    recopulsemap = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, name_of_pulsemap)
    for pmtkey, list_of_pulses in recopulsemap.items():
         for pulse in list_of_pulses:
              timing.append(pulse.time)
    timing = np.array(timing)
    frame[name_of_pulsemap + 'TimeRange'] = I3TimeWindow(timing.min()-100, timing.max()+100)
    return
tray.Add(readout, Streams=[icetray.I3Frame.Physics])


tray.AddModule('MuMillipede', 'millipede_highenergy',
    CascadePhotonicsService=cascade_service,
    Boundary=3000,
    PhotonsPerBin=15, MuonRegularization=0, ShowerRegularization=0,
    MuonSpacing=0, ShowerSpacing=10, SeedTrack='linefit',
    Output='MillipedeHighEnergy', Pulses='Gen2InIcePulsesHLC')

#tray.AddModule('MuMillipede', 'millipede_lowenergy',
#    MuonPhotonicsService=muon_service, CascadePhotonicsService=cascade_service,
#    PhotonsPerBin=10, MuonRegularization=2, ShowerRegularization=0,
#    MuonSpacing=15, ShowerSpacing=0, SeedTrack='linefit',
#    Output='MillipedeLowEnergy', Pulses='I3RecoPulseSeriesMapGen2')

tray.AddModule('I3Writer', 'writer', 
    filename=args.output,
    Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics]
    )

##print(sys.argv[1])
tray.Execute()



