#!/usr/bin/env python

"""
Example millipede muon energy loss fits.  The first fits the
loss pattern as stochastic losses (e.g. from a high-energy
muon), and the second as continuous losses.
input: offline reconstructed .i3 file(s)
"""

from I3Tray import *
import sys
import numpy as np
from icecube.STTools.seededRT.configuration_services import I3DOMLinkSeededRTConfigurationService

from icecube import icetray, dataio, dataclasses, photonics_service, finiteReco, lilliput, gulliver, rpdf, spline_reco
import icecube.lilliput.segments
from icecube.dataclasses import I3RecoPulseSeriesMap, I3TimeWindow
load('millipede')
from icecube import linefit
#from icecube import MPEFit
from icecube import phys_services



import argparse
parser = argparse.ArgumentParser()
GCD = '/home/manisha/t3store3/inice/src/millipede/resources/examples/IceCubeHEX_Sunflower_240m_v4.0beta_ExtendedDepthRange.GCDScintAnt-shift.i3.gz'
parser.add_argument('--gcd', type=str, metavar='FILE', default=GCD, help='GCD file')
parser.add_argument('--output', type=str, metavar='FILE', default='/home/manisha/t3store3/inice/src/millipede/resources/examples/hit_cleaning_COG_hits_millipede_Iron_4_Resampling_SplineMPE_400m_R_2000ns_T_lgE_16.8_sin2_0.0/Iron_199.i3.gz', help='Output')


parser.add_argument('result', type=str, nargs='+', metavar='FILE', help='Simulation Files')
args = parser.parse_args()


table_base = os.path.expandvars('/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/cascade_single_spice_bfr-v2_flat_z20_a5.%s.fits')
cascade_service = photonics_service.I3PhotoSplineService(table_base % 'abs', table_base % 'prob', timingSigma=0)

tray = I3Tray()


timingSplinePath = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_prob_z20a10_V2.fits'
amplitudeSplinePath = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_abs_z20a10_V2.fits'
stochTimingSplinePath = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfHighEStoch_mie_prob_z20a10.fits'
stochAmplitudeSplinePath = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfHighEStoch_mie_abs_z20a10.fits'

pulses = "I3RecoPulseSeriesMapGen2"

EnEstis = ["SplineMPETruncatedEnergy_SPICEMie_AllDOMS_Muon",
        "SplineMPETruncatedEnergy_SPICEMie_DOMS_Muon",
        "SplineMPETruncatedEnergy_SPICEMie_AllBINS_Muon",
        "SplineMPETruncatedEnergy_SPICEMie_BINS_Muon",
        "SplineMPETruncatedEnergy_SPICEMie_ORIG_Muon"]


tray.AddModule('I3Reader', 'reader', FilenameList = [args.gcd] + args.result)
tray.Add("I3NullSplitter")

def readout(frame, name_of_pulsemap = 'I3RecoPulseSeriesMapGen2'):
    timing = []
    recopulsemap = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, name_of_pulsemap)
    for pmtkey, list_of_pulses in recopulsemap.items():
         for pulse in list_of_pulses:
              timing.append(pulse.time)
    timing = np.array(timing)
    frame[name_of_pulsemap + 'TimeRange'] = I3TimeWindow(timing.min()-100, timing.max()+100)
    return
tray.Add(readout, Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics])

srt_config = I3DOMLinkSeededRTConfigurationService(
    ic_ic_RTRadius              = 400.0*I3Units.m,
    ic_ic_RTTime                = 2000.0*I3Units.ns,
    useDustlayerCorrection      = True,
    allowSelfCoincidence        = True,
    ic_strings               = ["1-35", "36", "36", "37-78", "1000-1120"],
    ic_oms                   = ["1-60", "1-19", "25-39", "1-60", "1-80"])    # Should this be true?


tray.AddModule("I3SeededRTCleaning_RecoPulseMask_Module", "SeededRT",
    InputHitSeriesMapName  = "I3RecoPulseSeriesMapGen2",
    OutputHitSeriesMapName = "CleanedPulses",
    STConfigService        = srt_config,
    SeedProcedure          = "HLCCOGSTHits",
    MaxNIterations         = -1,
    Streams                = [icetray.I3Frame.Physics],
    #HLCcoreThreshold       = 2
    )

tray.Add(linefit.simple, "linefit",
         inputResponse = "CleanedPulses", # or whatever pulses you have
         fitName = 'linefit')

tray.Add(lilliput.segments.I3SinglePandelFitter, "SPEFitSingle",
        pulses = "CleanedPulses",
        seeds = ['linefit',],
        fitname = "SPEFitPulses")

tray.Add(lilliput.segments.I3IterativePandelFitter, "SPEFit4",
        n_iterations=4,
        pulses = "CleanedPulses",
        seeds = ['SPEFitPulses',],
        fitname = "SPEFit4Out")

tray.AddSegment(spline_reco.SplineMPE, "SplineMPEmax",
        configuration="max", PulsesName="CleanedPulses", TrackSeedList=["SPEFit4Out"],
        fitname ="SplineMPEOut",
        BareMuTimingSpline=timingSplinePath,
        BareMuAmplitudeSpline=amplitudeSplinePath,
        StochTimingSpline=stochTimingSplinePath,
        StochAmplitudeSpline=stochAmplitudeSplinePath,
        EnergyEstimators=EnEstis)


def readout(frame, name_of_pulsemap = 'CleanedPulses'):
    timing = []
    cleanpulsemap = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, name_of_pulsemap)
    for pmtkey, list_of_pulses in cleanpulsemap.items():
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
    MuonSpacing=0, ShowerSpacing=10, SeedTrack='SplineMPEOut',
    Output='MillipedeHighEnergy', Pulses='CleanedPulses')



tray.AddModule('I3Writer', 'writer', 
    filename=args.output,
    Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics]
    )

##print(sys.argv[1])
tray.Execute()



