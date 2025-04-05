#!/usr/bin/env python

#---------------------------------------------------------------------
#---------------------------------------------------------------------
#    This is JUST an example how to run 3 STEP reconstruction  
#    with scintillators. Only ESTIMATION on scint pulses is used !!!     
#---------------------------------------------------------------------
#---------------------------------------------------------------------

from I3Tray import *
from icecube.icetray import I3Units
from icecube import icetray, dataclasses, dataio, recclasses
from icecube import rock_bottom, toprec, phys_services
from icecube import gulliver, gulliver_modules, lilliput
from icecube.rock_bottom.modules import *
from icecube.rock_bottom import *
from icecube.icetray.i3logging import log_info
import math
import numpy as np

import argparse
parser = argparse.ArgumentParser()
GCD = '/home/manisha/t3store3/IceCubeHEX_Sunflower_240m_v4.0beta_ExtendedDepthRange.GCDScintAnt-shift.i3.gz'
parser.add_argument('--gcd', type=str, metavar='FILE', default=GCD, help='GCD file')
parser.add_argument('--output', type=str, metavar='FILE', default='/home/manisha/t3store3/surfacearray/src/rock_bottom/resources/examples/Rockbottom_iron_4_times_lgE_17.0_sin2_0.0/iron_169.i3.gz', help='Output')
parser.add_argument('result', type=str, nargs='+', metavar='FILE', help='Simulation Files') 
args = parser.parse_args()
tray = I3Tray()

#---------------------------------------------------------------------
#                 Reader
#---------------------------------------------------------------------
tray.Add('I3Reader', 'TheReader',
    FilenameList = [args.gcd] +  args.result
    )
tray.Add('I3NullSplitter', 'TheSplitter',
    SubEventStreamName = 'NullSplit'
    )
    
#---------------------------------------------------------------------
#         Cut on 0.5VEM and
#---------------------------------------------------------------------
def Scint_VEM( frame, Pulses = '' ):
  ps = dataclasses.I3ScintRecoPulseSeriesMap()
  if Pulses in frame:
    charges_scint = frame[Pulses]
    strings   = np.asarray([i.station for i in charges_scint.keys()])
    oms       = np.asarray([ii.panel for ii in charges_scint.keys()])
    all_scint = [p for dom, pul in charges_scint for p in pul]
    pulse     = np.asarray([pi.charge for pi in all_scint])
    times     = np.asarray([pu.time for pu in all_scint])
    for u in range(0,(len(strings))):
      pulses = dataclasses.I3RecoPulse()
      if pulse[u]>=0.5:
##        if (len(pulse)>=5):	  
        charges       = pulse[u]
        pulses.time   = times[u]
        pulses.charge = charges
        ps[dataclasses.ScintKey(int(strings[u]),int(oms[u]))] = dataclasses.I3RecoPulseSeries()
        ps[dataclasses.ScintKey(int(strings[u]),int(oms[u]))].append(pulses)
  frame.Put(Pulses + '_VEM_05cut', ps)

tray.AddModule(Scint_VEM, Pulses='SiPMRecoPulses', streams = [icetray.I3Frame.DAQ])

#---------------------------------------------------------------------
#           Scintillator reconstruction 
#---------------------------------------------------------------------
ScintPulses = 'SiPMRecoPulses_VEM_05cut'
loglog      = rock_bottom.ldf.LogLog()
TimeFcn     = rock_bottom.showerfront.TimeGauss()


tray.AddModule('I3RbSeeds', "ScintSeed1",
  InputName     =  ScintPulses,
  GeoName       = 'I3ScintGeometry',
  ModelName     = 'ScintSignalModel',
  MapName       = 'ParameterSeed',
  OutputName    = 'ShowerSimpleSeed',
  FixSigmaPlane = True,
  Multiplicity  = True
  )

#---------------------------------------------------------------------
#          Minimizer preparation
#---------------------------------------------------------------------
tray.AddService('I3GulliverMinuitFactory', 'Minuit',
  MinuitPrintLevel = -2,
  FlatnessCheck    = True,
  Algorithm        = 'SIMPLEX',
  MaxIterations    = 2500,
  MinuitStrategy   = 2,
  Tolerance        = 0.001
  )  

#---------------------------------------------------------------------
#            1 STEP
#---------------------------------------------------------------------

tray.AddService('ScintSignalModel', 'ScintSignalModel',
  LDF = loglog,
  ParameterNames   = RbList([RbPair("slopeLdf", 0)]),
  ParameterValues  = [            2.7,        ],
  BoundNames       = RbList([RbPair("slopeLdf", 0), RbPair("lgSref", 0)]),
  BoundValues      = [         (1.2 , 3.6),      (-2.8, 2.8)],
  StepNames        = RbList([RbPair("slopeLdf", 0), RbPair("lgSref", 0) ]),
  StepValues       = [             1.,                1.,   ], 
  )

tray.AddService('I3RbLDFLikelihoodFactory', 'ScintLikelihood_LDF',
  DetectorType = rock_bottom.DetectorTypes.Scint,
  Model        ='ScintSignalModel',
  Pulses1      = ScintPulses,
  UseSilent    = True,
  MinSignal    = 0.5
  )

tray.AddService("I3MultiSurfaceSeedServiceFactory", "ScintSeeds_Step1",
  FirstGuesses    = ['ShowerSimpleSeed'],
  SignalModels    = ['ScintSignalModel'],
  SeedsMap        = 'ParameterSeed',
  ParticleXStep   = 20.,
  ParticleYStep   = 20.,
  ParticleXRelBounds = [-250., 250.],
  ParticleYRelBounds = [-250., 250.] 
  )


tray.AddService("I3MultiSurfaceParametrizationFactory","ParamsStep1",
  SeedService = 'ScintSeeds_Step1' )

tray.AddModule("I3SimpleFitter", "ScintReconstructionStep1",
  SeedService     = "ScintSeeds_Step1",
  Parametrization = "ParamsStep1",
  LogLikelihood   = "ScintLikelihood_LDF",
  OutputName      = "ScintReconstructionStep1a",
  Minimizer       = "Minuit"
  )

#---------------------------------------------------------------------
#            2 STEP
#---------------------------------------------------------------------

tray.AddService('ScintSignalModel', 'ScintSignalModel2',
  LDF = loglog,
  BoundNames       = RbList([RbPair("slopeLdf", 0), RbPair("lgSref", 0)]),
  BoundValues      = [         (0.6 , 6.),      (-2.8, 2.8)],
  StepNames        = RbList([RbPair("slopeLdf", 0), RbPair("lgSref", 0) ]),
  StepValues       = [             0.25,                0.1,   ], 
  )

tray.AddService('I3RbLDFLikelihoodFactory', 'ScintLikelihood_LDF2',
  DetectorType = rock_bottom.DetectorTypes.Scint,
  Model        ='ScintSignalModel2',
  Pulses1      = ScintPulses,
  UseSilent    = True,
  MinSignal    = 0.5
  )

tray.AddService('FrontScintModel', 'FrontScintModel2',
  TimeFcn     = TimeFcn,
  BoundNames  = RbList([RbPair("parTime", 2) ]),
  BoundValues = [ (5.e-6, 100.e-4) ], 
  StepNames   = RbList([RbPair("parTime", 2) ]),
  StepValues  = [2.e-4], 
  )

tray.AddService("I3RbTimingLikelihoodFactory", "ScintLikelihood_Curve2",
	DetectorType = rock_bottom.DetectorTypes.Scint,
	Model			   = "FrontScintModel2",
	Pulses1			 = ScintPulses,
	MinSignal 	 = 0.5   
	)

tray.AddService("I3MultiSurfaceSeedServiceFactory", "ScintSeeds_Step2",
  FirstGuesses    = ['ScintReconstructionStep1a'],
  SignalModels    = ['ScintSignalModel2', 'FrontScintModel2'],
  SeedsMap        = 'ScintReconstructionStep1aParams',
  ParticleTStep   = 100.,
  ParticleZenStep  = 0.01,
  ParticleAziStep  = 0.02)
  

tray.AddService("I3MultiSurfaceParametrizationFactory","ParamsStep2",
  SeedService = 'ScintSeeds_Step2' )


tray.AddService("I3EventLogLikelihoodCombinerFactory", "ScintLikelihood_Combo2",
	InputLogLikelihoods = ["ScintLikelihood_LDF2","ScintLikelihood_Curve2"],
	Multiplicity 		= "Sum",
	RelativeWeights 	= [1.,1.]
	)


tray.AddModule("I3SimpleFitter", "ScintReconstructionStep2",
  SeedService     = "ScintSeeds_Step2",
  Parametrization = "ParamsStep2",
  LogLikelihood   = "ScintLikelihood_Combo2",
  OutputName      = "ScintReconstructionStep2a",
  Minimizer       = "Minuit"
  )

#---------------------------------------------------------------------
#            3rd STEP
#---------------------------------------------------------------------

tray.AddService('ScintSignalModel', 'ScintSignalModel3',
  LDF = loglog,
  BoundNames       = RbList([RbPair("slopeLdf", 0), RbPair("lgSref", 0)]),
  BoundValues      = [         (0.6 , 6.),      (-2.8, 2.8)],
  StepNames        = RbList([RbPair("slopeLdf", 0), RbPair("lgSref", 0) ]),
  StepValues       = [           0.5,                0.1,   ], 
  )

tray.AddService('I3RbLDFLikelihoodFactory', 'ScintLikelihood_LDF3',
  DetectorType = rock_bottom.DetectorTypes.Scint,
  Model        ='ScintSignalModel3',
  Pulses1      = ScintPulses,
  UseSilent    = True,
  MinSignal    = 0.5
  )

tray.AddService("I3MultiSurfaceSeedServiceFactory", "ScintSeeds_Step3",
  FirstGuesses    = ['ScintReconstructionStep2a'],
  SignalModels    = ['ScintSignalModel3'],
  SeedsMap        = 'ScintReconstructionStep2aParams',
  ParticleXStep   = 20.,
  ParticleYStep   = 20.,
  ParticleXRelBounds = [-100., 100.],
  ParticleYRelBounds = [-100., 100.] )


tray.AddService("I3MultiSurfaceParametrizationFactory","ParamsStep3",
  SeedService = 'ScintSeeds_Step3' )


tray.AddModule("I3SimpleFitter", "ScintReconstructionStep3",
  SeedService     = "ScintSeeds_Step3",
  Parametrization = "ParamsStep3",
  LogLikelihood   = "ScintLikelihood_LDF3",
  OutputName      = "ScintReconstructionStep3a",
  Minimizer       = "Minuit"
  )

#---------------------------------------------------------------------
#           Writing
#---------------------------------------------------------------------

tray.AddModule("I3Writer","i3writer",
  Filename = args.output,
  Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics]
  )
  
tray.Execute()
tray.Finish()
