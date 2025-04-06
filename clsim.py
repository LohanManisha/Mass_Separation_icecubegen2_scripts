#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/gen2-optical-sim/software/icetray/build
## script for photon propagation
import os
import sys,getopt,string
from os.path import expandvars
import subprocess
import threading
import time
import copy
import math
import numpy as np
from functools import reduce

from os.path import expandvars, exists, isdir, isfile
from icecube import icetray, dataclasses, clsim, photonics_service, simclasses
from icecube.mDOM_WOM_simulation.traysegments import mdom, degg
from icecube.mDOM_WOM_simulation.traysegments.mDOMMakeHitsFromPhotons import mDOMMakeHitsFromPhotons
from icecube.mDOM_WOM_simulation.traysegments.DEggMakeHitsFromPhotons import DEggMakeHitsFromPhotons
from icecube.simprod import ipmodule
from icecube.icetray import I3Units, logging
logging.set_level_for_unit('I3PhotonToMCPEConverterForDOMs',
                           logging.I3LogLevel.LOG_TRACE)
from icecube.dataclasses import I3OMGeo
from icecube.gen2_sim.utils import iceprod_wrapper


def construct_mDOM_wavelength_acceptance():
    """
    The mDOM wavelentgth acceptance is implemented separately for 
    glass, gel and the pmt. Provide a wrapper so that it can be 
    called with a GetValue() like the wavelength acceptances of 
    degg and pdom, which are CLSimFunctionsFromTable.
    """
    import numpy as np
    qe = mdom.GetQuantumEfficiency()
    gel = mdom.GetGelAbsorptionLength()
    glass = mdom.GetGlassAbsorptionLength()

    def mDOMacc(wavelength):
        return np.exp(-(mdom.GlassThickness/(glass.GetValue(wavelength)) )- (mdom.GelThickness/(gel.GetValue(wavelength))))*qe.GetValue(wavelength)

    wlens = np.arange(qe.GetMinWlen(), qe.GetMaxWlen(), qe.GetWavelengthStepping())
    values = np.array([mDOMacc(k) for k in wlens]) 
    mDOMacceptance = clsim.I3CLSimFunctionFromTable(qe.GetMinWlen(),\
                                                    qe.GetWavelengthStepping(),\
                                                    values)
    return mDOMacceptance

def construct_supremum(acceptances, scales):
    """
    Combine the xvalues for the broadest range
    with smallest binsize
    """
    import numpy as np
    from icecube.clsim import I3CLSimFunctionFromTable
    # construct combined range
    assert len(acceptances) == len(scales)
    steplen = min([a.GetWavelengthStepping() for a in acceptances])
    minwlen = min([a.GetMinWlen() for a in acceptances])
    maxwlen = max([a.GetMaxWlen() for a in acceptances])
    def get_supremum(wavelength):
        return max([scale*a.GetValue(wavelength) for a,scale in zip(acceptances, scales)])

    wlen = np.arange(minwlen, maxwlen+steplen, steplen)
    values = [get_supremum(_) for _ in wlen]
    return I3CLSimFunctionFromTable(minwlen, steplen, values)


def GetWavelengthAcceptanceEnvelope(sensors, oversize=1., qe_scale=1.):
    # TODO: construct envelope of quantum efficiencies for all sensor types
    # TODO: scale envelope by qe_scale*(mDOM radius/DOM radius)**2
    # TODO: scale evelope back down by oversize**2

    """
    A function to get the wavelength acceptance envelope of all OMs to be used

    This code loops over all the OMTypes needed for photon propagation.
    For each, it calculates the wavelength acceptance, including a 15% safety margin.
    Then, it calculates an acceptance "envelope", which is the maximum acceptance
    of any OMType at a given frequency. In pseudo-code:

        acceptance_envelope = []
        for w in wavelengths:
            acceptance_envelope.append( max (acceptance_mdom(w), acceptance_degg(w), ...))

    It will return an I3CLSimFunction, which will give the acceptance of the 
    envelope as a function of frequency


    Paramters
    ---------

    sensors: array of OMTypes
        An array of I3OMGeo.OMTypes which specifies the list of OMs for which
        the user wants an acceptance enveloppe

    oversize: float or double
        The oversizing factor applied to the DOM during photon propagation

    qe_scale: float or double
        Quantum Efficiency scaling to apply


    Returns
    -------

    I3CLSimFunction object
        The aceptance envelope as a function of frequency

    """

    if sensors is None:
        icetray.logging.log_fatal('sensors is {} -- abort'.format(sensors))

    # FIXME import constants from appropriate location
    domRadius=0.16510*icetray.I3Units.m
    mDOM_radius = 0.190 * icetray.I3Units.m #FIXME is this correct

    acceptances = []
    scales = []
    scale = qe_scale/(oversize**2)
    for sensor in sensors:
        if sensor == I3OMGeo.OMType.IceCube:
            # Oversample the IceCube DOMs
            # Cover both the standard and high-QE DOMs
            # including a typical SPE compensation factor of 1.32

            # standard IceCube DOM with spe compensation factor
            acceptances.append(clsim.GetIceCubeDOMAcceptance(domRadius=0.16510*icetray.I3Units.m,
                                                             efficiency=1.33))
            scales.append(1.15*scale)

            # high-QE IceCube DOM with spe compensation factor
            acceptances.append(clsim.GetIceCubeDOMAcceptance(domRadius=0.16510*icetray.I3Units.m,
                                                             efficiency=1.33, highQE=True))
            scales.append(1.15*scale)
        elif sensor == I3OMGeo.OMType.PDOM: 
            acceptances.append(clsim.GetIceCubeDOMAcceptance(domRadius=0.16510*icetray.I3Units.m,
                                                             efficiency=1., highQE=True))
            scales.append(1.15*scale)
        elif sensor == I3OMGeo.OMType.mDOM:
            acceptances.append(construct_mDOM_wavelength_acceptance())
            rqe = sensor_area(I3OMGeo.OMType.mDOM)/sensor_area(I3OMGeo.OMType.PDOM)
            scales.append(1.15*scale/rqe*((mDOM_radius/domRadius)**2))
        elif sensor == I3OMGeo.OMType.DEgg:
            acceptances.append(degg.GetDEggAcceptance())
            rqe = sensor_area(I3OMGeo.OMType.DEgg)/sensor_area(I3OMGeo.OMType.PDOM)
            scales.append(1.15*scale/rqe)
        else:
            raise ValueError("Unknown sensor {}".format(sensor))
    return construct_supremum(acceptances, scales)



def taskset(pid,tt=None):
    # get/set the taskset affinity for pid
    # uses a binary number string for the core affinity
    l = ['/bin/taskset','-p']
    if tt:
        l.append(hex(int(tt,2))[2:])
    l.append(str(pid))
    p = subprocess.Popen(l,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output = p.communicate()[0].split(':')[-1].strip()
    if not tt:
        return bin(int(output,16))[2:]

def tasksetInUse():
    # check for cpu affinity (taskset)
    try:
        num_cpus = reduce(lambda b,a: b+int('processor' in a),open('/proc/cpuinfo').readlines(),0)
        affinity = taskset(os.getpid())
        if len(affinity) < num_cpus:
            return True
        for x in affinity[:num_cpus]:
            if x != '1':
                return True
        return False
    except Exception:
        return False

def resetTasksetThreads(main_pid):
    # reset thread taskset affinity
    time.sleep(60)
    num_cpus = reduce(lambda b,a: b+int('processor' in a),open('/proc/cpuinfo').readlines(),0)
    tt = '1'*num_cpus
    #tt = taskset(main_pid)
    p = subprocess.Popen(['/bin/ps','-Lo','tid','--no-headers','%d'%main_pid],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    for tid in p.communicate()[0].split():
        tid = tid.strip()
        if tid:
            taskset(tid,tt)

def LoadCascadeTables(IceModel = "SpiceMie", TablePath = expandvars("$I3_DATA/photon-tables/splines")):
    if IceModel == "Spice1":
        amplitudetable = TablePath+'/ems_spice1_z20_a10.abs.fits'
        timingtable = TablePath+'/ems_spice1_z20_a10.prob.fits'
    elif IceModel == "SpiceMie":
        amplitudetable = TablePath+'/ems_mie_z20_a10.abs.fits'
        timingtable = TablePath+'/ems_mie_z20_a10.prob.fits'
    elif IceModel == "SpiceMieNoHoleIce":
        amplitudetable = TablePath+'/NoHoleIceCascades_250_z20_a10.abs.fits'
        timingtable = TablePath+'/NoHoleIceCascades_250_z20_a10.prob.fits'
    elif IceModel == "SpiceLea":
        amplitudetable = TablePath+'/cascade_single_spice_lea_flat_z20_a10.abs.fits'
        timingtable = TablePath+'/cascade_single_spice_lea_flat_z20_a10.prob.fits' 
    else:
        raise RuntimeError("Unknown ice model: %s", IceModel)
   
    print("Loading cascade tables : ")
    print("  amp  = ", amplitudetable)
    print("  time = ", timingtable)
    cascade_service = photonics_service.I3PhotoSplineService(
        amplitudetable=amplitudetable,
        timingtable=timingtable,
        timingSigma=0.)
    return cascade_service

@icetray.traysegment
def CleanupSlicedMCTree(tray, name, MCPESeriesName='I3MCPESeriesMap', PhotonSeriesMapName=None, KeepSlicedMCTree=False, InputMCTree='I3MCTree'):
    sliceRemoverAdditionalParams = dict()
    if PhotonSeriesMapName is not None:
        sliceRemoverAdditionalParams["InputPhotonSeriesMapName"] = PhotonSeriesMapName
        sliceRemoverAdditionalParams["OutputPhotonSeriesMapName"] = PhotonSeriesMapName
    # re-assign the output hits from the sliced tree to the original tree
    tray.AddModule("I3MuonSliceRemoverAndPulseRelabeler",
        InputMCTreeName = InputMCTree+"_sliced",
        OldMCTreeName = InputMCTree,
        InputMCPESeriesMapName = MCPESeriesName,
        OutputMCPESeriesMapName = MCPESeriesName,
        **sliceRemoverAdditionalParams
        )
    if not KeepSlicedMCTree:
        tray.AddModule("Delete", name+"_cleanup_clsim_sliced_MCTree",
            Keys = [InputMCTree+"_sliced"])

@icetray.traysegment
def PropagatePhotons(tray, name, GCDFile,
    RandomService = None,
    KeepIndividualMaps = False,
    HybridMode = False,
    IgnoreMuons = False,
    IgnoreCascades = False,
    UseGPUs = True,
    UseAllCPUCores = False,
    KeepSlicedMCTree = False,
    IceModel = "SpiceBFRV1Complete",
    CascadeService = None,
    IceModelLocation = None,
    UseGeant4=False, 
    CrossoverEnergyEM=None, 
    CrossoverEnergyHadron=None, 
    UnshadowedFraction = 0.99, #changed 2014-10-16 to IC86 nominal preset, IC79 used 0.9
    DOMOversizeFactor=5.0,
    HoleIceParameterization=expandvars("$I3_SRC/ice-models/resources/models/angsens/as.h2-50cm"),
    #UseHoleIceParameterization=True,
    InputMCTree="I3MCTree",
    OutputPESeriesMapName="I3MCPESeriesMap",
    OutputPhotonSeriesName=None,
    SimulateMDOMs=[],
    IgnoreSubdetectors=[],
    DisableTilt=False,
    MaxPhotonsPerSource=-1):
    """ This traysegment offers multiple tweaks for adapted processing in different energy ranges,
    for example GEANT4 in conjunction with Parametrizations for the treatment for lowest energies
    and a HybridMode with the use of tables for the treatment of high energies.
    In any case, please refer to the documentation of clsim to find suitable settings for your
    simulation needs """
    from I3Tray import I3Units
    from icecube import icetray, dataclasses, dataio
    from icecube import phys_services, sim_services
    from icecube import clsim

    if RandomService is None:
        randomService = tray.context['I3RandomService']
    elif isinstance(RandomService, str):
        randomService = tray.context[RandomService]
    else:
        randomService = RandomService
   
    if IgnoreMuons and not HybridMode:
        raise RuntimeError("Can currently only ignore muons in hybrid mode")
    if IceModelLocation is None:
       IceModelLocation = expandvars("$I3_SRC/ice-models/resources/models")
    if IceModel == "SpiceMie":
        clsimIceModel = expandvars(IceModelLocation+"/spice_mie")
    elif IceModel == "SpiceLea":
        clsimIceModel = expandvars(IceModelLocation+"/spice_lea")
    elif IceModel == "Spice3":
        clsimIceModel = expandvars(IceModelLocation+"/spice_3")
    elif IceModel == "Spice3.1":
        clsimIceModel = expandvars(IceModelLocation+"/spice_3.1")
    elif IceModel == "Spice3.2":
        clsimIceModel = expandvars(IceModelLocation+"/spice_3.2")
    elif IceModel == "Spice3.2.1":
        clsimIceModel = expandvars(IceModelLocation+"/spice_3.2.1")
    elif IceModel == "SpiceBFRV1":
        clsimIceModel = expandvars(IceModelLocation+"/spice_bfr-dv1")
    elif IceModel == "SpiceBFRV1Complete":
        clsimIceModel = expandvars(IceModelLocation+"/spice_bfr-dv1_complete")
    elif IceModel == "SpiceBFRV2Complete":
        clsimIceModel = expandvars(IceModelLocation+"/spice_BFRv2_complete/ice-bfr-v2-x5")
    elif os.path.exists(os.path.join(IceModelLocation, IceModel)):
        clsimIceModel = expandvars(os.path.join(IceModelLocation, IceModel))
    else:
        raise RuntimeError("Unknown ice model: %s", IceModel)
    # if 'mie' not in IceModel.lower() and HybridMode:
    #     raise RuntimeError(f"Cannot use {IceModel} in hybrid mode.")

    if (not IgnoreCascades) and HybridMode:
        if CascadeService is None:
            print("*** no cascades tables provided. Loading tables for", IceModel)
           
            # If we can see CVMFS, we'll get the splines from there.
            TablePath=expandvars("$I3_DATA/photon-tables/splines")
            if not os.path.isdir(TablePath):
               TablePath="/data/sim/sim-new/spline-tables"
           
            print("Using splines from: ", TablePath)
            # Work out which splines to use based on ice model preferences
            #if(UseHoleIceParameterization):

            # FIXME: Where goes the hole ice parametrization?
            CascadeModel=IceModel
            #else:
            #    if IceModel=="SpiceMie":
            #        CascadeModel="SpiceMieNoHoleIce"
            #    else:
            #        raise RuntimeError("No no-hole-ice spline for %s", IceModel)
                                       
            cascade_service = LoadCascadeTables(IceModel=CascadeModel, TablePath=TablePath)
        else:
            cascade_service = CascadeService
    else:
        cascade_service = None
    if HybridMode:
        if OutputPhotonSeriesName is not None:
            raise RuntimeError("saving photons is not supported in hybrid mode")
        if UseGeant4:
            raise RuntimeError("Geant4 not supported in hybrid mode")
        if ((CrossoverEnergyEM is not None) or (CrossoverEnergyHadron is not None)):
            raise RuntimeError("CrossoverEnergyEM or CrossoverEnergyHadron not supported in hybrid mode")
        if SimulateMDOMs:
            raise ValueError("Can't simulate mDOMs in hybrid mode")
        # split the MCTree into a cascade-only and a track-only version
        tray.AddModule("I3MCTreeHybridSimulationSplitter", name+"_splitMCTree",
            InputMCTreeName=InputMCTree,
            OutputMCTreeNameTracks=InputMCTree+"Tracks",
            OutputMCTreeNameCascades=InputMCTree+"Cascades")
        tray.AddModule("I3TauSanitizer", name+"_sanitize_taus",
            InputMCTreeName = InputMCTree+"Tracks",
            OutputMCTreeName = InputMCTree+"Tracks") # overwrite the input
        if not IgnoreMuons:
            if UseGPUs:
                DoNotParallelize=False
            else:
                DoNotParallelize=not UseAllCPUCores
                threading.Thread(target=resetTasksetThreads,args=(os.getpid(),)).start()
            print('tasksetInUse = ',tasksetInUse())
            print('DoNotParallelize = ',DoNotParallelize)
            # simulate tracks (with clsim)
            tray.AddSegment(
                clsim.I3CLSimMakeHits, name+"_makeCLSimHits",
                GCDFile = GCDFile,
                PhotonSeriesName = None,
                MCTreeName = InputMCTree+"Tracks",
                OutputMCTreeName = InputMCTree+"Tracks_sliced",
                MCPESeriesName = OutputPESeriesMapName + "Tracks",
                RandomService = randomService,
                UnshadowedFraction=UnshadowedFraction,
                DoNotParallelize=DoNotParallelize,
                UseGeant4=False, # never use this with Geant4!
                UseGPUs=UseGPUs,
                UseCPUs=not UseGPUs,
                IceModelLocation=clsimIceModel,
                DOMOversizeFactor=DOMOversizeFactor,
                DisableTilt=True)

            if KeepSlicedMCTree:
                raise RuntimeError("cannot use KeepSlicedMCTree=True in hybrid simulation mode")
            tray.AddModule("Delete", name+"_cleanup_clsim_sliced_MCTree",
                Keys = [InputMCTree+"Tracks_sliced"])
        if not IgnoreCascades:
            tray.AddModule("I3PhotonicsHitMaker", name+"_hitsFromTheTable",
                CascadeService = cascade_service,
                TrackService = None, # tracks are handled by clsim
                UnshadowedFraction = UnshadowedFraction,
                Input = InputMCTree+"Cascades",
                Output = OutputPESeriesMapName + "Cascades",
                RandomService = randomService,
                MaxPhotonsPerSource = MaxPhotonsPerSource
                )
        MCPEsToCombine = []
        if not IgnoreMuons:
            MCPEsToCombine.append(OutputPESeriesMapName + "Tracks")
        if not IgnoreCascades:
            MCPEsToCombine.append(OutputPESeriesMapName + "Cascades")
        # combine the resulting I3MCHitSeriesMaps
        tray.AddModule("I3CombineMCPE", name+"_combine_pes",
            InputResponses = MCPEsToCombine,
            OutputResponse = OutputPESeriesMapName)
        if not KeepIndividualMaps:
            # delete the original maps and the split I3MCTrees
            tray.AddModule("Delete", name+"_cleanup_peseriesmaps",
                Keys = MCPEsToCombine)
            tray.AddModule("Delete", name+"_cleanup_MCTree",
                Keys=[InputMCTree+"Tracks", InputMCTree+"Cascades"])
    else:
        # non-hybrid clsim-only simulation
        tray.AddModule("I3TauSanitizer", name+"_sanitize_taus",
            InputMCTreeName = InputMCTree,
            OutputMCTreeName = InputMCTree) # overwrite the input
        if UseGPUs:
            DoNotParallelize=False
        else:
            DoNotParallelize=not UseAllCPUCores
            threading.Thread(target=resetTasksetThreads,args=(os.getpid(),)).start()
        print('tasksetInUse = ',tasksetInUse())
        print('DoNotParallelize = ',DoNotParallelize)
        # simulate tracks (with clsim)
        kwargs = dict(
            GCDFile = GCDFile,
            PhotonSeriesName = OutputPhotonSeriesName,
            MCTreeName = InputMCTree,
            OutputMCTreeName = InputMCTree+"_sliced",
            MCPESeriesName = OutputPESeriesMapName+"_IC",
            RandomService = randomService,
            UnshadowedFraction = UnshadowedFraction,
            DoNotParallelize = DoNotParallelize,
            UseGeant4=UseGeant4,
            CrossoverEnergyEM=CrossoverEnergyEM,
            CrossoverEnergyHadron=CrossoverEnergyHadron,
            UseGPUs=UseGPUs,
            UseCPUs=not UseGPUs,
            DOMOversizeFactor=DOMOversizeFactor,
            #UseHoleIceParameterization=UseHoleIceParameterization,
            IceModelLocation=clsimIceModel,
            DisableTilt=DisableTilt)
        
        tray.Add('I3GeometryDecomposer', If=lambda frame: "I3OMGeoMap" not in frame)
        def hex_subdetector(frame):
            subdetectors = copy.copy(frame['Subdetectors'])
            for k in subdetectors.keys():
                if k.string >= 1000:
                    subdetectors[k] = "HEX"
            del frame['Subdetectors']
            frame['Subdetectors'] = subdetectors
        tray.Add(hex_subdetector, Streams=[icetray.I3Frame.Geometry])
        
        if len(SimulateMDOMs) > 0:
            # simulate new strings with mDOMs
            ignored = set(IgnoreSubdetectors + ['IceTop'])
            assert len(ignored.intersection(SimulateMDOMs)) == 0
            
            from icecube.mDOM_WOM_simulation import mDOMMakePhotons, mDOMMakeHitsFromPhotons
            
            hex_photons = "HEXPhotons"
            hex_kwargs = dict(kwargs)
            hex_kwargs['PhotonSeriesName'] = hex_photons
            tray.AddSegment(mDOMMakePhotons, ExtraArgumentsToI3CLSimClientModule=dict(IgnoreSubdetectors=list(ignored)), **hex_kwargs)
            
            tray.AddSegment(mDOMMakeHitsFromPhotons,
                InputPhotonSeriesMapName=hex_photons,
                OutputMCHitSeriesMapName=OutputPESeriesMapName+"_HEX",
                DOMOversizeFactor=DOMOversizeFactor,
                MCTreeName="I3MCTree_sliced",
                RandomService=randomService)
            
            tray.Add(CleanupSlicedMCTree, MCPESeriesName=OutputPESeriesMapName+"_HEX")
            
            to_combine = [OutputPESeriesMapName+"_IC", OutputPESeriesMapName+"_HEX"]
            tray.AddModule("I3CombineMCPE",
                InputResponses = to_combine,
                OutputResponse = OutputPESeriesMapName)
            
            tray.Add("Delete", keys=[hex_photons] + to_combine)
        
        elif 'IceCube' in IgnoreSubdetectors:
            # simulate new strings with pDOMs
            hex_kwargs = dict(kwargs)
            hex_kwargs['MCPESeriesName'] = OutputPESeriesMapName + "_HEX"
            ignored = set(IgnoreSubdetectors + ['IceTop'])
            tray.AddSegment(clsim.I3CLSimMakeHits,
                            ExtraArgumentsToI3CLSimClientModule=dict(IgnoreSubdetectors=list(ignored)),
                            **hex_kwargs)
            
            tray.Add(CleanupSlicedMCTree, MCPESeriesName=hex_kwargs['MCPESeriesName'])
            
            to_combine = [OutputPESeriesMapName+"_IC", OutputPESeriesMapName+"_HEX"]
            tray.AddModule("I3CombineMCPE",
                InputResponses = to_combine,
                OutputResponse = OutputPESeriesMapName)
            
            tray.Add("Delete", keys=to_combine)
        else:
            # just simulate IceCube
            ignored = set(IgnoreSubdetectors + ['IceTop'])
            tray.AddSegment(clsim.I3CLSimMakeHits, ExtraArgumentsToI3CLSimClientModule=dict(IgnoreSubdetectors=list(ignored)), **kwargs)
            
            if len(IgnoreSubdetectors) > 0:
                tray.Add(CleanupSlicedMCTree, MCPESeriesName=OutputPESeriesMapName+"_IC")

def sensor_area(sensor):
    """
    Return the sensor's effective area averaged over all angles and wavelengths
    between 300 and 600 nm (Cherenkov-weighted)
    """
    if sensor == I3OMGeo.OMType.PDOM:
        return 29.2*I3Units.cm2
    elif sensor == I3OMGeo.OMType.DEgg:
        return 43.3*I3Units.cm2
    elif sensor == I3OMGeo.OMType.mDOM:
        return 65.4*I3Units.cm2
    else:
        raise ValueError("Unsupported sensor '%s'" % sensor)

@icetray.traysegment
def MakePhotonsMultiSensor(
    tray, name, GCDFile,
    Sensors=[I3OMGeo.OMType.IceCube, I3OMGeo.OMType.PDOM, I3OMGeo.OMType.mDOM, I3OMGeo.OMType.DEgg],
    InputMCTree="I3MCTree",
    # OutputPESeriesMapName="I3MCPESeriesMap",
    OutputPhotonSeriesName="I3PhotonSeriesMap",
    UseGPUs=True,
    UseCPUs=False,
    DOMOversizeFactor=1.,
    EfficiencyScale=1.,
    DisableTilt=False,
    RandomService=None,
    IceModelLocation="$I3_SRC/ice-models/resources/models/",
    IceModel="spice_bfr-dv1_complete",
    OutputMCTreeName=None,
    UseGeant4 = False, 
    UseI3PropagatorService=True,
    CrossoverEnergyEM = 0.1, # only when UseGeant4 = True
    CrossoverEnergyHadron = 30.0, # only when UseGeant4 = True
    HoleIceParameterization=expandvars("$I3_SRC/ice-models/resources/models/angsens/as.nominal"),
    StopDetectedPhotons = True, 
    UnshadowedFraction=0.99, # XXX: FIXME: this one is default from clsim
    IgnoreSubdetectors=['IceTop'],
    RunMPHitFilter=True, #Run polyplopia's mphitfilter (removes events that don't leave photons)
    ):

    """
    Propagate photons to 13 inch spheres, using an effective Cherenkov spectrum
    scaled so that it can be downsampled to the quantum efficiencies of each
    sensor, scaled up so that they are equivalent to `EfficiencyScale` times
    the the pDOM area
    """
    
    from icecube import clsim
   
    wavelengthAcceptance = GetWavelengthAcceptanceEnvelope(Sensors, oversize=DOMOversizeFactor, qe_scale=EfficiencyScale)
    if RandomService is None:
        randomService = tray.context['I3RandomService']
    elif isinstance(RandomService, str):
        randomService = tray.context[RandomService]
    else:
        randomService = RandomService
    
    kwargs = dict(
        GCDFile = GCDFile,
        PhotonSeriesName = OutputPhotonSeriesName,
        MCPESeriesName = None, # turn off conversion to MCPESeries at this step
        MCTreeName = InputMCTree,
        OutputMCTreeName = OutputMCTreeName,
        RandomService = randomService,
        UseGPUs=UseGPUs,
        UseCPUs=not UseGPUs,
        UnWeightedPhotons=False,
        UnWeightedPhotonsScalingFactor=None,
        DOMOversizeFactor=DOMOversizeFactor,
        WavelengthAcceptance=wavelengthAcceptance,
        IceModelLocation=expandvars(os.path.join(IceModelLocation, IceModel)),
        DisableTilt=DisableTilt, 
        UseGeant4 = UseGeant4,
        UseCascadeExtension=True,
        CrossoverEnergyEM = CrossoverEnergyEM,
        CrossoverEnergyHadron=CrossoverEnergyHadron,
        StopDetectedPhotons=StopDetectedPhotons,
        HoleIceParameterization=HoleIceParameterization,
        UnshadowedFraction = UnshadowedFraction,
        PhotonHistoryEntries=0,
        DoNotParallelize=False,
        IgnoreSubdetectors = IgnoreSubdetectors,
        UseI3PropagatorService=UseI3PropagatorService,
        ExtraArgumentsToI3CLSimClientModule={
                       "StatisticsName":"clsim_stats",
        }
        )

    icetray.logging.log_warn("env:")
    for k in sorted(os.environ.keys()):
        icetray.logging.log_warn("%20s: %s" % (k, os.environ[k]))
    icetray.logging.log_warn("OpenCL devices:")
    for d in clsim.I3CLSimOpenCLDevice.GetAllDevices():
        icetray.logging.log_warn(str(d))
    tray.Add(clsim.traysegments.I3CLSimMakePhotons,
             "multisensphot", **kwargs)

def trim_calibration(frame):
    """Slim down I3Calibration object by removing Gen2 entries not needed for photon->PE conversion"""
    calib = frame['I3Calibration']
    geo = frame['I3Geometry'].omgeo
    del frame['I3Calibration']
    for k in calib.dom_cal.keys():
        if k.string > 86:
            if geo[k].omtype == I3OMGeo.OMType.mDOM:
                del calib.dom_cal[k]
    frame['I3Calibration'] = calib

@icetray.traysegment
def MakePEFromPhotons(
    tray, name, GCDFile,
    Sensors=[I3OMGeo.OMType.IceCube,
             I3OMGeo.OMType.PDOM,
             I3OMGeo.OMType.mDOM,
             I3OMGeo.OMType.DEgg,],
    PhotonSeriesName="I3PhotonSeriesMap",
    DOMOversizeFactor=1.,
    EfficiencyScales=[1.0, 1.0, 2.2, 1.5],
    RandomService=None,
    IceModelLocation='$I3_SRC/ice-models/resources/models/',
    IceModel="spice_bfr-dv1_complete",
    HoleIceParameterization = expandvars("$I3_SRC/ice-models/resources/models/angsens/as.nominal"),
    HoleIceParameterizationGen2 = "",
    KeepPhotonSeries = True
    ):
    """
    :param EfficiencyScale: scale QE so that the average photon effective area
    of the sensor is equivalent to this multiple of the pDOM area
    """
    from icecube import simclasses
    from icecube.clsim.traysegments import I3CLSimMakeHitsFromPhotons
    from icecube.gen2_sim.utils import split_photons, ModuleToString

    tray.Add("I3GeometryDecomposer", If=lambda frame: "I3OMGeoMap" not in frame)

    # Break the photons up now to make life easier
    tray.Add(split_photons,
             Input = PhotonSeriesName,
             Sensors = Sensors,
             Streams = [icetray.I3Frame.DAQ],
             KeepPhotonSeries = KeepPhotonSeries)

    photonList = []

    randomService = tray.context['I3RandomService'] if RandomService is None else RandomService

    if randomService is None :
        assert tray.context['I3RandomService'] is not None, "Both `randomService` argument and random sercice in IceTray `context` are `None`"
        randomService = tray.context['I3RandomService']

    for i, Sensor in enumerate(Sensors): 
        currentPhotonSeriesName = PhotonSeriesName + ModuleToString(Sensor)
        currentMCPESeriesName = currentPhotonSeriesName.replace("Photon", "MCPE")
        photonList.append(currentPhotonSeriesName)
        
        EfficiencyScale, rqe = EfficiencyScales[i], 1.0
        if not Sensor == I3OMGeo.OMType.IceCube:
            rqe = EfficiencyScale*sensor_area(I3OMGeo.OMType.PDOM)/sensor_area(Sensor)

        # TODO: add a PE converter based on the value of `Sensor`
        # be sure to scale mDOM QE *up* by (mDOM radius/DOM radius)**2
        # to account for the fact that photons were captured on a smaller area
        assert Sensor in [I3OMGeo.OMType.IceCube, I3OMGeo.OMType.PDOM, I3OMGeo.OMType.mDOM, I3OMGeo.OMType.DEgg], "Unknown sensor! {}".format(sensor)
        icetray.logging.log_info("Will calculate PEs for {} with photocathode area equivalent to {}x pDOM ({}x {})".format(Sensor, EfficiencyScale, rqe, Sensor))

        HoleIceParameterization= expandvars(HoleIceParameterization)
        HoleIceParameterizationGen2 = expandvars(HoleIceParameterizationGen2)
        if HoleIceParameterizationGen2 == "":
            print("No separate hole ice for Gen2, using IceCube one. Warning applicable to pDOM/DEgg Gen2 only")
            HoleIceParameterizationGen2 = HoleIceParameterization

        if Sensor == I3OMGeo.OMType.PDOM: 
            tray.Add(I3CLSimMakeHitsFromPhotons, name + "_PEconverterForPDOMs",
                     PhotonSeriesName        = currentPhotonSeriesName,
                     MCPESeriesName          = currentMCPESeriesName,
                     RandomService           = randomService,
                     DOMOversizeFactor       = DOMOversizeFactor,
                     UnshadowedFraction      = rqe,
                     HoleIceParameterization = HoleIceParameterizationGen2,
                     GCDFile                 = GCDFile,
                     IceModelLocation        = expandvars(
                         os.path.join(IceModelLocation, IceModel)),
                     IgnoreSubdetectors      = ['IceTop', 'HEX']
            )
            tray.Add(CleanupSlicedMCTree,
                     MCPESeriesName = currentMCPESeriesName,
                     KeepSlicedMCTree=True, If=lambda frame:frame.Has('I3MCTree_sliced'))

        elif Sensor == I3OMGeo.OMType.mDOM:
            # in case of Gen2 mDOMs I3Calibration needs to be reduced, otherwise it's
            # too big to hold in memory
            tray.Add(trim_calibration, Streams=[icetray.I3Frame.Calibration])
            tray.Add(mDOMMakeHitsFromPhotons, name + "_PEconverterForMDOMs",
                     InputPhotonSeriesMapName = currentPhotonSeriesName,
                     DOMOversizeFactor        = DOMOversizeFactor,
                     EfficiencyScale          = rqe,
                     OutputMCHitSeriesMapName = currentMCPESeriesName,
                     RandomService            = randomService)
            tray.Add(CleanupSlicedMCTree,
                     MCPESeriesName = currentMCPESeriesName,
                     KeepSlicedMCTree=True, If=lambda frame:frame.Has('I3MCTree_sliced'))

        elif Sensor == I3OMGeo.OMType.DEgg:
            tray.Add(DEggMakeHitsFromPhotons, name + "_PEconverterForDEggs",
                     PhotonSeriesName        = currentPhotonSeriesName,
                     MCPESeriesName          = currentMCPESeriesName,
                     RandomService           = randomService,
                     DOMOversizeFactor       = DOMOversizeFactor,
                     UnshadowedFraction      = rqe,
                     HoleIceParameterization = HoleIceParameterizationGen2
            )
            tray.Add(CleanupSlicedMCTree,
                     MCPESeriesName = currentMCPESeriesName,
                     KeepSlicedMCTree=True, If=lambda frame:frame.Has('I3MCTree_sliced'))

        elif Sensor == I3OMGeo.OMType.IceCube:
            tray.Add(I3CLSimMakeHitsFromPhotons, name + "_PEconverterForDOMS",
                     PhotonSeriesName        = currentPhotonSeriesName,
                     MCPESeriesName          = currentMCPESeriesName,
                     RandomService           = randomService,
                     DOMOversizeFactor       = DOMOversizeFactor,
                     UnshadowedFraction      = 0.99,
                     HoleIceParameterization = HoleIceParameterization,
                     GCDFile                 = GCDFile,
                     IceModelLocation        = expandvars(
                         os.path.join(IceModelLocation, IceModel)),
                     IgnoreSubdetectors      = ['IceTop', 'HEX']
            )
            tray.Add(CleanupSlicedMCTree,
                     MCPESeriesName = currentMCPESeriesName,
                     KeepSlicedMCTree=True, If=lambda frame:frame.Has('I3MCTree_sliced'))

            # TODO: Need to switch to the CLSim ExtensionCylinder branch or merge it back into
            # trunk in order to handle WOMs

    # Clean up the photons before we move on
    tray.Add("Delete", Keys = photonList)
    tray.Add("Delete", Keys = ['I3MCTree_sliced',], If=lambda frame:frame.Has('I3MCTree_sliced'))


@icetray.traysegment
def CLSimTraySegment(tray,name, GCDFile,
        gpu=-1,  efficiency=0.99, oversize=5, 
        amplitude_fits=None, timing_fits=None,
        IceModelLocation=expandvars("$I3_SRC/ice-models/resources/models/"),
        IceModel="SpiceBFRV1Complete",
        UseGeant4=False,
        mcpeseries = "I3MCPESeriesMap",
        photonseries = None,
        DisableTilt = False,
        HybridMode = False,
        IgnoreMuons = False,
        SimulateMDOMs = [],
        IgnoreSubdetectors = [],):
        # Load libraries
        from icecube import icetray
        from icecube import dataclasses, dataio, phys_services, interfaces
        from icecube import interfaces,dataclasses,simclasses
        from icecube import sim_services,dataio
        UseGPU = gpu >= 0 
        # if UseGPU and (isinstance(gpu,str) and gpu.isdigit() or isinstance(gpu,int)):
        #    os.putenv("CUDA_VISIBLE_DEVICES",str(gpu))
        #    os.putenv("COMPUTE",":0."+str(gpu))
        #    os.putenv("GPU_DEVICE_ORDINAL",str(gpu))
        # Now we import the photon stuff
        from icecube import photonics_service
        randomService = tray.context["I3RandomService"]
        
        tray.AddSegment(
            PropagatePhotons, "normalpes",
            GCDFile = GCDFile,
            RandomService = randomService,
            KeepIndividualMaps = False,
            HybridMode = HybridMode,
            IgnoreMuons = IgnoreMuons,
            IgnoreCascades = False,
            UseGPUs = UseGPU,
            UseAllCPUCores = False,
            KeepSlicedMCTree = False,
            IceModel = IceModel,
            IceModelLocation = IceModelLocation,
            UnshadowedFraction = efficiency,
            DOMOversizeFactor = oversize,
            InputMCTree="I3MCTree",
            UseGeant4=UseGeant4,
            OutputPESeriesMapName=mcpeseries,
            OutputPhotonSeriesName=photonseries,
            DisableTilt=DisableTilt,
            SimulateMDOMs=SimulateMDOMs,
            IgnoreSubdetectors=IgnoreSubdetectors)

class ClSim(ipmodule.ParsingModule):
   """
    GPU Photon propagation
   """
   def __init__(self):
        ipmodule.ParsingModule.__init__(self)
        self.AddParameter('gcdfile','GeoCalibDetStatus filename','')
        self.AddParameter('inputfilelist','list of input filenames',[])
        self.AddParameter('outputfile','Output filename','')
        self.AddParameter('seed','RNG seed',0)
        self.AddParameter('procnum','job number (RNG stream number)',0)
        self.AddParameter('nproc','Number of jobs (Number of RNG streams)',1)
        self.AddParameter('summaryfile','XMLSummary filename','summary.xml')
        self.AddParameter('MJD','MJD (0=do not modify)',0)
        self.AddParameter("GPU", "Graphics Processing Unit number",0)
        self.AddParameter('weightsumname','Name of WeighSum in Frame','WeightSum')
        self.AddParameter('weightpartsum','Name of WeighPartialSum in Frame','WeightPartSum')
        self.AddParameter('NumberOfPrimaries',
                        'Polyplopia force multiple tracks (integer value of coincidences)',0)
        self.AddParameter('ParticleType',
                        'type of particle that we are simulating','corsika')
        self.AddParameter('RunMPHitFilter',"Run polyplopia's mphitfilter", True)
        self.AddParameter("oversize","over-R: DOM radius oversize scaling factor",5)
        self.AddParameter("efficiency","overall DOM efficiency correction",[0.99])
        self.AddParameter("volumecyl","set volume to regular cylinder (set to False for 300m spacing from the DOMs)",True)
        self.AddParameter("IceModelLocation","Location of ice model param files", expandvars("$I3_SRC/ice-models/resources/models"))
        self.AddParameter("IceModel","ice model subdirectory", "SpiceBFRV1Complete") 
        self.AddParameter("PhotonSeriesName","Photon Series Name","I3MCPESeriesMap") 
        self.AddParameter("RawPhotonSeriesName","Raw Photon Series Name",None) 
        self.AddParameter("UseGeant4","Enable Geant4 propagation",False) 
        self.AddParameter("DisableTilt","Turn off tilt in ice model",False)
        self.AddParameter("SimulateMDOMs","Simulate mDOMs for these subdetectors",[])
        self.AddParameter("IgnoreSubdetectors", "Ignore these subdetectors",[])

   def Execute(self,stats):
        if not ipmodule.ParsingModule.Execute(self,stats): return 0
        from icecube import icetray, dataclasses, dataio, phys_services, interfaces

        from I3Tray import I3Tray,I3Units
        # Instantiate a tray
        tray = I3Tray()
        inputfiles  = self.GetParameter('inputfilelist')
        summary_in  = self.summaryfile
        summary_out = self.summaryfile
        if not os.path.exists(self.summaryfile):
           summary_in  = ''
        # Configure IceTray services
        tray.AddService("I3XMLSummaryServiceFactory","summary",outputfilename=summary_out,inputfilename=summary_in)
        # RNG
        rngstate    = "rng.state"
        if not os.path.exists(rngstate): 
           rngstate = ''
           print("Warning: no RNG state found. Using seed instead.")
        tray.AddService("I3SPRNGRandomServiceFactory","sprngrandom",
            seed=self.seed, streamNum=self.procnum,nstreams=self.nproc,
            instatefile=rngstate,outstatefile="rng.state")
        # Configure IceTray modules
        tray.AddModule("I3Reader", "reader",filenamelist=[self.gcdfile]+inputfiles)
        efficiency = self.efficiency
        if type(self.efficiency) == list or type(self.efficiency) == tuple:
           if len(self.efficiency) == 1:
              efficiency=float(self.efficiency[0])
           elif len(self.efficiency) > 1:
              efficiency=map(float,self.efficiency)
           elif len(self.efficiency) > 1:
              raise Exception("Configured empty efficiency list")
        else:
           efficiency=float(self.efficiency)

        tray.AddSegment(CLSimTraySegment,"photons",
                        self.gcdfile,
                        gpu              = self.gpu,
                        efficiency       = efficiency,
                        oversize         = self.oversize,
                        IceModelLocation = self.icemodellocation,
                        IceModel         = self.icemodel,
                        mcpeseries       = self.photonseriesname,
                        UseGeant4        = self.usegeant4,
                        photonseries     = self.rawphotonseriesname,
                        DisableTilt      = self.disabletilt,
                        SimulateMDOMs    = self.simulatemdoms,
                        IgnoreSubdetectors = self.ignoresubdetectors,
                )
        if self.runmphitfilter:
            from icecube import polyplopia
            tray.AddModule("MPHitFilter","hitfilter",
                           HitOMThreshold=1,
                           I3MCPESeriesMapName=self.photonseriesname)

        high_eff_pe_series_name = self.photonseriesname
        if type(efficiency) == list or type(efficiency) == tuple:
              high_eff_pe_series_name = self.photonseriesname+"_"+str(max(efficiency))
        # Remove original
        if type(efficiency) == list or type(efficiency) == tuple:
           tray.AddModule("Delete","clean_"+self.photonseriesname,Keys=[self.photonseriesname])
        tray.AddModule("I3Writer","writer", filename=self.outputfile,streams=[icetray.I3Frame.DAQ])

        tray.AddModule("TrashCan","trashcan")

        # Execute the Tray
        tray.Execute()
        tray.Finish()
        del tray
        return 0
   def Finish(self,stats={}):
        self.logger.info("finish %s: %s" % (self.__class__.__name__,self.GetParameter("execute")))
        return 0


MultiSensorPhotons = iceprod_wrapper(MakePhotonsMultiSensor)
DownsamplePEs = iceprod_wrapper(MakePEFromPhotons)

    
if __name__ == "__main__":

    import inspect
    spec = inspect.getfullargspec(CLSimTraySegment)
    offset = len(spec.args) - len(spec.defaults)
    defaults = dict(((k, spec.defaults[i-offset]) for i, k in enumerate(spec.args) if i >= offset))

    from argparse import ArgumentParser
    gcddir = '/cvmfs/icecube.opensciencegrid.org/users/gen2-optical-sim/gcd'
    parser = ArgumentParser()
    parser.add_argument('--infile', dest='infile', nargs='+')
    parser.add_argument('--outfile', dest='outfile')
    parser.add_argument('--gcd', type=str, help='GCDFile to use for simulation',
                        default=expandvars(
                            os.path.join(gcddir,
                                         'IceCubeHEX_Sunflower_240m_v4.0beta_ExtendedDepthRange.GCDScintAnt-shift.i3.gz')))
    parser.add_argument('--gcd1size', type=str, help='GCDFile with equal-size DOMs for photon prop',
                        default=expandvars(
                            os.path.join(gcddir,
                                         'IceCubeHEX_Sunflower_240m_v4.0beta_ExtendedDepthRange.GCDScintAnt-shift.i3.gz')))
    parser.add_argument('--hybrid', action="store_true", default=False)
    parser.add_argument('--sensors', nargs="+", default=['IceCube', 'mDOM'])
    parser.add_argument('--effs', type=float, nargs="+", default=[1.,2.2])
    parser.add_argument('--no-gpu', dest='gpu', action="store_false", default=True)
    parser.add_argument('--dataset', type=int, default=1)
    parser.add_argument('--nfiles', type=int, default=1000)
    parser.add_argument('--fileno', type=int, default=0)
    opts = parser.parse_args()

    from icecube import icetray, dataclasses, dataio, phys_services
    from I3Tray import I3Tray

    tray = I3Tray()

    randomService = phys_services.I3SPRNGRandomService(opts.dataset, 2*opts.nfiles, 2*opts.fileno)
    tray.context['I3RandomService'] = randomService

    # evilly monkey-patch GridFTP stager to create parent directories if needed
    # gridftp_args = dataio.GridFTPStager.__init__.func_defaults[1]
    # if not '-cd' in gridftp_args:
    #     gridftp_args.append('-cd')
    tray.context['I3FileStager'] = dataio.get_stagers()

    tray.Add("I3Reader",
             filenamelist=[opts.gcd]+opts.infile)

    tray.Add("Rename", Keys=['I3MCTreeInIce', 'I3MCTree'])
    tray.Add("Rename", Keys=['I3MCTree', 'I3MCTree_preMuonProp'])
    from icecube.simprod.segments import PropagateMuons
    randomService = phys_services.I3GSLRandomService(0)
    #tray.context['I3RandomService'] = randomService
    tray.AddSegment(PropagateMuons, 'propagator', 
		   RandomService=randomService) 

    if opts.hybrid:
        tray.AddSegment(CLSimTraySegment,"photons",
                        GCDFile = opts.gcd1size,
			IceModel = 'SpiceMie',
                        HybridMode = opts.hybrid,
                        IgnoreMuons = opts.hybrid,
                        gpu              = -1,
                        DisableTilt      = opts.hybrid,
                )
    else:
        sensors = [getattr(I3OMGeo.OMType, _) for _ in opts.sensors]
        tray.AddSegment(MakePhotonsMultiSensor,
                        GCDFile=opts.gcd1size,
                        Sensors=sensors,
                        UseGPUs=opts.gpu,
                        EfficiencyScale=max(opts.effs))
        # make pDOM hits at maximum scale as a sanity check. if the QE envelope
        # is wrong, this will raise an error
        tray.AddSegment(MakePEFromPhotons,
                        GCDFile=opts.gcd,
                        Sensors=sensors,
                        EfficiencyScales=opts.effs)
        tray.Add("Delete", Keys=['I3MCPESeriesMap'])
    tray.Add("Dump")

    tray.Add("I3Writer", Streams=[icetray.I3Frame.Stream(_) for _ in 'SIQP'],
        # DropOrphanStreams=[icetray.I3Frame.DAQ],
        filename=opts.outfile)

    tray.Execute()
