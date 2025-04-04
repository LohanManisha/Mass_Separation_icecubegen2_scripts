#!/usr/bin/env python
# coding: utf-8
## Using this script, can raed the value of reconstructed beta as well as Sref variables. 
# In[ ]:


#!/usr/bin/env python
import os
from os import listdir
from optparse import OptionParser
from os.path import expandvars, isfile, join
import matplotlib.pyplot as plt
import pylab, math
import numpy as np
import argparse
import glob


import icecube
from I3Tray import I3Tray
from icecube import dataio, dataclasses, icetray, simclasses
from icecube import rock_bottom, toprec, phys_services
from icecube.icetray import I3Units
from icecube.icetray.i3logging import log_info, log_warn, log_fatal
icetray.I3Logger.global_logger.set_level(icetray.I3LogLevel.LOG_INFO)


# --------------------------------------------------------------------------------------------------



parser = argparse.ArgumentParser()
parser.add_argument('--path', type=str, metavar='FILE', help='Set path to read', \
                    default="/user/surface7January/lowest_beta_Proton_shower/") ## directory containing the simulation files
parser.add_argument('--gcd', type=str, metavar='FILE', help='GCD file', \
  default="/user/IceCubeHEX_Sunflower_240m_v4.0beta_ExtendedDepthRange.GCDScintAnt-shift.i3.gz") ## specify GCD file
args  = parser.parse_args()
files = [args.path + f for f in listdir(args.path) if isfile(join(args.path, f)) and (f.endswith(".i3.gz") or f.endswith(".i3"))] # reads in all filenames in input directory
# ------------------------------------------------------------------------------------------------- #

def dist_perp_to_axis(I3particle, I3position):
    """ calculate distance of detector perpendicular to shower axis
    """
    xc, yc, zc = I3particle.pos.x, I3particle.pos.y, I3particle.pos.z # core position from MCPrimary
    nx, ny = I3particle.dir.x, I3particle.dir.y # direction from MCPrimary

    XPos, YPos, ZPos = I3position.x, I3position.y, I3position.z # coordinates from scintillator panel

    abs_x_sq = (XPos - xc)*(XPos - xc) + (YPos - yc)*(YPos - yc) + (ZPos - zc)*(ZPos - zc) # distance to core position

    n_prod_x = nx*(XPos - xc) + ny*(YPos - yc) - np.sqrt(1. - nx*nx - ny*ny)*(ZPos - zc) # vector projection

    return np.sqrt(abs_x_sq - n_prod_x*n_prod_x)

class ReadScints(icetray.I3Module):
    def __init__(self,ctx):

        icetray.I3Module.__init__(self,ctx)

    def Configure(self):

        log_info("Configuring " + self.name)

    def get_data(self, SiPMMC, MC, geometry, frame):

        vems = np.array([pulses[0].charge for key, pulses in SiPMMC]) # measured charge in scintillator panels (in units of MIP)
        vems_cond = np.logical_not(vems<0.5) # trigger condition: charge > 0.5 MIP
        vems = vems[vems_cond]
        
        if (len(vems)>=5): # typical event trigger -> at least 5 hit scintillator panels

            time_reco = np.array([pulses[0].time for key, pulses in SiPMMC]) # read in 'raw' scintillator times (in ns)
            tmin_reco = min(time_reco) # t0
            timeRel_reco = np.array([ t-tmin_reco for t in time_reco ]) # time relative to t0

            scintkeyAll = [i for i in geometry.scintgeo.keys()]
            scintkey = np.asarray([dataclasses.ScintKey(key.station,key.panel) for key, pulses in SiPMMC])[vems_cond] # ScintKeys
            stations = np.array([key.station for key, pulses in SiPMMC])[vems_cond] # station number
            panels = np.array([key.panel for key, pulses in SiPMMC])[vems_cond] # panel number

            rad = np.array([dist_perp_to_axis(MC, geometry.scintgeo[sk].position) for sk in scintkey])
            res = np.array([item for item in scintkeyAll if item not in scintkey])
            
 
            sil_rad = np.array([dist_perp_to_axis(MC, geometry.scintgeo[i].position) for i in res]) 
            print(sil_rad.shape[0])
            if not "ScintReconstructionStep3aParams" in frame: return 
            Reco_pars = frame["ScintReconstructionStep3aParams"]

            lgSref = Reco_pars.GetParameterByName("lgSref", "ScintSignalModel", 0) # log_10 Sref (energy estimator)
            beta = Reco_pars.GetParameterByName("slopeLdf", "ScintSignalModel", 0) # beta (slope of the LDF)
            kappa = Reco_pars.GetParameterByName("slopeLdf", "ScintSignalModel", 1) # kappa parameter of the DLP function
            ##with open("output.txt", "a") as f:
            print(beta)
            print(lgSref)
            ##print(rad)

    def Physics(self, frame):

        log_info("Entering Physics...")
        
        MC = frame["MCPrimary"]
        SiPMMC = frame['SiPMRecoPulses']
        geometryScint = frame["I3ScintGeometry"]
    
        self.get_data(SiPMMC, MC, geometryScint, frame)
        
    def Finish(self):
        log_info("Finished!")



def FilterFrames(frame):

    return True


tray = I3Tray()



tray.Add('I3Reader', 'TheReader',
          FilenameList = [args.gcd] + files
          )


tray.AddModule(FilterFrames, "filter",
    streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics])




tray.AddModule(ReadScints, "ReadScints",
          )




tray.Execute()

