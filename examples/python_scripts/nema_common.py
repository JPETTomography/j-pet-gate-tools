#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os,  random
from numpy  import sqrt

# Physical constants:
c_scin = 12.6 # [cm/s]

class strip:

  def __init__(self,  name,  L,  crtPMT,  crtSi):

    self.name = name
    self.L = L

    self.crtPMT = crtPMT
    self.crtSi = crtSi
    self.crtWLS = crtSi

    self.fwhmPMT = crtPMT*c_scin/sqrt(2)
    self.sigmaPMT = self.fwhmPMT/2.35

    self.fwhmSi = crtSi*c_scin/sqrt(2)
    self.sigmaSi = self.fwhmSi/2.35

    self.sigmaWLS = 0.5
    self.fwhmWLS = self.sigmaWLS*2.35

strip020 = strip("L020",  20.,  0.248695,  0.178523)
strip050 = strip("L050",  50.,  0.314308,  0.235196)
strip100 = strip("L100",  100.,  0.435673,  0.365183)

strips = {
  20: strip020,
  50: strip050,
  100: strip100,
}

geometries_Sensitivity = ["D75_1lay_L020_7mm", "D75_1lay_L050_7mm", "D75_1lay_L100_7mm",
                          "D75_2lay_L020_7mm", "D75_2lay_L050_7mm", "D75_2lay_L100_7mm",
                          "D85_1lay_L020_7mm", "D85_1lay_L050_7mm", "D85_1lay_L100_7mm",
                          "D85_2lay_L020_7mm", "D85_2lay_L050_7mm", "D85_2lay_L100_7mm",
                          "D95_1lay_L020_7mm", "D95_1lay_L050_7mm", "D95_1lay_L100_7mm",
                          "D95_2lay_L020_7mm", "D95_2lay_L050_7mm", "D95_2lay_L100_7mm"]

geometries_NECR = ["D85_1lay_L020_7mm", "D85_1lay_L050_7mm", "D85_1lay_L100_7mm",
                   "D85_2lay_L020_7mm", "D85_2lay_L050_7mm", "D85_2lay_L100_7mm"]

# Workdir directories:
workdir_Results = "./Results/"
workdir_Sensitivity = workdir_Results + "Sensitivity/"
workdir_NECR = workdir_Results + "NECR/"

def prepare_directories():

  if (not os.path.isdir(workdir_Results)): os.system("mkdir " + workdir_Results)
  if (not os.path.isdir(workdir_Sensitivity)): os.system("mkdir " + workdir_Sensitivity)
  if (not os.path.isdir(workdir_NECR)): os.system("mkdir " + workdir_NECR)

#Gaussian distribution. sourcePosZ is the mean, and sigma is the standard deviation.
# L is the lenght of the scintillators.
def blur_sourcePosZ(sourcePosZ,  sigma,  L):

  result = random.gauss(sourcePosZ, sigma)
  #TODO:
  #while result < -L/2. or result > L/2.:
    #result = random.gauss(sourcePosZ, sigma)
  return result

if __name__ == "__main__":

  # Testing:
  print(strips[20.].sigmaPMT)
  print(strips[100].fwhmSi)
