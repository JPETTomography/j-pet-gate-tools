#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import random as rnd
from numpy import *

# Physical constants:

c_scin = 12.6 # [cm/s]

# Technical constants:

## NEMA_DISPLACEMENT_CUT
#
#  Cut from the NEMA norm: when sinograms are analyzed, all LORs with
#  displacement larger than this threshold may be removed from the dataset.
NEMA_DISPLACEMENT_CUT = 12 # [cm]

## ELLIPSE_PARAM
#
#  Parameter of the ellipse used as a threshold for the 2nd level of the
#  2-level selection method (in space Da vs. Dt).
ELLIPSE_PARAM = 2.2

def ellipsoid_threshold(time_difference):

  return 180-80*sqrt(1.-time_difference*time_difference/(ELLIPSE_PARAM*ELLIPSE_PARAM))

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
strip200 = strip("L200",  200.,  0.435673,  0.365183) #TODO calculate values of crtPMT and crtSi (set ones are for the 100 cm strip)

strips = {
  20: strip020,
  50: strip050,
  100: strip100,
  200: strip200
}

geometries_sensitivity = [
  "D75_1lay_L020_7mm", "D85_1lay_L020_7mm", "D95_1lay_L020_7mm",
  "D75_2lay_L020_7mm", "D85_2lay_L020_7mm", "D95_2lay_L020_7mm",
  "D75_1lay_L050_7mm", "D85_1lay_L050_7mm", "D95_1lay_L050_7mm",
  "D75_2lay_L050_7mm", "D85_2lay_L050_7mm", "D95_2lay_L050_7mm",
  "D75_1lay_L100_7mm", "D85_1lay_L100_7mm", "D95_1lay_L100_7mm",
  "D75_2lay_L100_7mm", "D85_2lay_L100_7mm", "D95_2lay_L100_7mm",
  "D75_1lay_L200_7mm", "D85_1lay_L200_7mm", "D95_1lay_L200_7mm",
  "D75_2lay_L200_7mm", "D85_2lay_L200_7mm", "D95_2lay_L200_7mm"
]

geometries_NECR = [
  "D85_1lay_L020_7mm", "D85_2lay_L020_7mm",
  "D85_1lay_L050_7mm", "D85_2lay_L050_7mm",
  "D85_1lay_L100_7mm", "D85_2lay_L100_7mm",
  "D85_1lay_L200_7mm", "D85_2lay_L200_7mm",
]

activities_NECR = [
  "0001", "0100", "0200", "0300", "0400", "0500", "0600", "0700",
  "0800", "0900", "1000", "1100", "1200", "1300", "1400", "1500",
  "1600", "1700", "1800", "1900", "2000"
]

# Workdir directories:
workdir_Results = "./Results/"
workdir_Sensitivity = workdir_Results + "Sensitivity/"
workdir_NECR = workdir_Results + "NECR/"

def create_work_directories():

  if (not os.path.isdir(workdir_Results)): os.system("mkdir " + workdir_Results)
  if (not os.path.isdir(workdir_Sensitivity)): os.system("mkdir " + workdir_Sensitivity)
  if (not os.path.isdir(workdir_NECR)): os.system("mkdir " + workdir_NECR)
  for g in geometries_NECR:
    if (not os.path.isdir(workdir_NECR + "/" + g)): os.system("mkdir " + workdir_NECR + "/" + g)

# Gaussian distribution. sourcePosZ is the mean, and sigma is the standard deviation.
# L is the lenght of the scintillators.
def blur_sourcePosZ(sourcePosZ,  sigma,  L):

  result = rnd.gauss(sourcePosZ, sigma)
  #TODO:
  #while result < -L/2. or result > L/2.:
    #result = rnd.gauss(sourcePosZ, sigma)
  return result

def calculate_counters(tim_diffs, ang_diffs):

  counter_above = 0
  counter_below = 0
  for i in xrange(len(tim_diffs)):
    t = tim_diffs[i]
    a = ang_diffs[i]
    try:
      newa = ellipsoid_threshold(t)
      if a>newa: counter_above += 1
      else: counter_below += 1
    except:
      pass

  return [counter_above, counter_below]

def calculate_differences(tim1, tim2, posx1, posy1, posx2, posy2):

  tim_diffs = absolute(tim1 - tim2)/1e3 # in ns

  vx = posx1
  vy = posy1
  ux = posx2
  uy = posy2

  vu = vx*ux + vy*uy
  modv = sqrt(vx*vx + vy*vy)
  modu = sqrt(ux*ux + uy*uy)

  ang_diffs = arccos(vu/(modv*modu))/pi*180.

  return [tim_diffs, ang_diffs]

if __name__ == "__main__":

  # Testing:
  print(strips[20.].sigmaPMT)
  print(strips[100].fwhmSi)
