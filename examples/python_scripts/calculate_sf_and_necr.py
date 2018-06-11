#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from numpy import *
from matplotlib import rcParams
from matplotlib.colors import LogNorm
from math import *
import argparse
import os
import sys

from nema_common import *

SLICE_WIDTH = 2. # in cm
BINS_DISPLACEMENTS = 100
BINS_ANGLES = 90

SUFFIX_COINCIDENCES = "_NECR_COINCIDENCES_short"
SUFFIX_REALTIME = "_NECR_REALTIME_short"
OUTPUT_FORMAT = ""

CALCULATE_SF = False

class sinogram:

  def __init__(self):
    self.displacements = []
    self.angles = []
    self.H = []
    self.dedges = []
    self.aedges = []

## The function discretizes the position of the hit along the strip.
#
#  @param    z   position of the hit along the strip [cm]
#  @param    L   legnth of the strip [cm]
#
#  @return   i   index of the virtual slice (along the strip)
def z2i (z, L):

  N = int(L/2.)
  SINOGRAM_WIDTH = L/N # in cm

  i = int(floor((z+L/2.)/SINOGRAM_WIDTH))
  if i==N:
    i=N-1 # TODO: correct the way pf calculating the index at the edge
  return i

## The function performs rebinning using the SSRB algorithm.
#
#  @param    sinograms   dictionary of all oblique sinograms
#  @param    N           number of virtual slices
#  @param    H_shape     shape of the hisotgrammed sinogram
#
#  @return   val         dictionary of all rebinned sinograms
def perform_rebinning(sinograms, N, H_shape):

  sinograms_rebinned = dict()

  indices = linspace(0,2*N-2,2*N-1)/2. # indices of the rebinned sinogram, they may be partial: 0, 0.5, 1, ...

  for i in indices:

    s = sinogram()
    Di = min(i, indices[-1]-i)
    H_sum = zeros(H_shape)

    for k in xrange(int(Di)+1):

      current_slice = int(indices[int(i)]*2.)

      couple = []
      if Di == int(Di):
        couple = sort([current_slice-k, current_slice+k])
      else:
        couple = sort([current_slice-k, current_slice+k+1])

      if ((couple[0],couple[1]) in sinograms.keys()): # checking if oblique sinogram exists in the set of all sinograms
        H_sum = H_sum + sinograms[(couple[0],couple[1])].H

    s.H = H_sum
    sinograms_rebinned[i] = s

  return sinograms_rebinned

def perform_analysis(activity, filepath, workdir):

  #===========================================
  # Matplotlib plotting parameters
  #===========================================

  rcParams['font.size'] = 24
  rcParams['legend.fontsize'] = 18
  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111)

  #===========================================
  # Length of the scintillators
  #===========================================

  L = 0.
  if "L020" in filepath:
      L = 20.
  elif "L050" in filepath:
      L = 50.
  elif "L100" in filepath:
      L = 100.
  elif "L200" in filepath:
      L = 200.
  else:
      L = 50.

  #===========================================
  # Number of slices
  #===========================================

  N = int(L/SLICE_WIDTH)

  #===========================================
  # Loading coincidences from the file
  #===========================================

  coincidences = loadtxt(filepath + SUFFIX_COINCIDENCES)
  posX1 = coincidences[:,0]
  posY1 = coincidences[:,1]
  posZ1 = coincidences[:,2]
  times1 = coincidences[:,3]
  posX2 = coincidences[:,4]
  posY2 = coincidences[:,5]
  posZ2 = coincidences[:,6]
  times2 = coincidences[:,7]
  toc = coincidences[:,12]

  #===========================================
  # Setting the time of the full simulation
  #===========================================

  time = 0.
  realtime_path = filepath + SUFFIX_REALTIME
  if os.path.exists(realtime_path):
    time = loadtxt(realtime_path)
  else:
    time = max(max(times1),max(times2))/1e12

  #===========================================
  # Creating and filling the dictionary of sinograms
  #===========================================

  sinograms = dict()

  # Counters initialization
  N_true = 0
  N_psca = 0
  N_dsca = 0
  N_acci = 0
  N_all = 0

  for i in xrange(len(coincidences)):

    # Displacement - distance between LOR and (0,0,0) in XY crossection
    displacement = 0.
    denominator = sqrt((posY2[i]-posY1[i])**2 + (posX2[i]-posX1[i])**2)
    if denominator != 0:
      displacement = (posX2[i]*posY1[i]-posY2[i]*posX1[i]) / denominator

    # Angle
    angle = 0.
    if (posX2[i]-posX1[i])!=0:
      angle = atan((posY1[i]-posY2[i])/(posX2[i]-posX1[i]))

    if displacement>0:
      angle=angle+pi/2.
    else:
      angle=angle+3.*pi/2.
    if angle > pi:
      angle = angle - pi
      displacement = -displacement

    cond1 = displacement != 0. or angle != 0.
    cond2 = (not isnan(displacement)) and (not isnan(angle))

    if cond1 and cond2 and displacement<NEMA_DISPLACEMENT_CUT:

      # Counters incrementation
      if toc[i]==1: N_true = N_true + 1
      elif toc[i]==2: N_psca = N_psca + 1
      elif toc[i]==3: N_dsca = N_dsca + 1
      elif toc[i]==4: N_acci = N_acci + 1
      N_all = N_all + 1

      indices = [z2i(posZ1[i],L), z2i(posZ2[i],L)]
      indices.sort()

      # the dictionary is filled using sorted pair of slice indexes, the event is attached to its sinogram
      if (indices[0],indices[1]) in sinograms:
        sinograms[(indices[0],indices[1])].displacements.append(displacement)
        sinograms[(indices[0],indices[1])].angles.append(angle)

      else:
        # if there is no such a sinogram, it is created
        s = sinogram()
        sinograms[(indices[0],indices[1])] = s
        sinograms[(indices[0],indices[1])].displacements.append(displacement)
        sinograms[(indices[0],indices[1])].angles.append(angle)

  #===========================================
  # Converting sinograms into 2D histograms
  #===========================================

  BINS = [BINS_DISPLACEMENTS, BINS_ANGLES]
  RANGE = [[0, 12], [0, pi]]

  for s in sinograms:

    H, dedges, aedges = histogram2d(sinograms[s].displacements, sinograms[s].angles, bins=BINS, range=RANGE)

    sinograms[s].H = H
    sinograms[s].dedges = dedges
    sinograms[s].aedges = aedges

  #===========================================
  # Rebinning (SSRB)
  #===========================================

  sinograms_rebinned = perform_rebinning(sinograms, N, H.shape)

  #===========================================
  # Summing rebinned sinograms
  #===========================================

  H_final = zeros(H.shape)
  for key in sinograms_rebinned:
    H_final = H_final + sinograms_rebinned[key].H

  #===========================================
  # Plotting rebinned sinograms
  #===========================================

  plt.subplots_adjust(left=0.15, right=0.96, top=0.97, bottom=0.15)
  plt.imshow(H_final.T, interpolation='nearest', origin='low', extent=[dedges[0], dedges[-1], aedges[0], aedges[-1]], aspect='auto')
  plt.colorbar()
  plt.xlabel("Displacement [cm]")
  plt.ylabel("Angle [rad.]")
  plt.savefig(workdir + "sinogram_final_" + activity + OUTPUT_FORMAT)
  plt.clf()

  #===========================================
  # Summing lines of rebinned final sinogram
  #===========================================

  step = dedges[1]-dedges[0]
  vecsum = zeros(2*BINS_DISPLACEMENTS-1)

  for i in xrange(BINS_ANGLES):
    vec = H_final[:,i].T
    vec = list(vec)
    maxind = vec.index(max(vec))

    startindex = BINS_DISPLACEMENTS-maxind-1

    for j in xrange(len(vec)):
      vecsum[startindex+j] = vecsum[startindex+j]+vec[j]

  veci = linspace(-BINS_DISPLACEMENTS+1, BINS_DISPLACEMENTS-1, 2*BINS_DISPLACEMENTS-1)
  vecd = step*veci

  #===========================================

  ind_m20mm = BINS_DISPLACEMENTS - int(2./step) - 1
  ind_p20mm = BINS_DISPLACEMENTS + int(2./step) - 1
  val_m20mm = vecsum[ind_m20mm]
  val_p20mm = vecsum[ind_p20mm]

  # y = a*x + b
  a = float(val_m20mm-val_p20mm)/float(ind_m20mm-ind_p20mm)
  b = val_m20mm - (float(val_m20mm-val_p20mm))/(float(ind_m20mm-ind_p20mm))*float(ind_m20mm)
  thres = a*(veci+BINS_DISPLACEMENTS-1)+b

  # T - above the threshold line (treated as true)
  # S - below the threshold line (treated as scattered but includes also the accidental coincidences)

  T = 0 # true
  S = 0 # scattered (and accidental)
  for i in xrange(len(vecsum)):
    if vecd[i]<-2 or vecd[i]>2:
      S = S + vecsum[i]
    else:
      if vecsum[i]<thres[i]:
        S = S + vecsum[i]
      else:
        S = S + thres[i]
        T = T + vecsum[i] - thres[i]

  T = T/(S+T)*N_all
  S = S/(S+T)*N_all

  #===========================================

  plt.subplots_adjust(left=0.15, right=0.96, top=0.97, bottom=0.15)
  plt.plot(vecd,vecsum/1000.)
  plt.plot([vecd[ind_m20mm],vecd[ind_p20mm]],[thres[ind_m20mm]/1000.,thres[ind_p20mm]/1000.],color='r')
  plt.plot([-2,-2],[0,2*thres[ind_m20mm]/1000.],color='k')
  plt.plot([2,2],[0,2*thres[ind_p20mm]/1000.],color='k')

  plt.xlabel("Displacement [cm]")
  plt.ylabel("Kilo counts")

  plt.xlim(-12,12) # from NEMA cut
  plt.ylim(0,1.1*max(vecsum/1000.))

  plt.savefig(workdir + "sumhist_" + activity + OUTPUT_FORMAT)
  plt.clf()
  plt.close()

  #===========================================

  if not CALCULATE_SF:
    float_activity = float(activity)/22000.*1000 # in kBq/cc (activity is in MBq, volume in cc)
  else:
    float_activity = 0.001/22000.*1000 # in kBq/cc (activity is in MBq, volume in cc)

  SF_sin = S/(S+T)*100. # sin - sinogram
  SF_ctr = float(N_dsca+N_psca)/(N_dsca+N_psca+N_true)*100. # ctr - counter

  NECR_sin = T*T/(S+T)/time
  NECR_ctr = N_true*N_true/N_all/time

  NECR_sin = NECR_sin/1000. # cps -> kcps
  NECR_ctr = NECR_ctr/1000. # cps -> kcps

  ratio_acci = N_acci/float(N_all)

  # Printing calculated values
  data_for_single_activity = str(float_activity) + "\t"
  data_for_single_activity += str(SF_sin) + "\t" + str(SF_ctr) + "\t"
  data_for_single_activity += str(NECR_sin) + "\t" + str(NECR_ctr) + "\t"
  data_for_single_activity += str(T) + "\t" + str(S) + "\t"
  data_for_single_activity += str(N_true) + "\t" + str(N_dsca) + "\t" + str(N_psca) + "\t" + str(N_acci) + "\t"
  data_for_single_activity += str(time) + "\t"
  data_for_single_activity += str(ratio_acci)

  print data_for_single_activity

  with open(workdir + "necr_dependency.txt", "a") as necr_dependency:
    necr_dependency.write(data_for_single_activity + '\n')

  #===========================================

  tim_diffs = []
  ang_diffs = []

  for i in xrange(len(coincidences)):
    tdiff = abs(times1[i]-times2[i])/1e3 # in ns TODO abs ??????????????
    tim_diffs.append(tdiff)
    #v=[vx,vy] = posX1, posY1
    #u=[ux,uy] = posX2, posY2
    vx = posX1[i]
    vy = posY1[i]
    ux = posX2[i]
    uy = posY2[i]
    vu = posX1[i]*posX2[i] + posY1[i]*posY2[i]
    modv = sqrt(posX1[i]*posX1[i] + posY1[i]*posY1[i])
    modu = sqrt(posX2[i]*posX2[i] + posY2[i]*posY2[i])
    try:
      a = acos(vu/(modv*modu))/pi*180.
      adiff = a
    except:
      pass
    ang_diffs.append(adiff)

  #===========================================

  # Saving 2D histogram into the image

  H, xedges, yedges = histogram2d(tim_diffs, ang_diffs, bins=(BINS_DISPLACEMENTS,BINS_ANGLES), range=[[0, 3],[0, 180]])

  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111)
  plt.subplots_adjust(left=0.20, right=0.9, top=0.9, bottom=0.1)

  VMAX = H.max()
  plt.imshow(H.T, interpolation='None', origin='low',
             extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
             aspect='auto', norm=LogNorm(vmin=1, vmax=VMAX))
  plt.colorbar()

  param = 2.2
  xxx = linspace(0,param,100)
  yyy = []
  # Ellipse curve
  for i in xrange(len(xxx)):
    yyy.append(180-80*sqrt(1.-xxx[i]*xxx[i]/(param*param)))
  plt.plot(xxx,yyy,color='red')
  plt.xlabel("Time difference [ns]")
  plt.ylabel("Angle difference [deg.]")
  #plt.show()
  plt.ylim(90,180)
  plt.savefig(workdir + "2D_differences_" + activity + OUTPUT_FORMAT, bbox_inches='tight')
  plt.clf()
  plt.close()

def is_coincidences_directory_valid(coincidences_directory):
  #TODO
  return True

if __name__ == "__main__":

  #===========================================
  # The script assumes that the files with next activities are in the "directory"
  # and they have name matching the pattern: geometry_NECR_activity, for example
  # D85_1lay_L050_7mm_NECR_1000 (GOJA format)
  #===========================================

  parser = argparse.ArgumentParser(description='Calculate sensitivity and sensitivity profiles using the GOJA results.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-cd', '--coincidences-directory',
                      dest='coincidences_directory',
                      type=str,
                      default="/media/pkowalski/TOSHIBA EXT/NCBJ/GATE/NEMA/4_NECR/N0_1000_short/",
                      help='path to dir with the GOJA sensitivity results')

  parser.add_argument('-sw', '--slice-width',
                      dest='slice_width',
                      type=float,
                      default=SLICE_WIDTH,
                      help='width of the virtual slice')

  parser.add_argument('-bd', '--bins-displacements',
                      dest='bins_displacements',
                      type=float,
                      default=BINS_DISPLACEMENTS,
                      help='nr of bins for displacements')

  parser.add_argument('-ba' ,'--bins-angles',
                      dest='bins_angles',
                      type=float,
                      default=BINS_ANGLES,
                      help='nr of bins for angles')

  parser.add_argument('--scatter-fraction',
                      dest='scatter_fraction',
                      action='store_true',
                      help='use to calculate the SF using the set of 1 kBq simulations')

  parser.add_argument('-of', '--outputformat',
                      dest='outputformat',
                      type=str,
                      default="png",
                      help='output format of images')

  args = parser.parse_args()

  OUTPUT_FORMAT = "." + args.outputformat

  if not args.coincidences_directory:
    print "No directory with coincidences provided. Analysis cannot be performed. Check --help option."
    sys.exit()
  elif not os.path.isdir(args.coincidences_directory):
    print "Directory " + args.coincidences_directory + " is unavailable. Check --help option."
    sys.exit()
  elif not is_coincidences_directory_valid(args.coincidences_directory):
    print "Directory " + args.coincidences_directory + " is not valid. It should contain coincidences files with proper names. Check --help option."
    sys.exit()

  if args.scatter_fraction:
    CALCULATE_SF = True
    SUFFIX_COINCIDENCES = "SF_COINCIDENCES_short"
    SUFFIX_REALTIME = "SF_REALTIME_short"
    activities_NECR = [""]

  SLICE_WIDTH = args.slice_width
  BINS_DISPLACEMENTS = int(args.bins_displacements)
  BINS_ANGLES = int(args.bins_angles)

  directory = args.coincidences_directory

  create_work_directories()

  for geometry in geometries_NECR:

    print geometry

    workdir = workdir_NECR + geometry + "/"
    if (not os.path.isdir(workdir)):
      os.system("mkdir " + workdir)
    if (os.path.isfile(workdir + "necr_dependency.txt")):
      os.system("rm " + workdir + "necr_dependency.txt")

    for i in xrange(len(activities_NECR)):

      activity = activities_NECR[i]
      filepath = directory + geometry + "_" + activity
      perform_analysis(activity, filepath, workdir)
