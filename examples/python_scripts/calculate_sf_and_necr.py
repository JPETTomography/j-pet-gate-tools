#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from numpy import *
from matplotlib import rcParams, rcdefaults, rc
from matplotlib.colors import LogNorm
from math import *
import os
import sys

import matplotlib.image as mpimg
from collections import OrderedDict

SLICE_WIDTH = 2. # in cm
BINS_DISPLACEMENTS = 100
BINS_ANGLES = 90

class sinogram:

  def __init__(self):
    self.displacements = []
    self.angles = []
    self.H = []
    self.dedges = []
    self.aedges = []

#===========================================

# TODO: correct the way pf calculating the index at the edge

def z2i (z, L):

  N = int(L/2.)
  SINOGRAM_WIDTH = L/N # in cm

  i = int(floor((z+L/2.)/SINOGRAM_WIDTH))
  if i==N: i=N-1
  return i

def mainfunction(geometry, activity, filepath, workdir):

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
  else:
      L = 50.

  #===========================================
  # Number of slices
  #===========================================

  N = int(L/SLICE_WIDTH)

  #===========================================
  # Loading coincidences from the file
  #===========================================

  coincidences = loadtxt(filepath + "_coincidences")
  #coincidences = genfromtxt(filepath + "_coincidences")
  posX1 = coincidences[:,0]
  posY1 = coincidences[:,1]
  posZ1 = coincidences[:,2]
  times1 = coincidences[:,3]
  posX2 = coincidences[:,4]
  posY2 = coincidences[:,5]
  posZ2 = coincidences[:,6]
  times2 = coincidences[:,7]
  IDs1 = coincidences[:,8]
  IDs2 = coincidences[:,9]
  toc = coincidences[:,12]

  #===========================================
  # Setting the time of the full simulation
  #===========================================

  time = 0.
  if os.path.exists(filepath + "_realtime"):
    time = loadtxt(filepath + "_realtime")
  else:
    time = max(max(times1),max(times2))/1e12

  #===========================================
  # Creating and filling the dictionary of sinograms
  #===========================================

  sinograms = dict() # dictionary of sinograms is created

  # Counters initialization
  N_true = 0
  N_psca = 0
  N_dsca = 0
  N_acci = 0
  N_all = 0

  for i in range(len(coincidences)):

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
    if cond1 and cond2 and displacement<12: # cut from the NEMA norm

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

  sinograms_rebinned = dict()

  indices = linspace(0,2*N-2,2*N-1)/2. # indices of the rebinned sinogram, they may be partial: 0, 0.5, 1, ...

  for i in indices:

    s = sinogram()
    Di = min(i, indices[-1]-i)
    H_sum = zeros(H.shape)

    for k in range(int(Di)+1):

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
  plt.savefig(workdir + "sinogram_final_" + activity + ".png")
  plt.clf()

  #===========================================
  # Summing lines of rebinned final sinogram
  #===========================================

  step = dedges[1]-dedges[0]
  vecsum = zeros(2*BINS_DISPLACEMENTS-1)

  for i in range(BINS_ANGLES):
    vec = H_final[:,i].T
    vec = list(vec)
    maxind = vec.index(max(vec))

    startindex = BINS_DISPLACEMENTS-maxind-1

    for j in range(len(vec)):
      vecsum[startindex+j] = vecsum[startindex+j]+vec[j]

  vecd = step*linspace(-BINS_DISPLACEMENTS+1, BINS_DISPLACEMENTS-1, 2*BINS_DISPLACEMENTS-1)

  ind_m20mm = BINS_DISPLACEMENTS - int(2./step)
  ind_p20mm = BINS_DISPLACEMENTS + int(2./step)
  val_m20mm = vecsum[ind_m20mm-1]
  val_p20mm = vecsum[ind_p20mm-1]

  thres = (val_m20mm+val_p20mm)/2.

  plt.subplots_adjust(left=0.15, right=0.96, top=0.97, bottom=0.15)
  plt.plot(vecd,vecsum/1000)
  plt.plot([-L,L],[thres/1000,thres/1000],color='r')
  plt.plot([-2,-2],[0,2*thres/1000],color='k')
  plt.plot([2,2],[0,2*thres/1000],color='k')
  plt.xlabel("Displacement [cm]")
  plt.ylabel("Kilo counts")

  plt.xlim(-12,12) # from NEMA cut
  plt.ylim(0,1.1*max(vecsum/1000))

  plt.savefig(workdir + "sumhist_"+activity+".png")
  plt.clf()
  plt.close()

  #===========================================

  # T - above the threshold line (treated as true)
  # S - below the threshold line (treated as scattered but includes also the accidental coincidences)

  T = 0 # true
  S = 0 # true (and accidental)
  for i in range(len(vecsum)):
    if vecsum[i]<thres and (vecd[i]<-2 or vecd[i]>2):
      S = S + vecsum[i]
    else:
      S = S + thres
      T = T + vecsum[i] - thres

  newT = T/(S+T)*N_all
  newS = S/(S+T)*N_all

  T = newT
  S = newS

  float_activity = float(activity)/22000.*1000 # in kBq/cc (activity is in MBq, volume in cc)

  SF_sin = S/(S+T)*100. # sin - sinogram
  SF_ctr = float(N_dsca+N_psca)/(N_dsca+N_psca+N_true)*100. # ctr - counter

  NECR_sin = T*T/(S+T)/time
  NECR_ctr = N_true*N_true/N_all/time

  NECR_sin = NECR_sin/1000. # cps -> kcps
  NECR_ctr = NECR_ctr/1000. # cps -> kcps

  # Printing calculated values
  data_for_single_activity = str(float_activity) + "\t"
  data_for_single_activity += str(SF_sin) + "\t" + str(SF_ctr) + "\t"
  data_for_single_activity += str(NECR_sin) + "\t" + str(NECR_ctr) + "\t"
  data_for_single_activity += str(T) + "\t" + str(S) + "\t"
  data_for_single_activity += str(N_true) + "\t" + str(N_dsca) + "\t" + str(N_psca) + "\t" + str(N_acci) + "\t"
  data_for_single_activity += str(time)

  print data_for_single_activity

  with open(workdir + "necr_dependency.txt", "a") as necr_dependency:
    necr_dependency.write(data_for_single_activity + '\n')

  #===========================================

  tim_diffs = []
  ang_diffs = []

  for i in range(len(coincidences)):
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
    #print tdiff, adiff

  #===========================================

  # Saving 2D histogram into the image

  H, xedges, yedges = histogram2d(tim_diffs, ang_diffs, bins=(BINS_DISPLACEMENTS,BINS_ANGLES), range=[[0, 3],[0, 180]])

  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111)
  #plt.subplots_adjust(left=0.20, right=0.9, top=0.9, bottom=0.1)
  ##plt.xlabel(r'$z_S$')
  ##  plt.xticks([0, 30, 60, 90, 120, 150, 180], ["0", "30", "60", "90", "120", "150", "180"])
  ##plt.ylabel(r'$x_S$')
  ##  plt.yticks([0, 30, 60, 90, 120, 150, 180], ["0", "30", "60", "90", "120", "150", "180"])
  ##  plt.text(100, 10, r'$E_{th}='+ene.split('_')[0] + '$' + '\n' + r'$\theta_{12} < \theta_{23} < \theta_{31}$')
  ##  VMAX_2g = 10**ceil(log10(H_12_23.max()))
  ##  if (VMAX_2g==1): VMAX_2g=10
  ##  VMAX_3g = 10**ceil(log10(H_12_23_3g.max()))
  ##  if (VMAX_3g==1): VMAX_3g=10
  ##plt.title("Map of efficiency\nlongitudinal (y=0)")
  VMAX = H.max()
  plt.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', norm=LogNorm(vmin=1, vmax=VMAX))
  plt.colorbar()
  ##  plt.imshow(H_12_23_3g.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', norm=LogNorm(vmin=1, vmax=VMAX_3g), cmap='summer_r', alpha=0.5)
  ##  plt.colorbar()
  #plt.savefig("./2_histogram_2D.png")
  ##plt.savefig("./"+name+"_longitudinal.pdf")

  param = 2.2
  xxx = linspace(0,param,100)
  yyy = []
  # Ellipse curve
  for i in range(len(xxx)):
    yyy.append(180-80*sqrt(1.-xxx[i]*xxx[i]/(param*param)))
  plt.plot(xxx,yyy,color='red')
  plt.xlabel("Time difference [ns]")
  plt.ylabel("Angle difference [deg.]")
  #plt.show()
  plt.savefig(workdir + "2D_differences_"+activity+".png")
  plt.clf()
  plt.close()

if __name__ == "__main__":

  #===========================================
  # The script assumes that the files with next activities are in the "directory"
  # and they have name matching the pattern: geometry_NECR_activity, for example
  # D85_1lay_L050_7mm_NECR_1000 (GOJA format)
  #===========================================

  activities = ["0001","0100","0200","0300","0400","0500","0600","0700","0800","0900","1000","1100","1200","1300","1400","1500","1600","1700","1800","1900","2000"]

  if (not os.path.isdir("./Results")): os.system("mkdir ./Results")
  if (not os.path.isdir("./Results/NECR")): os.system("mkdir ./Results/NECR")

  geometries = ["D85_1lay_L020_7mm", "D85_1lay_L050_7mm", "D85_1lay_L100_7mm", "D85_2lay_L020_7mm", "D85_2lay_L050_7mm", "D85_2lay_L100_7mm"]

  for geometry in geometries:

    print geometry

    workdir = "./Results/NECR/" + geometry + "/"
    if (not os.path.isdir(workdir)): os.system("mkdir " + workdir)
    if (not os.path.isdir(workdir + "rebinned")): os.system("mkdir " + workdir + "rebinned")
    if (os.path.isfile(workdir + "necr_dependency.txt")):
      os.system("rm " + workdir + "necr_dependency.txt")

    for i in range(len(activities)):

      activity = activities[i]

      directory = "/media/pkowalski/TOSHIBA EXT/NCBJ/GATE/NEMA/4_NECR/PMB_realtime/"

      filepath = directory + geometry + "_NECR_" + activity

      mainfunction(geometry, activity, filepath, workdir)
