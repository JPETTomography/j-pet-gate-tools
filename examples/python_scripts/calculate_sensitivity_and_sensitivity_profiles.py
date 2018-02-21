#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from numpy import *
from matplotlib import rcParams, rcdefaults, rc
from matplotlib.colors import LogNorm
import matplotlib.image as mpimg
import argparse, os, sys

def is_coincidences_directory_valid(coincidences_directory):
  #TODO
  return 1

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Calculate sensitivity and sensitivity profiles using the GOJA results.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('--coincidences-directory',
                      dest='coincidences_directory',
                      type=str,
                      default="/media/pkowalski/TOSHIBA EXT/NCBJ/GATE/NEMA/1_Sensitivity/PMB_realtime/",
                      help='path to dir with the GOJA sensitivity results')

  args = parser.parse_args()

  if not args.coincidences_directory:
    print "No directory with coincidences provided. Analysis cannot be performed. Check --help option."
    sys.exit()
  elif not os.path.isdir(args.coincidences_directory):
    print "Directory " + args.coincidences_directory + " is unavailable. Check --help option."
    sys.exit()
  elif not is_coincidences_directory_valid(args.coincidences_directory):
    print "Directory " + args.coincidences_directory + " is not valid. It should contain coincidences files with proper names. Check --help option."
    sys.exit()

  rcParams['font.size'] = 24
  rcParams['legend.fontsize'] = 24

  activity = 1000. # in kBq

  coincidences_directory = args.coincidences_directory

  geometries = ["D75_1lay_L020_7mm", "D75_1lay_L050_7mm", "D75_1lay_L100_7mm", "D75_2lay_L020_7mm", "D75_2lay_L050_7mm", "D75_2lay_L100_7mm", "D85_1lay_L020_7mm", "D85_1lay_L050_7mm", "D85_1lay_L100_7mm", "D85_2lay_L020_7mm", "D85_2lay_L050_7mm", "D85_2lay_L100_7mm", "D95_1lay_L020_7mm", "D95_1lay_L050_7mm", "D95_1lay_L100_7mm", "D95_2lay_L020_7mm", "D95_2lay_L050_7mm", "D95_2lay_L100_7mm"]

  for geometry in geometries:

    N = 0
    L = 0
    if "L020" in geometry:
        N = 20.
        L = 20.
    elif "L050" in geometry:
        N = 50.
        L = 50.
    elif "L100" in geometry:
        N = 100.
        L = 100.

    file_to_load = coincidences_directory + geometry + "_sensitivity_coincidences"
    print file_to_load
    tmp = loadtxt(file_to_load)

    toc = tmp[:,12]
    sourcePosZ = tmp[:,15]

    posX1 = tmp[:,0]
    posY1 = tmp[:,1]
    posZ1 = tmp[:,2]
    time1 = tmp[:,3]

    posY2 = tmp[:,4]
    posZ2 = tmp[:,5]
    posZ2 = tmp[:,6]
    time2 = tmp[:,7]

    #time = max([max(time1),max(time2)])/1e12 # time in seconds
    time = loadtxt(coincidences_directory + geometry + "_sensitivity_realtime")

    true = 0
    dsca = 0
    psca = 0
    acci = 0

    for t in toc:
     if t==1: true+=1
     elif t==2: dsca+=1
     elif t==3: psca+=1
     elif t==4: acci+=1

    sensitivity = true/time/activity

    #print geometry, sensitivity, time, float(acci)/(true+psca+dsca+acci)*100.

    #vec_true = zeros(int(N))
    #vec_all = zeros(int(N))

    #hist / ilosc emisji w hist w jednym binie (= aktywnosc/szerokosc binu)

    #tu mozna tez wstawic rozmycie

    #######################################
    #         BLURRING
    #######################################

    # przesuniecie jednego z hitow o wektor wyznaczony przez roznice czasow

    #######################################

    sourcePosZ_true = []
    sourcePosZ_true_PMT = []
    for i in range(len(sourcePosZ)):
      if toc[i]==1:
        sourcePosZ_true.append(sourcePosZ[i])

    b = linspace(-N/2,N/2,N+1)

    hist, bin_edges = histogram(sourcePosZ_true, bins=b)

    norm_factor = activity/N # activity per slice

    sensitivity_profile = hist/norm_factor/time

    #print sensitivity_profile
    #print (sensitivity_profile[int(N/2-1)]+sensitivity_profile[int(N/2)])/2.

    #slice_width = L/N
    #xxx = linspace(-N/2+slice_width/2,N/2-slice_width/2,N)
    #plt.plot(xxx,sensitivity_profile)

    #for i in range(len(sourcePosZ)):
      #index = int(sourcePosZ[i]+N/2.)
      #if index>=0 and index<=N-1:
        #if toc[i]==1: vec_true[index]+=1
        #vec_all[index]+=1

    #positions = linspace(-N/2+0.5,N/2-0.5,N)

    #fig = plt.figure(figsize=(8, 6))
    #ax = fig.add_subplot(111)
    #plt.subplots_adjust(left=0.2, right=0.95, top=0.9, bottom=0.15)
    #plt.xlabel('Axial position z [cm]')
    #plt.ylabel('Sensitivity [cps/kBq]')
    #plt.title(geometry)

    #sensitivities = vec_true/time/activity*N

    workdir_Results = "./Results/"
    if (not os.path.isdir(workdir_Results)): os.system("mkdir " + workdir_Results)
    workdir = workdir_Results + "Sensitivity/"
    if (not os.path.isdir(workdir)): os.system("mkdir " + workdir)

    savetxt(workdir + geometry + "_sensitivity.txt", [sensitivity])
    savetxt(workdir + geometry + "_sensitivity_profile.txt", sensitivity_profile)
