#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from numpy import *

from nema_common import *

activity = 1000. # in kBq

def is_coincidences_directory_valid(coincidences_directory):

  result = True
  for geometry in geometries_Sensitivity:
    coincidences_file = coincidences_directory + geometry + "_sensitivity_coincidences";
    if not os.path.isfile(coincidences_file):
      print("File " + coincidences_file + " is missing.")
      result = False
    realtime_file = coincidences_directory + geometry + "_sensitivity_realtime";
    if not os.path.isfile(realtime_file):
      print("File " + realtime_file + " is missing.")
      result = False
  return result

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Calculate sensitivity and sensitivity profiles using the GOJA results.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('--coincidences-directory',
                      dest='coincidences_directory',
                      type=str,
                      default="/home/pkowalski/Pulpit/PMB_realtime_sensitivity/",
                      help='path to dir with the GOJA sensitivity results')

  args = parser.parse_args()

  if not args.coincidences_directory:
    print("No directory with coincidences provided. Analysis cannot be performed. Check --help option.")
    sys.exit()
  elif not os.path.isdir(args.coincidences_directory):
    print("Directory " + args.coincidences_directory + " is unavailable. Check --help option.")
    sys.exit()
  elif not is_coincidences_directory_valid(args.coincidences_directory):
    print("Directory " + args.coincidences_directory + " is not valid. It should contain coincidences files with proper names. Check --help option.")
    sys.exit()

  coincidences_directory = args.coincidences_directory

  for geometry in geometries_Sensitivity:

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

    norm_factor = activity/N # activity per slice

    file_to_load = coincidences_directory + geometry + "_sensitivity_coincidences"
    print(file_to_load)
    tmp = loadtxt(file_to_load)

    toc = tmp[:,12]
    sourcePosZ = tmp[:,15]
    time = loadtxt(coincidences_directory + geometry + "_sensitivity_realtime")

    true = 0
    for t in toc:
      if t==1: true+=1

    sensitivity = true/time/activity

    sourcePosZ_true = []
    sourcePosZ_true_PMT = []

    for i in range(len(sourcePosZ)):
      if toc[i]==1:
        sourcePosZ_true.append(sourcePosZ[i])
        posPMT = blur_sourcePosZ(sourcePosZ[i], strips[L].sigmaPMT, L)
        sourcePosZ_true_PMT.append(posPMT)

    pos_bins = linspace(-N/2,N/2,N+1)

    hist = histogram(sourcePosZ_true, bins=pos_bins)
    hist_PMT = histogram(sourcePosZ_true_PMT, bins=pos_bins)

    sensitivity_profile = hist[0]/norm_factor/time
    sensitivity_profile_PMT = hist_PMT[0]/norm_factor/time

    prepare_directories()

    savetxt(workdir_Sensitivity + geometry + "_sensitivity.txt", [sensitivity])
    savetxt(workdir_Sensitivity + geometry + "_sensitivity_profiles.txt",
            array([sensitivity_profile, sensitivity_profile_PMT]).T)
