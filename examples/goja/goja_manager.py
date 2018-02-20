#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from numpy import *
from matplotlib import rcParams, rcdefaults, rc
from matplotlib.colors import LogNorm

import os
import argparse

import matplotlib.image as mpimg

import re
import sys

def verify_goja_output(gate_path, goja_path):
    nr_of_missing_files = 0
    for fname in os.listdir(gate_path):
      if ".root" in fname:
        path_coincidences = goja_path + fname[:-5] + "_coincidences"
        if not os.path.isfile(path_coincidences):
          print "File ", path_coincidences, " is missing."
          nr_of_missing_files += 1
        path_realtime = goja_path + fname[:-5] + "_realtime"
        if not os.path.isfile(path_realtime):
          print "File ", path_realtime, " is missing."
          nr_of_missing_files += 1
    print "Number of missing files: ", nr_of_missing_files
    return nr_of_missing_files

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Analyze, verify and concatenate the GOJA results.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('--gate-output',
                      dest='path_gate_output',
                      type=str,
                      default="../output/",
                      help='path to dir with the splitted GATE results')

  parser.add_argument('--goja-output',
                      dest='path_goja_output',
                      type=str,
                      default="./",
                      help='path to dir with the GOJA results')

  parser.add_argument('--mode',
                      dest='mode',
                      type=str,
                      default="analyze",
                      help='analyze, verify or concatenate')

  parser.add_argument('--simulation-name',
                      dest='simulation_name',
                      type=str,
                      default="simulation",
                      help='name of the simulation')

  args = parser.parse_args()

  if not os.path.isdir(args.path_gate_output):
    print "Directory " + args.path_gate_output + " does not exist. Type goja_manager.py --help."
    sys.exit()

  if not os.path.isdir(args.path_goja_output):
    print "Directory " + args.path_goja_output + " does not exist. Type goja_manager.py --help."
    sys.exit()

  if args.mode=="analyze":

    print "Analysis:"

    nr_of_root_files = 0

    fnames = os.listdir(args.path_gate_output)
    fnames = [fname for fname in fnames if ".root" in fname]
    fnames = sorted(fnames, key=lambda x: (int(re.sub('\D','',x)),x))

    for fname in fnames:
      nr_of_root_files += 1
      basename = fname[:-5]
      basepath = args.path_goja_output + basename
      goja_command = "goja --root " + args.path_gate_output + fname \
                     + " --N0 1000" \
                     + " --save-real-time-to " + basepath + "_realtime" \
                     + " > " + basepath + "_coincidences &"
      print goja_command
      os.system(goja_command)

  elif args.mode=="verify":

    print "Verification:"

    if not verify_goja_output(args.path_gate_output, args.path_goja_output):
      print "Goja output is ok."

  elif args.mode=="concatenate":

    print "Concatenation:"

    nr_of_missing_files = verify_goja_output(args.path_gate_output, args.path_goja_output)

    if nr_of_missing_files:
      print "Concatenation cannot be performed."

    else:
      path_coincidences = args.path_goja_output + args.simulation_name + "_COINCIDENCES"
      if os.path.isfile(path_coincidences):
        os.system('rm ' + path_coincidences)
      path_realtime = args.path_goja_output + args.simulation_name + "_REALTIME"
      if os.path.isfile(path_realtime):
        os.system('rm ' + path_realtime)

      realtime = 0.

      fnames = os.listdir(args.path_goja_output)
      fnames = [fname for fname in fnames if "_coincidences" in fname]
      fnames = sorted(fnames, key=lambda x: (int(re.sub('\D','',x)),x))

      for fname in fnames:

        basename = fname[0:-13]
        basepath = args.path_goja_output + basename

        realtime += loadtxt(basepath + "_realtime")
        os.system("cat " + basepath + "_coincidences >> " + path_coincidences)

      savetxt(path_realtime, [realtime])

      print "Goja output succesfully concatenated."

  else:

    print 'Inproper mode. Type goja_manager.py --help.'
