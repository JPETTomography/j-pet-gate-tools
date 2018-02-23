#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import re
from numpy import *
import string

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

  help_message = "Type goja_manager.py --help."

  parser = argparse.ArgumentParser(description='Analyze, verify and concatenate the GOJA results.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-gt', '--gate-output',
                      dest='path_gate_output',
                      type=str,
                      default="../output/",
                      help='path to dir with the splitted GATE results')

  parser.add_argument('-gj', '--goja-output',
                      dest='path_goja_output',
                      type=str,
                      default="./",
                      help='path to dir with the GOJA results')

  parser.add_argument('-m', '--mode',
                      dest='mode',
                      type=str,
                      default="analyze",
                      help='analyze, analyze-missing, verify or concatenate')

  parser.add_argument('--N0',
                      dest='N0',
                      type=int,
                      default=1000,
                      help='maximum number of events above the noise energy threshold in the coincidence window [for mode \'analyze\']')

  parser.add_argument('-r', '--run',
                      dest='type_of_run',
                      type=str,
                      default='locally',
                      help='run \'locally\' or \'on-cluster\' [for mode \'analyze\']')

  parser.add_argument('-sn', '--simulation-name',
                      dest='simulation_name',
                      type=str,
                      default="simulation",
                      help='name of the simulation [for mode \'concatenate\']')

  parser.add_argument('-c', '--clean',
                      action='store_true',
                      help='remove partial files aftre concatenation [for mode \'concatenate\']')

  args = parser.parse_args()

  if not os.path.isdir(args.path_gate_output):
    print "Directory " + args.path_gate_output + " does not exist. " + help_message
    sys.exit()

  if not os.path.isdir(args.path_goja_output):
    print "Directory " + args.path_goja_output + " does not exist. " + help_message
    sys.exit()

  if args.mode=="analyze":

    print "Analysis:"

    fnames = os.listdir(args.path_gate_output)
    fnames = [fname for fname in fnames if ".root" in fname]
    fnames = sorted(fnames, key=lambda x: (int(re.sub('\D','',x)),x))

    if args.type_of_run == 'locally':
      for fname in fnames:
        goja_command = "goja --root " + args.path_gate_output + fname \
                     + " --N0 " + str(args.N0) \
                     + " --save-real-time-to " + args.path_goja_output + fname[0:-5] + "_realtime" \
                     + " > " + args.path_goja_output + fname[0:-5] + "_coincidences &"
        os.system(goja_command)
        print goja_command

    elif args.type_of_run == 'on-cluster':
      basename = fnames[0][0:-5]
      basepath_goja = (args.path_goja_output + basename).rstrip(string.digits)
      basepath_gate = (args.path_gate_output + basename).rstrip(string.digits)
      # generate array.pbs:
      os.system('echo \'#!/bin/sh\' > array.pbs')
      os.system('echo \'#PBS -q i3d\' >> array.pbs')
      os.system('echo \'#PBS -l nodes=1:ppn=1\' >> array.pbs')
      os.system('echo \'#PBS -N GOJA\' >> array.pbs')
      os.system('echo \'#PBS -V\' >> array.pbs')
      os.system('echo \'cd ${PBS_O_WORKDIR}\' >> array.pbs')
      goja_command = "goja --root " + basepath_gate + "${PBS_ARRAYID}" + ".root" \
                   + " --N0 " + str(args.N0) \
                   + " --save-real-time-to " + basepath_goja + "${PBS_ARRAYID}" + "_realtime" \
                   + " > " + basepath_goja + "${PBS_ARRAYID}" + "_coincidences"
      os.system('echo \'' + goja_command + '\' >> array.pbs')
      os.system('echo \'exit 0;\' >> array.pbs')
      # push into queue:
      qsub_command = 'qsub -t 1-' + str(len(fnames)) + ' array.pbs'
      os.system(qsub_command)
      # remove array.pbs:
      os.system('rm array.pbs')

    else:
      print "Inproper type of run. " + help_message

  elif args.mode=="analyze-missing":

    print "Analysis of missing files:"

    fnames = os.listdir(args.path_gate_output)
    fnames = [fname for fname in fnames if ".root" in fname]
    fnames = sorted(fnames, key=lambda x: (int(re.sub('\D','',x)),x))

    for fname in fnames:
      path_coincidences = args.path_goja_output + fname[:-5] + "_coincidences"
      path_realtime = args.path_goja_output + fname[:-5] + "_realtime"
      if not os.path.isfile(path_coincidences) or not os.path.isfile(path_realtime):
        basepath_goja = args.path_goja_output + fname[:-5]
        goja_command = "goja --root " + args.path_gate_output + fname \
                       + " --N0 " + str(args.N0) \
                       + " --save-real-time-to " + basepath_goja + "_realtime" \
                       + " > " + basepath_goja + "_coincidences &"
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
        basepath_goja = args.path_goja_output + basename

        realtime += loadtxt(basepath_goja + "_realtime")
        os.system("cat " + basepath_goja + "_coincidences >> " + path_coincidences)

      savetxt(path_realtime, [realtime])

      if args.clean:
        for fname in fnames:
          basename = fname[0:-13]
          basepath_goja = args.path_goja_output + basename
          os.system("rm " + basepath_goja + "_coincidences")
          os.system("rm " + basepath_goja + "_realtime")

      print "Goja output succesfully concatenated."

  else:

    print 'Inproper mode. ' + help_message
