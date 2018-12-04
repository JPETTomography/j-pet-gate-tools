#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import errno
import argparse
import re
from numpy import *
import string

def verify_gate_output(gate_path):

  nr_of_missing_files = 0
  missing_files = []
  nr_of_splits = 0

  with open(gate_path + '../Gate_parallel.sh', 'r') as gate_parallel_sh:
    for line in gate_parallel_sh:
      if line.split('=')[0] == 'NR_OF_SPLITS':
        nr_of_splits = int((line.split('=')[1]).replace('\'', '').replace('\n', ''))

  for s in xrange(nr_of_splits):
    output_root = gate_path + 'output' + str(s+1) + '.root'
    if not os.path.isfile(output_root):
      print "\tFile ", output_root, " is missing."
      nr_of_missing_files += 1
      missing_files.append(output_root)

  if nr_of_missing_files == 0:
    print "\tGATE output is ok." #TODO check if files where closed properly

  return nr_of_missing_files, missing_files

def verify_goja_output(gate_path, goja_path):

    nr_of_missing_files, missing_files = verify_gate_output(gate_path)

    #TODO add filling missing goja files

    for fname in os.listdir(gate_path):
      if ".root" in fname:
        path_coincidences = goja_path + fname[:-5] + "_coincidences"
        if not os.path.isfile(path_coincidences):
          print "\tFile ", path_coincidences, " is missing."
          nr_of_missing_files += 1
        path_realtime = goja_path + fname[:-5] + "_realtime"
        if not os.path.isfile(path_realtime):
          print "\tFile ", path_realtime, " is missing."
          nr_of_missing_files += 1
        path_statistics = goja_path + fname[:-5] + "_statistics"
        if not os.path.isfile(path_statistics):
          print "\tFile ", path_statistics, " is missing."
          nr_of_missing_files += 1
    if nr_of_missing_files>0:
      print "\tNumber of missing files: ", nr_of_missing_files
    return nr_of_missing_files

def concatenate_files(fnames):

      path_coincidences = args.path_goja_output + args.simulation_name + "_COINCIDENCES"
      try:
        os.unlink(path_coincidences)
      except OSError as e:
        if e.errno == errno.ENOENT:
          pass
        else:
          raise e
      path_realtime = args.path_goja_output + args.simulation_name + "_REALTIME"
      try:
        os.unlink(path_realtime)
      except OSError as e:
        if e.errno == errno.ENOENT:
          pass
        else:
          raise e
      path_statistics = args.path_goja_output + args.simulation_name + "_STATISTICS"
      try:
        os.unlink(path_statistics)
      except OSError as e:
        if e.errno == errno.ENOENT:
          pass
        else:
          raise e

      realtime = 0.
      counter_all_compton_hits = 0
      counter_compton_hits_over_the_ETH0 = 0
      counter_compton_hits_over_the_ETH = 0

      with open(path_coincidences, 'w') as outfile:
        for fname in fnames:
          basename = fname[0:-13]
          basepath_goja = args.path_goja_output + basename

          realtime += loadtxt(basepath_goja + "_realtime")

          with open(basepath_goja + "_coincidences") as infile:
            for line in infile:
              outfile.write(line)

          counters = genfromtxt(basepath_goja + "_statistics", usecols=(0))
          counter_all_compton_hits += counters[0]
          counter_compton_hits_over_the_ETH0 += counters[1]
          counter_compton_hits_over_the_ETH += counters[2]

      savetxt(path_realtime, [realtime])

      with open(path_statistics, 'w') as stats:
        stats.write(str(int(counter_all_compton_hits)) + " # all Compton hits\n")
        stats.write(str(int(counter_compton_hits_over_the_ETH0)) + " # compton hits with edep over the ETH0\n")
        stats.write(str(int(counter_compton_hits_over_the_ETH)) + " # compton hits with edep over the ETH\n")

      if args.clean:
        for fname in fnames:
          basename = fname[0:-13]
          basepath_goja = args.path_goja_output + basename
          os.unlink(basepath_goja + "_coincidences")
          os.unlink(basepath_goja + "_realtime")
          os.unlink(basepath_goja + "_statistics")
      print "Goja output succesfully concatenated."

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
                      help='analyze, analyze-missing, verify, concatenate or concatenate-force')

  parser.add_argument('--eth0',
                      dest='eth0',
                      type=float,
                      default=0.01,
                      help='noise energy threshold in MeV [for mode \'analyze\']')

  parser.add_argument('--tw',
                      dest='tw',
                      type=float,
                      default=3,
                      help='time window in ns [for mode \'analyze\']')

  parser.add_argument('--N0',
                      dest='N0',
                      type=int,
                      default=1000,
                      help='maximum number of events above the noise energy threshold in the coincidence window [for mode \'analyze\']')

  parser.add_argument('-r', '--run',
                      dest='type_of_run',
                      type=str,
                      default='on-cluster',
                      help='run \'locally\' or \'on-cluster\' [for mode \'analyze\']')

  parser.add_argument('-sn', '--simulation-name',
                      dest='simulation_name',
                      type=str,
                      default="simulation",
                      help='name of the simulation [for modes \'concatenate\' and \'concatenate-force\']')

  parser.add_argument('-c', '--clean',
                      action='store_true',
                      help='remove partial files after concatenation [for modes \'concatenate\' and \'concatenate-force\']')

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
    print fnames
    if len(fnames)>1:
      fnames = sorted(fnames, key=lambda x: (int(re.sub('\D','',x)),x))

    if args.type_of_run == 'locally':
      for fname in fnames:
        goja_command = "goja --root " + args.path_gate_output + fname \
                     + " --eth0 " + str(args.eth0) \
                     + " --tw " + str(args.tw) \
                     + " --N0 " + str(args.N0) \
                     + " --save-real-time-to " + args.path_goja_output + fname[0:-5] + "_realtime" \
                     + " --save-statistics-to " + args.path_goja_output + fname[0:-5] + "_statistics" \
                     + " > " + args.path_goja_output + fname[0:-5] + "_coincidences &"
        os.system(goja_command)
        print goja_command

    elif args.type_of_run == 'on-cluster':
      basename = fnames[0][0:-5]
      basepath_goja = (args.path_goja_output + basename).rstrip(string.digits)
      basepath_gate = (args.path_gate_output + basename).rstrip(string.digits)
      # generate array.pbs:
      with open('array.pbs', 'a') as array_pbs:
        array_pbs.write('#!/bin/sh\n')
        array_pbs.write('#PBS -q i3d\n')
        array_pbs.write('#PBS -l nodes=1:ppn=1\n')
        array_pbs.write('#PBS -N GOJA\n')
        array_pbs.write('#PBS -V\n')
        array_pbs.write('cd ${PBS_O_WORKDIR}\n')
        goja_command = "goja --root " + basepath_gate + "${PBS_ARRAYID}" + ".root" \
                     + " --eth0 " + str(args.eth0) \
                     + " --tw " + str(args.tw) \
                     + " --N0 " + str(args.N0) \
                     + " --save-real-time-to " + basepath_goja + "${PBS_ARRAYID}" + "_realtime" \
                     + " --save-statistics-to " + basepath_goja + "${PBS_ARRAYID}" + "_statistics" \
                     + " > " + basepath_goja + "${PBS_ARRAYID}" + "_coincidences"
        array_pbs.write(goja_command + '\n')
        array_pbs.write('exit 0;\n')
      # push into queue:
      qsub_command = 'qsub -t 1-' + str(len(fnames)) + ' array.pbs'
      os.system(qsub_command)
      # remove array.pbs:
      os.unlink('array.pbs')

    else:
      print "Improper type of run. " + help_message

  elif args.mode=="analyze-missing":

    print "Analysis of missing files:"

    fnames = os.listdir(args.path_gate_output)
    fnames = [fname for fname in fnames if ".root" in fname]
    fnames = sorted(fnames, key=lambda x: (int(re.sub('\D','',x)),x))

    for fname in fnames:

      path_coincidences = args.path_goja_output + fname[:-5] + "_coincidences"
      path_realtime = args.path_goja_output + fname[:-5] + "_realtime"
      path_statistics = args.path_goja_output + fname[:-5] + "_statistics"

      if not os.path.isfile(path_coincidences) or not os.path.isfile(path_realtime) or not os.path.isfile(path_statistics):

        if args.type_of_run == 'locally':
          basepath_goja = args.path_goja_output + fname[:-5]
          goja_command = "goja --root " + args.path_gate_output + fname \
                        + " --eth0 " + str(args.eth0) \
                        + " --tw " + str(args.tw) \
                        + " --N0 " + str(args.N0) \
                        + " --save-real-time-to " + basepath_goja + "_realtime" \
                        + " --save-statistics-to " + basepath_goja + "_statistics" \
                        + " > " + basepath_goja + "_coincidences &"
          os.system(goja_command)
          print goja_command

        elif args.type_of_run == 'on-cluster':
          basename = fname[0:-5]
          basepath_goja = args.path_goja_output + basename
          basepath_gate = args.path_gate_output + basename
          # generate array.pbs:
          with open('array.pbs', 'a') as array_pbs:
            array_pbs.write('#!/bin/sh\n')
            array_pbs.write('#PBS -q i3d\n')
            array_pbs.write('#PBS -l nodes=1:ppn=1\n')
            array_pbs.write('#PBS -N GOJA\n')
            array_pbs.write('#PBS -V\n')
            array_pbs.write('cd ${PBS_O_WORKDIR}\n')
            goja_command = "goja --root " + basepath_gate + ".root" \
                        + " --eth0 " + str(args.eth0) \
                        + " --tw " + str(args.tw) \
                        + " --N0 " + str(args.N0) \
                        + " --save-real-time-to " + basepath_goja + "_realtime" \
                        + " --save-statistics-to " + basepath_goja + "_statistics" \
                        + " > " + basepath_goja + "_coincidences"
            array_pbs.write(goja_command + '\n')
            array_pbs.write('exit 0;\n')
          # push into queue:
          qsub_command = 'qsub array.pbs'
          os.system(qsub_command)
          # remove array.pbs:
          os.unlink('array.pbs')

        else:
          print "Improper type of run. " + help_message

  elif args.mode=="verify":

    print "Verification:"

    if not verify_goja_output(args.path_gate_output, args.path_goja_output):
      print "\tGOJA output is ok."

  elif args.mode=="concatenate":

    print "Concatenation:"

    nr_of_missing_files = verify_goja_output(args.path_gate_output, args.path_goja_output)

    if nr_of_missing_files:
      print "Concatenation cannot be performed because of missing files. " + \
            "Obtain missing files using the analyze-missing mode or " + \
            "force the concatenation (concatenate-force) mode."
    else:
      fnames = os.listdir(args.path_goja_output)
      fnames = [fname for fname in fnames if "_coincidences" in fname]
      if len(fnames)>1:
        fnames = sorted(fnames, key=lambda x: (int(re.sub('\D','',x)),x))
      concatenate_files(fnames)

  elif args.mode=="concatenate-force":

    print "Concatenation (force):"

    fnames_tmp = os.listdir(args.path_goja_output)
    fnames_tmp = [fname for fname in fnames_tmp if "_coincidences" in fname]
    fnames = []
    for fname in fnames_tmp:
      path_coincidences = args.path_goja_output + fname
      path_realtime = args.path_goja_output + fname.replace("_coincidences", "") + "_realtime"
      path_statistics = args.path_goja_output + fname.replace("_coincidences", "") + "_statistics"
      if os.path.isfile(path_coincidences) and os.path.isfile(path_realtime) and os.path.isfile(path_statistics):
        fnames.append(fname)
    if len(fnames)>1:
      fnames = sorted(fnames, key=lambda x: (int(re.sub('\D','',x)),x))
    concatenate_files(fnames)

  else:

    print 'Improper mode. ' + help_message
