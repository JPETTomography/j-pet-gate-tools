#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import errno
from numpy import *
import os
import re
import string
import subprocess
import sys
from getpass import getuser

QUEUE_RUN = 'o14d'
QUEUE_ANALYZE = 'o12h'

# Paths:
ARRAY_PBS_GOJA = "array_goja.pbs"
ARRAY_PBS_MISSING = "array.pbs.missing"

def run_simulation(type_of_run):

  command_run = ''
  if type_of_run == "locally":
    command_run += 'Gate main.mac'
    print('\t' + command_run)
    p = subprocess.Popen(command_run, stdout=subprocess.PIPE, shell=True)
    p.wait()
  else:
    command_run += './Gate_parallel.sh'
    print('\t' + command_run)
    os.system(command_run)

def run_missing_simulations_on_cluster():

  with open("./missing_gate_results.txt", 'r') as missing_gate_results_txt:
    missing_gate_results = missing_gate_results_txt.readlines()
    if len(missing_gate_results)>0:
      for m in missing_gate_results:
        m = m.replace('\n', '')
        with open("array.pbs", "r") as array_pbs:
          with open(ARRAY_PBS_MISSING, "w") as array_pbs_missing:
            for line in array_pbs:
              array_pbs_missing.write(line.replace('${PBS_ARRAYID}', m))
        os.system('qsub ' + ARRAY_PBS_MISSING)

def get_goja_command(gate_result, goja_result, eth, eth0, tw, N0):

  goja_command = "goja --root " + gate_result \
               + " --eth " + str(eth) \
               + " --eth0 " + str(eth0) \
               + " --tw " + str(tw) \
               + " --N0 " + str(N0) \
               + " --save-real-time-to " + goja_result + "_realtime" \
               + " --save-statistics-to " + goja_result + "_statistics" \
               + " > " + goja_result + "_coincidences"
  return goja_command

def analyze_simulations_on_cluster(path_gate_output, single_file_prefix, path_goja_output, splits, eth, eth0, tw, N0):

  for s in splits:
    ss = str(int(s))
    gate_result = path_gate_output + single_file_prefix + ss + ".root"
    goja_result = path_goja_output + single_file_prefix + ss
    goja_command = get_goja_command(gate_result, goja_result, eth, eth0, tw, N0)
    # generate ARRAY_PBS_GOJA:
    with open(ARRAY_PBS_GOJA, 'w') as array_pbs:
      array_pbs.write('#!/bin/sh\n')
      array_pbs.write('#PBS -q ' + QUEUE_ANALYZE + '\n')
      array_pbs.write('#PBS -l nodes=1:ppn=1\n')
      array_pbs.write('#PBS -N GOJA' + ss + '\n')
      array_pbs.write('#PBS -V\n')
      array_pbs.write('cd ${PBS_O_WORKDIR}\n')
      array_pbs.write(goja_command + '\n')
      array_pbs.write('exit 0;\n')
    # push into queue:
    qsub_command = 'qsub ' + ARRAY_PBS_GOJA
    os.system(qsub_command)

def get_nr_of_splits(simulation_path):

  nr_of_splits = 0
  with open('./Gate_parallel.sh', 'r') as gate_parallel_sh:
    for line in gate_parallel_sh:
      if line.split('=')[0] == 'NR_OF_SPLITS':
        nr_of_splits = int((line.split('=')[1]).replace('\"', '').replace('\'', '').replace('\n', ''))
  return nr_of_splits

def verify_gate_output(path_gate_output, type_of_run):

  nr_of_missing_files = 0

  if type_of_run == "locally":
    missing_files = []
    output_root = path_gate_output + 'output.root'
    if not os.path.isfile(output_root):
      print("\tFile ", output_root, " is missing.")
      missing_files.append(output_root)
    nr_of_missing_files = len(missing_files)

  elif type_of_run == "on-cluster":
    VERIFY_GATE_RESULTS_PBS = "verify_gate_results.pbs"
    with open(VERIFY_GATE_RESULTS_PBS, 'w') as file_pbs:
      file_pbs.write('#!/bin/sh\n')
      file_pbs.write('#PBS -q o12h\n')
      file_pbs.write('#PBS -l nodes=1:ppn=1\n')
      file_pbs.write('#PBS -N verify_gate_results.py\n')
      file_pbs.write('#PBS -V\n')
      file_pbs.write('cd ${PBS_O_WORKDIR}\n')
      file_pbs.write('verify_gate_results.py\n')
      file_pbs.write('exit 0;\n')
    # push into queue:
    qsub_command = 'qsub ' + VERIFY_GATE_RESULTS_PBS
    os.system(qsub_command)

  return nr_of_missing_files

def verify_goja_output(path_gate_output, path_goja_output, type_of_run):

  nr_of_missing_files = 0
  missing_files = []

  if type_of_run == "locally":
    for fname in os.listdir(path_gate_output):
      if ".root" in fname:
        path_coincidences = path_goja_output + fname[:-5] + "_coincidences"
        if not os.path.isfile(path_coincidences):
          print("\tFile ", path_coincidences, " is missing.")
          missing_files.append(path_coincidences)
        path_realtime = path_goja_output + fname[:-5] + "_realtime"
        if not os.path.isfile(path_realtime):
          print("\tFile ", path_realtime, " is missing.")
          missing_files.append(path_realtime)
        path_statistics = path_goja_output + fname[:-5] + "_statistics"
        if not os.path.isfile(path_statistics):
          print("\tFile ", path_statistics, " is missing.")
          missing_files.append(path_statistics)
    nr_of_missing_files = len(missing_files)

  elif type_of_run == "on-cluster":
    VERIFY_GOJA_RESULTS_PBS = "verify_goja_results.pbs"
    with open(VERIFY_GOJA_RESULTS_PBS, 'w') as file_pbs:
      file_pbs.write('#!/bin/sh\n')
      file_pbs.write('#PBS -q o12h\n')
      file_pbs.write('#PBS -l nodes=1:ppn=1\n')
      file_pbs.write('#PBS -N verify_goja_results.py\n')
      file_pbs.write('#PBS -V\n')
      file_pbs.write('cd ${PBS_O_WORKDIR}\n')
      file_pbs.write('verify_goja_results.py\n')
      file_pbs.write('exit 0;\n')
    # push into queue:
    qsub_command = 'qsub ' + VERIFY_GOJA_RESULTS_PBS
    os.system(qsub_command)

  return nr_of_missing_files, missing_files

def concatenate_files(fnames):

      path_coincidences = path_goja_output + args.simulation_name + "_COINCIDENCES"
      try:
        os.unlink(path_coincidences)
      except OSError as e:
        if e.errno == errno.ENOENT:
          pass
        else:
          raise e
      path_realtime = path_goja_output + args.simulation_name + "_REALTIME"
      try:
        os.unlink(path_realtime)
      except OSError as e:
        if e.errno == errno.ENOENT:
          pass
        else:
          raise e
      path_statistics = path_goja_output + args.simulation_name + "_STATISTICS"
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
          basepath_goja = path_goja_output + basename

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
          basepath_goja = path_goja_output + basename
          os.unlink(basepath_goja + "_coincidences")
          os.unlink(basepath_goja + "_realtime")
          os.unlink(basepath_goja + "_statistics")
      print("Goja output succesfully concatenated.")

if __name__ == "__main__":

  help_message = "Type goja_manager.py --help."

  parser = argparse.ArgumentParser(description='Analyze, verify and concatenate the GOJA results.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-gt', '--gate-output',
                      dest='path_gate_output',
                      type=str,
                      default="",
                      help='path to dir with the splitted GATE results')

  parser.add_argument('-sfp', '--single-file-prefix',
                      dest='single_file_prefix',
                      type=str,
                      default="output",
                      help='single file prefix')

  parser.add_argument('-gj', '--goja-output',
                      dest='path_goja_output',
                      type=str,
                      default="",
                      help='path to dir with the GOJA results')

  parser.add_argument('-sp', '--simulation-path',
                      dest='simulation_path',
                      type=str,
                      default=".",
                      help='path to dir with the simulation (for the purpose of the Simulations Manager)')

  parser.add_argument('-m', '--mode',
                      dest='mode',
                      type=str,
                      default="analyze",
                      help='mode of the script: run, run-missing, analyze, analyze-missing, verify, verify-gate, verify-goja, concatenate, clear-gate, clear-goja, clear-cluster-artifacts')

  parser.add_argument('--eth',
                      dest='eth',
                      type=float,
                      default=0.2,
                      help='fixed energy threshold in MeV [for mode \'analyze\']')

  parser.add_argument('--eth0',
                      dest='eth0',
                      type=float,
                      default=0.0,
                      help='noise energy threshold in MeV [for mode \'analyze\']')

  parser.add_argument('-tw', '--time-window',
                      dest='tw',
                      type=float,
                      default=3,
                      help='time window in ns [for mode \'analyze\']')

  parser.add_argument('--N0',
                      dest='N0',
                      type=int,
                      default=1000,
                      help='maximum number of events above the noise energy threshold in the coincidence window [for mode \'analyze\']')

  parser.add_argument('-nos', '--nr-of-splits',
                      dest='nr_of_splits',
                      type=int,
                      default=0,
                      help='number of splits')

  parser.add_argument('-r', '--run',
                      dest='type_of_run',
                      type=str,
                      default='on-cluster',
                      help='run \'locally\' or \'on-cluster\' [for modes: run, run-missing, analyze, analyze-missing]')

  parser.add_argument('-sn', '--simulation-name',
                      dest='simulation_name',
                      type=str,
                      default="simulation",
                      help='name of the simulation [for modes \'concatenate\' and \'concatenate-force\']')

  parser.add_argument('-c', '--clean',
                      action='store_true',
                      help='remove partial files after concatenation [for modes \'concatenate\' and \'concatenate-force\']')

  parser.add_argument('--lustre',
                      action='store_true',
                      help='use lustre file system (if not nfs is used)')

  args = parser.parse_args()

  current_path = os.getcwd()

  path_gate_output = ""
  if args.path_gate_output == "":
    if args.type_of_run == "locally":
      path_gate_output = current_path + "/output/"
    elif args.type_of_run == "on-cluster":
      if args.lustre:
        path_gate_output = '/mnt/lustre/home/' + getuser() + '/' + args.simulation_path + "/output/"
      else:
        path_gate_output = current_path + '/output/'
  else:
    path_gate_output = args.path_gate_output

  path_goja_output = ""
  if args.path_goja_output == "":
    path_goja_output = current_path + "/goja/"
  else:
    path_goja_output = args.path_goja_output

  if args.type_of_run=="locally" and not ('run' in args.mode) and not os.path.isdir(path_gate_output):
    print("Directory " + path_gate_output + " does not exist. " + help_message)
    sys.exit()

  if not os.path.isdir(path_goja_output):
    try:
      os.system('mkdir -p ' + path_goja_output)
    except:
      pass
    if not os.path.isdir(path_goja_output):
      print("Directory " + path_goja_output + " does not exist. " + help_message)
      sys.exit()

  if args.mode == "run":

    print("Run:")

    run_simulation(args.type_of_run)

  elif args.mode == "run-missing":

    print("Run missing:")

    if args.type_of_run == "locally":
      if verify_gate_output(path_gate_output, args.type_of_run)>0:
        run_simulation(args.type_of_run)
    else:
      if os.path.isfile("./missing_gate_results.txt"):
        run_missing_simulations_on_cluster()

  elif args.mode == "analyze":

    print("Analyze:")

    if args.type_of_run == 'locally':
      gate_result = path_gate_output + "output.root"
      goja_result = path_goja_output + args.single_file_prefix
      goja_command = get_goja_command(gate_result, goja_result, args.eth, args.eth0, args.tw, args.N0)
      print(goja_command)
      p = subprocess.Popen(goja_command, shell=True)
      p.wait()

    elif args.type_of_run == 'on-cluster':
      if args.nr_of_splits == 0:
        nr_of_splits = get_nr_of_splits(args.simulation_path)
      else:
          nr_of_splits = args.nr_of_splits
      splits = linspace(1, nr_of_splits, nr_of_splits)
      analyze_simulations_on_cluster(path_gate_output, args.single_file_prefix, path_goja_output, splits, args.eth, args.eth0, args.tw, args.N0)

  elif args.mode == "analyze-missing":

    print("Analyze missing:")

    if args.type_of_run == 'locally':
      gate_result = path_gate_output + "output.root"
      goja_result = path_goja_output + args.single_file_prefix
      if not os.path.isfile(goja_result + "coincidences") or \
         not os.path.isfile(goja_result + "realtime") or \
         not os.path.isfile(goja_result + "statistics"):
        goja_command = get_goja_command(gate_result, goja_result, args.eth, args.eth0, args.tw, args.N0)
        print(goja_command)
        p = subprocess.Popen(goja_command, shell=True)
        p.wait()

    elif args.type_of_run == 'on-cluster':
      if os.path.isfile("./missing_goja_results.txt"):
        missing_goja_results = loadtxt("./missing_goja_results.txt")
        if len(missing_goja_results)>0:
          analyze_simulations_on_cluster(path_gate_output, path_goja_output, missing_goja_results, args.eth, args.eth0, args.tw, args.N0)

  elif args.mode == "verify":

    print("Verify:")

    verify_gate_output(path_gate_output, args.type_of_run)
    verify_goja_output(path_gate_output, path_goja_output, args.type_of_run)

  elif args.mode == "verify-gate":

    print("Verify (GATE):")

    verify_gate_output(path_gate_output, args.type_of_run)

  elif args.mode == "verify-goja":

    print("Verify (GOJA):")

    verify_goja_output(path_gate_output, path_goja_output, args.type_of_run)

  elif args.mode == "concatenate":

    print("Concatenate:")

    fnames_tmp = os.listdir(path_goja_output)
    fnames_tmp = [fname for fname in fnames_tmp if "_coincidences" in fname]
    fnames = []
    for fname in fnames_tmp:
      path_coincidences = path_goja_output + fname
      path_realtime = path_goja_output + fname.replace("_coincidences", "") + "_realtime"
      path_statistics = path_goja_output + fname.replace("_coincidences", "") + "_statistics"
      if os.path.isfile(path_coincidences) and os.path.isfile(path_realtime) and os.path.isfile(path_statistics):
        fnames.append(fname)
    if len(fnames)>1:
      fnames = sorted(fnames, key=lambda x: (int(re.sub('\D','',x)),x))
    concatenate_files(fnames)

  elif args.mode == "clear-gate":

    print("Clear (GATE):")
    command = 'rm -f ' + path_gate_output + '*'
    print('\t' + command)
    os.system(command)

  elif args.mode == "clear-goja":

    print("Clear (GOJA):")
    command = 'rm -f ' + path_goja_output + '*'
    print('\t' + command)
    os.system(command)

  elif args.mode == "clear-cluster-artifacts":

    print("Clear (cluster artifacts):")
    command = 'rm -f *.o* *.e* ' + ARRAY_PBS_GOJA + ' ' + ARRAY_PBS_MISSING
    print('\t' + command)
    os.system(command)

  else:

    print('Improper mode. ' + help_message)
