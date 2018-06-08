#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from numpy import *
import numpy as np
from matplotlib import rcParams
from matplotlib.colors import LogNorm
from math import *
import argparse
from scipy import interpolate

from nema_common import *

OUTPUT_FORMAT = ".png"
ELLIPSE_PARAM = 2.2

def calculate_differences(tim1, tim2, posx1, posy1, posx2, posy2):

  tim_diffs = absolute(tim1 - tim2)/1e3 # in ns

  vx = posx1
  vy = posy1
  ux = posx2
  uy = posy2

  vu = vx*ux + vy*uy
  modv = np.sqrt(vx*vx + vy*vy)
  modu = np.sqrt(ux*ux + uy*uy)

  ang_diffs = np.arccos(vu/(modv*modu))/pi*180.

  return [tim_diffs, ang_diffs]

def calculate_counters(tim_diffs, ang_diffs, param):

  counter_above = 0
  counter_below = 0
  for i in xrange(len(tim_diffs)):
    t = tim_diffs[i]
    a = ang_diffs[i]
    try:
      newa = 180-80*sqrt(1.-t*t/(param*param))
      if a>newa: counter_above += 1
      else: counter_below += 1
    except:
      pass

  return [counter_above, counter_below]

# Plot using data from GATE simulations
def plot_Da_vs_Dt(goja_output_file, result_figure_path, show_cut, ylim=[0,180], toc=0):

  # Load data from file

  tmp = loadtxt(goja_output_file)

  posX1 = tmp[:,0]
  posY1 = tmp[:,1]
  times1 = tmp[:,3]
  posX2 = tmp[:,4]
  posY2 = tmp[:,5]
  times2 = tmp[:,7]
  type_of_coincidence = tmp[:,12]

  if toc != 0:
    tmp_posX1 = []
    tmp_posY1 = []
    tmp_times1 = []
    tmp_posX2 = []
    tmp_posY2 = []
    tmp_times2 = []
    for i in xrange(len(type_of_coincidence)):
      if type_of_coincidence[i]==toc:
        tmp_posX1.append(posX1[i])
        tmp_posY1.append(posY1[i])
        tmp_times1.append(times1[i])
        tmp_posX2.append(posX2[i])
        tmp_posY2.append(posY2[i])
        tmp_times2.append(times2[i])
    posX1 = array(tmp_posX1)
    posY1 = array(tmp_posY1)
    times1 = array(tmp_times1)
    posX2 = array(tmp_posX2)
    posY2 = array(tmp_posY2)
    times2 = array(tmp_times2)

  label = ""
  if toc==1: label = "_true"
  elif toc==2: label = "_psca"
  elif toc==3: label = "_dsca"
  elif toc==4: label = "_acci"

  print len(posX1)

  # Calculate vectors of times and angles differences

  [tim_diffs, ang_diffs] = calculate_differences(times1, times2, posX1, posY1, posX2, posY2)

  # Plot 2D histogram

  rcParams['font.size'] = 24
  rcParams['legend.fontsize'] = 18

  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111)
  plt.subplots_adjust(left=0.20, right=0.9, top=0.9, bottom=0.1)

  t_bins = 100
  a_bins = 180

  H, xedges, yedges = histogram2d(tim_diffs, ang_diffs, bins=(t_bins,a_bins), range=[[0, 3],[0, 180]])
  VMAX = H.max()
  plt.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', norm=LogNorm(vmin=1, vmax=VMAX))
  plt.colorbar()

  # Ellipse curve
  xxx = linspace(0,ELLIPSE_PARAM,100)
  yyy = []
  for i in xrange(len(xxx)):
    yyy.append(180-80*sqrt(1.-xxx[i]*xxx[i]/(ELLIPSE_PARAM*ELLIPSE_PARAM)))

  if show_cut:
    plt.plot(xxx, yyy, color='r', linewidth=2)

  [counter_above, counter_below] = calculate_counters(tim_diffs, ang_diffs, ELLIPSE_PARAM)
  print "Number of events above the cut: ", counter_above
  print "Number of events below the cut: ", counter_below
  print "Percentage of cut events: ", float(counter_below)/(counter_above+counter_below)*100.

  plt.ylim(ylim)

  plt.xlabel("Time difference [ns]")
  plt.ylabel("Angle difference [deg.]")
  plt.savefig(result_figure_path + label + OUTPUT_FORMAT, bbox_inches='tight')
  plt.clf()
  plt.close()

# Plot using data from experiment
def plot_Da_vs_Dt_exp(datafile, result_figure_path):

  # Load data from file

  tmp = loadtxt(datafile)
  print size(tmp)

  tim_diffs = tmp[:,6]
  ang_diffs = tmp[:,7]

  # Plot 2D histogram

  rcParams['font.size'] = 24
  rcParams['legend.fontsize'] = 18

  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111)
  plt.subplots_adjust(left=0.20, right=0.9, top=0.9, bottom=0.1)

  t_bins = 25
  a_bins = 45

  H, xedges, yedges = histogram2d(tim_diffs, ang_diffs, bins=(t_bins,a_bins)) #, range=[[0, 3],[0, 180]])
  VMAX = H.max()
  plt.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', norm=LogNorm(vmin=1, vmax=VMAX))
  plt.colorbar()

  plt.xlabel("Time difference [ps]")
  plt.ylabel("Angle difference [deg.]")
  plt.savefig(result_figure_path + OUTPUT_FORMAT, bbox_inches='tight')
  plt.clf()
  plt.close()

def plot_sourcePosX_vs_sourcePosY(goja_output_file, result_figure_path):

  # Load data from file

  tmp = loadtxt(goja_output_file)

  sourcePosX = tmp[:,13]
  sourcePosY = tmp[:,14]

  # Plot 2D histogram

  rcParams['font.size'] = 24
  rcParams['legend.fontsize'] = 18

  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111)
  plt.subplots_adjust(left=0.20, right=0.9, top=0.9, bottom=0.1)

  H, xedges, yedges = histogram2d(sourcePosX, sourcePosY, bins=100)
  VMAX = H.max()
  plt.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', norm=LogNorm(vmin=1, vmax=VMAX))
  plt.colorbar()

  plt.xlabel("sourcePosX [cm]")
  plt.ylabel("sourcePosY [cm]")
  plt.savefig(result_figure_path + OUTPUT_FORMAT, bbox_inches='tight')
  plt.clf()
  plt.close()

def calculate_ratios(goja_output_file):

  filename = goja_output_file.split("/")[-1]

  tmp = loadtxt(goja_output_file)

  type_of_coincidence = tmp[:,12]

  N_true = 0
  N_acci = 0
  N_all = len(type_of_coincidence)

  for i in xrange(N_all):
    if type_of_coincidence[i]==1: N_true += 1
    elif type_of_coincidence[i]==4: N_acci += 1

  ratio_acci = float(N_acci)/float(N_all)*100.
  ratio_acci_to_true = float(N_acci)/float(N_true)*100.

  print("filename={0}, ratio_acci={1:.2f}%, ratio_acci_to_true={2:.2f}%".format(filename, ratio_acci, ratio_acci_to_true))

def calculate_reduction_for_necr_simulations(necr_simulations):

  for g in geometries_NECR:
    sls_file = workdir_NECR + g + "/second_lvl_selection.txt"
    if os.path.exists(sls_file):
      os.system('rm ' + sls_file)
    for a in activities_NECR:
      coincidences_file = necr_simulations + "/" + g + "_" + a + "_NECR_COINCIDENCES_short"
      tmp = loadtxt(coincidences_file)
      posX1 = tmp[:,0]
      posY1 = tmp[:,1]
      times1 = tmp[:,3]
      posX2 = tmp[:,4]
      posY2 = tmp[:,5]
      times2 = tmp[:,7]
      [tim_diffs, ang_diffs] = calculate_differences(times1, times2, posX1, posY1, posX2, posY2)
      [counter_above, counter_below] = calculate_counters(tim_diffs, ang_diffs, ELLIPSE_PARAM)
      with open(sls_file, "a") as myfile:
        myfile.write("{0}\t{1}\t{2}\n".format(counter_above, counter_below, counter_above+counter_below))
      print g + "\t" + a + "\t" + str(counter_above) + "\t" + str(counter_below) + "\t" + str(counter_above+counter_below)

def plot_reduction_for_necr_simulations():

  rcParams['font.size'] = 24
  rcParams['legend.fontsize'] = 18

  activities = []
  for a in activities_NECR:
    activities.append(float(a)/22000.*1000) # in kBq/cc
  new_activities = linspace(activities[0],activities[-1],100)

  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111)

  plt.ylim(ymin=0,ymax=80)
  plt.xlim(xmin=0,xmax=90)

  for g in geometries_NECR:

    lab = ""
    c = ""
    l = ""
    if "1lay" in g:
      lab += "1 layer, "
      c = 'k'
    elif "2lay" in g:
      lab += "2 layers, "
      c = 'r'
    if "L020" in g:
      lab += "L = 20 cm"
      l = '--'
    elif "L050" in g:
      lab += "L = 50 cm"
      l = '-'
    elif "L100" in g:
      lab += "L = 100 cm"
      l = '-.'

    sls_file = workdir_NECR + g + "/second_lvl_selection.txt"
    if os.path.exists(sls_file):
      tmp = loadtxt(sls_file)
      counter_above = tmp[:,0]
      counter_below = tmp[:,1]
      reduction = counter_below/(counter_above+counter_below)*100.
      new_reduction = interpolate.splev(new_activities, interpolate.splrep(activities, reduction, s=5), der=0)
      plt.plot(new_activities, new_reduction, linestyle=l, color=c, label=lab)

  plt.legend(loc=4)
  plt.xlabel("Activity concentration [kBq/cc]")
  plt.ylabel("Reduction [%]")
  plt.savefig(workdir_NECR + "second_lvl_selection" + OUTPUT_FORMAT, bbox_inches='tight')

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Plot DA (angles differences) vs. DT (times differences) using the coincidecnes file.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-cf', '--coincidences-file',
                      dest='path_coincidences_file',
                      type=str,
                      help='path to the coincidences file obtained using the GOJA tool')

  parser.add_argument('-ns', '--necr-simulations',
                      dest='necr_simulations',
                      type=str,
                      help='path to the base directory of the NECR simulations')

  parser.add_argument('-m', '--mode',
                      dest='mode',
                      type=str,
                      default="plot",
                      help='mode of the script: plot or stats')

  parser.add_argument('-oat', '--output-da-dt',
                      dest='path_output_da_dt',
                      type=str,
                      default="./DA_vs_DT",
                      help='path to the output figure')

  parser.add_argument('-osxsy', '--output-sposx-sposy',
                      dest='path_output_sposx_sposy',
                      type=str,
                      default="./sourcePosX_vs_sourcePosY",
                      help='path to the output figure sourcePosX vs. sourcePosY')

  parser.add_argument('-sc', '--show-cut',
                      dest='show_cut',
                      action='store_true',
                      help='set if you want to show cut')

  parser.add_argument('-of', '--outputformat',
                      dest='outputformat',
                      type=str,
                      default="png",
                      help='output format of images')

  args = parser.parse_args()

  OUTPUT_FORMAT = "." + args.outputformat

  if args.mode == "plot":
    plot_Da_vs_Dt(args.path_coincidences_file, args.path_output_da_dt, args.show_cut, ylim=[90,180])
    #plot_Da_vs_Dt(args.path_coincidences_file, args.path_output_da_dt, args.show_cut, toc=1)
    #plot_Da_vs_Dt(args.path_coincidences_file, args.path_output_da_dt, args.show_cut, toc=2)
    #plot_Da_vs_Dt(args.path_coincidences_file, args.path_output_da_dt, args.show_cut, toc=3)
    #plot_Da_vs_Dt(args.path_coincidences_file, args.path_output_da_dt, args.show_cut, toc=4)
    plot_sourcePosX_vs_sourcePosY(args.path_coincidences_file, args.path_output_sposx_sposy)

  elif args.mode == "stats":
    calculate_ratios(args.path_coincidences_file)

  elif args.mode == "necr":
    calculate_reduction_for_necr_simulations(args.necr_simulations)
    plot_reduction_for_necr_simulations()
