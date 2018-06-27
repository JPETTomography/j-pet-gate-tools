#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from numpy import *
from matplotlib import rcParams
from matplotlib.colors import LogNorm
from math import *
import argparse

from nema_common import *

OUTPUT_FORMAT = ".png"

## Plot Da vs. Dt using data from GATE simulations.
#
#  Da - differences of central angles between hits in the coincidence
#  Dt - differences of times of interactions of hits in the coincidence
#
#  Parameters
#  ----------
#  coincidences : ndarray
#      Data read from the listmode goja output file using the numpy.loadtxt.
#  result_figure_path : str
#      Path to the output image file.
#  show_cut : bool
#      Switch to enable or disable the ellipsoidal cut line on the image.
#  ylim : tuple
#      Limits of the angles differences axis.
#  toc : int
#      Type of the coincidence. When set to 0, all coincidences are plotted.
#      When set to 1 (true), 2 (psca), 3 (psca) or 4 (acci), only coincidences
#      with a chosen type of the coincidence are plotted.
def plot_Da_vs_Dt(coincidences, result_figure_path, show_cut, ylim=[0,180], toc=0):

  posX1 = coincidences[:,0]
  posY1 = coincidences[:,1]
  times1 = coincidences[:,3]
  posX2 = coincidences[:,4]
  posY2 = coincidences[:,5]
  times2 = coincidences[:,7]
  type_of_coincidence = coincidences[:,12]

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
  plt.imshow(H.T, interpolation='none', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', norm=LogNorm(vmin=1, vmax=VMAX))
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

## Plot Da vs. Dt using data from experiment.
#
#  Da - differences of central angles between hits in the coincidence
#  Dt - differences of times of interactions of hits in the coincidence
#
#  Parameters
#  ----------
#  coincidences : ndarray
#      Data read from the listmode goja output file using the numpy.loadtxt.
#  result_figure_path : str
#      Path to the output image file.
def plot_Da_vs_Dt_exp(coincidences, result_figure_path):

  tim_diffs = coincidences[:,6]
  ang_diffs = coincidences[:,7]

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
  plt.imshow(H.T, interpolation='none', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', norm=LogNorm(vmin=1, vmax=VMAX))
  plt.colorbar()

  plt.xlabel("Time difference [ps]")
  plt.ylabel("Angle difference [deg.]")
  plt.savefig(result_figure_path + OUTPUT_FORMAT, bbox_inches='tight')
  plt.clf()
  plt.close()

## Plot source position x vs. y using data from GATE simulations.
#
#  source position - position of the annihilation
#
#  Parameters
#  ----------
#  coincidences : ndarray
#      Data read from the listmode goja output file using the numpy.loadtxt.
#  result_figure_path : str
#      Path to the output image file.
def plot_sourcePosX_vs_sourcePosY(coincidences, result_figure_path):

  sourcePosX = coincidences[:,13]
  sourcePosY = coincidences[:,14]

  # Plot 2D histogram

  rcParams['font.size'] = 24
  rcParams['legend.fontsize'] = 18

  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111)
  plt.subplots_adjust(left=0.20, right=0.9, top=0.9, bottom=0.1)

  H, xedges, yedges = histogram2d(sourcePosX, sourcePosY, bins=100)
  VMAX = H.max()
  plt.imshow(H.T, interpolation='none', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect='auto', norm=LogNorm(vmin=1, vmax=VMAX))
  plt.colorbar()

  plt.xlabel("sourcePosX [cm]")
  plt.ylabel("sourcePosY [cm]")
  plt.savefig(result_figure_path + OUTPUT_FORMAT, bbox_inches='tight')
  plt.clf()
  plt.close()

## Calculate (and print) ratios of accidental coincidences using data from GATE simulations.
#
#  2 ratios are calculated:
#   - ratio_acci - ratio between accidental and all coincidences
#   - ratio_acci_to_true - ratio between accidental and true coincidences
#
#  Parameters
#  ----------
#  coincidences : ndarray
#      Data read from the listmode goja output file using the numpy.loadtxt.
#  filename : str
#      Name of the goja output file used in the analysis.
def calculate_ratios(coincidences, filename):

  type_of_coincidence = coincidences[:,12]

  N_true = 0
  N_acci = 0
  N_all = len(type_of_coincidence)

  for i in xrange(N_all):
    if type_of_coincidence[i]==1: N_true += 1
    elif type_of_coincidence[i]==4: N_acci += 1

  ratio_acci = float(N_acci)/float(N_all)*100.
  ratio_acci_to_true = float(N_acci)/float(N_true)*100.

  print("filename={0}, ratio_acci={1:.2f}%, ratio_acci_to_true={2:.2f}%".format(filename, ratio_acci, ratio_acci_to_true))

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Analyze coincidences file and plot results of the analysis',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-cf', '--coincidences-file',
                      dest='path_coincidences_file',
                      type=str,
                      help='path to the coincidences file obtained using the GOJA tool')

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

  # Load data from file
  coincidences = loadtxt(args.path_coincidences_file)

  if args.mode == "plot":
    plot_Da_vs_Dt(coincidences, args.path_output_da_dt, args.show_cut, ylim=[90,180])
    #plot_Da_vs_Dt(coincidences, args.path_output_da_dt, args.show_cut, toc=1)
    #plot_Da_vs_Dt(coincidences, args.path_output_da_dt, args.show_cut, toc=2)
    #plot_Da_vs_Dt(coincidences, args.path_output_da_dt, args.show_cut, toc=3)
    #plot_Da_vs_Dt(coincidences, args.path_output_da_dt, args.show_cut, toc=4)
    plot_sourcePosX_vs_sourcePosY(coincidences, args.path_output_sposx_sposy)

  elif args.mode == "stats":
    calculate_ratios(coincidences, args.path_coincidences_file.split("/")[-1])
