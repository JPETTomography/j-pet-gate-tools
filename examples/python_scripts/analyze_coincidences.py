#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from numpy import *
import numpy as np
from matplotlib import rcParams
from matplotlib.colors import LogNorm
from math import *
import argparse

OUTPUT_FORMAT = ".png"

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

  tim_diffs = absolute(times1 - times2)/1e3 # in ns

  vx = posX1
  vy = posY1
  ux = posX2
  uy = posY2

  vu = vx*ux + vy*uy
  modv = np.sqrt(vx*vx + vy*vy)
  modu = np.sqrt(ux*ux + uy*uy)

  ang_diffs = np.arccos(vu/(modv*modu))/pi*180.

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
  param = 2.2
  xxx = linspace(0,param,100)
  yyy = []
  for i in xrange(len(xxx)):
    yyy.append(180-80*sqrt(1.-xxx[i]*xxx[i]/(param*param)))

  if show_cut:
    plt.plot(xxx, yyy, color='r', linewidth=2)

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

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Plot DA (angles differences) vs. DT (times differences) using the coincidecnes file.',
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

  args = parser.parse_args()

  if args.mode == "plot":
    plot_Da_vs_Dt(args.path_coincidences_file, args.path_output_da_dt, args.show_cut)
    #plot_Da_vs_Dt(args.path_coincidences_file, args.path_output_da_dt, args.show_cut, toc=1)
    #plot_Da_vs_Dt(args.path_coincidences_file, args.path_output_da_dt, args.show_cut, toc=2)
    #plot_Da_vs_Dt(args.path_coincidences_file, args.path_output_da_dt, args.show_cut, toc=3)
    #plot_Da_vs_Dt(args.path_coincidences_file, args.path_output_da_dt, args.show_cut, toc=4)
    plot_sourcePosX_vs_sourcePosY(args.path_coincidences_file, args.path_output_sposx_sposy)

  elif args.mode == "stats":
    calculate_ratios(args.path_coincidences_file)
