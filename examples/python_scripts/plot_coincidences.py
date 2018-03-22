#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from numpy import *
from matplotlib import rcParams, rcdefaults, rc
from matplotlib.colors import LogNorm
from math import *
import os
import argparse

import matplotlib.image as mpimg
from collections import OrderedDict

# Plot using data from GATE simulations
def plot_Da_vs_Dt(goja_output_file, result_figure_path, show_cut, ylim=[0,180]):

  # Load data from file

  tmp = loadtxt(goja_output_file)

  posX1 = tmp[:,0]
  posY1 = tmp[:,1]
  posZ1 = tmp[:,2]
  times1 = tmp[:,3]
  posX2 = tmp[:,4]
  posY2 = tmp[:,5]
  posZ2 = tmp[:,6]
  times2 = tmp[:,7]
  IDs1 = tmp[:,8]
  IDs2 = tmp[:,9]
  toc = tmp[:,12]

  # Approx time of a whole simulation taken into account in seconds as a max value of registered hits

  time = max(max(times1),max(times2))/1e12
  I = len(tmp)

  # Calculate vectors of times and angles differences

  tim_diffs = []
  ang_diffs = []

  for i in range(len(tmp)):

    tdiff = abs(times1[i]-times2[i])/1e3 # in ns

    # v = [vx,vy] = posX1, posY1
    # u = [ux,uy] = posX2, posY2

    vx = posX1[i]
    vy = posY1[i]
    ux = posX2[i]
    uy = posY2[i]

    vu = vx*ux + vy*uy
    modv = sqrt(vx*vx + vy*vy)
    modu = sqrt(ux*ux + uy*uy)

    try:
        a = arccos(vu/(modv*modu))/pi*180.
        adiff = a
    except:
        pass

    tim_diffs.append(tdiff)
    ang_diffs.append(adiff)

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
  for i in range(len(xxx)):
    yyy.append(180-80*sqrt(1.-xxx[i]*xxx[i]/(param*param)))

  if show_cut:
    plt.plot(xxx, yyy, color='r', linewidth=2)

  counter_above = 0
  counter_below = 0
  for i in range(len(tim_diffs)):
    t = tim_diffs[i]
    a = ang_diffs[i]
    newt = 0
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
  plt.savefig(result_figure_path, bbox_inches='tight')
  plt.clf()
  plt.close()

# Plot using data from experiment
def plot_Da_vs_Dt_exp(datafile, result_figure_path):

  # Load data from file

  tmp = loadtxt(datafile)
  print size(tmp)

  posX1 = tmp[:,0]
  posY1 = tmp[:,1]
  posZ1 = tmp[:,2]
  posX2 = tmp[:,3]
  posY2 = tmp[:,4]
  posZ2 = tmp[:,5]

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
  plt.savefig(result_figure_path, bbox_inches='tight')
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
  plt.savefig(result_figure_path, bbox_inches='tight')
  plt.clf()
  plt.close()

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Plot DA (angles differences) vs. DT (times differences) using the coincidecnes file.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-cf', '--coincidences-file',
                      dest='path_coincidences_file',
                      type=str,
                      help='path to the coincidences file obtained using the GOJA tool')

  parser.add_argument('-oat', '--output-da-dt',
                      dest='path_output_da_dt',
                      type=str,
                      default="./DA_vs_DT.png",
                      help='path to the output figure')

  parser.add_argument('-osxsy', '--output-sposx-sposy',
                      dest='path_output_sposx_sposy',
                      type=str,
                      default="./sourcePosX_vs_sourcePosY.png",
                      help='path to the output figure sourcePosX vs. sourcePosY')

  parser.add_argument('-sc', '--show-cut',
                      dest='show_cut',
                      action='store_true',
                      help='set if you want to show cut')

  args = parser.parse_args()

  plot_Da_vs_Dt(args.path_coincidences_file, args.path_output_da_dt, args.show_cut)
  plot_sourcePosX_vs_sourcePosY(args.path_coincidences_file, args.path_output_sposx_sposy)
