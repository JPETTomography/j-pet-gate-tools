#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from math import *
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter
from numpy import *
import sys

from nema_common import *

OUTPUT_FORMAT = ".png"

def get_strips_centers(geometry):

  xs = []
  ys = []

  if geometry == "lab192": #TODO add other geometries

    r1 = 42.5
    r2 = 46.75
    r3 = 57.5

    n1 = 48
    n3 = 96

    a_step1 = 2*pi/n1

    angles1 = linspace(0, 2*pi-a_step1, n1)
    xs1 = r1*cos(angles1)
    ys1 = r1*sin(angles1)

    angles2 = linspace(0, 2*pi-a_step1, n1)+a_step1/2.
    xs2 = r2*cos(angles2)
    ys2 = r2*sin(angles2)

    a_step3 = 2*pi/n3

    angles3 = linspace(0, 2*pi-a_step3, n3)+a_step3/2.
    xs3 = r3*cos(angles3)
    ys3 = r3*sin(angles3)

    xs = concatenate((xs1, xs2, xs3), axis=None)
    ys = concatenate((ys1, ys2, ys3), axis=None)

  strips_centers = []
  for i in range(len(xs)):
    strips_centers.append((xs[i], ys[i]))

  return strips_centers

def get_closest_strip_center(strip_center, strips_centers):

  strips_centers = asarray(strips_centers)
  dist_2 = sum((strips_centers - strip_center)**2, axis=1)
  result = strips_centers[argmin(dist_2)]
  return result[0], result[1]

## Plot Da vs. Dt using data from GATE simulations.
#
#  Da - differences of central angles between hits in the coincidence
#  Dt - differences of times of interactions of hits in the coincidence
#
#  Parameters
#  ----------
#  coincidences : ndarray
#      Data read from the listmode goja output file using the numpy.loadtxt.
#  tw : float
#      Time window, for which the coincidence set was produced [ns].
#  result_figure_path : str
#      Path to the output image file.
#  show_cut : bool
#      Switch to enable or disable the ellipsoidal cut line on the image.
#  t_bins : int
#      Number of bins for times differences axis.
#  a_bins : int
#      Number of bins for angles differences axis for continous values or number of strips for discrete values.
#  ylim : tuple
#      Limits of the angles differences axis.
#  toc : int
#      Type of the coincidence. When set to 0, all coincidences are plotted.
#      When set to 1 (true), 2 (psca), 3 (psca) or 4 (acci), only coincidences
#      with a chosen type of the coincidence are plotted.
#  discrete : bool
#      If True, then align x and y of the hits to strips centers.
#  recalculate : bool
#      If True, then recalculate strips_centers.txt file needed for discretization.
def plot_Da_vs_Dt(coincidences, tw, result_figure_path, show_cut, t_bins, a_bins, ylim=[0,180], toc=0, discrete=False, recalculate=False):

  posX1 = coincidences[:,0]
  posY1 = coincidences[:,1]
  times1 = coincidences[:,3]
  posX2 = coincidences[:,4]
  posY2 = coincidences[:,5]
  times2 = coincidences[:,7]
  type_of_coincidence = coincidences[:,12]

  vol1 = coincidences[:,8]
  vol2 = coincidences[:,9]

  centers_x = zeros(192)
  centers_y = zeros(192)
  xxx = linspace(1,192,192)

  if discrete:

    if recalculate or not os.path.exists("./strips_centers.txt"):

      for i in range(len(posX1)):

        strips_centers = get_strips_centers('lab192')
        x, y = get_closest_strip_center((posX1[i], posY1[i]), strips_centers)
        ind = int(vol1[i])-1
        centers_x[ind] = x
        centers_y[ind] = y

        if i == 10000:
          savetxt('./strips_centers.txt', array([xxx, centers_x, centers_y]).T, fmt='%d\t%.18e\t%.18e')
          print 'File strips_centers.txt generated. Try again.'
          sys.exit(1)

    else:

      strips_centers = loadtxt("./strips_centers.txt")

      for i in range(len(posX1)):
        ind1 = int(vol1[i])-1
        posX1[i] = strips_centers[:,1][ind1]
        posY1[i] = strips_centers[:,2][ind1]
        ind2 = int(vol2[i])-1
        posX2[i] = strips_centers[:,1][ind2]
        posY2[i] = strips_centers[:,2][ind2]

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

  # Calculate vectors of times and angles differences:

  [tim_diffs, ang_diffs] = calculate_differences(times1, times2, posX1, posY1, posX2, posY2)

  # Plot 2D histogram:

  rcParams['font.size'] = 24
  rcParams['legend.fontsize'] = 18
  ASPECT = tw/(ylim[1]-ylim[0]) # force square pixels

  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111)
  plt.subplots_adjust(left=0.20, right=0.9, top=0.9, bottom=0.1)

  a_bin_size = 180./a_bins
  a_bins_vec = linspace(0-a_bin_size/2., 180.+a_bin_size/2., a_bins+2)
  t_bins_vec = linspace(0, tw, t_bins+1)

  H, edges_tim_diffs, edges_ang_diffs = histogram2d(tim_diffs, ang_diffs, bins=(t_bins_vec,a_bins_vec), range=[[0, tw],ylim])
  VMAX = H.max()
  print "VMAX=", VMAX
  VMAX = 100000 #TODO
  plt.imshow(H.T, interpolation='none', origin='low', extent=[edges_tim_diffs[0], edges_tim_diffs[-1], edges_ang_diffs[0], edges_ang_diffs[-1]], aspect=ASPECT, norm=LogNorm(vmin=1, vmax=VMAX))
  plt.colorbar()

  # Save 2D histogram to txt file:

  savetxt(result_figure_path + label + ".txt", H)

  # Ellipsoidal cut obtained using simulations with SF phantom and rod source
  # displaced by 25 cm from the axis of the scanner:

  xxx = linspace(0,ELLIPSE_PARAM,100)
  yyy = []
  for i in xrange(len(xxx)):
    yyy.append(ellipsoid_threshold(xxx[i]))

  if show_cut:
    plt.plot(xxx, yyy, color='r', linewidth=2)

  [counter_above, counter_below] = calculate_counters(tim_diffs, ang_diffs)
  print "Number of events above the cut: ", counter_above
  print "Number of events below the cut: ", counter_below
  print "Percentage of cut events: ", float(counter_below)/(counter_above+counter_below)*100.

  plt.ylim(ylim)

  plt.xlabel("Time difference [ns]")
  plt.ylabel("Angle difference [deg.]")
  plt.savefig(result_figure_path + label + OUTPUT_FORMAT, bbox_inches='tight')
  plt.clf()
  plt.close()

  H_tim_diffs = zeros([t_bins])
  H_ang_diffs = zeros([a_bins+1])

  T, A = shape(H)
  for t in range(T):
    for a in range(A):
      H_tim_diffs[t] += H[t,a]
      H_ang_diffs[a] += H[t,a]

  YLIM = 1e5

  times = []
  for i in range(len(edges_tim_diffs)-1):
    times.append((edges_tim_diffs[i]+edges_tim_diffs[i+1])/2.)

  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111)
  plt.subplots_adjust(left=0.1, right=0.95, bottom=0.15, top=0.92)
  ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  bar_width = (times[-1]-times[0])/(len(times)-1)
  plt.bar(times, H_tim_diffs, bar_width)
  plt.xlabel("Time difference [ns]")
  plt.xlim(0,5)
  plt.ylim(0,YLIM)
  plt.savefig(result_figure_path + label + "_1D_tim_diffs" + OUTPUT_FORMAT)
  plt.clf()
  plt.close()

  angles = []
  for i in range(len(edges_ang_diffs)-1):
    angles.append((edges_ang_diffs[i]+edges_ang_diffs[i+1])/2.)

  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111)
  plt.subplots_adjust(left=0.1, right=0.95, bottom=0.15, top=0.92)
  ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  bar_width = (angles[-1]-angles[0])/(len(angles)-1)
  plt.bar(angles, H_ang_diffs, width=bar_width)
  plt.xlabel("Angle difference [deg.]")
  plt.xlim(-bar_width/2.,180+bar_width/2.)
  plt.ylim(0,YLIM)
  plt.savefig(result_figure_path + label + "_1D_ang_diffs" + OUTPUT_FORMAT)
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

  # Plot 2D histogram:

  rcParams['font.size'] = 24
  rcParams['legend.fontsize'] = 18

  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111)
  plt.subplots_adjust(left=0.20, right=0.9, top=0.9, bottom=0.1)

  H, edges_tim_diffs, edges_ang_diffs = histogram2d(sourcePosX, sourcePosY, bins=100)
  VMAX = H.max()
  plt.imshow(H.T, interpolation='none', origin='low', extent=[edges_tim_diffs[0], edges_tim_diffs[-1], edges_ang_diffs[0], edges_ang_diffs[-1]], aspect='auto', norm=LogNorm(vmin=1, vmax=VMAX))
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

def plot_diff(hist1, hist2):

  hist = abs(hist2 - hist1)

  TLIM = 5.
  ALIM = 180.
  TBIN = TLIM/shape(hist)[0]
  ABIN = ALIM/shape(hist)[1]

  hist = delete(hist, (0), axis=1)

  rcParams['font.size'] = 24
  rcParams['legend.fontsize'] = 18

  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111)
  plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)

  ASPECT = TLIM/ALIM # force square pixels

  plt.imshow(hist.T, interpolation='none', origin='low', extent=[0, 5, 0, 180], aspect=ASPECT, norm=LogNorm(vmin=1, vmax=hist.max()))
  plt.colorbar()

  plt.xlabel("Time difference [ns]")
  plt.ylabel("Angle difference [deg.]")
  plt.savefig("./diff" + OUTPUT_FORMAT)
  plt.clf()
  plt.close()

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Analyze coincidences file and plot results of the analysis',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-cf', '--coincidences-file',
                      dest='path_coincidences_file',
                      default='',
                      type=str,
                      help='path to the coincidences file obtained using the GOJA tool')

  parser.add_argument('-tw', '--time-window',
                      dest='tw',
                      type=float,
                      default=3,
                      help='time window, for which the coincidence set was produced [ns]')

  parser.add_argument('-m', '--mode',
                      dest='mode',
                      type=str,
                      default="plot",
                      help='mode of the script: plot, plot_by_types, plot_diff, stats')

  parser.add_argument('-oat', '--output-da-dt',
                      dest='path_output_da_dt',
                      type=str,
                      default="./DA_vs_DT",
                      help='path to the output figure')

  parser.add_argument('-sc', '--show-cut',
                      dest='show_cut',
                      action='store_true',
                      help='set if you want to show cut')

  parser.add_argument('-osxsy', '--output-sposx-sposy',
                      dest='path_output_sposx_sposy',
                      type=str,
                      default="./sourcePosX_vs_sourcePosY",
                      help='path to the output figure sourcePosX vs. sourcePosY')

  parser.add_argument('-of', '--outputformat',
                      dest='outputformat',
                      type=str,
                      default="png",
                      help='output format of images')

  parser.add_argument('--t_bins',
                      dest='t_bins',
                      type=int,
                      default=100,
                      help='number of bins for times differences axis')

  parser.add_argument('--a_bins',
                      dest='a_bins',
                      type=int,
                      default=180,
                      help='number of bins for angles differences axis')

  parser.add_argument('-d', '--discrete',
                      dest='discrete',
                      action='store_true',
                      help='if set, then align x and y of the hits to strips centers')

  parser.add_argument('--hist1',
                      dest='hist1',
                      type=str,
                      help='path to hist1 (txt file)')

  parser.add_argument('--hist2',
                      dest='hist2',
                      type=str,
                      help='path to hist2 (txt file)')

  args = parser.parse_args()

  OUTPUT_FORMAT = "." + args.outputformat

  # Load data from file
  if args.path_coincidences_file != '':
    coincidences = loadtxt(args.path_coincidences_file)

  a_bin_size = 180./args.a_bins

  if args.mode == "plot":
    plot_Da_vs_Dt(coincidences, args.tw, args.path_output_da_dt, args.show_cut,
      args.t_bins, args.a_bins, ylim=[-a_bin_size/2.,180.+a_bin_size/2.], discrete=args.discrete)
    #plot_sourcePosX_vs_sourcePosY(coincidences, args.path_output_sposx_sposy)

  elif args.mode == "plot_by_types":
    plot_Da_vs_Dt(coincidences, args.tw, args.path_output_da_dt, args.show_cut,
      args.t_bins, args.a_bins, toc=1)
    try:
      plot_Da_vs_Dt(coincidences, args.tw, args.path_output_da_dt, args.show_cut,
        args.t_bins, args.a_bins, toc=2)
    except:
      print "There are no phantom coincidences in the goja output file."
    plot_Da_vs_Dt(coincidences, args.tw, args.path_output_da_dt, args.show_cut,
      args.t_bins, args.a_bins, toc=3)
    plot_Da_vs_Dt(coincidences, args.tw, args.path_output_da_dt, args.show_cut,
      args.t_bins, args.a_bins, toc=4)

  elif args.mode == "plot_diff":

    print str(args.hist1)

    hist1 = loadtxt(args.hist1)
    hist2 = loadtxt(args.hist2)

    plot_diff(hist1, hist2)

  elif args.mode == "stats":
    calculate_ratios(coincidences, args.path_coincidences_file.split("/")[-1])
