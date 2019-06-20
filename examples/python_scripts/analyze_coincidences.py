#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from math import *
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt2
from matplotlib import rcParams
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter
from numpy import *
import sys

from nema_common import *

OUTPUT_FORMAT = ".png"
NR_OF_STRIPS_LAB192 = 192 # nr of strips in the 3-layer prototype of the J-PET scanner

def get_strips_centers(geometries_set):

  xs = []
  ys = []

  if geometries_set == "lab192":

    r1 = 42.5*10
    r2 = 46.75*10
    r3 = 57.5*10

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

  elif geometries_set == "djpet-total-body":

    r1 = 408.08027535 # in mm, radius of inside layer
    r2 = 443.08027535 # in mm, radius of outside layer
    n = 24 # number of modules in each layer
    alpha = 2*pi/n

    strips_per_module = 16
    thickness = 6. # in mm, thickness of a single scintillator
    distance = 0.5 # in mm, distance between the strips in a module
    module_width = strips_per_module*thickness + (strips_per_module-1)*distance
    strips_ys = linspace(-module_width/2+thickness/2., module_width/2.-thickness/2., strips_per_module)

    rs = [r1, r2]

    for r in rs:

      for y in strips_ys:
        xs.append(r)
        ys.append(y)

      tmp_xs = xs
      tmp_ys = ys

      for i in range(n-1):
        xprims = []
        yprims = []
        for t in range(len(tmp_xs)):
          xprim = tmp_xs[t]*cos(alpha) - tmp_ys[t]*sin(alpha)
          yprim = tmp_xs[t]*sin(alpha) + tmp_ys[t]*cos(alpha)
          xprims.append(xprim)
          yprims.append(yprim)
        tmp_xs = xprims
        tmp_ys = yprims
        for t in range(len(tmp_xs)):
          xs.append(xprims[t])
          ys.append(yprims[t])

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
def plot_Da_vs_Dt(coincidences,
                  geometries_set,
                  tw,
                  result_figure_path,
                  show_cut,
                  t_bins,
                  a_bins,
                  ylim=[0,180],
                  toc=0,
                  discrete=False):

  posX1 = coincidences[:,0]
  posY1 = coincidences[:,1]
  times1 = coincidences[:,3]
  posX2 = coincidences[:,4]
  posY2 = coincidences[:,5]
  times2 = coincidences[:,7]
  type_of_coincidence = coincidences[:,12]

  volhits = coincidences[:,8]
  vol2 = coincidences[:,9]

  if discrete:

    if not os.path.exists("./strips_centers_" + geometries_set + ".txt"):

      print "You need to run calculate_strips_centers mode firstly."

    else:

      strips_centers = loadtxt("./strips_centers_" + geometries_set + ".txt")

      # Overwrite x and y hits positions with hits centers
      for i in range(len(coincidences)):
        ind1 = int(volhits[i])-1
        posX1[i] = strips_centers[:,1][ind1]/10.
        posY1[i] = strips_centers[:,2][ind1]/10.
        ind2 = int(vol2[i])-1
        posX2[i] = strips_centers[:,1][ind2]/10.
        posY2[i] = strips_centers[:,2][ind2]/10.

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
  ASPECT = tw/float(ylim[1]-ylim[0]) # force square pixels

  fig = plt.figure(figsize=(8, 6))
  ax = fig.add_subplot(111)
  plt.subplots_adjust(left=0.20, right=0.9, top=0.9, bottom=0.1)

  a_bin_size = 180./a_bins
  a_bins_vec = linspace(0-a_bin_size/2., 180.+a_bin_size/2., a_bins+2)
  t_bins_vec = linspace(0, tw, t_bins+1)

  H, edges_tim_diffs, edges_ang_diffs = histogram2d(tim_diffs, ang_diffs, bins=(t_bins_vec,a_bins_vec), range=[[0, tw],ylim])
  VMAX = H.max()
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

  YLIM = 1.1 * max(H_tim_diffs)

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
  #plt.xlim(0,5)
  plt.ylim(0,YLIM)
  plt.savefig(result_figure_path + label + "_1D_tim_diffs" + OUTPUT_FORMAT)
  plt.clf()
  plt.close()

  YLIM = 1.1 * max(H_ang_diffs)

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

  ASPECT = TLIM/float(ALIM) # force square pixels

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
                      help='mode of the script: plot, plot_by_types, plot_diff, stats, plot_strips_centers, calculate_strips_centers')

  parser.add_argument('-gs', '--geometries-set',
                      dest='geometries_set',
                      type=str,
                      default='djpet-total-body',
                      help='geometries set: djpet-total-body, lab')

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

  rcParams['font.size'] = 24
  rcParams['legend.fontsize'] = 18

  OUTPUT_FORMAT = "." + args.outputformat

  # Load data from file
  if args.path_coincidences_file != '':
    coincidences = loadtxt(args.path_coincidences_file)

  a_bin_size = 180./args.a_bins

  if args.mode == "plot":

    "Mode to plot all coincidences from the coincidence file"

    plot_Da_vs_Dt(coincidences, args.geometries_set, args.tw, args.path_output_da_dt, args.show_cut,
      args.t_bins, args.a_bins, ylim=[-a_bin_size/2.,180.+a_bin_size/2.], discrete=args.discrete)
    plot_sourcePosX_vs_sourcePosY(coincidences, args.path_output_sposx_sposy)

  elif args.mode == "plot_by_types":

    "Mode to plot all coincidences from the coincidence file by type"

    plot_Da_vs_Dt(coincidences, args.geometries_set, args.tw, args.path_output_da_dt, args.show_cut,
      args.t_bins, args.a_bins, toc=1)
    try:
      plot_Da_vs_Dt(coincidences, args.geometries_set, args.tw, args.path_output_da_dt, args.show_cut,
        args.t_bins, args.a_bins, toc=2)
    except:
      print "There are no phantom coincidences in the goja output file."
    plot_Da_vs_Dt(coincidences, args.geometries_set, args.tw, args.path_output_da_dt, args.show_cut,
      args.t_bins, args.a_bins, toc=3)
    plot_Da_vs_Dt(coincidences, args.geometries_set, args.tw, args.path_output_da_dt, args.show_cut,
      args.t_bins, args.a_bins, toc=4)

  elif args.mode == "plot_diff":

    "Mode to plot differential image of two DA vs. DT histograms given with commandline arguments hist1 and hist2"

    hist1 = loadtxt(args.hist1)
    hist2 = loadtxt(args.hist2)
    plot_diff(hist1, hist2)

  elif args.mode == "stats":

    "Mode to calculate some statistics from the coincidence file"

    calculate_ratios(coincidences, args.path_coincidences_file.split("/")[-1])

  elif args.mode == "plot_strips_centers":

    geometries_sets = ["lab192", "djpet-total-body"]

    for geometries_set in geometries_sets:

      strips_centers = get_strips_centers(geometries_set)

      fig = plt.figure(figsize=(8, 6))
      ax = fig.add_subplot(111)
      ax.set_aspect(aspect=1.)
      plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
      plt.plot(array(strips_centers)[:,0], array(strips_centers)[:,1], 'o', markersize=0.1)
      plt.xlabel('x [mm]')
      plt.ylabel('y [mm]')
      plt.savefig("./strips_centers_" + geometries_set + OUTPUT_FORMAT)

  elif args.mode == "calculate_strips_centers":

    if args.path_coincidences_file == '':
      print "You need to provide the coincidences file for the " + args.geometries_set \
            + " geometries set (the best would be from the point source simulation)."
      sys.exit(1)

    thickness = 0.
    nr_of_strips = 0

    if args.geometries_set == "lab192":
      thickness = 7. # in mm
      nr_of_strips = NR_OF_STRIPS_LAB192

    elif args.geometries_set == "djpet-total-body":
      thickness = 6. # in mm
      nr_of_strips = 24*2*16

    else:
     print "Not supported geometries set."
     sys.exit(1)

    nrs_of_hits = zeros(nr_of_strips)
    xs = zeros(nr_of_strips)
    ys = zeros(nr_of_strips)

    xhits = concatenate([coincidences[:,0], coincidences[:,4]])*10. # in mm
    yhits = concatenate([coincidences[:,1], coincidences[:,5]])*10. # in mm
    volhits = concatenate([coincidences[:,8], coincidences[:,9]])

    for i in xrange(len(coincidences)):
      index = int(volhits[i])-1
      nrs_of_hits[index] += 1
      xs[index] += xhits[i]
      ys[index] += yhits[i]
    xs = xs/nrs_of_hits
    ys = ys/nrs_of_hits

    strips_centers = get_strips_centers(args.geometries_set)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    ax.set_aspect(aspect=1.)
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
    plt.plot(xhits, yhits, '.', markersize=0.001)
    plt.plot(xs, ys, '.', markersize=0.5, color='green', label='average hits positions')
    plt.plot(array(strips_centers)[:,0], array(strips_centers)[:,1], '.', markersize=0.5, color='red', label='strips centers')
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')
    plt.savefig("./strips_centers_vs_averages_vs_hits_" + args.geometries_set + OUTPUT_FORMAT)

    identifiers = {}
    for i in xrange(len(coincidences)):
      x, y = get_closest_strip_center((xhits[i], yhits[i]), strips_centers)
      d = sqrt((x-xhits[i])**2 + (y-yhits[i])**2)
      if d<thickness/2.:
        identifiers[int(volhits[i])] = [x, y]
        if len(identifiers) == len(xs):
          break

    ids = []
    strips_centers_x = []
    strips_centers_y = []
    averages_x = []
    averages_y = []

    for i in sorted (identifiers.keys()):
      ids.append(i)
      strips_centers_x.append(identifiers[i][0])
      strips_centers_y.append(identifiers[i][1])
      averages_x.append(xs[i-1])
      averages_y.append(ys[i-1])

    strips_centers_txt = "./strips_centers_" + args.geometries_set + ".txt"
    savetxt(strips_centers_txt, array([ids, strips_centers_x, strips_centers_y]).T, fmt='%d\t%.2f\t%.2f')

    fig2 = plt2.figure(figsize=(8, 6))
    ax2 = fig2.add_subplot(111)
    plt2.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
    ax2.set_aspect(aspect=1.)
    for i in range(len(ids)):
      plt2.plot([strips_centers_x[i], averages_x[i]], [strips_centers_y[i], averages_y[i]], '-', color='k')
    plt2.xlabel('x [mm]')
    plt2.ylabel('y [mm]')
    plt2.savefig("./strips_centers_vs_averages_differences_" + args.geometries_set + OUTPUT_FORMAT)
