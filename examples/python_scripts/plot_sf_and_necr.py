#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.pyplot as plt2
from numpy import *
from matplotlib import rcParams, rcParamsDefault
import argparse

from nema_common import *

outputformat = ""

def plot_rates(geometry,float_activity,N_true,N_dsca,N_psca,N_acci,time):

  rcParams['font.size'] = 20
  rcParams['legend.fontsize'] = 16

  fig, axs = plt2.subplots(nrows=1, ncols=1, sharex=True)
  plt2.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.17)

  plt2.plot(float_activity,N_true/time/1000.,'o-',label="true", markersize=4)
  plt2.plot(float_activity,N_dsca/time/1000.,'o-',label="dsca", markersize=4)
  plt2.plot(float_activity,N_psca/time/1000.,'o-',label="psca", markersize=4)
  plt2.plot(float_activity,N_acci/time/1000.,'o-',label="acci", markersize=4)

  plt2.legend(loc=2)
  plt2.xlim(0,90)
  plt2.ylim(ymin=0)
  plt2.xlabel("Activity concentration [kBq/cc]")
  plt2.ylabel("Rate [kcps]")
  plt2.savefig(workdir_NECR + geometry + "_rates." + outputformat)
  plt2.clf()
  plt2.close()

  rcParams.update(rcParamsDefault)

def plot_necrs(float_activities, NECRs, colors, labels, necr_type, lstyles):

  fig, axs = plt.subplots(nrows=1, ncols=1, sharex=True)

  for i in xrange(len(NECRs)):
    plt.plot(float_activities[i], NECRs[i], lstyles[i], color=colors[i], label=labels[i], markersize=4)

  rcParams.update(rcParamsDefault)
  rcParams['legend.fontsize'] = 11
  rcParams['font.size'] = 20

  FONTSIZE = 20

  plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.17)
  plt.legend(loc=1)
  plt.xlim(0,90)
  plt.ylim(0,1.1*NECR_sin_max)
  plt.xticks(fontsize=FONTSIZE)
  plt.yticks(fontsize=FONTSIZE)
  plt.xlabel("Activity concentration [kBq/cc]", fontsize=FONTSIZE)
  plt.ylabel("NECR [kcps]", fontsize=FONTSIZE)
  plt.savefig(workdir_NECR + "NECR_all_geometries_" + necr_type + '.' + outputformat)
  plt.clf()
  plt.close()

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
    elif "L200" in g:
      lab += "L = 200 cm"
      l = '..'

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
  plt.savefig(workdir_NECR + "second_lvl_selection" + outputformat, bbox_inches='tight')

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Plot NECR.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-of', '--outputformat',
                      type=str,
                      default="png",
                      help='output format of images')

  parser.add_argument('-r', '--reduction',
                      dest='reduction',
                      action='store_true',
                      help='set if you want to calculate and plot the reduction ' +
                           'given by the 2nd level of the reduction method in case ' +
                           'of the NECR simulations (must be used with -ns option)')

  parser.add_argument('-ns', '--necr-simulations',
                      dest='necr_simulations',
                      type=str,
                      help='path to the base directory of the NECR simulations')

  args = parser.parse_args()
  outputformat = args.outputformat

  if args.reduction:
    create_work_directories()
    calculate_reduction_for_necr_simulations(args.necr_simulations)
    plot_reduction_for_necr_simulations()

  else:
    plt.subplots(nrows=1, ncols=1, sharex=True)

    float_activities = []
    NECRs_sin = []
    NECRs_ctr = []
    colors = []
    labels = []
    lstyles = []

    for geometry in geometries_NECR:

      tmp = loadtxt(workdir_NECR + geometry + "/necr_dependency.txt")
      float_activity = tmp[:,0]
      SF_sin = tmp[:,1]
      SF_crt = tmp[:,2]
      NECR_sin = tmp[:,3]
      NECR_sin_max = max(NECR_sin)
      NECR_ctr = tmp[:,4]
      T = tmp[:,5]
      S = tmp[:,6]
      N_true = tmp[:,7]
      N_dsca = tmp[:,8]
      N_psca = tmp[:,9]
      N_acci = tmp[:,10]
      time = tmp[:,11]

      plot_rates(geometry,float_activity,N_true,N_dsca,N_psca,N_acci,time)

      new_label = ""
      if "1lay" in geometry:
          linestyle='o-'
          new_label += "1 layer"
      else:
          linestyle = 'o--'
          new_label += "2 layers"
      if "L020" in geometry:
          datacolor = 'r'
          new_label += ", L = 20 cm"
      elif "L050" in geometry:
          datacolor = 'b'
          new_label += ", L = 50 cm"
      elif "L100" in geometry:
          datacolor = 'y'
          new_label += ", L = 100 cm"
      elif "L200" in geometry:
          datacolor = 'g'
          new_label += ", L = 200 cm"

      float_activities.append(float_activity)
      NECRs_sin.append(NECR_sin)
      NECRs_ctr.append(NECR_ctr)
      colors.append(datacolor)
      labels.append(new_label)
      lstyles.append(linestyle)

    plot_necrs(float_activities, NECRs_sin, colors, labels, "sin", lstyles)
    plot_necrs(float_activities, NECRs_ctr, colors, labels, "ctr", lstyles)
