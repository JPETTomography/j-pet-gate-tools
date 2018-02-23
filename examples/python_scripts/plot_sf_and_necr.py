#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.pyplot as plt2
from numpy import *
from matplotlib import rcParams, rcParamsDefault
from matplotlib.colors import LogNorm

from nema_common import *

outputformat = ""

def plot_rates(geometry,float_activity,N_true,N_dsca,N_psca,N_acci,time):

  rcParams['font.size'] = 20
  rcParams['legend.fontsize'] = 16

  fig, axs = plt2.subplots(nrows=1, ncols=1, sharex=True)
  plt2.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.17)

  new_activity=linspace(float_activity[0],float_activity[-1],100)

  plt2.plot(float_activity,N_true/time/1000.,'o-',label="true", markersize=4)
  plt2.plot(float_activity,N_dsca/time/1000.,'o-',label="dsca", markersize=4)
  plt2.plot(float_activity,N_psca/time/1000.,'o-',label="psca", markersize=4)
  plt2.plot(float_activity,N_acci/time/1000.,'o-',label="acci", markersize=4)

  plt2.legend(loc=2)
  plt2.xlim(0,90)
  plt2.ylim(ymin=0)
  plt2.xlabel("Activity concentration [kBq/cc]")
  plt2.ylabel("Rate [kcps]")
  plt2.savefig(workdir_Results + geometry + "_rates" + outputformat)
  plt2.clf()
  plt2.close()

  rcParams.update(rcParamsDefault)

def plot_necrs(float_activities, NECRs, colors, labels, necr_type):

  fig, axs = plt.subplots(nrows=1, ncols=1, sharex=True)

  for i in range(len(NECRs)):
    plt.plot(float_activities[i], NECRs[i], 'o-', color=colors[i], label=labels[i], markersize=4)

  rcParams.update(rcParamsDefault)
  rcParams['legend.fontsize'] = 8
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
  plt.savefig(workdir_Results + "NECR_all_geometries_" + necr_type + outputformat)
  plt.clf()
  plt.close()

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Plot NECR.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-of', '--outputformat',
                      type=str,
                      default="png",
                      help='output format of images')

  args = parser.parse_args()
  outputformat = args.outputformat

  float_activity = []
  NECRs_sin = []
  NECRs_crt = []

  fig_ctr, axs_ctr = plt.subplots(nrows=1, ncols=1, sharex=True)
  #plt.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.17)

  float_activities = []
  NECRs_sin = []
  NECRs_ctr = []
  colors = []
  labels = []

  for geometry in geometries_NECR:

    tmp = loadtxt(workdir_Results + geometry + "/necr_dependency.txt")
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
        linestyle='-'
        new_label += "1 layer"
    else:
        linestyle = '--'
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

    float_activities.append(float_activity)
    NECRs_sin.append(NECR_sin)
    NECRs_ctr.append(NECR_ctr)
    colors.append(datacolor)
    labels.append(new_label)

  plot_necrs(float_activities, NECRs_sin, colors, labels, "sin")
  plot_necrs(float_activities, NECRs_ctr, colors, labels, "ctr")
