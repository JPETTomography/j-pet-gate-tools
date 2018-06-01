#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from numpy import *
from matplotlib import rcParams
import operator
from matplotlib.ticker import FormatStrFormatter
import argparse

from nema_common import *

def create_label (template, label_old):

    label_new = ""
    tmp = label_old.split('_')

    flag = 0
    for i in xrange(len(tmp)):
        if tmp[i]!=template and tmp[i]!="7mm":
            if tmp[i]=="D75":
                label_new += "D=75cm"
            elif tmp[i]=="D85":
                label_new += "D=85cm"
            elif tmp[i]=="D95":
                label_new += "D=95cm"
            elif tmp[i]=="1lay":
                label_new += "1 layer"
            elif tmp[i]=="2lay":
                label_new += "2 layers"
            elif tmp[i]=="L020":
                label_new += "L=20cm"
            elif tmp[i]=="L050":
                label_new += "L=50cm"
            elif tmp[i]=="L100":
                label_new += "L=100cm"
            elif tmp[i]=="L200":
                label_new += "L=200cm"
            else:
                label_new += tmp[i]
            if (flag==0):
                label_new += ', '
                flag = 1

    return label_new

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Plot sensitivity and sensitivity profiles.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-b', '--blurred',
                      type=int,
                      default=0,
                      help='for plotting blurred sensitivity profiles set to 1')

  parser.add_argument('-of', '--outputformat',
                      type=str,
                      default="png",
                      help='output format of images')

  args = parser.parse_args()

  rcParams['font.size'] = 24
  rcParams['legend.fontsize'] = 10

  data = {}

  for geometry in geometries_sensitivity:

    sensitivity = float(loadtxt(workdir_Sensitivity + geometry + "_sensitivity.txt"))
    sensitivityProfiles = loadtxt(workdir_Sensitivity + geometry + "_sensitivity_profiles.txt")
    sensitivityProfile = sensitivityProfiles[:,0]
    sensitivityProfile_PMT = sensitivityProfiles[:,1]

    data[geometry] = [sensitivity, sensitivityProfile, sensitivityProfile_PMT]

  N = 0

  templates = ["D75", "D85", "D95", "L020", "L050", "L100", "L200", "1lay", "2lay"]

  titles = ["D = 75 cm", "D = 85 cm", "D = 95 cm",
            "L = 20 cm", "L = 50 cm", "L = 100 cm", "L = 200 cm",
            "1 layer", "2 layers"]

  for i in xrange(len(templates)):

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    plt.subplots_adjust(left=0.15, right=0.99, top=0.91, bottom=0.17)

    t = templates[i]
    data_tmp = {}

    for d in data:

      if t in d:
        data_tmp[d] = data[d]

    data_sorted = sorted(data_tmp.items(), key=operator.itemgetter(0))

    plt.clf()

    for d in data_sorted:

      if "1lay" in d[0]:
        lw = 1
      elif "2lay" in d[0]:
        lw = 3

      if "D75" in d[0]:
          ls = '--'
      elif "D85" in d[0]:
          ls = '-'
      elif "D95" in d[0]:
          ls = '-.'

      if "L020" in d[0]:
          N = 20
          c = 'red'
      elif "L050" in d[0]:
          N = 50
          c = 'black'
      elif "L100" in d[0]:
          N = 100
          c = 'blue'
      elif "L200" in d[0]:
          N = 200
          c = 'green'

      arguments = linspace(-N/2.+0.5,N/2-0.5, N)

      rcParams['font.size'] = 30
      rcParams['legend.fontsize'] = 18

      # Group specific settings
      if t in ["D75", "D85", "D95"]:
        plt.ylim([0,25])

      label_blurred = ""
      if not args.blurred:
        plt.plot(arguments, d[1][1], label=create_label(t,d[0]), linewidth=lw, linestyle=ls, color=c)
      else:
        plt.plot(arguments, d[1][2], label=create_label(t,d[0]), linewidth=lw, linestyle=ls, color=c)
        label_blurred = "_blurred"

    if t=="L020":
      plt.gca().get_yaxis().set_major_formatter(FormatStrFormatter('%d'))
      plt.gca().get_yaxis().set_ticks([0,1])
    elif t=="L100":
      plt.gca().get_yaxis().set_major_formatter(FormatStrFormatter('%d'))
      plt.gca().get_yaxis().set_ticks([0,5,10,15,20])

    plt.title(titles[i] + ", 7 mm")
    plt.xlabel("Position along z axis [cm]")
    plt.ylabel("Sensitivity [cps/kBq]")
    #plt.legend(loc=1)
    plt.xlim(-35,35)
    plt.ylim(ymin=0)
    plt.savefig(workdir_Sensitivity + "sensitivityProfile_" + t + label_blurred + "." + args.outputformat)

    # create a second figure for the legend
    #figlegend = pylab.figure(figsize = (8,6))
    # produce a legend for the objects in the other figure
    #pylab.figlegend(*ax.get_legend_handles_labels(), loc = 'upper left')
    #figlegend.savefig(workdir + "sensitivityProfile_" + t + '_legend' + outputformat)
    #figlegend.clf()

    plt.clf()

  # Plot sensitivities:

  MARKERSIZE = 10

  L = [20,50,100,200]
  D75_1lay = [data["D75_1lay_L020_7mm"][0], data["D75_1lay_L050_7mm"][0], data["D75_1lay_L100_7mm"][0], data["D75_1lay_L200_7mm"][0]]
  D85_1lay = [data["D85_1lay_L020_7mm"][0], data["D85_1lay_L050_7mm"][0], data["D85_1lay_L100_7mm"][0], data["D85_1lay_L200_7mm"][0]]
  D95_1lay = [data["D95_1lay_L020_7mm"][0], data["D95_1lay_L050_7mm"][0], data["D95_1lay_L100_7mm"][0], data["D95_1lay_L200_7mm"][0]]
  D75_2lay = [data["D75_2lay_L020_7mm"][0], data["D75_2lay_L050_7mm"][0], data["D75_2lay_L100_7mm"][0], data["D75_2lay_L200_7mm"][0]]
  D85_2lay = [data["D85_2lay_L020_7mm"][0], data["D85_2lay_L050_7mm"][0], data["D85_2lay_L100_7mm"][0], data["D85_2lay_L200_7mm"][0]]
  D95_2lay = [data["D95_2lay_L020_7mm"][0], data["D95_2lay_L050_7mm"][0], data["D95_2lay_L100_7mm"][0], data["D95_2lay_L200_7mm"][0]]

  plt.subplots_adjust(left=0.19, right=0.99, top=0.97, bottom=0.17)
  plt.plot(L, D75_1lay, '<', markersize=MARKERSIZE, color='r', label="D=75cm, 1 layer")
  plt.plot(L, D85_1lay, 's', markersize=MARKERSIZE, color='r', label="D=85cm, 1 layer")
  plt.plot(L, D95_1lay, '+', markersize=MARKERSIZE, color='r', label="D=95cm, 1 layer")
  plt.plot(L, D75_2lay, '>', markersize=MARKERSIZE, color='k', label="D=75cm, 2 layers")
  plt.plot(L, D85_2lay, 'D', markersize=MARKERSIZE, color='k', label="D=85cm, 2 layers")
  plt.plot(L, D95_2lay, 'P', markersize=MARKERSIZE, color='k', label="D=95cm, 2 layers")
  plt.xlabel("Scintillator length [cm]")
  plt.ylabel("Sensitivity [cps/kBq]")
  rcParams['legend.fontsize'] = 16
  plt.legend(loc=2)
  plt.ylim(ymin=0)
  plt.savefig(workdir_Sensitivity + "Sensitivities." + args.outputformat)
  plt.clf()

  # Plot sensitivities ratios (between 1- and 2-layer geometries):

  D75_ratio = array(D75_2lay)/array(D75_1lay)
  D85_ratio = array(D85_2lay)/array(D85_1lay)
  D95_ratio = array(D95_2lay)/array(D95_1lay)

  plt.subplots_adjust(left=0.19, right=0.99, top=0.97, bottom=0.17)
  plt.plot(L, D75_ratio, '<', markersize=MARKERSIZE, color='r', label = "D=75cm")
  plt.plot(L, D85_ratio, 's', markersize=MARKERSIZE, color='r', label = "D=85cm")
  plt.plot(L, D95_ratio, '+', markersize=MARKERSIZE, color='r', label = "D=95cm")
  plt.xlabel("Scintillator length [cm]")
  plt.ylabel("Ratio")
  rcParams['legend.fontsize'] = 18
  plt.legend(loc=1)
  plt.savefig(workdir_Sensitivity + "Sensitivities_ratios." + args.outputformat)
  plt.clf()

  # Print mean and max values:

  print "GEOMETRY         \tMEAN\tMAX"
  for geometry in geometries_sensitivity:
    print "%s\t%.2f\t%.2f" % (geometry, mean(data[geometry][1]), max(data[geometry][1]))
