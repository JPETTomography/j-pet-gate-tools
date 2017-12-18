#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from numpy import *
from matplotlib import rcParams, rcdefaults, rc
from matplotlib.colors import LogNorm

import os

import matplotlib.image as mpimg

if __name__ == "__main__":

  root_directory = "../output/"

  os.system("rm ./coincidences.txt")

  nr_of_files = 100 # number of output*.root files (NR_OF_SPLITS from Gate_parallel.sh macro)

  for i in range(1,nr_of_files+1):

      root_file = root_directory + "output" + str(i) + ".root"

      command = "./goja --root " + root_file + " >> ./coincidences.txt"

      os.system(command)
