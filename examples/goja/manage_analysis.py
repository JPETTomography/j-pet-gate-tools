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

  for i in range(100):

      i = i+1

      root_file = root_directory + "output" + str(i) + ".root"

      command = "./goja --root " + root_file + " >> ./coincidences.txt"
      os.system(command)