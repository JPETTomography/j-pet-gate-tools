#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

if __name__ == "__main__":

  # Analyze Sensitivity simulation files:
  os.system("./calculate_sensitivity_and_sensitivity_profiles.py")
  os.system("./plot_sensitivity_and_sensitivity_profiles.py")

  # Analyze Scatter Fraction simulation files:

  # Analyze NECR simulation files:
  os.system("./calculate_sf_and_necr.py")
  os.system("./plot_sf_and_necr.py")
