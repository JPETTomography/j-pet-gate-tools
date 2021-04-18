#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from numpy import array, savetxt

if __name__ == "__main__":

  missing_goja_results = []

  with open(".Gate/main/main.split") as f:
    lines = f.readlines()
    for line in lines:
      if line[0:14] == "Root filename:":
        output = line.split(' ')[2]
        path = output.replace('\n','') + '.root'
        if os.path.isfile(path):
          gate_result = int(output.split('output/output')[-1])
          goja_result = "./goja/output" + str(gate_result) + "_"
          if not os.path.isfile(goja_result + "coincidences") or \
             not os.path.isfile(goja_result + "realtime") or \
             not os.path.isfile(goja_result + "statistics"):
            missing_goja_result = gate_result
            missing_goja_results.append(missing_goja_result)

  savetxt("./missing_goja_results.txt", array(missing_goja_results).T, fmt="%d")
