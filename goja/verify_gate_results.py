#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from numpy import array, savetxt

if __name__ == "__main__":

  missing_gate_results = []

  with open(".Gate/main/main.split") as f:
    lines = f.readlines()
    for line in lines:
      if line[0:14] == "Root filename:":
        output = line.split(' ')[2]
        path = output.replace('\n','') + '.root'
        if not os.path.isfile(path):
          missing_gate_result = int(output.split('output/output')[-1])
          missing_gate_results.append(missing_gate_result)

  savetxt("./missing_gate_results.txt", array(missing_gate_results).T, fmt="%d")
