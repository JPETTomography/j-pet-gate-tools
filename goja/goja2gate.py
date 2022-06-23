#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from ctypes import c_int32, c_float
import numpy as np
import uproot

if __name__ == "__main__":

  help_message = "Type goja2gate.py --help."

  parser = argparse.ArgumentParser(description='Convert GOJA coincidences from cylindricalPET geometry singles to GATE coincidences tree.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-in', '--input-path',
                      dest='input_path',
                      type=str,
                      help='path to input GOJA file')

  parser.add_argument('-out', '--output-path',
                      dest='output_path',
                      type=str,
                      help='path to output GATE file')

  args = parser.parse_args()

  goja = np.loadtxt(args.input_path)
  if goja.shape[1]<39:
    print("GOJA input should contain 39 columns with data. Use newer GOJA version to generate an input.")
    exit(1)

  gate = uproot.recreate(args.output_path, compression=uproot.ZLIB(4))

  gate["Coincidences"] = {

    "globalPosX1": goja[:,0]*10,
    "globalPosY1": goja[:,1]*10,
    "globalPosZ1": goja[:,2]*10,
    "time1": goja[:,3]*1e-12,
    "globalPosX2": goja[:,4]*10,
    "globalPosY2": goja[:,5]*10,
    "globalPosZ2": goja[:,6]*10,
    "time2": goja[:,7]*1e-12,

    # Needed for classification:
    "eventID1": goja[:,19].astype(c_int32),
    "eventID2": goja[:,20].astype(c_int32),
    "comptonPhantom1": goja[:,21].astype(c_int32),
    "comptonPhantom2": goja[:,22].astype(c_int32),
    "RayleighPhantom1": goja[:,23].astype(c_int32),
    "RayleighPhantom2": goja[:,24].astype(c_int32),
    "comptonCrystal1": goja[:,25].astype(c_int32),
    "comptonCrystal2": goja[:,26].astype(c_int32),
    "RayleighCrystal1": goja[:,27].astype(c_int32),
    "RayleighCrystal2": goja[:,28].astype(c_int32),

    # Needed for obtaining IDs for CASToR:
    "rsectorID1": goja[:,29].astype(c_int32),
    "rsectorID2": goja[:,30].astype(c_int32),
    "moduleID1": goja[:,31].astype(c_int32),
    "moduleID2": goja[:,32].astype(c_int32),
    "submoduleID1": goja[:,33].astype(c_int32),
    "submoduleID2": goja[:,34].astype(c_int32),
    "crystalID1": goja[:,35].astype(c_int32),
    "crystalID2": goja[:,36].astype(c_int32),
    "layerID1": goja[:,37].astype(c_int32),
    "layerID2": goja[:,38].astype(c_int32),

    # Fake branches:
    "sourceID1": np.zeros(len(goja)).astype(c_int32),
    "sourceID2": np.zeros(len(goja)).astype(c_int32),
    "sourcePosX1": np.zeros(len(goja)).astype(c_float),
    "sourcePosY1": np.zeros(len(goja)).astype(c_float),
    "sourcePosZ1": np.zeros(len(goja)).astype(c_float)

  }
