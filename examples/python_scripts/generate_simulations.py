#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

from nema_common import *

workdir_Generated = './Generated'
workdir_Simulations = '../simulations/'

NR_OF_SPLITS = 100
NR_OF_ANNIHILATIONS = 10000
SOURCE_LENGTH = 70.
ENABLE_VISUALISATION = False

def copy_to_destination(what, where):

  command_cp = 'cp ' + workdir_Simulations + what + ' ' + where
  os.system(command_cp)

def generate_sf_necr_phantom(length):

  phantom = ""

  phantom += "/gate/world/daughters/name phantom\n"
  phantom += "/gate/world/daughters/insert cylinder\n"
  phantom += "/gate/phantom/geometry/setRmax 101.5 mm\n"
  phantom += "/gate/phantom/geometry/setHeight " + str(length) + " cm\n"
  phantom += "/gate/phantom/setMaterial HDPE\n"
  phantom += "/gate/phantom/vis/setVisible 1\n"
  phantom += "/gate/phantom/vis/forceWireframe\n"
  phantom += "/gate/phantom/vis/setColor blue\n"

  phantom += "/gate/phantom/daughters/name hole\n"
  phantom += "/gate/phantom/daughters/insert cylinder\n"
  phantom += "/gate/hole/setMaterial Air\n"
  phantom += "/gate/hole/geometry/setRmax 3.2 mm\n"
  phantom += "/gate/hole/geometry/setHeight " + str(length) + " cm\n"
  phantom += "/gate/hole/placement/setTranslation 0 -45 0 mm\n"
  phantom += "/gate/hole/vis/setVisible 1\n"
  phantom += "/gate/hole/vis/forceWireframe\n"
  phantom += "/gate/hole/vis/setColor blue\n"

  phantom += "/gate/phantom/daughters/name sourceinsert_tube\n"
  phantom += "/gate/phantom/daughters/insert cylinder\n"
  phantom += "/gate/sourceinsert_tube/setMaterial HDPE\n"
  phantom += "/gate/sourceinsert_tube/geometry/setRmin 1.6 mm\n"
  phantom += "/gate/sourceinsert_tube/geometry/setRmax 2.4 mm\n"
  phantom += "/gate/sourceinsert_tube/geometry/setHeight " + str(length) + " cm\n"
  phantom += "/gate/sourceinsert_tube/placement/setTranslation 0 -45 0 mm\n"
  phantom += "/gate/sourceinsert_tube/vis/setVisible 1\n"
  phantom += "/gate/sourceinsert_tube/vis/forceWireframe\n"
  phantom += "/gate/sourceinsert_tube/vis/setColor white\n"

  phantom += "/gate/phantom/daughters/name sourceinsert_water\n"
  phantom += "/gate/phantom/daughters/insert cylinder\n"
  phantom += "/gate/sourceinsert_water/setMaterial Water\n"
  phantom += "/gate/sourceinsert_water/geometry/setRmax 1.6 mm\n"
  phantom += "/gate/sourceinsert_water/geometry/setHeight " + str(length) + " cm\n"
  phantom += "/gate/sourceinsert_water/placement/setTranslation 0 -45 0 mm\n"
  phantom += "/gate/sourceinsert_water/vis/setVisible 1\n"
  phantom += "/gate/sourceinsert_water/vis/forceWireframe\n"
  phantom += "/gate/sourceinsert_water/vis/setColor blue\n"

  phantom += "/gate/phantom/attachPhantomSD\n"

  return phantom

def generate_sf_necr_source(length, activity):

  source = ""

  source += '/gate/source/addSource linear gps\n'
  source += '/gate/source/linear/gps/type Volume\n'
  source += '/gate/source/linear/gps/shape Cylinder\n'
  source += '/gate/source/linear/gps/radius 1.6 mm\n'
  source += '/gate/source/linear/gps/halfz ' + str(length/2.) + ' cm\n'
  source += '/gate/source/linear/gps/pos/centre 0. -45. 0. mm\n'
  source += '/gate/source/linear/setActivity ' + str(activity) + ' Bq\n'
  source += '/gate/source/linear/setType backtoback\n'
  source += '/gate/source/linear/gps/particle gamma\n'
  source += '/gate/source/linear/gps/ene/mono 511. keV\n'
  source += '/gate/source/linear/gps/angtype iso\n'
  source += '/gate/source/list\n'

  return source

def generate_sensitivity_source(length, activity):

  source = ""

  source += '/gate/source/addSource linear gps\n'
  source += '/gate/source/linear/gps/type Volume\n'
  source += '/gate/source/linear/gps/shape Cylinder\n'
  source += '/gate/source/linear/gps/radius 0.5 mm\n'
  source += '/gate/source/linear/gps/halfz ' + str(length/2.) + ' cm\n'
  source += '/gate/source/linear/gps/pos/centre 0. 0. 0. mm\n'
  source += '/gate/source/linear/setActivity ' + str(activity) + ' Bq\n'
  source += '/gate/source/linear/setType backtoback\n'
  source += '/gate/source/linear/gps/particle gamma\n'
  source += '/gate/source/linear/gps/ene/mono 511. keV\n'
  source += '/gate/source/linear/gps/angtype iso\n'
  source += '/gate/source/list\n'

  return source

def generate_gate_parallel_sh(nr_of_splits):

  gate_parallel_sh = ""

  gate_parallel_sh += '#!/bin/bash\n'
  gate_parallel_sh += 'NR_OF_SPLITS=\"' + str(nr_of_splits) + '\"\n'
  gate_parallel_sh += 'module load gate/7.2\n'
  gate_parallel_sh += 'export GC_DOT_GATE_DIR=$PWD\n'
  gate_parallel_sh += 'export GC_GATE_EXE_DIR=$GATE_DIR/bin/\n'
  gate_parallel_sh += 'if [[ $PWD/ != $GC_DOT_GATE_DIR ]]; then\n'
  gate_parallel_sh += '\tcd $GC_DOT_GATE_DIR\n'
  gate_parallel_sh += 'fi\n'
  gate_parallel_sh += 'rm -rf .Gate *.submit\n'
  gate_parallel_sh += 'qstat -a\n'
  gate_parallel_sh += 'gjs -clusterplatform openPBS -openPBSscript ./array.pbs -numberofsplits $NR_OF_SPLITS main.mac\n'
  gate_parallel_sh += 'qsub -t 1-$NR_OF_SPLITS array.pbs\n'
  gate_parallel_sh += 'qstat -t -u $USER\n'

  return gate_parallel_sh

def generate_array_pbs(geometry, activity):

  array_pbs = ""

  array_pbs += '#!/bin/sh\n'
  array_pbs += '#PBS -q i3d\n'
  array_pbs += '#PBS -l nodes=1:ppn=1\n'
  array_pbs += '#PBS -N ' + geometry + '_' + activity + '\n'
  array_pbs += '#PBS -V\n'
  array_pbs += 'cd ${PBS_O_WORKDIR}\n'
  array_pbs += 'Gate ./.Gate/main/main${PBS_ARRAYID}.mac\n'
  array_pbs += 'exit 0;\n'

  return array_pbs

def geometry_to_strip_length(geometry):

  return float((geometry.split('_')[2]).split('L')[1])

def generate_necr_simulations():

  workdir_Generated_NECR = workdir_Generated + '/NECR'
  if (os.path.isdir(workdir_Generated_NECR)):
    os.system('rm -rf ' + workdir_Generated_NECR)
  os.system('mkdir ' + workdir_Generated_NECR)

  for g in geometries_NECR:
    os.system('mkdir ' + workdir_Generated_NECR + '/' + g)

    for a in activities_NECR:
      os.system('mkdir ' + workdir_Generated_NECR + '/' + g + '/' + a)
      os.system('mkdir ' + workdir_Generated_NECR + '/' + g + '/' + a + '/output')

      where = workdir_Generated_NECR + '/' + g + '/' + a + '/'
      copy_to_destination('materials/GateMaterials.db', where)
      copy_to_destination('materials/Materials.xml', where)
      copy_to_destination('materials/Surfaces.xml', where)
      copy_to_destination('geometries/' + g.split('_')[0] + '/Geometry_' + g + '.mac', where)
      copy_to_destination('Visualisation.mac', where)
      copy_to_destination('physics/PhysicsEMLivermorePolar.mac', where)

      # nr of annihilations in the axial field of view
      time_simulated = NR_OF_ANNIHILATIONS/(float(a)*1e6)
      STRIP_LENGTH = geometry_to_strip_length(g)
      if STRIP_LENGTH<SOURCE_LENGTH:
        time_simulated = time_simulated*float(SOURCE_LENGTH)/float(STRIP_LENGTH)
      time_slice = time_simulated/NR_OF_SPLITS

      with open(where + '/NEMA_IEC_2001_SF_CL_CR_Phantom.mac', 'a') as phantom_mac:
        phantom_mac.write(generate_sf_necr_phantom(SOURCE_LENGTH))

      with open(where + '/NEMA_IEC_2001_SF_CL_CR_Source_Gamma.mac', 'a') as source_mac:
        source_mac.write(generate_sf_necr_source(SOURCE_LENGTH, int(float(a)*1e6)))

      with open(where + '/main.mac', 'a') as main_mac:
        if (ENABLE_VISUALISATION):
          main_mac.write('/vis/enable\n')
          main_mac.write('/control/execute ./Visualisation.mac\n')
        else:
          main_mac.write('/vis/disable\n')
        main_mac.write('/gate/geometry/setMaterialDatabase ./GateMaterials.db\n')
        main_mac.write('/control/execute ./Geometry_' + g + '.mac\n')
        main_mac.write('/control/execute ./NEMA_IEC_2001_SF_CL_CR_Phantom.mac\n')
        main_mac.write('/control/execute ./PhysicsEMLivermorePolar.mac\n')
        main_mac.write('/gate/run/initialize\n')
        main_mac.write('/control/execute ./NEMA_IEC_2001_SF_CL_CR_Source_Gamma.mac\n')
        main_mac.write('/gate/output/root/enable\n')
        main_mac.write('/gate/output/root/setFileName output/output\n')
        main_mac.write('/gate/output/root/setRootHitFlag 1\n')
        main_mac.write('/gate/random/setEngineName MersenneTwister\n')
        main_mac.write('/gate/random/setEngineSeed auto\n')
        main_mac.write('/gate/application/setTimeSlice ' + str('{0:.10f}'.format(time_slice)) + ' s\n')
        main_mac.write('/gate/application/setTimeStart 0. s\n')
        main_mac.write('/gate/application/setTimeStop ' + str('{0:.10f}'.format(time_simulated)) + ' s\n')
        main_mac.write('/gate/application/startDAQ\n')

      with open(where + '/array.pbs', 'a') as array_pbs:
        array_pbs.write(generate_array_pbs(g, a))

      with open(where + '/Gate_parallel.sh', 'a') as gate_parallel_sh:
        gate_parallel_sh.write(generate_gate_parallel_sh(NR_OF_SPLITS))

def generate_sf_simulations():

  workdir_Generated_SF = workdir_Generated + '/SF'
  if (os.path.isdir(workdir_Generated_SF)):
    os.system('rm -rf ' + workdir_Generated_SF)
  os.system('mkdir ' + workdir_Generated_SF)

  for g in geometries_NECR:
    os.system('mkdir ' + workdir_Generated_SF + '/' + g)
    os.system('mkdir ' + workdir_Generated_SF + '/' + g + '/output')

    where = workdir_Generated_SF + '/' + g + '/'
    copy_to_destination('materials/GateMaterials.db', where)
    copy_to_destination('materials/Materials.xml', where)
    copy_to_destination('materials/Surfaces.xml', where)
    copy_to_destination('geometries/' + g.split('_')[0] + '/Geometry_' + g + '.mac', where)
    copy_to_destination('Visualisation.mac', where)
    copy_to_destination('physics/PhysicsEMLivermorePolar.mac', where)

    # nr of annihilations in the axial field of view
    a = 1000. # in Bq
    time_simulated = NR_OF_ANNIHILATIONS/a
    STRIP_LENGTH = geometry_to_strip_length(g)
    if STRIP_LENGTH<SOURCE_LENGTH:
      time_simulated = time_simulated*float(SOURCE_LENGTH)/float(STRIP_LENGTH)
    time_slice = time_simulated/NR_OF_SPLITS

    with open(where + '/NEMA_IEC_2001_SF_CL_CR_Phantom.mac', 'a') as phantom_mac:
      phantom_mac.write(generate_sf_necr_phantom(SOURCE_LENGTH))

    with open(where + '/NEMA_IEC_2001_SF_CL_CR_Source_Gamma.mac', 'a') as source_mac:
      source_mac.write(generate_sf_necr_source(SOURCE_LENGTH, int(a)))

    with open(where + '/main.mac', 'a') as main_mac:
      if (ENABLE_VISUALISATION):
        main_mac.write('/vis/enable\n')
        main_mac.write('/control/execute ./Visualisation.mac\n')
      else:
        main_mac.write('/vis/disable\n')
      main_mac.write('/gate/geometry/setMaterialDatabase ./GateMaterials.db\n')
      main_mac.write('/control/execute ./Geometry_' + g + '.mac\n')
      main_mac.write('/control/execute ./NEMA_IEC_2001_SF_CL_CR_Phantom.mac\n')
      main_mac.write('/control/execute ./PhysicsEMLivermorePolar.mac\n')
      main_mac.write('/gate/run/initialize\n')
      main_mac.write('/control/execute ./NEMA_IEC_2001_SF_CL_CR_Source_Gamma.mac\n')
      main_mac.write('/gate/output/root/enable\n')
      main_mac.write('/gate/output/root/setFileName output/output\n')
      main_mac.write('/gate/output/root/setRootHitFlag 1\n')
      main_mac.write('/gate/random/setEngineName MersenneTwister\n')
      main_mac.write('/gate/random/setEngineSeed auto\n')
      main_mac.write('/gate/application/setTimeSlice ' + str('{0:.10f}'.format(time_slice)) + ' s\n')
      main_mac.write('/gate/application/setTimeStart 0. s\n')
      main_mac.write('/gate/application/setTimeStop ' + str('{0:.10f}'.format(time_simulated)) + ' s\n')
      main_mac.write('/gate/application/startDAQ\n')

    with open(where + '/array.pbs', 'a') as array_pbs:
      array_pbs.write(generate_array_pbs(g, "1kBq"))

    with open(where + '/Gate_parallel.sh', 'a') as gate_parallel_sh:
      gate_parallel_sh.write(generate_gate_parallel_sh(NR_OF_SPLITS))

def generate_sensitivity_simulations():

  workdir_Generated_Sensitivity = workdir_Generated + '/Sensitivity'
  if (os.path.isdir(workdir_Generated_Sensitivity)):
    os.system('rm -rf ' + workdir_Generated_Sensitivity)
  os.system('mkdir ' + workdir_Generated_Sensitivity)

  for g in geometries_sensitivity:
    os.system('mkdir ' + workdir_Generated_Sensitivity + '/' + g)
    os.system('mkdir ' + workdir_Generated_Sensitivity + '/' + g + '/output')

    where = workdir_Generated_Sensitivity + '/' + g + '/'
    copy_to_destination('materials/GateMaterials.db', where)
    copy_to_destination('materials/Materials.xml', where)
    copy_to_destination('materials/Surfaces.xml', where)
    copy_to_destination('geometries/' + g.split('_')[0] + '/Geometry_' + g + '.mac', where)
    copy_to_destination('Visualisation.mac', where)
    copy_to_destination('physics/PhysicsEMLivermorePolar.mac', where)

    # nr of annihilations in the axial field of view
    a = 1000000. # in Bq
    time_simulated = NR_OF_ANNIHILATIONS/a
    STRIP_LENGTH = geometry_to_strip_length(g)
    if STRIP_LENGTH<SOURCE_LENGTH:
      time_simulated = time_simulated*float(SOURCE_LENGTH)/float(STRIP_LENGTH)
    time_slice = time_simulated/NR_OF_SPLITS

    with open(where + '/NEMA_Sensitivity_Source_Gamma.mac', 'a') as source_mac:
      source_mac.write(generate_sensitivity_source(SOURCE_LENGTH, int(a)))

    with open(where + '/main.mac', 'a') as main_mac:
      if (ENABLE_VISUALISATION):
        main_mac.write('/vis/enable\n')
        main_mac.write('/control/execute ./Visualisation.mac\n')
      else:
        main_mac.write('/vis/disable\n')
      main_mac.write('/gate/geometry/setMaterialDatabase ./GateMaterials.db\n')
      main_mac.write('/control/execute ./Geometry_' + g + '.mac\n')
      main_mac.write('/control/execute ./PhysicsEMLivermorePolar.mac\n')
      main_mac.write('/gate/run/initialize\n')
      main_mac.write('/control/execute ./NEMA_Sensitivity_Source_Gamma.mac\n')
      main_mac.write('/gate/output/root/enable\n')
      main_mac.write('/gate/output/root/setFileName output/output\n')
      main_mac.write('/gate/output/root/setRootHitFlag 1\n')
      main_mac.write('/gate/random/setEngineName MersenneTwister\n')
      main_mac.write('/gate/random/setEngineSeed auto\n')
      main_mac.write('/gate/application/setTimeSlice ' + str('{0:.10f}'.format(time_slice)) + ' s\n')
      main_mac.write('/gate/application/setTimeStart 0. s\n')
      main_mac.write('/gate/application/setTimeStop ' + str('{0:.10f}'.format(time_simulated)) + ' s\n')
      main_mac.write('/gate/application/startDAQ\n')

    with open(where + '/array.pbs', 'a') as array_pbs:
      array_pbs.write(generate_array_pbs(g, "1kBq"))

    with open(where + '/Gate_parallel.sh', 'a') as gate_parallel_sh:
      gate_parallel_sh.write(generate_gate_parallel_sh(NR_OF_SPLITS))

if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='Analyze coincidences file and plot results of the analysis',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-m', '--mode',
                      dest='mode',
                      type=str,
                      default='all',
                      help='mode of the script: necr, sf, sensitvity')

  parser.add_argument('--nr-of-splits',
                      dest='nr_of_splits',
                      type=int,
                      default=NR_OF_SPLITS,
                      help='nr of splits')

  parser.add_argument('--nr-of-annihilations',
                      dest='nr_of_annihilations',
                      type=int,
                      default=NR_OF_ANNIHILATIONS,
                      help='nr of annihilations in the AFOV')

  parser.add_argument('--source-length',
                      dest='source_length',
                      type=float,
                      default=SOURCE_LENGTH,
                      help='length of the source/phantom [cm]')

  parser.add_argument('-w', '--workdir',
                      dest='workdir',
                      type=str,
                      default=workdir_Generated,
                      help='main workdir with generated macros')

  parser.add_argument('--enable-visualisation',
                      dest='enable_visualisation',
                      action='store_true',
                      help='enable visualisation in simulations')

  args = parser.parse_args()

  workdir_Generated = args.workdir
  if (not os.path.isdir(workdir_Generated)):
    os.system('mkdir ' + workdir_Generated)

  NR_OF_SPLITS = args.nr_of_splits
  NR_OF_ANNIHILATIONS = args.nr_of_annihilations
  SOURCE_LENGTH = args.source_length
  ENABLE_VISUALISATION = args.enable_visualisation

  if args.mode == 'necr':
    generate_necr_simulations()

  elif args.mode == 'sf':
    generate_sf_simulations()

  elif args.mode == 'sensitivity':
    generate_sensitivity_simulations()

  elif args.mode == 'all':
    generate_necr_simulations()
    generate_sf_simulations()
    generate_sensitivity_simulations()
