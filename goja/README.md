# GOJA

GOJA is an acronym of GATE Output J-PET Analyzer. It may be used to
parse the GATE *.root output files and extract sets of coincidences
from them.

Available processing options
----------------------------

1) GOJA uses both scanner (default) or cylindricalPET GATE geometries. It assumes that the hits tree is called Hits, while the singles tree (only in cylindricalPET geometry) is called HESingles.

2) Coincidences may be formed with raw hits (based on tree Hits), from GOJA-merged singles (based on tree Hits) or GATE-marged singles (based on tree HEsingles).

3) There are currently 3 modes of merging hits into singles in GOJA (may be set with a switch in the code):
   a) kCentroidWinnerNaivelyWeighted: time and position of the single is set to mean values for merged hits,
   b) kCentroidWinnerEnergyWeighted: time and position of the single is set to mean values for merged hits weighted with the energies of these hits,
   c) kCentroidWinnerEnergyWeightedFirstTime: time of the single is set to time of the first merged hit, while the position of the single is set to mean value for merged hits weighted with the energies of these hits (default option),
   d) kEnergyWinner: single takes time and position from the hit with highest energy.
   In all modes, energy of the single is the sum of the energies of merged hits.

4) The definition of the coincidence may be adjusted with the lenght of time window (--tw) and energy threshold (--eth). Default: 3 ns, 200 keV.

Dependencies
------------

ROOT v. 5.X
Boost v. 1.50.X

Tested for following configurations:
 ROOT v. 5.34/26, Boost v. 1.58.0.1ubuntu1 (Desktop Ubuntu 16.04 LTS)
 ROOT v. 5.34-x86_64-gcc48-python27, Boost v. 1.58.0 (cis.gov.pl cluster)

Building and installation
-------------------------

mkdir build (in source directory)
cd build
cmake ..
make
sudo make install

Files 'goja' and 'goja_manager.py' will be installed by default in /usr/local/bin directory.

GOJA output
-----------

GOJA output consists of 19 columns. Each line of the output contains a single coincidence (LOR).
Columns 1-4 contain information about the 1st hit:
1) posX1 [cm]
2) posY1 [cm]
3) posZ1 [cm]
4) time1 [ns]
Columns 5-8 contain information about the 2nd hit:
5) posX2 [cm]
6) posY2 [cm]
7) posZ2 [cm]
8) time2 [ns]
Next columns contain following information:
9) volumeID of the 1st hit
10) volumeID of the 2nd hit
11) edep of the 1st hit [keV]
12) edep of the 2nd hit [keV]
13) type of coincidence (1 - true, 2 - phantom-scattered, 3 - detector-scattered, 4 - accidental)
14) sourcePosX1 [cm]
15) sourcePosY1 [cm]
16) sourcePosZ1 [cm]
17) sourcePosX2 [cm]
18) sourcePosY2 [cm]
19) sourcePosZ2 [cm]

Usage
-----

Examples:

Example 1:
goja --eth 0.2 --eth0 0.01 --tw 3 --root output.root > coincidences.txt
If the user won't redirect the output to a file, it will be written into the
standard output (console).

Example 2:
goja_manager.py --mode analyze --gate-output ../output/ --goja-output ./

Example 3: Analyze the NEMA sensitivity GATE output

If one wants to use goja and goja_manager.py on cis.gov.pl cluster, he needs to load gate/7.2 module:

  module load gate/7.2

Let's assume that calcualtions were performed for some geomtries, which names start with letter 'D'.
In each geometry directory, there is an 'output' directory with splitted root files. Then analysis
may be performed in following way:

  for i in D*; do mkdir $i/goja_N0_1000; done
  for i in D*; do cd $i/goja_N0_1000 && goja_manager.py --N0 1000 --run on-cluster && cd ../..; done

After calcualtions one may verify if there are no missing files:

  for i in D*; do cd $i/goja_N0_1000 && goja_manager.py --mode verify && cd ../..; done

If there are no missing files one can clean goja output irectory form cluster logs and concatenate results:

  for i in D*; do rm $i/goja_N0_1000/*.e* $i/goja_N0_1000/*.o*; done
  for i in D*; do cd $i/goja_N0_1000 && goja_manager.py -m concatenate -sn ${i}_sensitivity -c && cd ../..; done

Example 4: NECR

Let's assume that in directory D85_1lay_L020_7mm there are subudirectories with names 0001, 01000, 0200, ... 2000.
In this case, name of subdirectory relates to the subsequent activity.
The analysis may be performed as in Example 3, using follwoing commands:

  cd D85_1lay_L020_7mm
  for i in *0*; do mkdir $i/goja_N0_1000; done
  for i in *0*; do cd $i/goja_N0_1000 && goja_manager.py --N0 1000 --run on-cluster && cd ../..; done
  for i in *0*; do cd $i/goja_N0_1000 && goja_manager.py -m verify && cd ../..; done
  for i in *0*; do cd $i/goja_N0_1000 && goja_manager.py -m concatenate -sn D85_1lay_L020_7mm_${i}_NECR -c && cd ../..; done

Hints
-----

1) How to check number of coincidences of a given type?

  analyze_coincidences.py -cf ./output_coincidences -m stats

2) How to kill series of matrix tasks on cluster?

  for i in {2680660..2680671}; do qdel ${i}[]; done

Help
----

Type goja --help (or just goja without arguments) for details.
Type goja_manager.py --help for details fo the GOJA manager.
