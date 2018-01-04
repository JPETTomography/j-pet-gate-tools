# GOJA

GOJA is an acronym of GATE Output J-PET Analyzer. It may be used to
parse the GATE *.root output files and extract sets of coincidences
from them.

Dependencies
------------
ROOT v. 5.X
Boost v. 1.50.X

Tested for following configurations:
 ROOT v. 5.34/26, Boost v. 1.58.0.1ubuntu1 (Desktop Ubuntu 16.04 LTS)
 ROOT v. 5.34-x86_64-gcc48-python27, Boost v. 1.58.0 (cis.gov.pl cluster)

Building
--------
mkdir build
cd build
cmake ..
make
sudo make install

Usage
-----
Example:
goja --eth 0.2 --eth0 0.01 --tw 3 --root output.root > coincidences.txt
If the user won't redirect the output to a file, it will be written into the
standard output (console).

GOJA output
-----------
GOJA output consists of 16 columns. Each line of the output contains a single coincidences (LOR).
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

Help
----
Type goja --help (or just goja without arguments) for details.
