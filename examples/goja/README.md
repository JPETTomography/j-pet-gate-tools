# GOJA

GOJA is an acronym of GATE Output J-PET Analyzer. It may be used to
parse the GATE *.root output files and extract sets of coincidences
from them.

Dependencies
------------
ROOT
Boost

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
./goja --eth 0.2 --eth0 0.01 --tw 3 --root output.root > coincidences.txt
Type ./goja --help for details