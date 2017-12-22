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

Help
----
Type goja --help (or just goja without arguments) for details.
