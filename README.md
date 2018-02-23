# j-pet-gate-tools
Extra tools for GATE simulations

detector splitter
-----------------

  Create geometry macro with split strpis on regions.

  Template pattern: J-PET cylindrical system 348 strips 500 mm L

  python splitter.py  -c <crystal_number> -o <outpu_tfile> -s <number_strips>

examples
--------

  Macros for GATE simulations of the J-PET scanner and
  python codes for analysis of the NEMA simulations.

examples/goja
-------------

  Code of the c++ app that reads .root output from GATE software and
  converts hits to set of coincidences saved in plain text files.

jplot
-----

  Draw line of response. Input lm file generate aplication analysis
  made by Pawel Kowalski.

  python jplot.py  <input_file>
