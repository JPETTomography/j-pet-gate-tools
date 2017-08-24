# j-pet-gate-tools
Extra tools for GATE simulations

detector splitter
-----------------
Create geometry macro with split strpis on regions.

Template pattern: J-PET cylindrical system 348 strips 500 mm L

python splitter.py  -c <crystal_number> -o <outpu_tfile> -s <number_strips>

jplot
-----
Draw line of response. Input lm file  generate aplication analysis made by Pawel Kowalski.

python jplot.py  <input_file>
