# j-pet-gate-tools
Extra tools for GATE simulations

detector splitter
-----------------
Create geometry macro with split strpis on regions.

Template pattern: Geometry_384_strips_L50.mac made by Pawel Kowalski.

python splitter.py -l <strip_lenght> -n <part_number> -o <outpu_tfile> -s <number_strips>

jplot
-----
Draw line of response. Input file generate aplication analysis made by Pawel Kowalski.

python jplot.py  <input_file>