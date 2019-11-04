# J-PET layers controller

Transforms list-mode ASCII data from 3-layer scanner geometry (see [GOJA output](https://github.com/JPETTomography/j-pet-gate-tools/tree/master/goja#goja-output "GOJA output") description) into 1-layer.

The default parameters implies [big barrell](http://koza.if.uj.edu.pl/petwiki/index.php/Simulated_geometries#Laboratory_geometry_-_3_layers.2C_192_strips_.28big_barell.29 "JPET Wiki") geometry, with the options to perform TOF adjustment (since LORs, redefined to 1-layer, are of different lengths). An arbitrary layer ID which the data will be remapped to could be chosen, with the 1st one as the default.

## Prerequisites

The application is written in Python 2.7 (exact version Python 2.7.13), with the addition of the following packages:

* numpy (tested for 1.11.3)
* pandas (tested for 0.19.2)

## Usage

Both files (```LayersController.py``` and ```reduce_multilayer.py```) must be located in the same directory. The generic usage is as follows:
```
$ python reduce_multilayer.py <input_file> [-o <output_file>] [-l <layer ID>] [-s <no of strips>] [-with_tof]
```
Here, ```<layer ID>``` could be set as

* 0 - 'zero' layer, where 7-mm wide strips would be composed tightly (R=22.34 cm for 192 strips).
* 1/2/3 - layers of J-PET scanner, R=42.5/46.75/57.5 cm, respectively.
* 4 - ideal geometry: R=43.73 cm, 384 stips if tightly composed.

The program will map the hits onto the fixed number of strips (default is 192), assuming they are adjacent to each other regardless of width. One can set this number  manually as ```<no of strips>``` by using ```-s``` flag.

The flag ```-with_tof``` forces all hit times to be recalculated in order to match new LORs lengths. It is also better to turn off warnings during the execution:
```
$ python -W ignore reduce_multilayer.py SOME_INPUT_ASCII_DATA -l 0 -with_tof
```
(!) IMPORTANT: the application implies that hit times are in ps, and not ns, as described in [GOJA description](https://github.com/JPETTomography/j-pet-gate-tools/tree/master/goja "see GOJA output").
