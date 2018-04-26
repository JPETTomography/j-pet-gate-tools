# Smear events

A tool which allows smearing both hits for each event/coincidence of back-to-back gamma quanta emission from the source 
inside J-PET scanner, according to the type of the readout (photomultiplier). It proceeds 
[list mode ASCII data file](https://github.com/JPETTomography/j-pet-gate-tools/tree/master/examples/goja 
"see the description here") sequentially, line by line.

Compilation and usage
---------------------

As the source files utilize ```unordered_map```, they ought to be compiled with ```-std=c++11``` option:

```
g++ -Wall -std=c++11 -o SmearEvents SmearEvents.cpp
```
For the executable above, the usage will be as follows:

```
./SmearEvents -i <input_data_file> -o <output_file_name> -r PMT [-l 2]
```

The key for the readout ```-r``` is restricted to three options: ```PMT```, ```SI``` (Silicon PM or SiPM) 
and ```WLS```. ```-l``` defines the scanner length (1/2/3 for 20/50/100 [cm], with 2 as the default).
