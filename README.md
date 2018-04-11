[![DOI](https://zenodo.org/badge/110665320.svg)](https://zenodo.org/badge/latestdoi/110665320) [![Build Status](https://travis-ci.org/gujans/THERMAID.svg?branch=master)](https://travis-ci.org/gujans/THERMAID)

# THERMAID
A numerical code to solve flow and heat transport in fractured porous media by the embedded discrete fracture method

THERMAID is a thermo-hydraulic code for fractured media that accounts for mechanical stability of the fractures, slip-induced transmissivity increase and thermally induced stresses.

**Important note:** *Our documentation is currently under construction. Please don't hesistate to contact the author's directly should any questions arise.*

What is Thermaid?
------------------

THERMAID is short for *Thermo-Hydraulic Energy Resource Modelling for Application and Development* and is a MATLAB program solving flow and heat transport in fractured porous media using the embedded discrete fracture method. It is targeted at researchers in the field of subsurface simulations who want to investigate energy resources in the subsurface and understand the embedded discrete fracture model.



For the impatient:
------------------

Download or clone the THERMAID repository. If you downloaded the repository as .zip or .tar.gz unpack it. If you have downloaded the repository, please rename it to ```THERMAID```

Let's say you've stored everything in a directory /path/to/THERMAID. 
Then you can run and test the program by opening matlab and running:

    >> cd /path/to/THERMAID
    >> open examples/README.md
    >> ex1

This will run first example of THERMAID. The other examples can be used similarly:

    >> ex2
    >> ex3
    >> ex4
    >> ex5

If you recieve an error between running different examples or simulations complaining about **matrix dimensions that do not agree** 
it is advised to use 

```
close all; clear all;
```

in order to remove the persistent matrices which are resolution and DFN dependent. Afterwards the simulation should be able to start!


Detailed instructions can be found in the quick start guide.

License:
--------

We wan't to provide you with the freedom to use and modify the code, but care about sharing with the community. Therefore the EDFM source code is released under the GNU General Public License (GPL) v3. The GNU GPL v3 is a copyleft license that requires anyone who distributes the source code or a derivative work to make the source available under the same terms. Please see the file ./LICENSE for details

Further information:
--------------------

Please see the THERMAID publication, the quick start guide or direct questions directly to the author!
