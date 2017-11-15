# Examples in Thermaid

The following example cases are currenlty included in Thermaid

 - ex1 : Crossed shaped fracture geometry. Steady-state pressure validation against COMSOL.
 - ex2 : Crossed shaped fracture geometry. Steady-state heat transport test validation against COMSOL.
 - ex3 : Complex DFN with 13 fractures. Dynamic heat transport example.
 - ex4 : Loads a fracSim DFN (see also DFN/frac_load_fracsun.m for details) and runs transient fluid pressure diffusion.
 - ex5 : Loads a fracman DFN (see also DFN/frac_load_fracman.m for details) and runs transient fluid pressure diffusion.

Note that between calling different examples (or simulations in general), it might be requiered to use 

```
close all; clear all;
```

in order to remove the persistent matrices which are resolution and DFN dependent.
