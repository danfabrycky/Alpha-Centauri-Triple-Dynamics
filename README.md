# Alpha-Centauri-Triple-Dynamics

These codes and results are in support of 
Jones & Fabrycky 2017.

Section 2 of the paper derives the orbit, and IDL code that gives these results can be found in mcorbitker2.pro. 
The resulting chain of parameters is found in mcorbitker2.txt.

Section 3 of the paper gives dynamical results. 
These are converted into dynamical simulation input files by mkin2ker2.pro; 
those files are the 100 in insker2.tar -- one extra one is to test the GR precession without Proxima. 
The dynamical code that is used for the investigation is called tripledyn.c; its header describes how it is compiled.
The input files were sent to this code on the Midway supercomputer at University of Chicago, resulting in the
output files within outsker2.tar.gz.  
Those output files were processed into the main results of the paper using plotkepfbplbker2gr.pro. 
