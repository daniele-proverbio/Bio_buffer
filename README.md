# Bio_buffer
Code for the project "Mechanisms to buffer variability in cell regulation motifs next to criticality".  
Author: Daniele Proverbio, LCSB, 2022

## Analytical results
Analytical results are in the Mathematica notebooks, ordered by topic.  
Calculation and characterizatoin of the bifurcation diagram and its focal width is performed in ``equilibria_coop.m''.

## Simulations
Computer simulations of the biological system are performed by sequential run of Matlab files (Simulate -> Analysis -> Plotting). 
Note that, due to heavy files generated by multiple simulatoins, only the generatin code "simulate.m" is provided.  

Fiigures are mostly generated by "plotting.m". However, a few more are generated by intermediate files to decrease computaiional cost associated to savings. Final manuscrpt figures are collated using Inkscape.


## Notes
Due to limited computing resources, neither Matlab file is optimised to run sequentially. Hence, changing parameter should be performed manuallly at each highlighted line of the code, for the desired runs. Results are then saved and analysed later.
