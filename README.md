# HaploDynamics
HaploDynamics, or HaploDX for short, is a python library providing functions that can be used to population-specific simulate genomic data.

You can find a tutorial (and presentation) of this library on my personal webpage (click <a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">here</a>).

The following lines show how to use the library to generate genomic data with LD-blocks of length 20kb, 5kb, 20kb, 35kb, 30kb and 15kb. 
```python
matrix = genmatrix([20,5,20,35,30,15],strength=1,population=0.1)
create_vcfgz("genomic-data.simulation.v1",matrix)
```
The equation ```stregnth=1``` forces a high amount of linkage disequilibrium and the equation ```population=0.1``` increases the likelyhood of the simulated population to have rare mutations (e.g. to simulate a population profile close to African and South-Asian populations). 
