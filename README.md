# HaploDynamics
The python library **HaploDynamics**, or **HaploDX** for short, provides a collection of functions that can be used to simulate population-specific genomic data. This package is part of the Genetic Simulator Resources (GSR) catalog (click the link below).

<div style="width: 180px; margin: 0 auto;"><a href="https://surveillance.cancer.gov/genetic-simulation-resources/"><img src="https://surveillance.cancer.gov/gsr/static/img/gsr_tile.jpg" alt="Catalogued on GSR" width="180" height="60" /></a></div>

## Installation

```bash
$ git clone https://github.com/remytuyeras/HaploDynamics.git
```
Then, use your script where the folder ```HaploDynamics``` is located.
```bash
$ ls
HaploDynamics
$ touch myscript.py
```

## Quick start

The following script generates a VCF file containing simulated diploid genotypes for a population of 1000 individuals with LD-blocks of length 20kb, 5kb, 20kb, 35kb, 30kb and 15kb. 
```python
from HaploDynamics.HaploDX import genmatrix , create_vcfgz

simulated_data = genmatrix([20,5,20,35,30,15],strength=1,population=0.1,Npop=1000)
create_vcfgz("genomic-data.simulation.v1",*simulated_data)
```
The equation ```stregnth=1``` forces a high amount of linkage disequilibrium and the equation ```population=0.1``` increases the likelyhood of the simulated population to have rare mutations (e.g. to simulate a population profile close to African and South-Asian populations). 

More generally, the function ```genmatrix()``` takes the following types of parameters:
Parameters | Type | Values
| :--- | :--- | :---
```blocks```  | ```list[int]``` | List of positive integers, ideally between 1 and 40.
```strength```  | ```float``` | From -1 (little linkage) to 1 (high linkage)
```population```  | ```float``` | From 0 (for more rare mutations) to 1 (for less rare mutations)
```Npop```  | ```int```  | Positive integer specifying the number of individuals in the genomic matrix

The generation of each locus in a VCF file tend to be linear in the parameter ```Npop```. For example, one genetic variant can take from 0.3 to 0.6 seconds to be generated/simulated when we set ```Npop=100000``` (this may vary depending on your machine). The estimated time complexity for an average machine is shown below.

![GitHub Logo](/time_complexity.png)

The following script shows how to display the linkage disequilibirum correlations associated with the simulated data.
```python
import matplotlib.pyplot as plt
from HaploDynamics.HaploDX import genmatrix, create_vcfgz, display, LD_corr_matrix

simulated_data = genmatrix([20,5,20,35,30,15],strength=1,population=0.1,Npop=1000)
create_vcfgz("genomic-data.simulation.v1",*simulated_data)

rel, m, _ = LD_corr_matrix(simulated_data[0])
plt.imshow(display(rel,m))
plt.show()
```
The following plot is an example of the output that can be returned by the previous script when using 6 LD-blocks of 20kb each.

![alt text](http://www.normalesup.org/~tuyeras/node_diss/blg/blg_stat/img/LD_block_corr_strength_high.png)


## HaploDX Functions

You can find a complete presentation (or in fact a thorough tutorial) of the **HaploDX** library on my personal webpage (<a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">here</a>). Below is the list of all functions accessible from the library. It is recommended to first experiment with the functions presented in the [Data Generation](#data-generation) section.

### Population and Allele Frequency Spectrum Modeling

```python
def stochastic_line(a: float,b: float,sigma: float) -> Callable[float,float]
```
```python
population_mut = stochastic_line(0.08,0.17,0.01)
```
```python
def afs_distribution(index: int,alpha: float = 4/30) -> float
```
```python
def afs_intervals(pick: float,alpha: float = 4/30) -> list[float]
```
```python
def afs_sample(alpha: float = 4/30) -> float
```
```python
def genotype_schema(alpha: float = 4/30) -> tuple[float,list[float]]
```
```python
def genotype(hwp: list[float],minor: int) -> tuple[int,int]
```
```python
gref = lambda g: g[0] if g[0] in [2,0] else g[1]
```
```python
def population_mld(t: float) -> tuple[float,float,float]
```

### LD and Hardy–Weinberg Principle Modeling

```python
def decay(initial: float,halfwidth: float,shift: float) -> Callable[float,float]
```
```python
def ref_alt_function(y: float,x: float) -> float
```
```python
def alt_alt_function(y: float,z: float,x: float) -> float
```
```python
def amplifier(beta: float,p: float,q: float,s: float = 1) -> float
```
```python
def lb_freq(beta: float,gamma: float,previous_freq: float,distance: float,shift: float) -> float
```
```python
def ub_freq(beta: float,gamma: float,previous_freq: float,distance: float,shift: float) -> float
```
```python
def linkage_disequilibrium(alpha: float,beta: float,gamma: float,strength: float = -1) -> Callable[float,Callabel[float,tuple[float,float]]]
```
```python
def cond_genotype_schema(previous_maf: float,distance: float,alpha: float,beta: float,gamma: float,strength: float = -1) -> tuple[float,list[float],float]
```

### Data Generation

```python
def SNP_distribution(reference: float,length: float) -> list[float]
```
**Description**
- generates the list of positions for the VCF file
```python
def initiate_block(reference: float,alpha: float,Npop: int = 1000) -> tuple[float,list[list],list[list]]
```
**Description**
- initializes the first LD-block of the simulation
- ```reference``` refers to the first locus position at which the generation starts
```python
def continue_block(maf0: float,pre_matrix: list[list],matrix: list[list],positions: list[float],alpha: float,beta: float,gamma: float,strength: int = -1,Npop: int = 1000) -> tuple[float,list[list],list[list]]
```
**Description**
- generates an LD-block from a given position with a given minor allele frequency
- augments the genomic matrix with further genetic variants as specified by the arguments
```python
def genmatrix(blocks: list[int],strength: float,population: float,Npop: int)
```
**Description**
- implements a basic genomic matrix generator using ```SNP_distribution```, ```initiate_block``` and ```continue_block```
- _note:_ many possible improvements and variation of this function are possible (see tutorial <a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">here</a> for more detail)
  
```python
def gt_vcf(value: int)-> str
```
```python
def create_vcfgz(vcf_name: str,matrix: list[list],alpha: float,beta: float,gamma: float,system: str = "unix") -> None
```
**Description**
- generates a VCF file from a ```matrix``` generated by ```initiate_block```, ```continue_block``` or ```genmatrix```
  
### LD Analytics

```python
def LD_corr_matrix(matrix: list[list]) -> tuple[list[list],float,list[float]]
```
**Description**
- generates a (non-normalized) correlation matrix from a ```matri```x generated by ```initiate_block```, ```continue_block``` or ```genmatrix```
```python
def LD_r2_matrix(pre_matrix: list[list]) -> tuple[list[list],float,list[float]]
```
**Description**
- generates a (non-normalized) LD-r2 matrix from a ```pre_matrix``` generated by ```initiate_block```, ```continue_block```
- _note:_ this function cannot be used with ```genmatrix``` since the ```pre_matrix``` is not returned
```python
def display(rel: list[list],m: float) -> list[list[tuple[float,float,float]]]
```
```python
def minor_haplotype(sub_pre_matrix: list[list]) -> float
```


