# HaploDynamics
The python library **HaploDynamics**, or **HaploDX** for short, provides a collection of functions that can be used to simulate population-specific genomic data. This package is part of the Genetic Simulator Resources (GSR) catalog (click the link below for more information).

<div style="width: 180px; margin: 0 auto;"><a href="https://surveillance.cancer.gov/genetic-simulation-resources/"><img src="https://surveillance.cancer.gov/gsr/static/img/gsr_tile.jpg" alt="Catalogued on GSR" width="180" height="60" /></a></div>

## New features recently added

* Installation via ```pip```;
* Module ```Framework``` -- this module will serve as a class-based library for development and experimentation purposes.

## Installation

### Manual installation
Download the github package by using the following command in a terminal.
```bash
$ git clone https://github.com/remytuyeras/HaploDynamics.git
```
Then, from your current directory, record the absolute path leading to the main module of the package, as shown below.
```bash
$ ls
HaploDynamics
$ cd HaploDynamics/
$ pwd
absolute/path/to/HaploDynamics
```
To import the modules of the library to your script, use the following syntax where you must replace ```absolute/path/to/HaploDynamics``` with the path that you obtained above.
```python
import sys
sys.path.insert (0,"absolute/path/to/HaploDynamics")
import HaploDynamics.HaploDX as hdx
import HaploDynamics.Framework as hdx_frm
```
### Installation via Pip
Install the package by using ```pip install``` as follows.
```bash
$ pip install HaploDynamics
```
Then, you can import the modules of the library to your script as follows.
```python
import HaploDynamics.HaploDX as hdx
import HaploDynamics.Framework as hdx_frm
```
## Quick start

The following script generates a VCF file containing simulated diploid genotypes for a population of 1000 individuals with LD-blocks of length 20kb, 5kb, 20kb, 35kb, 30kb and 15kb. 
```python
import HaploDynamics.HaploDX as hdx

simulated_data = hdx.genmatrix([20,5,20,35,30,15],strength=1,population=0.1,Npop=1000)
hdx.create_vcfgz("genomic-data.simulation.v1",*simulated_data)
```
The equation ```strength=1``` forces a high amount of linkage disequilibrium and the equation ```population=0.1``` increases the likelyhood of the simulated population to have rare mutations (e.g. to simulate a population profile close to African and South-Asian populations). 

More generally, the function ```genmatrix()``` takes the following types of parameters:
Parameters | Type | Values
| :--- | :--- | :---
```blocks```  | ```list[int]``` | List of positive integers, ideally between 1 and 40.
```strength```  | ```float``` | From -1 (little linkage) to 1 (high linkage)
```population```  | ```float``` | From 0 (for more rare mutations) to 1 (for less rare mutations)
```Npop```  | ```int```  | Positive integer specifying the number of individuals in the genomic matrix

The generation of each locus in a VCF file tends to be linear in the parameter ```Npop```. On average, a genetic variant can take from 0.3 to 0.8 seconds to be generated/simulated when ```Npop=100000``` (this may vary depending on your machine). The estimated time complexity for an average machine is shown below.

![](/img/time_complexity.png) 

## Use cases
The following script shows how to display linkage disequilibirum correlations for the simulated data.
```python
import matplotlib.pyplot as plt
import HaploDynamics.HaploDX as hdx

simulated_data = hdx.genmatrix([20,20,20,20,20,20],strength=1,population=0.1,Npop=1000)
hdx.create_vcfgz("genomic-data.simulation.v1",*simulated_data)

rel, m, _ = hdx.LD_corr_matrix(simulated_data[0])
plt.imshow(hdx.display(rel,m))
plt.show()
```
A typical output for the previous script should look as follows.
![](/img/simulation_LD_0.png) 

The following script shows that you can control linkage disequilibrium quite easily by using sequences of small LD-blocks. You can display the graph relating _distance between pairs of SNPS_ to _average correlation scores_ by using the last output of the function ```LD_corr_matrix()```.

```python
import matplotlib.pyplot as plt
import HaploDynamics.HaploDX as hdx

ld_blocks = [5,5,5,10,20,5,5,5,5,5,5,1,1,1,2,2,10,20,40]
strength=1
population=0.1
Npop = 1000
simulated_data = hdx.genmatrix(ld_blocks,strength,population,Npop)
hdx.create_vcfgz("genomic-data.simulation.v1",*simulated_data)
#Correlations
rel, m, dist = hdx.LD_corr_matrix(simulated_data[0])
plt.imshow(hdx.display(rel,m))
plt.show()
#from SNP-distance to average correlaions
plt.plot([i for i in range(len(dist)-1)],dist[1:])
plt.ylim([0, 1])
plt.show()
```
Typical outputs for the previous script should look as follows, where the right graph shows the graph linking **distance between SNPS** and **average correlations**.

Correlations            |  SNP-distance to average correlations
:-------------------------:|:-------------------------:
![](/img/simulation_LD_1.png)  |  ![](/img/simulation_dist_1.png)

Finally, the following script shows how you can generate large regions of linkage.

```python
import matplotlib.pyplot as plt
import HaploDynamics.HaploDX as hdx

ld_blocks = [1] * 250
strength=1
population=0.1
Npop = 1000
simulated_data = hdx.genmatrix(ld_blocks,strength,population,Npop)
hdx.create_vcfgz("genomic-data.simulation.v1",*simulated_data)
#Correlations
rel, m, dist = hdx.LD_corr_matrix(simulated_data[0])
plt.imshow(hdx.display(rel,m))
plt.show()
#from SNP-distance to average correlaions
plt.plot([i for i in range(len(dist)-1)],dist[1:])
plt.ylim([0, 1])
plt.show()
```
Typical outputs for the previous script should look as follows.

Correlations            |  SNP-distance to average correlations
:-------------------------:|:-------------------------:
![](/img/simulation_LD_2.png)  |  ![](/img/simulation_dist_2.png)

## Functions from the HaploDX module

You can find a complete presentation (or in fact a thorough tutorial) of the **HaploDX** module on my personal webpage (<a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">here</a>). Below is the list of all functions accessible from the library. It is recommended to first experiment with the functions presented in the [Data Generation](#data-generation) section.

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
- the parameter ```reference``` refers to the first locus position at which the generation starts
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
- implements a basic genomic matrix generator using ```SNP_distribution()```, ```initiate_block()``` and ```continue_block()```
- _note:_ many possible improvements or variations of this function are possible (see tutorial <a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">here</a> for more detail)
  
```python
def gt_vcf(value: int)-> str
```
```python
def create_vcfgz(vcf_name: str,matrix: list[list],alpha: float,beta: float,gamma: float,system: str = "unix") -> None
```
**Description**
- generates a VCF file from a ```matrix``` generated by ```initiate_block()```, ```continue_block()``` or ```genmatrix()```
  
### LD Analytics

```python
def LD_corr_matrix(matrix: list[list]) -> tuple[list[list],float,list[float]]
```
**Description**
- generates a (non-normalized) correlation matrix from a ```matrix``` generated by ```initiate_block()```, ```continue_block()``` or ```genmatrix()```
```python
def LD_r2_matrix(pre_matrix: list[list]) -> tuple[list[list],float,list[float]]
```
**Description**
- generates a (non-normalized) LD-r2 matrix from a ```pre_matrix``` generated by ```initiate_block()```, ```continue_block()```
- _note:_ this function cannot be used with ```genmatrix()``` since the ```pre_matrix``` is not returned
```python
def display(rel: list[list],m: float) -> list[list[tuple[float,float,float]]]
```
```python
def minor_haplotype(sub_pre_matrix: list[list]) -> float
```
