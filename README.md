# Presentation of HaploDynamics 
**HaploDynamics** (**HaploDX**) is a Python 3+ library that provides a collection of functions for simulating population-specific genomic data. It is part of the Genetic Simulator Resources (GSR) catalog. You can access the GSR catalog by clicking on the image below.

<div style="width: 180px; margin: auto;"><a href="https://surveillance.cancer.gov/genetic-simulation-resources/"><img src="https://surveillance.cancer.gov/gsr/static/img/gsr_tile.jpg" alt="Catalogued on GSR" width="180" height="60" /></a></div>

## New features added

* Installation via ```pip```;
* A module ```Framework``` for software development and experimentations. Planned features:
  - classes to model varied populations;
  - multiprocessing functionalities for simulation on the cloud.

## Installation

### Installation via ```pip```
Install the HaploDynamics package by using the following command.
```bash
$ pip install HaploDynamics
```
After this, you can import the modules of the library to your script as follows.
```python
import HaploDynamics.HaploDX as hdx
import HaploDynamics.Framework as hdx_frm
```

### Manual installation
HaploDynamics uses the [SciPy](https://docs.scipy.org/doc/scipy/reference/stats.html) library for certain calculations. To install SciPy, run the following command, or see SciPy's [installation instructions](https://scipy.org/install/) for more options.
```bash
$ python -m pip install scipy
```
You can install the HaploDynamics GitHub package by using the following command in a terminal.
```bash
$ git clone https://github.com/remytuyeras/HaploDynamics.git
```
Then, use the ```pwd``` command to get the absolute path leading to the downloaded package.
```bash
$ ls
HaploDynamics
$ cd HaploDynamics/
$ pwd
absolute/path/to/HaploDynamics
```
To import the modules of the library to your script, you can use the following syntax where the path ```absolute/path/to/HaploDynamics``` should be replaced with the path obtained earlier.
```python
import sys
sys.path.insert (0,"absolute/path/to/HaploDynamics")
import HaploDynamics.HaploDX as hdx
import HaploDynamics.Framework as hdx_frm
```
## Quickstart

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

The generation of each locus in a VCF file tends to be linear in the parameter ```Npop```. On average, a genetic variant can take from 0.3 to 0.8 seconds to be generated when ```Npop=100000``` (this may vary depending on your machine). The estimated time complexity for an average machine is shown below.

![](img/time_complexity.png) 

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

![](img/simulation_LD_0.png) 

The following script shows how you can control linkage disequilibrium by using LD-blocks of varying legnths. You can display the graph relating distances between pairs of SNPs to average correlation scores by using the last output of the function ```LD_corr_matrix()```.

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
Typical outputs for the previous script should look as follows.

Correlations            |  SNP-distance to average correlations
:-------------------------:|:-------------------------:
![](img/simulation_LD_1.png)  |  ![](img/simulation_dist_1.png)

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
![](img/simulation_LD_2.png)  |  ![](img/simulation_dist_2.png)

## To cite this work

Tuyeras, R. (2023). _HaploDynamics: A python library to develop genomic data simulators_ (Version 0.2-beta.4) [Computer software]. [![DOI](https://zenodo.org/badge/609227235.svg)](https://zenodo.org/badge/latestdoi/609227235)

<br/>

# Documentation

## Functions available in the ```HaploDX``` module

A comprehensive tutorial of the **HaploDX** module can be found on my personal webpage (<a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">here</a>). This tutorial provides a holistic overview of the module, including its functions, features, and applications.

The present documentation, on the other hand, focuses on each function of the HaploDX module individually. It is recommended to first experiment with the functions presented in the [Data Generation](#data-generation) section, which will give you a good foundation for understanding the rest of the module.

### Population and Allele Frequency Spectrum Modeling
The functions presented here are based on research results borrowed from the literature. For more information on the functions' designs, please refer to the <a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">tutorial</a>.

\
**stochastic_line**
```python
def stochastic_line(a: float,b: float,sigma: float) -> Callable[[float],float]
```
> <u>Description</u>
> - **inputs:**
>   - ```a```: extreme end of a line;
>   - ```b```: other extreme end of the line;
>   - ```sigma```: standard deviation.
> - **output:**
>   - a function taking an input ```t: float``` an returning a floating-point number sampled from a Gaussian distribution whose mean is ```a*t+b*(1-t)``` and whose standard deviation is ```sigma```.
>
\
**population_mut**
```python
population_mut = stochastic_line(0.08,0.17,0.01)
```
> <u>Description</u>
> - **inputs:**
>   -  ```t```: floating-point number in the interval $[0,1]$.
> - **output:**
>   - floating-point number sampled from a Gaussian distribution whose mean is ```a*t+b*(1-t)``` and whose standard deviation is ```sigma```.
>
\
**afs_distribution**
```python
def afs_distribution(index: int,alpha: float = 4/30) -> float
```
> <u>Description</u>
> - **inputs:**
>   -  ```index```: integer from 0 to 15. See the section on allele frequency modeling in the <a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">tutorial</a>.
>   - ```alpha```: bottleneck parameter defining a population-specific mutation rate profile (e.g. YRI is ```0.16```, CEU is ```0.9``` and CHB+JPT is ```0.8```). See the section on allele frequency modeling in the <a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">tutorial</a>.
> - **output:**
>   - probability that a given allele has an allele frequency equal to ```(0.5+index)/15```.
>
\
**afs_intervals**
```python
def afs_intervals(pick: float,alpha: float = 4/30) -> list[float]
```
> <u>Description</u>
> - **inputs:**
>   -  ```pick```: allele frequency;
>   - ```alpha```: bottleneck parameter defining a population-specific mutation rate profile (e.g. YRI is ```0.16```, CEU is ```0.9``` and CHB+JPT is ```0.8```). See the section on allele frequency modeling in the <a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">tutorial</a>.
> - **output:**
>   - group bracket```[i/15,(i+1)/15]``` for the given allele frequency distribution ```pick```. See the section on allele frequency modeling in the <a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">tutorial</a>.
>
\
**afs_sample**
```python
def afs_sample(alpha: float = 4/30) -> float
```
> <u>Description</u>
> - **inputs:**
>   - ```alpha```: bottleneck parameter defining a population-specific mutation rate profile (e.g. YRI is ```0.16```, CEU is ```0.9``` and CHB+JPT is ```0.8```). See the section on allele frequency modeling in the <a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">tutorial</a>.
> - **output:**
>   - floating-point number encoding an (stochastic) allele frequency for a population characterized by the input value ```alpha```.
>
\
**genotype_schema**
```python
def genotype_schema(alpha: float = 4/30) -> tuple[float,list[float]]
```
> <u>Description</u>
> - **inputs:**
>   - ```alpha```: bottleneck parameter defining a population-specific mutation rate profile (e.g. YRI is ```0.16```, CEU is ```0.9``` and CHB+JPT is ```0.8```). See the section on allele frequency modeling in the <a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">tutorial</a>.
> - **output:**
>   - ```maf```: floating-point number encoding a (stochastic) minor allele frequency $p$;
>   - ```hwp```: list of four floating-point numbers defining a tiling of the interval $[0,1]$. The tiling is based on Hardy-Weinberg distribution given by the four values $0$, $p$, $p+p(1-p)$, $p+2p(1-p)$ and $1$.
>
\
**genotype**
```python
def genotype(hwp: list[float],minor: int) -> tuple[int,int]
```
> <u>Description</u>:
> - This function simulates information about a hypothetical genotype.
> - **inputs:**
>   - ```hwp```: a Hardy-Weinberg distribution defined by four values $0$, $p$, $p+p(1-p)$, $p+2p(1-p)$ and $1$, where the variable $p$ represents a minor allele frequency. The input ```hwp``` is typically returned by the function ```genotype_schema()```;
>   - ```minor```: an integer value from the set $\{0,2\}$ to indicate whether a given minor allele corresponds to the allele of an imaginary reference genome sequence. The value $0$ should be used to simulate a situation where the minor allele is the _reference allele_.
> - **output:**
>   - ```minor_counts```: integer value counting the number of minor alleles present in a simulated genotype (either $0$, $1$ or $2$);
>   - ```genotype_code```: integer value representing the simulated genotype. This value is either equal to ```minor``` if the simulated genotype is homozygous, or equal to $1$ or $-1$ if the simulated genotype is heterozygous. See the function ```gref()``` for an interpretation of these values.
>
\
**gref**
```python
gref = lambda g: g[0] if g[0] in [2,0] else g[1]
```
> <u>Description</u>
> - **inputs:**
>   - ```g```: an output of the function ```genotype()```;
> - **output:**
>   - integer value from the set $\{-1,0,1,2\}$. Specificaly, the value $0$ represents the major homozygous genotype, the value $2$ represents the minor homozygous genotype, the value $1$ represents the minor-major genotype and the value $-1$ represents the major-minor genotype.
>
\
**population_mld**
```python
def population_mld(t: float) -> tuple[float,float,float]
```
> <u>Description</u>
> - **inputs:**
>   - ```t```: floating-point number in the interval $[0,1]$ to model human population characterics through a linear representation. See the section on Hardy-Weinberg principle and linkage disequilibrium modeling in the <a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">tutorial</a>.
> - **output:**
>   - a triple of floating-point numbers that represent a point near the representative line for human population characteristics, specifically the bottleneck and LD-decay parameters.
>

### LD and Hardy–Weinberg Principle Modeling
The functions presented here use calculations of probabilities that generalize the Hardy-Weinberg principle to haplotypes. For more information, please refer to the <a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">tutorial</a>.

\
**decay**
```python
def decay(initial: float,halfwidth: float,shift: float) -> Callable[[float],float]
```
> <u>Description</u>
> - **inputs:**
>   - ```initial```: .
>   - ```halfwidth```: .
>   - ```shift```: .
> - **output:**
>   - funtion
>
\
**ref_alt_function**
```python
def ref_alt_function(y: float,x: float) -> float
```
> <u>Description</u>
> - **inputs:**
>   - ```y```: .
>   - ```x```: .
> - **output:**
>   - floating-point number
>
\
**alt_alt_function**
```python
def alt_alt_function(y: float,z: float,x: float) -> float
```
> <u>Description</u>
> - **inputs:**
>   - ```y```: .
>   - ```z```: .
>   - ```x```: .
> - **output:**
>   - floating-point number
>
\
**amplifier**
```python
def amplifier(beta: float,p: float,q: float,s: float = 1) -> float
```
> <u>Description</u>
> - **inputs:**
>   - ```beta```: .
>   - ```p```: .
>   - ```q```: .
>   - ```s```: .
> - **output:**
>   - floating-point number
>
\
**lb_freq**
```python
def lb_freq(beta: float,gamma: float,previous_freq: float,distance: float,shift: float) -> float
```
> <u>Description</u>
> - **inputs:**
>   - ```beta```: .
>   - ```gamma```: .
>   - ```previous_freq```: .
>   - ```distance```: .
>   - ```shift```: .
> - **output:**
>   - floating-point number
>
\
**ub_freq**
```python
def ub_freq(beta: float,gamma: float,previous_freq: float,distance: float,shift: float) -> float
```
> <u>Description</u>
> - **inputs:**
>   - ```beta```: .
>   - ```gamma```: .
>   - ```previous_freq```: .
>   - ```distance```: .
>   - ```shift```: .
> - **output:**
>   - floating-point number
>
\
**linkage_disequilibrium**
```python
def linkage_disequilibrium(alpha: float,beta: float,gamma: float,strength: float = -1) -> Callable[[float],Callable[[float],tuple[float,float]]]
```
> <u>Description</u>
> - **inputs:**
>   - ```alpha```: .
>   - ```beta```: .
>   - ```gamma```: .
>   - ```strength```: .
> - **output:**
>   - funtion
>
\
**cond_genotype_schema**
```python
def cond_genotype_schema(previous_maf: float,distance: float,alpha: float,beta: float,gamma: float,strength: float = -1) -> tuple[float,list[float],float]
```
> <u>Description</u>
> - **inputs:**
>   - ```previous_maf```: .
>   - ```distance```: .
>   - ```alpha```: .
>   - ```beta```: .
>   - ```gamma```: .
>   - ```strength```: .
> - **output:**
>   - ```maf```: .
>   - ```hwp```: .
>   - ```ld```: .
>

### Data Generation
The functions presented here can be used to implement simulators of variant call data. The function ```genmatrix()``` is a typical example of this. For more information about the implementation of ```genmatrix()```, please refer to the <a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">tutorial</a>.

\
**SNP_distribution**
```python
def SNP_distribution(reference: float,length: float) -> list[float]
```
> <u>Description</u>
> - generates the list of positions for the VCF file
> - **inputs:**
>   - ```reference```: .
>   - ```length```: .
> - **output:**
>   - list
>
\
**initiate_block**
```python
def initiate_block(reference: float,alpha: float,Npop: int = 1000) -> tuple[float,list[list],list[list]]
```
> <u>Description</u>
> - initializes the first LD-block of the simulation
> - the parameter ```reference``` refers to the first locus position at which the generation starts
> - **inputs:**
>   - ```reference```: .
>   - ```alpha```: .
>   - ```Npop```: .
> - **output:**
>   - tuple
>
\
**continue_block**
```python
def continue_block(maf0: float,pre_matrix: list[list],matrix: list[list],positions: list[float],alpha: float,beta: float,gamma: float,strength: int = -1,Npop: int = 1000) -> tuple[float,list[list],list[list]]
```
> <u>Description</u>
> - generates an LD-block from a given position with a given minor allele frequency
> - augments the genomic matrix with further genetic variants as specified by the arguments
> - **inputs:**
>   - ```reference```: .
>   - ```alpha```: .
>   - ```Npop```: .
> - **output:**
>   - tuple
>
\
**genmatrix**
```python
def genmatrix(blocks: list[int],strength: float,population: float,Npop: int)
```
> <u>Description</u>
> - implements a basic genomic matrix generator using ```SNP_distribution()```, ```initiate_block()``` and ```continue_block()```
> - _note:_ many possible improvements or variations of this function are possible (see tutorial <a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">here</a> for more detail)
> - **inputs:**
>   - ```reference```: .
>   - ```alpha```: .
>   - ```Npop```: .
> - **output:**
>   - tuple
>
\
**gt_vcf**
```python
def gt_vcf(value: int)-> str
```
> <u>Description</u>
> - **inputs:**
>   - ```value```: .
> - **output:**
>   - string
>
\
**create_vcfgz**
```python
def create_vcfgz(vcf_name: str,matrix: list[list],alpha: float,beta: float,gamma: float,system: str = "unix") -> None
```
> <u>Description</u>
> - generates a VCF file from a ```matrix``` generated by ```initiate_block()```, ```continue_block()``` or ```genmatrix()```
> - **inputs:**
>   - ```reference```: .
>   - ```alpha```: .
>   - ```Npop```: .
> - **output:**
>   - tuple
>

### LD Analytics
The functions presented in this section can be used to visualize linkage disequilibrium (LD)-related information that characterizes the data generated by the functions presented in the [Data Generation](#data-generation) section.

\
**LD_corr_matrix**
```python
def LD_corr_matrix(matrix: list[list]) -> tuple[list[list],float,list[float]]
```
> <u>Description</u>
> - generates a (non-normalized) correlation matrix from a ```matrix``` generated by ```initiate_block()```, ```continue_block()``` or ```genmatrix()```
> - **inputs:**
>   - ```matrix```: .
>   - ```alpha```: .
> - **output:**
>   - tuple
>
\
**LD_r2_matrix**
```python
def LD_r2_matrix(pre_matrix: list[list]) -> tuple[list[list],float,list[float]]
```
> <u>Description</u>
>  - generates a (non-normalized) LD-r2 matrix from a ```pre_matrix``` generated by ```initiate_block()```, ```continue_block()```
>  - _note:_ this function cannot be used with ```genmatrix()``` since the ```pre_matrix``` is not returned
> - **inputs:**
>   - ```matrix```: .
>   - ```alpha```: .
> - **output:**
>   - tuple
>
\
**display**
```python
def display(rel: list[list],m: float) -> list[list[tuple[float,float,float]]]
```
> <u>Description</u>
> - **inputs:**
>   - ```rel```: .
>   - ```m```: .
> - **output:**
>   - list
>
\
**minor_haplotype**
```python
def minor_haplotype(sub_pre_matrix: list[list]) -> float
```
> <u>Description</u>
> - **inputs:**
>   - ```sub_pre_matrix```: .
> - **output:**
>   - float
>


## Classes available in the ```Framework``` module

To be added.
