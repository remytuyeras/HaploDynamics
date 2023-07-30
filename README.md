# HaploDynamics
HaploDynamics, or HaploDX for short, is a python library providing functions that can be used to simulate population-specific genomic data.

### Quick start

The following lines of code illustrate how one can use the library to generate genomic data for 1000 individuals with LD-blocks of length 20kb, 5kb, 20kb, 35kb, 30kb and 15kb. 
```python
from HaploDynamics.HaploDX import genmatrix , create_vcfgz

simulated_data = genmatrix([20,5,20,35,30,15],strength=1,population=0.1,Npop=1000)
create_vcfgz("genomic-data.simulation.v1",*simulated_data)
```
The equation ```stregnth=1``` forces a high amount of linkage disequilibrium and the equation ```population=0.1``` increases the likelyhood of the simulated population to have rare mutations (e.g. to simulate a population profile close to African and South-Asian populations). 

The following script shows how to display the linkage disequilibirum correlations associated with the simulated data.
```python
import matplotlib.pyplot as plt
from HaploDynamics.HaploDX import genmatrix, create_vcfgz, display, LD_corr_matrix

simulated_data = genmatrix([20,5,20,35,30,15],strength=1,population=0.1,Npop=1000)
create_vcfgz("genomic-data.simulation.v1",*simulated_data)

rel, m, dist = LD_corr_matrix(simulated_data[0])
plt.imshow(display(rel,m))
plt.show()
```
![alt text](http://www.normalesup.org/~tuyeras/node_diss/blg/blg_stat/img/LD_block_corr_strength_high.png)


### HaploDX Functions

You can find a tutorial (and presentation) of this library on my personal webpage (click <a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">here</a>).

#### Population and Allele Frequency Spectrum Modeling

```python
def stochastic_line(a: float,b: float,sigma: float) -> Callable[float,float]
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

#### LD and Hardy–Weinberg Principle Modeling

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

#### Data Generation (in VCF files)

```python
def SNP_distribution(reference: float,length: float) -> list[float]
```
```python
def initiate_block(reference: float,alpha: float,Npop: int = 1000) -> tuple[float,list[list],list[list]]
```
```python
def continue_block(maf0: float,pre_matrix: list[list],matrix: list[list],positions: list[float],alpha: float,beta: float,gamma: float,strength: int = -1,Npop: int = 1000) -> tuple[float,list[list],list[list]]
```
```python
def genmatrix(blocks: list[int],strength: float,population: float,Npop: int)
```
```python
def gt_vcf(value: int)-> str
```
```python
def create_vcfgz(vcf_name: str,matrix: list[list],alpha: float,beta: float,gamma: float,system: str = "unix") -> None
```

#### LD Analytics

```python
def LD_corr_matrix(matrix: list[list]) -> tuple[list[list],float,list[float]]
```
```python
def LD_r2_matrix(pre_matrix: list[list]) -> tuple[list[list],float,list[float]]
```
```python
def display(rel: list[list],m: float) -> list[list[tuple[float,float,float]]]
```
```python
def minor_haplotype(sub_pre_matrix: list[list]) -> float
```


