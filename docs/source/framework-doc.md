# Documentation for the  ```Framework``` module
This module is a cloud-native tool for software development and experimentation. It provides enhanced versions of the functions contained in the module ```HaploDX``` for more efficient computations and memory usage. It is designed to be deployed and scaled on cloud computing platforms, making it ideal for large-scale projects. 

* [HaploDynamics.Framework.Model](#haplodynamicsframeworkmodel)
  * [HaploDynamics.Framework.Model.\_\_init\_\_](#haplodynamicsframeworkmodel__init__)
  * [HaploDynamics.Framework.Model.initiate_vcf](#haplodynamicsframeworkmodelinitiate_vcf)
  * [HaploDynamics.Framework.Model.generate_vcf](#haplodynamicsframeworkmodelgenerate_vcf)

> **Warning** (test)
> This module is currently under active development
>

&nbsp; 

---

## HaploDynamics.Framework.Model

* [HaploDynamics.Framework.Model.\_\_init\_\_](#haplodynamicsframeworkmodel__init__)
* [HaploDynamics.Framework.Model.initiate_vcf](#haplodynamicsframeworkmodelinitiate_vcf)
* [HaploDynamics.Framework.Model.generate_vcf](#haplodynamicsframeworkmodelgenerate_vcf)

&nbsp; 

### HaploDynamics.Framework.Model.\_\_init\_\_

This is the initializer method for the class ```Model```.

&nbsp; 

####  **Presentation**
```python
def __init__(self,
             fname: str,
             cpus: int = 1,
             system: str = "unix"
             ) -> None
```
 <ins>Description</ins>
 - **inputs:**
   - .
 - **output:**
   - .


&nbsp; 

####  **Tutorial**
```python
import HaploDynamics.Framework as fmx

model = fmx.Model("tutorial")
```

&nbsp; 

### HaploDynamics.Framework.Model.initiate_vcf

The function ```initiate_vcf()```  function creates a VCF file with only a header. This means that the file contains metadata and a row of attributes, but no genotypes. To be able to initialize the metadata and the labels for the population used in the file, the function ```initiate_vcf()``` requires the number of variant rows and the number of individuals to be contained in the file.

&nbsp; 

####  **Presentation**
```python
def initiate_vcf(self,
                 length: float,
                 Npop: int
                 ) -> None
```
 <ins>Description</ins>
 - **inputs:**
   - .
 - **output:**
   - .


&nbsp; 

####  **Tutorial**
```python
import HaploDynamics.Framework as fmx

model = fmx.Model("tutorial")

nb_variants = 120
nb_indiv = 1000

model.initiate_vcf(nb_variants,nb_indiv).close()
```

```shell
$ ls
tutorial.vcf.gz
$ zless tutorial.vcf.gz
...
```

&nbsp; 

### HaploDynamics.Framework.Model.generate_vcf

The process ```generate_vcf()``` combines the two processes ```genmatrix()``` and ```create_vcfgz()``` of the  module ```HaploDX``` into a single process. This is done to optimize execution time by reducing the amount of RAM space used during execution. The ```generate_vcf()``` function would ideally be used on machines with limited memory bandwidth due to cost or equipment constraints.

&nbsp; 

####  **Presentation**
```python
def generate_vcf(self,
                 blocks: list[int],
                 strength: float,
                 population: float,
                 Npop: int,
                 chrom: str = "23"
                 ) -> tuple[float,float,float]
```
 <ins>Description</ins>
 - **inputs:**
   - ```blocks```: a list of integer values used to encode blocks of genetic positions linked together by an amount of linkage disequilibrium determined by a decay function. The genetic positions generated for each block are used to index the data contained in the output VCF file;
   - ```strength```: a floating-point number in the interval $[-1,1]$ that represents the strength of the linkage disequilibrium. For example, a strength of $1$ refers to the maximum value that the linkage disequilibrium measure can take given the values of the parameters generated during the siumlation. See the <a href="https://www.normalesup.org/~tuyeras/node_diss/blg/home.php?page=blg_stat/stat_1/home.php">tutorial</a> to learn more about the parameter ```strength```;
   - ```population```: a floating-point number in the interval $[0,1]$, which is given as an input to the function ```population_mld()``` to generate a stochastic population profile;
   - ```Npop```: an integer $N$ that represents the number of individuals in the simulated dataset;
   - ```chrom```: a string that represents the chromosome number used to annotate the VCF file.
 - **output:**
   - ```speed```: execution time for the run of the process;
   - ```max_mem```: maximal memory usage during the run of the process;
   - ```cur_mem```: memory usage at the end the process;


&nbsp; 

####  **Tutorial**

The following snippet shows a redimentary example of how the user can use the method ```generate_vcf()```.

```python
import HaploDynamics.Framework as fmx

model = fmx.Model("tutorial")

blocks = [20]*7
strength = 0.1
population = 1
Npop = 1000
chrom = "1"

model.generate_vcf(blocks,strength,population,Npop,chrom)
```

```shell
$ python myscript.py
Model.generate_vcf: >||||||||||||||||||||< 100%
time (sec.): 0.7510931491851807
max. mem (MB): 0.11163139343261719
cur. mem (MB): 0.0834970474243164
$ ls
tutorial.vcf.gz
$ zless tutorial.vcf.gz
...
```
The following script compares the performance of the function ```generate_vcf()```, which is optimized for execution time and memory usage on a single processor, to the performance of the C++ software [MaCS](https://github.com/gchen98/macs/tree/master). Both MaCS and HaploDX rely on Bayesian networks to generate variants. This allows these two frameworks to output rows of genotypes sequentially, which contributes to their speed. This comparison is intended to demonstrate the performance of HaploDX relative to a similar tool.
```python
import matplotlib.pyplot as plt
import HaploDynamics.Framework as fmx

model = fmx.Model("tutorial")

#Region of 10 Mb
blocks = [20] * 500
strength = 0.1
population = 1
Npop = 1000
chrom = "1"

x = [100,300,1000,3000,10000]
s = []
m = []
for i in range(len(x)):
  Npop = x[i]
  ave_s = []
  ave_m = []
  for j in range(10):
    print("Npop:",Npop,j)
    speed, max_mem, cur_mem = model.generate_vcf(blocks,strength,population,Npop,chrom)
    ave_s.append(speed)
    ave_m.append(max_mem)
  s.append(sum(ave_s)/10)
  m.append(sum(ave_m)/10)

print(s)
print(m)
plt.plot(x,[13,43,3*60,13*60+5,50*60+50],c="r",linestyle='dashed',alpha = 0.5)
plt.plot(x,s,c="r")
plt.legend(["MaCS", "HaploDX"])
plt.ylabel("Average execution times for a region of "+str(sum(blocks))+"kb")
plt.xlabel("Npop")
plt.show()

plt.plot(x,[9.6,10,11.7,13.8,20.4],c="b",linestyle='dashed',alpha = 0.5)
plt.plot(x,m,c="b")
plt.legend(["MaCS", "HaploDX"])
plt.ylabel("Average memory usage for a region of "+str(sum(blocks))+"kb")
plt.xlabel("Npop")
plt.show()
```
The outputs for the previous script should typically look as shown below. Note that the execution time and memory usage statistics for MaCS were reportedly conducted on a 2 GHz processor with 32GB of RAM (see the supporting [research paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2612967/)). The same performance tests for HaploDX were conducted on a 2.4 GHz processor with 16GB of RAM, which is slightly slower.

Memory usage            |  Execution times
:-------------------------:|:-------------------------:
![](perform_space.png)  |  ![](perform_time.png)

&nbsp; 

####  **Dependencies from HaploDX**

The function ```generate_vcf()``` relies on the following processes from the module ```HaploDX```:
  - [HaploDynamics.HaploDX.population_mld](haplodx-doc.md#haplodynamicshaplodxpopulation_mld)
  - [HaploDynamics.HaploDX.gt_vcf](haplodx-doc.md#haplodynamicshaplodxgt_vcf)
  - [HaploDynamics.HaploDX.genotype_schema](haplodx-doc.md#haplodynamicshaplodxgenotype_schema)
  - [HaploDynamics.HaploDX.genotype](haplodx-doc.md#haplodynamicshaplodxgenotype)
  - [HaploDynamics.HaploDX.SNP_distribution](haplodx-doc.md#haplodynamicshaplodxsnp_distribution)
  - [HaploDynamics.HaploDX.cond_genotype_schema](haplodx-doc.md#haplodynamicshaplodxcond_genotype_schema)
  - [HaploDynamics.HaploDX.gref](haplodx-doc.md#haplodynamicshaplodxgref)
