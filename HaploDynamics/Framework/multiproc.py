from HaploDynamics.utils.load import Loading
from HaploDynamics.HaploDX import *
import multiprocessing as mp
import gzip
import random
import math
import datetime
#import os

# os.cpu_count()
# mp.cpu_count()

class Model(object):

  def __init__(self,fname,cpus = 1,system = "unix"):
    self.pool = mp.Pool(cpus)
    self.fname = fname
    self.eol = "\n" if system == "unix" else "\r\n"


  def initiate_vcf(self,length,Npop,chrom = "23"):
    f = gzip.open(self.fname+".vcf.gz","wt")
    order = int(math.log10(Npop))
    write = lambda s : f.write(s+self.eol)
    idnum = lambda s: "ID"+"0"*(order-len(str(s)))+str(s)
    write("##fileformat=VCFv4.2")
    write("##fileDate="+"".join(str(datetime.date.today()).split("-")))
    write("##source=HaploDX")
    write("##reference=https://github.com/remytuyeras/HaploDynamics/blob/main/README.md")
    write("##contig=<ID="+chrom+",length="+str(length)+",species=\"simulated Homo sapiens\">")
    write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">")
    write("##INFO=<ID=AP,Number=1,Type=Float,Description=\"Alpha Parameter for Simulation\">")
    write("##INFO=<ID=BP,Number=1,Type=Float,Description=\"Beta Parameter for Simulation\">")
    write("##INFO=<ID=CP,Number=1,Type=Float,Description=\"Gamma Parameter for Simulation\">")
    write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    INDIV = "\t".join(map(idnum,range(1,Npop)))
    write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+INDIV)
    f.flush()
    return f

  def generate_vcf(self,blocks,strength,population,Npop,chrom = "23"):
    #Values needed to make the VCF file
    alpha, beta, gamma = population_mld(population)
    nuc   = ["A","C","G","T"]
    #Function to write a variant's row in the VCF file
    def write_genotype(f,vector):
      A = random.choice(nuc)
      B = random.choice(nuc[:nuc.index(A)]+nuc[nuc.index(A)+1:])
      info = [chrom, str(int(1000*vector[0])), ".",A,B, ".", "PASS", "NS="+str(Npop)+";AP="+str(alpha)+";BP="+str(beta)+";CP="+str(gamma), "GT"]
      genotypes = [gt_vcf(vector[j]) for j in range(1,len(vector))]
      f.write("\t".join(info+genotypes)+self.eol)
      f.flush()
    #Important parameters to control the memory space
    pre_vector  = None
    vector      = None
    length      = sum(blocks)+1
    #For the loop
    reference   = 0
    maf         = None
    #For the output
    speed = None
    max_mem = None
    cur_mem = None
    with self.initiate_vcf(length,Npop,chrom) as f:
      ####### initiate_block
      maf, hwp = genotype_schema(alpha)
      minor = random.choice([0,2])
      vector = [reference]
      pre_vector = [genotype(hwp,minor) for i in range(Npop)]
      vector.extend([pre_vector[i][1] for i in range(Npop)])
      write_genotype(f,vector)
      ####### end of initiate_block

      #Loop: Loading bar on the standard output
      loading_frame = Loading("Model.generate_vcf: ",len(blocks))
      for i, block in enumerate(blocks):
        #Loading bar progress
        speed, max_mem, cur_mem = loading_frame.show(i+1)
        positions = SNP_distribution(reference,block)
        ####### continue_block
        for k in range(len(positions)):
          distance = abs(positions[k]-reference)
          maf, hwp, ld = cond_genotype_schema(maf,distance,alpha,beta,gamma,strength)
          minor = random.choice([0,2])
          vector = [positions[k]]
          pre_vector = [genotype(hwp(gref(pre_vector[i])),minor) for i in range(Npop)]
          vector.extend([pre_vector[i][1] for i in range(Npop)])
          write_genotype(f,vector)
        ####### end of continue_block

        reference = positions[-1]
      #end of loop
    return speed, max_mem, cur_mem

#if __name__ == '__main__':

  #model = Model("hello")
  ##model.pool.map(a_function,a_list)
  #model.close()
