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
    self.initial_schema = None
    self.landscape = []
    self.initiate_landscape()

  def initiate_landscape(self,afs=afs_sample,reference=0):
    self.initial_schema = (afs,reference)

  def extend_landscape(self,*schemas):
    self.landscape.extend(list(schemas))
    return len(self.landscape)

  @staticmethod
  def standard_schema(npos):
    pos_dist = lambda reference: SNP_distribution(reference,npos)
    return npos, afs_sample, pos_dist

  def genotype_schema(self,alpha=4/30,afs=afs_sample):
    p =  afs(alpha)
    #We use the minor allele frequency as a reference
    maf = min(p,1-p)
    #Probability for minor homogeneous genotypes
    hom = maf**2
    #Probability for heterogeneous genotypes
    het = maf*(1-maf)
    #Interval decomposition simulating the Hardy-Weinberg principle
    hwp = [0, hom, hom+het, hom+2*het, 1]
    return maf, hwp

  def linkage_disequilibrium(self,alpha,beta,gamma,strength=-1,afs=afs_sample):
    ub = gamma*math.log2(1/beta)
    shift = random.uniform(min(ub,strength*ub),ub)
    def LD(p):
      def ld(t):
        q = afs(alpha)
        while not(lb_freq(beta,gamma,p,t,shift) <= q <= ub_freq(beta,gamma,p,t,shift)):
          q = afs(alpha)
        l = decay(amplifier(beta,p,q),2*gamma,shift)(t)+1
        return l, q
      return ld
    return LD

  def cond_genotype_schema(self,previous_maf,distance,alpha,beta,gamma,strength=-1,afs=afs_sample):
    #Linkage disequilibrium
    ld, maf = self.linkage_disequilibrium(alpha,beta,gamma,strength,afs)(previous_maf)(distance)
    while maf > 0.5:
      ld, maf = self.linkage_disequilibrium(alpha,beta,gamma,strength,afs)(previous_maf)(distance)
    #Conditional Hardy-Weinberg principle encoded as a function
    def hwp(previous_genotype):
      #Previous genotype has minor alleles
      if previous_genotype == 2:
        #Probability for minor homogeneous genotypes
        hom = (ld**2) * (maf**2)
        #Probability for heterogeneous genotypes
        het = ld * ref_alt_function(maf,ld) * maf * (1-maf)
        #Interval decomposition simulating the conditional Hardy-Weinberg principle
        return [0, hom, hom+het,hom+2*het, 1]
      #Previous genotype has one minor allele
      if abs(previous_genotype) == 1:
        #Probability for minor homogeneous genotypes
        hom = ld * ref_alt_function(previous_maf,ld) * (maf**2)
        #Probability for heterogeneous genotypes
        het1 = ld * alt_alt_function(previous_maf,maf,ld) * maf * (1-maf)
        het2 = ref_alt_function(previous_maf,ld) * ref_alt_function(maf,ld) * maf * (1-maf)
        #Previous genotype consists of a minor-major pairing
        if previous_genotype == -1:
          return [0, hom, hom+het1,hom+het1+het2, 1]
        #Previous genotype consists of a major-minor pairing
        elif previous_genotype == 1:
          return [0, hom, hom+het2,hom+het2+het1, 1]
      #Previous genotype has no minor alleles
      if previous_genotype == 0:
        #Probability for minor homogeneous genotypes
        hom = (ref_alt_function(previous_maf,ld)**2) * (maf**2)
        #Probability for heterogeneous genotypes
        het = ref_alt_function(previous_maf,ld) * alt_alt_function(previous_maf,maf,ld) * maf * (1-maf)
        #Interval decomposition simulating the conditional Hardy-Weinberg principle
        return [0, hom, hom+het,hom+2*het, 1]
    return maf, hwp, ld

  def initiate_vcf(self,Npos,Npop,chrom = "23"):
    f = gzip.open(self.fname+".vcf.gz","wt")
    order = int(math.log10(Npop))+1
    write = lambda s : f.write(s+self.eol)
    idnum = lambda s: "ID"+"0"*(order-len(str(s)))+str(s)
    write("##fileformat=VCFv4.2")
    write("##fileDate="+"".join(str(datetime.date.today()).split("-")))
    write("##source=HaploDX")
    write("##reference=https://github.com/remytuyeras/HaploDynamics/blob/main/README.md")
    write("##contig=<ID="+chrom+",length="+str(Npos)+",species=\"simulated Homo sapiens\">")
    write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">")
    write("##INFO=<ID=AP,Number=1,Type=Float,Description=\"Alpha Parameter for Simulation\">")
    write("##INFO=<ID=BP,Number=1,Type=Float,Description=\"Beta Parameter for Simulation\">")
    write("##INFO=<ID=CP,Number=1,Type=Float,Description=\"Gamma Parameter for Simulation\">")
    write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    INDIV = "\t".join(map(idnum,range(1,Npop+1)))
    write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+INDIV)
    f.flush()
    return f

  def generate_vcf(self,strength,population,Npop,chrom = "23"):
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
    Npos        = sum([npos for npos,_,_ in self.landscape])+1
    #For the loop
    maf         = None
    #For the output
    speed = None
    max_mem = None
    cur_mem = None
    with self.initiate_vcf(Npos,Npop,chrom) as f:
      ####### initiate_block
      afs, reference = self.initial_schema
      maf, hwp = self.genotype_schema(alpha,afs)
      minor = random.choice([0,2])
      vector = [reference]
      pre_vector = [genotype(hwp,minor) for i in range(Npop)]
      vector.extend([pre_vector[i][1] for i in range(Npop)])
      write_genotype(f,vector)
      ####### end of initiate_block

      #Loop: Loading bar on the standard output
      loading_frame = Loading("Model.generate_vcf: ",len(self.landscape))
      for i, schema in enumerate(self.landscape):
        #Loading bar progress
        speed, max_mem, cur_mem = loading_frame.show(i+1)
        positions = schema[2](reference)
        ####### continue_block
        for k in range(len(positions)):
          distance = abs(positions[k]-reference)
          maf, hwp, ld = self.cond_genotype_schema(maf,distance,alpha,beta,gamma,strength,schema[1])
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
