from HaploDynamics.utils.load import Loading
import gzip
import random
import math
#import matplotlib.pyplot as plt
import scipy.stats as stat
import datetime

def stochastic_line(a,b,sigma):
  return lambda t: random.gauss(a*t+b*(1-t),sigma)

population_mut = stochastic_line(0.08,0.17,0.01)

#print(population_mut(0.5))

def afs_distribution(index,alpha=4/30):
  f = (0.5+index)/15
  p = 0
  if 0 < f < 0.25:
    p = 480*(0.1-alpha)*f + 120*alpha
  if 0.25 < f < 1:
    p = -16*(f-1)
  return p/(15*alpha+6)

#x = [(0.5+i)/15 for i in range(15)]
#plt.plot(x,[afs_distribution(i,alpha=0.16) for i in range(15)],color='green')
#plt.plot(x,[afs_distribution(i,alpha=0.133) for i in range(15)],color='grey',linestyle="--")
#plt.plot(x,[afs_distribution(i,alpha=0.09) for i in range(15)],color='red')
#plt.plot(x,[afs_distribution(i,alpha=0.08) for i in range(15)],color='blue')
#plt.grid(True)
#plt.show()

def afs_intervals(pick,alpha=4/30):
  intervals = list()
  beg = 0
  for i in range(15):
    end = beg + afs_distribution(i,alpha)
    if max(beg/15,0) <= pick < min(end/15,1):
      return [i/15,(i+1)/15]
    beg = end
  #in case afs_distribution(14,alpha)/15 <= pick < 1
  return [14/15,1]

def afs_sample(alpha=4/30):
  return random.uniform(*afs_intervals(random.random(),alpha))

#plt.hist([afs_sample(0.16) for i in range(15000)],15,facecolor='green',alpha=0.5)
#plt.hist([afs_sample(0.09) for i in range(15000)],15,facecolor='red',alpha=0.5)
#plt.hist([afs_sample(0.08) for i in range(15000)],15,facecolor='blue',alpha=0.5)
#plt.grid(True)
#plt.show()

def genotype_schema(alpha=4/30):
  p =  afs_sample(alpha)
  #We use the minor allele frequency as a reference
  maf = min(p,1-p)
  #Probability for minor homogeneous genotypes
  hom = maf**2
  #Probability for heterogeneous genotypes
  het = maf*(1-maf)
  #Interval decomposition simulating the Hardy-Weinberg principle
  hwp = [0, hom, hom+het, hom+2*het, 1]
  return maf, hwp

def genotype(hwp,minor):
  t = random.random()
  #print("generator: ",t,hwp)
  for i in range(4):
    #We sample an element from the Hardy-Weinberg distribution
    if hwp[i] <= t < hwp[i+1]:
      sign = 1-minor
      genotype_code = minor*(i==0) + sign*(i==1) - sign*(i==2) + (2-minor)*(i==3)
      return minor, genotype_code

#For reference, this function maps the homogeneous genotypes consisting of two major alleles to 0
gref= lambda g: g[1] if g[0] == 2 else (2-g[1])*(g[1] in [0,2]) - g[1]*(g[1] in [1,-1])

#alpha_mut = population_mut(0.5)
#sch = genotype_schema(alpha_mut)
#minor = random.choice([0,2])
#g = [genotype(sch[1],minor) for i in range(1000)]
#plt.hist([gref(g[i]) for i in range(1000)],facecolor='green',alpha=0.5)
#plt.grid(True)
#plt.show()

def population_mld(t):
  beta = stochastic_line(0.62,0.45,0.01)
  gamma = stochastic_line(19,9,0.1)
  return (population_mut(t), beta(t), gamma(t))

#fig = plt.figure()
#ax = plt.axes(projection='3d')
#ax.scatter3D([0.08,0.09,0.17], [0.62,0.6,0.45], [19,18,9], marker="o",color="red")
#p = list()
#for i in range(1000):
  #p.append(population_mld(0.01*int(0.1*i)))
  #print(0.01*int(0.1*i),p[i])
#x,y,z = zip(*p)
#ax.scatter3D(x,y,z,color="blue",alpha=0.05)
#plt.show()


def decay(initial,halfwidth,shift):
  ln_2 = math.log(2)/math.log(math.exp(1))
  l = ln_2/halfwidth
  return lambda t: initial * math.exp(-l*(t-shift))

#x = [0.2*i for i in range(1000)]
#y = [decay(0.5,20,random.gauss(0,0.02)*t)(t) for t in x]
#plt.plot(x,y,color='blue',alpha=0.5)
#plt.ylim([0, 1])
#plt.grid(True)
#plt.show()

def ref_alt_function(y,x):
  return (1-y*x)/(1-y)

#J = [(0,550),(550,850),(850,955)]
#for j in range(len(J)):
  #I = [0.001*i for i in range(*J[j])]
  #ld = [0.9,1.05,1.15,1.35,1.55,1.80]
  #colors = ["grey","blue","green","red","orange","purple"]
  #for i in range(len(ld)):
    #y = [ref_alt_function(p,ld[i]) for p in I]
    #plt.plot(I,y,color=colors[i],alpha=0.5)
  #plt.legend(["L = "+str(i) for i in ld])
  #plt.grid(True)
  #plt.show()

def alt_alt_function(y,z,x):
  return (1-y-z+y*z*x)/(1-y-z+y*z)

#J = [(1,1000),(1,130),(130,300),(300,950),(950,1000)]
#for j in range(len(J)):
  #I = [0.001*i for i in range(*J[j])]
  #ld = [0.9,1.05,1.15,1.35,1.55,1.80]
  #colors = ["grey","blue","green","red","orange","purple"]
  #for i in range(len(ld)):
    #y = [alt_alt_function(p,0.25/(p*ld[i]),ld[i]) for p in I]
    #plt.plot(I,y,color=colors[i],alpha=0.5)
  #plt.legend(["L = "+str(i) for i in ld])
  #plt.grid(True)
  #plt.show()

def amplifier(beta,p,q,s=1):
  return s*math.sqrt(beta*(1/p-1)*(1/q-1))

#print(amplifier(0.5,0.1,0.15,1))

def lb_freq(beta,gamma,previous_freq,distance,shift):
  top = math.sqrt(beta)
  e = decay(1,-2*gamma,shift)(distance)*(1/previous_freq-1)
  return top/(e+top)

def ub_freq(beta,gamma,previous_freq,distance,shift):
  top = decay(1,-2*gamma,shift)(distance)
  e = math.sqrt(beta)*(1/previous_freq-1)
  return top/(e+top)

#_, beta_mld, gamma_mld = population_mld(0.5)
#freqs = [0.001*i for i in range(1,1000)]
#style = ["-","--",":"]
#t = [3,25,50]
#for i in range(len(t)):
  #plt.plot(freqs,[ub_freq(beta_mld,gamma_mld,p,t[i],0) for p in freqs],color="green",linestyle=style[i])
  #plt.plot(freqs,[lb_freq(beta_mld,gamma_mld,p,t[i],0) for p in freqs],color="red",linestyle=style[i])
#plt.plot(freqs,freqs,color='grey',linestyle="--")
#plt.ylim([0, 1])
#plt.legend(["upperbound"*(i%2==0) + "lowerbound"*(i%2==1) + " for t="+ str(t[int(i/2)]) for i in range(2*len(t))])
#plt.grid(True)
#plt.show()

def linkage_disequilibrium(alpha,beta,gamma,strength=-1):
  ub = gamma*math.log2(1/beta)
  shift = random.uniform(min(ub,strength*ub),ub)
  def LD(p):
    def ld(t):
      q = afs_sample(alpha)
      while not(lb_freq(beta,gamma,p,t,shift) <= q <= ub_freq(beta,gamma,p,t,shift)):
        q = afs_sample(alpha)
      l = decay(amplifier(beta,p,q),2*gamma,shift)(t)+1
      return l, q
    return ld
  return LD

#I = [0.2*i for i in range(1,1000)]
#colors = ["green","red","blue"]
#parameters = [(0.17,0.45,9),(0.09,0.6,18),(0.08,0.62,19)]

#q_graph, l_graph = [], []
#for i in range(3):
  #l,q = zip(*[linkage_disequilibrium(*parameters[i])(0.1)(t) for t in I])
  #l_graph.append(l)
  #q_graph.append(q)

#for i in range(3):
  #plt.scatter(I,q_graph[i],color=colors[i],alpha=0.5,s=10)
#plt.legend(["African","European","Asian"])
#plt.grid(True)
#plt.show()

#for i in range(3):
  #plt.scatter(I,l_graph[i],color=colors[i],alpha=0.5,s=10)
#plt.legend(["African","European","Asian"])
#plt.grid(True)
#plt.show()

def cond_genotype_schema(previous_maf,distance,alpha,beta,gamma,strength=-1):
  #Linkage disequilibrium
  ld, maf = linkage_disequilibrium(alpha,beta,gamma,strength)(previous_maf)(distance)
  while maf > 0.5:
    ld, maf = linkage_disequilibrium(alpha,beta,gamma,strength)(previous_maf)(distance)
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

#alpha_mld, beta_mld, gamma_mld = population_mld(0.5)

#maf1, hwp1 = genotype_schema(alpha_mld)
#minor1 = random.choice([0,2])
#g1 = [gref(genotype(hwp1,minor1)) for i in range(1000)]
#plt.hist([g1[i] for i in range(1000)],facecolor='green',alpha=0.5)
#plt.title("maf ~" + str(int(maf1*100)/100),fontsize = 40)
#plt.grid(True)
#plt.show()

##We let the distance from SNP1 to SNP2 be 1kb
#distance = 1
##We let the distance from SNP1 to SNP2 be 200kb to model independence
##distance = 200

#maf2, hwp2, ld2 = cond_genotype_schema(maf1,distance,alpha_mld,beta_mld,gamma_mld)
#minor2 = random.choice([0,2])
#g2 = [gref(genotype(hwp2(g1[i]),minor2)) for i in range(1000)]
#plt.hist([g2[i] for i in range(1000)],facecolor="green",alpha=0.5)
#plt.title("maf ~" + str(int(maf2*100)/100),fontsize = 40)
#plt.grid(True)
#plt.show()

#plt.scatter([g1[i]+random.gauss(0,.05) for i in range(1000)],[g2[i]+random.gauss(0,.05) for i in range(1000)],color = "red", s=12,alpha=0.1)
#plt.title("ld ~" + str(int(ld2*100)/100),fontsize = 40)
#plt.grid(True)
#plt.show()
##Not in tutorial
##print("hwp1: ",hwp1)
##for x in [0,1,-1,2]:
  ##print("hwp2: "+str(x)+":",hwp2(x))

def SNP_distribution(reference,length):
  positions = list()
  d_saved = list()
  for i in range(int(length)):
    d = int(1000*random.uniform(0.002,length))/1000
    while d in d_saved:
      d = int(1000*random.uniform(0.002,length))/1000
    d_saved.append(d)
    positions.append(reference+d)
  return sorted(positions)

#positions = SNP_distribution(0,200)
#plt.scatter(range(200),positions,s=2)
#plt.show()

def initiate_block(reference,alpha,Npop=1000):
  #Two types of tables generated
  pre_matrix = []
  matrix = []
  #The genomic data is indexed by genomic positions
  matrix.append([reference])

  maf, hwp = genotype_schema(alpha)
  minor = random.choice([0,2])
  pre_matrix.append([genotype(hwp,minor) for i in range(Npop)])
  #Here, we do not use the function gref, but the genotype code
  matrix[-1].extend([pre_matrix[-1][i][1] for i in range(Npop)])
  return maf, pre_matrix, matrix

#alpha_mut = population_mut(0.5)
#maf, pre_matrix, matrix = initiate_block(0,alpha_mut)
#plt.hist(matrix[-1][1:],facecolor="green",alpha=0.5)
#plt.title("maf ~" + str(int(maf*100)/100),fontsize = 40)
#plt.grid(True)
#plt.show()

def continue_block(maf0,pre_matrix,matrix,positions,alpha,beta,gamma,strength=-1,Npop=1000):
  maf = None
  reference = matrix[-1][0]
  for k in range(len(positions)):
    distance = abs(positions[k]-reference)
    matrix.append([positions[k]])

    maf, hwp, ld = cond_genotype_schema(maf0,distance,alpha,beta,gamma,strength)
    minor = random.choice([0,2])
    pre_matrix.append([genotype(hwp(gref(pre_matrix[-1][i])),minor) for i in range(Npop)])
    #Here, we do not use the function gref, but the genotype code
    matrix[-1].extend([pre_matrix[-1][i][1] for i in range(Npop)])

  return maf, pre_matrix, matrix

#alpha_mld, beta_mld, gamma_mld = population_mld(0.5)
#maf, pre_matrix, matrix = initiate_block(0,alpha_mld)
#plt.hist(matrix[-1][1:],facecolor="green",alpha=0.5)
#plt.title("maf ~" + str(int(maf*100)/100),fontsize = 40)
#plt.grid(True)
#plt.show()

#positions = SNP_distribution(0,200)
#print(positions[:10])

#maf, pre_matrix, matrix = continue_block(maf,pre_matrix,matrix,positions,alpha_mld,beta_mld,gamma_mld)
#print(matrix[-1][:20])
#plt.hist(matrix[-1][1:],facecolor="green",alpha=0.5)
#plt.title("maf ~" + str(int(maf*100)/100),fontsize = 40)
#plt.grid(True)
#plt.show()

def LD_corr_matrix(matrix):
  rel, m = list(), 0
  dist = [[] for _ in range(len(matrix))]
  for i in range(len(matrix)):
    rel.append([])
    for j in range(len(matrix)):
      c = stat.pearsonr(matrix[i][1:],matrix[j][1:])[0] if j>i else 0
      dist[int(abs(matrix[i][0]-matrix[j][0]))].append(abs(c))
      m = max(abs(c),m) if not(math.isnan(c)) else m
      rel[i].append(c)
  dist = [sum(dist[k])/len(dist[k]) if len(dist[k]) !=0 else 0 for k in range(len(dist))]
  return rel, m, dist

def display(rel,m):
  mod = lambda x: int((1-abs(x))*255)
  fcolor = lambda x: [mod(0),mod(x),mod(x)] if x>0 else [mod(x),mod(x),mod(0)]
  image = list()
  for i in range(len(rel)):
    image.append([])
    for j in range(len(rel[i])):
      image[i].append(fcolor(rel[i][j]/m) if not(math.isnan(rel[i][j])) else [200,200,200])
  return image

#for strength in [-1,0,1]:
  #reference = 0
  #alpha_mld, beta_mld, gamma_mld = population_mld(0.5)
  #maf, pre_matrix, matrix = initiate_block(reference,alpha_mld)
  #for _ in range(6):
    #positions = SNP_distribution(reference,20)
    #maf, pre_matrix, matrix = continue_block(maf,pre_matrix,matrix,positions,alpha_mld,beta_mld,gamma_mld,strength)
    #reference = positions[-1]

  #rel, m, dist = LD_corr_matrix(matrix)
  #plt.imshow(display(rel,m))
  #plt.show()
  #plt.plot([i for i in range(len(dist)-1)],dist[1:])
  #plt.ylim([0, 1])
  #plt.show()

def minor_haplotype(sub_pre_matrix):
  count = 0
  n,m = len(sub_pre_matrix), (len(sub_pre_matrix[0]) if len(sub_pre_matrix) != 0 else 0)
  for i in range(m):
    chr1 = True
    chr2 = True
    for j in range(n):
      chr1 = chr1 and (gref(sub_pre_matrix[j][i]) == 2 or gref(sub_pre_matrix[j][i]) == 1)
      chr2 = chr2 and (gref(sub_pre_matrix[j][i]) == 2 or gref(sub_pre_matrix[j][i]) == -1)
    count = count + (1 if chr1 else 0) + (1 if chr2 else 0)
  return (count/(2*m) if m > 0 else 0)

def LD_r2_matrix(pre_matrix):
  rel, m = list(), 0
  dist = [[] for _ in range(len(pre_matrix))]
  for i in range(len(pre_matrix)):
    #Part of the code differing from LD_corr_matrix
    rel.append([])
    p1 = minor_haplotype([pre_matrix[i]])
    for j in range(len(pre_matrix)):
      p2, p12 = (minor_haplotype([pre_matrix[j]]), minor_haplotype([pre_matrix[i],pre_matrix[j]])) if j>i else (0, 0)
      r = (p12-p1*p2)**2/(p1*p2*(1-p2)*(1-p1)) if 0<p1<1 and 0<p2<1 else 0
      dist[abs(i-j)].append(r)
      m = max(abs(r),m)
      rel[i].append(r)
  #Same ending as in LD_corr_matrix
  dist = [sum(dist[k])/len(dist[k]) for k in range(len(dist))]
  return rel, m, dist

#for strength in [-1,0,1]:
  #reference = 0
  #alpha_mld, beta_mld, gamma_mld = population_mld(0.5)
  #maf, pre_matrix, matrix = initiate_block(reference,alpha_mld)
  #for _ in range(6):
    #positions = SNP_distribution(reference,20)
    #maf, pre_matrix, matrix = continue_block(maf,pre_matrix,matrix,positions,alpha_mld,beta_mld,gamma_mld,strength)
    #reference = positions[-1]

  #rel, m, dist = LD_r2_matrix(pre_matrix)
  #plt.imshow(display(rel,m))
  #plt.show()
  #plt.plot([i for i in range(len(dist)-1)],dist[1:])
  #plt.ylim([0, 1])
  #plt.show()

def genmatrix(blocks,strength,population,Npop):
  reference = 0
  alpha_mld, beta_mld, gamma_mld = population_mld(population)
  maf, pre_matrix, matrix = initiate_block(reference,alpha_mld,Npop)
  loading_frame = Loading("genmatrix: ",len(blocks))
  for i in range(len(blocks)):
    loading_frame.show(i+1)
    positions = SNP_distribution(reference,blocks[i])
    maf, pre_matrix, matrix = continue_block(maf,pre_matrix,matrix,positions,alpha_mld,beta_mld,gamma_mld,strength,Npop)
    reference = positions[-1]
  return matrix, alpha_mld, beta_mld, gamma_mld

def gt_vcf(value):
  if value == 0:
    return "0|0"
  elif value == 1:
    return "0|1"
  elif value == -1:
    return "1|0"
  elif value == 2:
    return "1|1"
  else:
    return ".|."

def create_vcfgz(vcf_name,matrix,alpha,beta,gamma,system="unix"):
  f = gzip.open(vcf_name+".vcf.gz","wt")
  eol = "\n" if system=="unix" else "\r\n"
  f.write("##fileformat=VCFv4.2"+eol)
  f.write("##fileDate="+"".join(str(datetime.date.today()).split("-"))+eol)
  f.write("##source=HaploDX"+eol)
  f.write("##reference=https://github.com/remytuyeras/HaploDynamics/blob/main/README.md"+eol)
  f.write("##contig=<ID=23,length="+str(len(matrix))+",species=\"simulated Homo sapiens\">"+eol)
  f.write("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"+eol)
  f.write("##INFO=<ID=AP,Number=1,Type=Float,Description=\"Alpha Parameter for Simulation\">"+eol)
  f.write("##INFO=<ID=BP,Number=1,Type=Float,Description=\"Beta Parameter for Simulation\">"+eol)
  f.write("##INFO=<ID=CP,Number=1,Type=Float,Description=\"Gamma Parameter for Simulation\">"+eol)
  f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"+eol)
  INDIV = "\t".join(map(lambda s: "ID"+"0"*(6-len(str(s)))+str(s),range(1,len(matrix[0]))))
  f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+INDIV+eol)
  f.flush()
  nuc = ["A","C","G","T"]
  loading_frame = Loading("create_vcfgz: ",len(matrix))
  for i in range(len(matrix)):
    loading_frame.show(i+1)
    A = random.choice(nuc)
    B = random.choice(nuc[:nuc.index(A)]+nuc[nuc.index(A)+1:])
    info = ["23", str(int(1000*matrix[i][0])), ".", A, B, ".", "PASS", "NS="+str(len(matrix[i])-1)+";AP="+str(alpha)+";BP="+str(beta)+";CP="+str(gamma), "GT"]
    genotypes = [gt_vcf(matrix[i][j]) for j in range(1,len(matrix[i]))]
    f.write("\t".join(info+genotypes)+eol)
    f.flush()
  f.close()

#simulated_data = genmatrix([20,5,20,35,30,15],strength=1,population=0.1,Npop=1000)
#create_vcfgz("genomic-data.simulation.v1",*simulated_data)
