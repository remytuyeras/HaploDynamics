import sys
import time
import tracemalloc

class Loading(object):

  def __init__(self,name: str,final: int):
    self.t0 = time.time()
    self.final = final
    self.previous = 0
    self.current = 0
    self.progress = 0
    self.name = name
    tracemalloc.start()
    sys.stdout.write(name+self.text(self.current))
    sys.stdout.flush()
    self.speed = None
    self.max_mem = None
    self.cur_mem = None

  def assess(self,state: int):
    return int(20*state/self.final)

  def text(self,d: int,color: bool = False):
    c1, c2 = ("", "") if not color else ("\033[93m", "\033[0m")
    info = ""
    if d == 20:
      self.speed = time.time() - self.t0
      t= "time (sec.): "+str(self.speed)
      mem = tracemalloc.get_traced_memory()
      self.max_mem = mem[0]/(1024*1024)
      self.cur_mem = mem[1]/(1024*1024)
      m1 = "max. mem (MB): " + str(self.max_mem)
      m0 = "cur. mem (MB): " + str(self.cur_mem)
      info = "\n".join(["",t,m1,m0,""])
    return ">"+c1+"|"*d+" "*(20-d)+c2+"< "+str(int(100*d/20))+"%"+info

  def show(self,state: int):
    self.current = state
    d = self.assess(self.current)
    if self.progress != d:
      self.progress = d
      sys.stdout.write("\b"*len(self.text(self.assess(self.previous))))
      sys.stdout.write(self.text(d,color = True))
      sys.stdout.flush()
    self.previous = self.current
    return self.speed, self.max_mem, self.cur_mem

#x = Loading("loading: ",1000)

#for i in range(1000):
  #x.show(i+1)
  #time.sleep(0.005)
