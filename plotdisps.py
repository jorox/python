import numpy as np
import sys
import matplotlib.pyplot as plt

def main():
  print "hello world"
  dispfname = sys.argv[1]
  infofname = sys.argv[2]
  
  print "... opening file "+dispfname
  dispfl = open(dispfname)
  setfl = open(infofname)
  dt = float(setfl.readline())
  [pxls,tags] = readsetfile(setfl)
  print "pxls = ", 
  print pxls
  print "tags = ",
  print tags
  dispdata = np.loadtxt(dispfl)
  n_omega = len(dispdata)
  n_k = len(dispdata[0])
  dw = 1/dt/n_omega
  print "    >>> w size = " + str(n_omega)
  print "    >>> k size = " + str(n_k)
  
  dispdata = np.log2(dispdata)
  dispdata.transpose()  #flip array for gnuplot
  dispdata = np.flipud(dispdata)
  
  edges = [0, n_k,0, n_omega*0.01677777/3.0e-2 ]
  
  plt.clf()
  
  plt.imshow(dispdata,extent=edges,aspect='auto',interpolation='bicubic')
  plt.colorbar()
  plt.show()
  
  
def readsetfile(fin):
    pxl = [0]
    tags = []
    while True:
        line1 = fin.readline().strip().split()
        line2 = fin.readline().strip().split()
        line3 = fin.readline().strip().split()
        
        if not line1: break

        pxl.append(pxl[-1] + int(line3[0])-1)
        pxl.append(pxl[-1]+1)
        
        tags.append(line1[-1])
        tags.append(line2[-1])
    
    pxl = pxl[:-1]
    return pxl,tags  #remove last entry since there are   no more sections
    
  

if __name__ == "__main__":
  main()