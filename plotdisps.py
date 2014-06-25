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
  dw = 1.0/2/dt/n_omega
  omega = range(n_omega)
  omega = [dw*x for x in omega] # in THz
  print "    >>> w size = " + str(n_omega)
  print "    >>> k size = " + str(n_k)
  
  dispdata = np.log2(dispdata)
  dispdata.transpose()  #flip array for gnuplot
  #dispdata = np.flipud(dispdata)
  
  edges = [0, n_k,0, omega[-1] ]
  
  plt.clf()
  
  plt.imshow(dispdata,extent=edges,aspect='auto',interpolation='bicubic')
  plt.xticks(pxls,tags)
  plt.yticks(omega[::8])
  plt.colorbar()
  plt.show()
  
  
def readsetfile(fin):
    pxl = []
    tags = []
    mrkr2 = -1
    
    while True:
        line1 = fin.readline().strip().split()
        if not line1:
            break
        line2 = fin.readline().strip().split()
        line3 = fin.readline().strip().split()
        mrkr1 = mrkr2+1
        mrkr2 += int(line3[-1])
                
        if len(tags)>0:
            if tags[-1] == line1[-1]: 
                pxl.append(mrkr2)
                tags.append(line2[-1])
                continue
        
        pxl.append(mrkr1)
        pxl.append(mrkr2)
        
        tags.append(line1[-1])
        tags.append(line2[-1])
        
    pxl[-1] -=1
    
    return pxl,tags  #remove last entry since there are   no more sections
    
  

if __name__ == "__main__":
  main()
