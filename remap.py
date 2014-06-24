"""
This file is meant to rempa atoms from one basis to the other
after running lammps in order to unfold the Brillouin zone.
Needed: old map file + old basis + new basis + some knowledge of the 
        atoms in the new unit cells
"""
import sys
import numpy as np

def main():
    
    base1_e = readbasefile(sys.argv[1])
    base2 = readbasefile(sys.argv[2])
    
### Claculate the old vectors in the new basis
    cbm = np.array(base2[:3])
    cbm = np.transpose(cbm)
    
    INVcbm = np.linalg.inv(cbm)
    
    base1_b2 = [np.dot(INVcbm,np.transpose(np.array(base1_e[i]))) for i in range(3)]
    base1_b2 = [x.tolist() for x in base1_b2]
    
    print " ... original = "
    print ("      %1.4f  %1.4f  %1.4f\n      %1.4f  %1.4f  %1.4f\n      %1.4f  %1.4f  %1.4f" % 
           (base1_e[0][0], base1_e[1][0], base1_e[2][0],
            base1_e[0][1], base1_e[1][1], base1_e[2][1],
            base1_e[0][2], base1_e[1][2], base1_e[2][2])) 
    
    print " ... new = "
    print ("      %1.4f  %1.4f  %1.4f\n      %1.4f  %1.4f  %1.4f\n      %1.4f  %1.4f  %1.4f" % 
           (base1_b2[0][0], base1_b2[1][0], base1_b2[2][0],
            base1_b2[0][1], base1_b2[1][1], base1_b2[2][1],
            base1_b2[0][2], base1_b2[1][2], base1_b2[2][2]))
    
    print " ... basis = "
    print base2[:3]
    
### Read map file to build atoms
    oldmap = sys.argv[3]
    outmap = sys.argv[4]
    map_file_in  = open(oldmap)
    map_file_out = open(outmap, 'w')
    atm = map_file_in.readline()
    atm = map_file_in.readline()
    
    newmap = []
    minx = miny = minz = 0
    atm = map_file_in.readline()
    
    while atm:
        atm = atm.strip().split()
        
        atmindx = int(atm[3])
        atmid = int(atm[5])
        atmtyp = int(atm[4])
        
        atm = atm[:3]
        atm = [float(x) for x in atm]
        
        for i in range(3):
            atm[i] += base1_e[atmindx+3][i]
        
        newatm = [0,0,0]
        for i in range(3):
            newatm[i] = atm[0]*base1_b2[0][i] + atm[1]*base1_b2[1][i] + atm[2]*base1_b2[2][i]
        
        newatm.append(float(atmtyp))
        newatm.append(float(atmid))
        
        newmap.append(newatm)
        
        if minx > newatm[0]: minx = newatm[0]
        if miny > newatm[1]: miny = newatm[1]
        if minz > newatm[2]: minz = newatm[2]
        atm = map_file_in.readline()
        
    minx = np.around(minx)
    miny = np.around(miny)
    minz = np.around(minz)
    
    print " ... min unit cell is = %1.4f %1.4f %1.4f " % (minx,miny,minz)
    print " ... shifting "
    
    for i in range(len(newmap)):
        newmap[i][0] -= minx
        newmap[i][1] -= miny
        newmap[i][2] -= minz
        
    minx=miny=minz=0
    
    
    maxs = np.amax(np.array(newmap),axis=0)
    print maxs
    
    map_file_out.write("%1.0f %1.0f %1.0f %1.0f\n#l1 l2 l3 k type id" % (maxs[0]+1,maxs[1]+1,maxs[2]+1,len(base2[3:]) ))
    # the +1 is because there is always a unit cell 0,0,0
    
    for newatm in newmap:
        temp = [x-int(x) for x in newatm]
        newatmindx = match(temp[:3],base2[3:])
        
        if newatmindx == -1:
            print "!!! Error unmatched atom id = " + str(atmid)
            print "!!! basis atoms = "
            print base2[3:]
            return 0
        
        map_file_out.write("\n%1.0f %1.0f %1.0f %1.0f %1.0f %1.0f" % 
                               (newatm[0], newatm[1], newatm[2],
                                newatmindx, newatm[3], newatm[4])
                          )
    
    map_file_out.close()
    print " ... done writing to new map " + str(outmap)
        
def readbasefile(fname):
    fin = open(fname)
    
    base = []
    for line in fin:
        if line[0] == "#": continue
        line = line.strip().split()
        line = line[:3]
        line = [float(x) for x in line]
        base.append(line)
        
    return base
        
def match(a,B):
    """
     match a to one of the elements in B
    """
    
    indx = 0
    for b in B:
        if a[0] == b[0] and a[1] == b[1] and a[2] == b[2]:
            return indx
        else:
            indx += 1
    print " Error %1.4f %1.4f %1.4f" % (a[0],a[1],a[2])
    return -1

def transUC(atms, oldbasis, newbasis):
    """
    change from old basis which is the input file basis to the new basis i.e. the lammps basis
    """
    
    newbasis = np.transpose(np.array(newbasis)) #form the V matix

    from numpy.linalg import inv
    inv_newbasis = inv(newbasis) # invert
    
    CBM = []
    for u in oldbasis:
        u = np.dot(inv_newbasis,np.transpose(u))
        CBM.append(u) # collect the vectors
    CBM = np.transpose(np.array(CBM)) # turn into matrix, transpose so that you get column vectors
    # CBM is the representation of the old basis vectors in the new vector space
    
    newatms = []
    for elem in atms:
        type = elem[3]
        elem = elem[:3]
        elem = np.transpose(np.array(elem))
        newatms.append(np.dot(CBM,elem).tolist())
        newatms[-1].append(type)  # restore type
    return newatms

if __name__ == "__main__":
    main()