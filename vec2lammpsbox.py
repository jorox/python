import sys
import numpy as np
def main():
    # Parse command line arguments
    n = 1
    douc = False  # flag input file contains atoms
    lammps = False # flag output file in lammps format
    domap = False # flag build map file
    flname = sys.argv[1]
    print "... number of args = " + str(len(sys.argv))
    if len(sys.argv)>2:
        if len(sys.argv[2]) == 1:
            n = int(sys.argv[2])*3
        else:
            n = sys.argv[2].strip().split(",")
            n = [int(x) for x in n]
    if len(sys.argv)>3:
        if sys.argv[3].strip() == "-lammps":
            lammps = True
    if len(sys.argv)>4:
        if sys.argv[4] == "-map":
            domap = True
    
    # load text from input file, ignore lines with # 
    ABC = np.loadtxt(flname)
    A = ABC[0]
    B = ABC[1]
    C = ABC[2]
    # build the lammps basis a,b,c
    a = [norm(A),0.0,0.0]
    Ahat = unit(A)
    b = [dot(B,Ahat), norm(np.cross(Ahat,B)), 0]
    ABhat = unit(np.cross(A,B))
    c = [dot(C,Ahat), dot(C,np.cross(ABhat,Ahat)), 0.]
    c[2] = np.sqrt(norm(C)**2-c[0]**2-c[1]**2)
    # Test for atoms in text file
    if len(ABC) > 3:
        nuc = len(ABC[3:])
        print "... found atoms (" + str(nuc) + ") in text file, changing basis"
        douc = True
    # if found count number of atoms in unit cell and transform. Build more if needed
    if douc:
        ucatms = ABC[3:]
        atoms = transUC(ucatms,[A,B,C], [a,b,c])
        if n[0]>1 or n[1]>1 or n[2]>1:
            print "... need more atoms adding " + str(n) 
            atoms = buildmore(n,[a,b,c],atoms)
    print "... " + str(len(atoms)) + " atoms"
    flname = flname + ".lammps"
    # write to file in either lammps format or in same format as input text
    fout = open(flname,'w')
    
    if lammps:
        fout.write("Cell with dimensions " + str(n) + "\n\n")
        fout.write(str(len(atoms)) + " atoms\n1 atom type\n\n")
    
    fout.write("%10.8f %10.8f xlo xhi\n" % (0.0, n[0]*a[0]))
    fout.write("%10.8f %10.8f ylo yhi\n" % (0.0, n[1]*b[1]))
    fout.write("%10.8f %10.8f zlo zhi\n" % (0.0, n[2]*c[2]))
    fout.write("%10.8f %10.8f %10.8f xy xz yz" % (n[1]*b[0], n[2]*c[0], n[2]*c[1]))
    if douc:
        ia = 0
        fout.write("\n\nAtoms\n")
        for elem in atoms:
            ia +=1
            fout.write("\n%1.0f 1 %10.8f %10.8f %10.8f" % (ia, elem[0], elem[1], elem[2]))
    fout.close()
    print "... Done writing to file " + flname
    # write to map file if needed
    if domap:
        flname = "map." + flname
        fout = open(flname, 'w')
        
        fout.write("%1.0f %1.0f %1.0f %1.0f" % (n[0], n[1], n[2], nuc)) #header
        fout.write("\n#l1 l2 l3 k type atom_id")
        
        inv_abc = np.linalg.inv(np.transpose(np.array([a,b,c])))
        ia = 0
        for elem in atoms:
            ia +=1
            elem = np.dot(inv_abc, np.transpose(np.array(elem)))
            elem = [int(x) for x in elem]
            fout.write("\n%1.0f %1.0f %1.0f %1.0f 1 %1.0f" % 
                       (elem[0], elem[1], elem[2], (ia-1)%nuc, ia))
        print "... Done writing map file: " + flname
        fout.close()
        
def dot(a,b):
    res = 0
    for x,y in zip(a,b):
        res += x*y
    return res

def norm(a):
    res = 0;
    for x in a:
        res += x**2
    res = res**0.5
    return res

def unit(a):
    return a/norm(a)
  
def transUC(atms, oldbasis, newbasis):
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
        elem = np.transpose(np.array(elem))
        newatms.append(np.dot(CBM,elem))
    return newatms

def buildmore(n,basis,uc):
    numofatoms = n[0]*n[1]*n[2]*len(uc)  #total number of atoms
    atoms = []
    map = []
    for ix in range(n[0]):
        for seed in uc:
            atoms.append([seed[0]+ix, seed[1], seed[2]])
    
    seedmax = len(atoms)
    for iy in range(1,n[1]):
        for ia in range(seedmax):
            atoms.append([atoms[ia][0], atoms[ia][1]+iy, atoms[ia][2]])
            
    seedmax = len(atoms)
    for iz in range(1,n[2]):
        for ia in range(seedmax):
            atoms.append([atoms[ia][0], atoms[ia][1], atoms[ia][2]+iz])
    
    basis = np.array(basis)
    for ia in range(len(atoms)):
        temp = (atoms[ia][0]*basis[0] + 
                atoms[ia][1]*basis[1] + 
                atoms[ia][2]*basis[2]).tolist()
        atoms[ia] = temp 
        
    return atoms
            

if __name__=="__main__":
    main()