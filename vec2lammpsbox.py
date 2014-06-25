import sys
import numpy as np
def main():
    # Parse command line arguments
    #   arg[1] = input file name
    #   arg[2] = output file name
    #   -n nx,ny,nz  = extra unit cells
    #   -lammps      = make out file into lammps format
    #   -map         = create an map file map.outname (PREREQUISATE: coordinates of atoms)
    #   -x t1,t2 xt2 = create an alloy from type t1 to type t2 with concentration xt2 
    
    n = [1]*3
    douc = False  # flag input file contains atoms
    lammps = False # flag output file in lammps format
    domap = False # flag build map file
    doalloy = False # alloy flag
    flname = sys.argv[1]
    print "... number of args = " + str(len(sys.argv))
    outname = sys.argv[2]
    
    args = sys.argv[3:]
    nargs = len(args)
    if nargs>0:
        for iarg in range(nargs):
            if args[iarg] == "-n":
                n = args[iarg+1].strip().split(",")
                n = [int(x) for x in n]
                continue
            
            if args[iarg] == "-lammps":
                lammps = True
                continue
            
            if args[iarg] == "-map":
                domap = True
                continue
            
            if args[iarg] == "-x":
                doalloy = True
                base = float(args[iarg+1].split(",")[0])
                alloy = float(args[iarg+1].split(",")[1])
                conc = float(args[iarg+2])
                continue
    
#### load text from input file, ignore lines with 
    ABC = np.loadtxt(flname)
    # unit cell vectors will be the first three lines 
    #   take only the first three entries 
    #   (trailing zeros for book-keeping)
    A = ABC[0][:3]
    B = ABC[1][:3]
    C = ABC[2][:3]
    # build the lammps basis a,b,c using the relations 
    #   specified on the website
    a = [norm(A),0.0,0.0]
    Ahat = unit(A)
    b = [dot(B,Ahat), norm(np.cross(Ahat,B)), 0]
    ABhat = unit(np.cross(A,B))
    c = [dot(C,Ahat), dot(C,np.cross(ABhat,Ahat)), 0.]
    c[2] = np.sqrt(norm(C)**2-c[0]**2-c[1]**2)
    
#### If there are more lines these will be the atoms in the unit cell
    if len(ABC) > 3:
        nuc = len(ABC[3:])
        ntypes = 1  # at least one atom type
        ## find the number of atom types by taking the maximum
        for elem in ABC[3:]:
            if elem[3] > ntypes: ntypes = elem[3] # number of atom types
        print "... found atoms (" + str(nuc) + ") in text file, changing basis"
        douc = True
    # if found count number of atoms in unit cell 
    #    and transform. Build more if needed
    # ucatms = [ [x1,y1,z1,t1], [x2,y2,z2,t2], .... , [xn,yn,zn,tn]]
    if douc:
        ucatms = ABC[3:]
        atoms = transUC(ucatms,[A,B,C], [a,b,c])
        if n[0]>1 or n[1]>1 or n[2]>1:
            print "... need more atoms adding " + str(n) 
            atoms = buildmore(n,[a,b,c],atoms)
    print "     >>> " + str(len(atoms)) + " atoms"
    
#### Alloy if requested
    if doalloy:
        print "... alloying type __" + str(base) + "__ to type __" + str(alloy) + "__"
        Nalloy = conc*len(atoms)
        atoms,changeList = changetype(atoms, base,alloy, Nalloy)
        print "     >>> changed atoms percentage: %.4f at." % (float(len(changeList))/len(atoms))
        if alloy > ntypes: ntypes +=1 # added a new type to the crystal
    
#### write to file in either lammps format or in same format as input text
    fout = open(outname,'w')
    
    if lammps:
        fout.write("Cell with dimensions " + str(n))
        if doalloy: fout.write(" changed " + str(base) + " to " + str(alloy) +  
                               " with conc. " + str(float(len(changeList))/len(atoms)) + " at.")
        fout.write("\n\n")
        fout.write(str(len(atoms)) + " atoms\n%1.0f atom types\n\n"%(ntypes))
    
    fout.write("%10.8f %10.8f xlo xhi" % (0.0, n[0]*a[0]))
    fout.write("\n%10.8f %10.8f ylo yhi" % (0.0, n[1]*b[1]))
    fout.write("\n%10.8f %10.8f zlo zhi" % (0.0, n[2]*c[2]))
    if b[0] > 0.0 and c[0] > 0.0 and c[1] > 0.0: 
        fout.write("\n%10.8f %10.8f %10.8f xy xz yz" % (n[1]*b[0], n[2]*c[0], n[2]*c[1]))
    if douc:
        ia = 0
        fout.write("\n\nAtoms\n")
        for elem in atoms:
            ia +=1
            fout.write("\n%1.0f %1.0f %10.8f %10.8f %10.8f" % (ia, elem[3], elem[0], elem[1], elem[2]))
    fout.close()
    print "... Done writing to file " + outname
#### Done writing to file
#### Write to map file
    if douc and domap:
        outname = "map." + outname
        fout = open(outname, 'w')
        
        fout.write("%1.0f %1.0f %1.0f %1.0f" % (n[0], n[1], n[2], nuc)) #header
        fout.write("\n#l1 l2 l3 k type atom_id")
        
        inv_abc = np.linalg.inv(np.transpose(np.array([a,b,c])))
        ia = 0
        for elem in atoms:
            ia +=1
            type = elem[3]
            elem = elem[:3]
            elem = np.dot(inv_abc, np.transpose(np.array(elem))) # transform to fractional coordinates
            elem = [int(x) for x in elem]                        # take the integer part of the coordinate
            
            fout.write("\n%1.0f %1.0f %1.0f %1.0f %1.0f %1.0f" % 
                       (elem[0], elem[1], elem[2], (ia-1)%nuc, type, ia))
        print "... Done writing map file: " + outname
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

def buildmore(n,basis,uc):
    numofatoms = n[0]*n[1]*n[2]*len(uc)  #total number of atoms
    atoms = []

    for ix in range(n[0]):
        for seed in uc:
            atoms.append([seed[0]+ix, seed[1], seed[2], seed[3]])
    
    seedmax = len(atoms)
    for iy in range(1,n[1]):
        for ia in range(seedmax):
            atoms.append([atoms[ia][0], atoms[ia][1]+iy, atoms[ia][2], atoms[ia][3]])
            
    seedmax = len(atoms)
    for iz in range(1,n[2]):
        for ia in range(seedmax):
            atoms.append([atoms[ia][0], atoms[ia][1], atoms[ia][2]+iz, atoms[ia][3]])
    
    basis = np.array(basis)
    for ia in range(len(atoms)):
        temp = (atoms[ia][0]*basis[0] + 
                atoms[ia][1]*basis[1] + 
                atoms[ia][2]*basis[2]).tolist()
        temp.append(atoms[ia][3])    # add the type
        atoms[ia] = temp 
        
    return atoms
            

def changetype(atms,t1,t2, nt2):
    """
    change 'nt2' atoms from type 't1' to type 't2'
    """

    # find total number of atoms
    ntotal = len(atms)
    t1ids = []
    nt1 = 0;
    
    # find total number of t1 atoms and their ids
    for ia in range(ntotal):
        if atms[ia][3] == t1:
            t1ids.append(ia)
            nt1 +=1
    
    nt2 = int(nt2)
    if nt2 < 1:
        print "$$$ ERROR: not enough atoms of type ## " + str(t1) + " ## " + str(nt1) + " are available to alloy $$$"
        return

    np.random.shuffle(t1ids)
    rnd = t1ids[:nt2]
    for zombie in rnd:
        atms[zombie][3] = t2

    return (atms,rnd)
                
    
    
if __name__=="__main__":
    main()