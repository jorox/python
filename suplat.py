# This file uses code from vec2lammpsbox - consider merging
import numpy as np
import sys

def main():
    ABC1,n1,ABC2,n2,shft = readparam(sys.argv[1]).strip()

    abc1 = doABC(ABC1)
    abc2 = doABC(ABC2)
    noAtoms1 = len(ABC1[3:]); noAtoms2 = len(ABC2[3:])
    super,nAB = buildInPlane(ABC1, ABC2)

    print "... superlattice vectors are: %1.5f %1.5f" % (super[0], super[1]])
    print "    --> %1.3f, %1.3f from A ¬¬ %1.3f, 1.3f from B" (nAB[0][0], nAB[1][0], nAB[0][1], nAB[1][1])

    apprv = False
    while not apprv:
        nb = raw_input(" How many parent super-cells are needed: x,y")
        np = [int(x.strip()) for x in np.split(",")]
        N = countTotalNumberOfAtoms(nAB)
        lx = nb[0]*norm(super[0]); ly = nb[1]*norm(super[1])
        print "   --> system size = %1.4f -by- %1.4f" % (lx, ly)
        print "   --> number of atoms = %1.0f" %

    atoms1,ntypes1 = build(ABC1,abc1,n1)
    atoms2,ntypes2 = build(ABC2,abc2,n2)


def buildInPlance(A,B):
    super = []
### Find the larger one in x
    if norm(A[0]) < norm(B[0]): super.append(B[0])
    else: super.append(A[0])

    if norm(A[1]) < norm(B[1]): super.append(B[1])
    else: super.append(A[1])

    n = [ [norm(super[0])/norm(A[0]), norm(super[0])/norm(B[0])],
          [norm(super[1])/norm(A[1]), norm(super[1])/norm(B[1])] ]
    return (super,n)



def readparam(fname):
#### Input
    fin = open(fname)
    for line in fin:
        line= line.strip().strip()
        if line[0] == "#":
            continue

        if line[0] == "A": A = np.loadtxt(line[2].strip()); continue
        if line[0] == "B": B = np.loadtxt(line[2].strip()); continue

        if line[0] == "nA":
            nA = [int(x.strip()) for x in line[2].split(",")]
            continue
        if line[0] == "nB":
            nB = [int(x.strip()) for x in line[2].split(",")]
            continue

        if line[0] == "shift":
            shft = [float(x.strip()) for x in line[2].strip().split(",")]
            continue

    return (A,nA,B,nB,shft)

def doABC(ABC):

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

    return [a,b,c]

def build(ABC,abc,n):

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
        atoms = transUC(ucatms,ABC[:3], abc])
        if n[0]>1 or n[1]>1 or n[2]>1:
            print "... need more atoms adding " + str(n)
            atoms = buildmore(n,abc,atoms)
    print "     >>> " + str(len(atoms)) + " atoms"

    return (atoms,ntypes)


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

def norm(a):
    res = 0.0
    res += a[0]**2 + a[1]**2 + a[2]**2
    res = res**0.5
    return res
