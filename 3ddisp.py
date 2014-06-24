import numpy as np
import sys
import os

class vdata:

    def __init__(self,fname,mname,maxsteps=0):
	print "... reading " + str(maxsteps) + " timesteps"
        self.readmap(mname)
        self.readfile(fname,maxsteps)

    def readmap(self,mname):
        m_in = open(mname,'r')
        line = m_in.readline().strip()
        line = line.split()
        line = [int(x) for x in line]

        self.n = line[:3]      #3D system
        self.atomspercell = line[3]
        self.Natoms = np.prod(line)
        Nlines = self.Natoms

        line = m_in.readline() # usually this is a comment line

        b = self.atomspercell
        n = self.n
##### initialize n-by-m-by-pxb array filled with -1
        amap = np.ones((n[0],n[1],n[2],b),'int')*-1
##### initialize a list of types for atoms
        atyp = np.zeros(self.Natoms)
        cntr = 0  # loop counter
        for line in m_in:
            cntr +=1
            if cntr>Nlines: print "too many lines in map file"
            line = line.strip().split()
            line = [int(x) for x in line]
##### format for map file: l0 l1 l2 index type id
            nx = line[0]; ny = line[1]; nz = line[2]
            b = line[3]
            typ = line[4]
            aid = line[5]
##### store the (python) id of the atom
            amap[nx,ny,nz,b] = aid-1   #lammps id and python id differ by 1, check how atm_indx is given in readfile()
##### store the type of atom
            atyp[aid-1] = typ
        self.amap = amap
        self.atype = atyp
        return 0

    def getAtomData(self,nxyz,alpha,b,t=-1):
        '''
        returns the data corresponding to dimension 'alpha' for the atom in
        unit cell 'nxyz[0],nxyz[1],nxyz[2]' with index 'b' for time 't'
        if time t=-1 return the entire temporal data

        return type (np.array)
        '''

        aid = self.amap[nxyz[0],nxyz[1],nxyz[2],b]
        if t>-1:
            return self.data[aid][alpha][t]
        else:
            return np.array(self.data[aid][alpha])

    def getrmap(self,navg=-1):
        '''
        return the position of each unit cell from the average position of
        the first atom in the unit cell.
        navg: number of samples when averaging

        returns : np.array
        '''
        n = self.n
        rmap = np.zeros((n[0],n[1],n[2],3))
        for nx in range(n[0]):
            for ny in range(n[1]):
                for nz in range(n[2]):
                    for ix in range(3):
                        temp = self.getAtomData([nx,ny,nz],ix,0)
                        temp = np.median(temp[:navg]) #use median because average can get into trouble with atoms crossing periodic boundaries
                        rmap[nx,ny,nz,ix] = temp
        return rmap


    def readfile(self,flname, maxSteps=0):
        '''
        parse input file for data
        format of input file:
        {   id type   vx vy vz   x y z   } (gratutous whitespaces)
        note that 'id' is critical but 'mass' is not
        format of data array:
        N atoms
        M timesteps
           data = [ atom0{[x0,x1,...,xM],[y0,...,yM],[z],[u],[v],[w]}
                    atom1{[x], [y], [z], [u], [v], [w]} ,
    		        ...
    		        ...
    		        ...,
    		        atomN{...} ]

        data[aid][alpha][t]
           aid    : atom id
           alpha  : dimension
           t      : timestep
        '''

        # Load data from file
        f_in = open(flname,'r')
        mrkr = 0
        yummy = ["|","/","-","\\"]
        for i in range(0,3):
            f_in.readline()
        Natoms = int(f_in.readline())

        data = []  #preallocate so that atom index can be used to reference

        for i in range(Natoms):
            data.append([[],[],[],[],[],[]])

        f_in.close()
        f_in = open(flname)
        Nsteps = 0
        lout = " "

        for line in f_in:
            if maxSteps and Nsteps == maxSteps:
                break
            if mrkr == 1:
                onstep = line.strip()
            if mrkr > 8:  #8 header lines per section
                if mrkr == 9+Natoms:
                    lout = "\b"*(len(lout)) + yummy[Nsteps%4] + ' reading file ... ' + onstep + " (" + "%1.2f"%((float(Nsteps)/maxSteps)*100) + "%)"
                    print lout,
                    Nsteps += 1
                    mrkr = 1
                    continue

                # get the id
                tmp = line.strip().split(" ")
                tmp = [float(x) for x in tmp]

                atm_indx = int(tmp[0]) - 1
                data[atm_indx][3].append(tmp[2]) #vx
                data[atm_indx][4].append(tmp[3]) #vx
                data[atm_indx][5].append(tmp[4]) #vx
                data[atm_indx][0].append(tmp[5]) #x
                data[atm_indx][1].append(tmp[6]) #y
                data[atm_indx][2].append(tmp[7]) #z

            mrkr += 1
        if not maxSteps: Nsteps +=1  # needed when reading all steps in file
        f_in.close()
        print("\nDone loading data ... "+str(Nsteps)+" timesteps, "+str(Natoms)+" atoms")
        self.data = data
        return 0

def setupkvecs(kvfname,avctr=None):
    fin = open(kvfname)
    kvecs = []  #initialize kvec array
    while True:
        line = fin.readline()
        if not line: break
        v1 = line.strip().split()
        v2 = fin.readline().strip().split()
        N =  fin.readline().strip()
        N = int(N)

        v1 = [float(x) for x in v1[:3]]
        v2 = [float(x) for x in v2[:3]]
        v1 = np.array(v1)
        v2 = np.array(v2)
        v = v2-v1

        for k in range(N+1):
            temp = v1 + v*float(k)/N
            if (len(kvecs)>1 and np.all(temp==kvecs[-1])):
                continue
            kvecs.append(temp)
    kvecs = np.array(kvecs)
    #build reciprocal vectors and transform kvecs to x,y,z system
    if avctr is not None:
        bvcs = buildrecp(avctr)

        for ik in range(len(kvecs)):
            kvecs[ik] = ( kvecs[ik,0]*bvcs[0] +
                          kvecs[ik,1]*bvcs[1] +
                          kvecs[ik,2]*bvcs[2] )
    return kvecs

def buildrecp(a):
    """
    gives the reciprocal lattice vectors defined by a
    a is a 3x3 array representing the lattice vectors
    """

    a = np.array(a)
    a12 = np.cross(a[1],a[2])
    a01 = np.cross(a[0],a[1])
    a20 = np.cross(a[2],a[0])

    V = np.sum(a[0]*a12)
    b = 2*np.pi*np.array([a12,a20,a01])/V

    return b

def parsepfile(fname):
    fin = open(fname)
    line = fin.readline()

    while line:
        line = line.strip().split()

        if line[0] == "ddir":
            ddir = line[2]

        if line[0] == "flname":
            vflname = line[2]
            outname = vflname+".phi"

        if line[0] == "mname":
            mname = line[2]


        if line[0] == "kname":
            kname = line[2]


        if line [0] == "dims":
            dims = line[2].split(",")
            dims = [int(x.strip()) for x in dims]

        if line[0] == "a":
            i = 0
            a = []
            while i<3:
                i +=1
                line = fin.readline().strip()
                at = line.split(",")
                at = [float(x.strip()) for x in at]
                a.append(at)


        if line[0] == "dt":
            dt = float(line[2].strip())

        if line[0] == "mtyp":
            mtyp = line[2].split(",")
            mtyp = [float(x.strip()) for x in mtyp]
        
        if line[0] == "maxsteps":
            maxstep = int(line[2])
        
        if line[0] == "outname":
            outname = line[2]

        line = fin.readline()

    return (ddir+vflname, ddir+mname, ddir+kname, dims, a, dt, mtyp, maxstep, outname)


def main():
    pfile = sys.argv[1]
    os.getcwd()
    flname, mname, kname, dims, a, dt, mtyp, maxsteps, outname = parsepfile(pfile)
    print "   --> file = " + flname
    print "   --> out = " + outname
    print "   --> map = " + mname
    print "   --> lattice = " + str(a)
    print "   --> dt = " + str(dt)
    print "   --> masses = " + str(mtyp)

    #ddir = "/Users/kassemw/Documents/LiClipse_Workspace/dispersion/"
    #flname = ddir+"xv_Cu.dat"
    #mname  = ddir+"map.in"
    #kname = ddir+"in.disp"
    #dims = [3,4,5]   #vx,vy,vz
    #a = [[2.55618,   0,        0],
    #     [1.2780875, 2.213726, 0],
    #     [1.2780875, 0.737909, 2.08712122]]       #unit cell length
    #a = [[3.1118, 0.0, 0.0],
    #     [0.0 , 5.3898 ,0.0],
    #     [0.0 , 0.0 ,5.0737899]]
    #dt = 0.03  #timestep ps
    #mtyp = [63.546]
    #mtyp = [26.9815386,14.007,15.9994] #masses of different atoms


    ###### read file into memory ######
    krnl = vdata(flname,mname,maxsteps)

    ###### calculated from input ######
    N_T = np.prod(krnl.n)    #number of unit cells
    N = len(krnl.data[0][0])    #number of time steps
    atomspercell = krnl.atomspercell       #number of atoms per unit cell
    r = krnl.getrmap(500) #coordinates of each unit cell
    cmap = buildmap(krnl.n)

    dw = 1./N/dt
    w_max = 1/dt/2  #maximum frequency
    tau = N*dt #total run time

    omega = np.arange(0,w_max,dw)
    kappa = setupkvecs(kname,a)
    time  = np.arange(0,N)*dt
    nw = len(omega)
    steps = float(len(kappa))/10

    phi = np.zeros((len(kappa),len(omega)))
    ik = 0
    print "kvec progress = " + '[           ]',
    print '\b'*12,
    for k in kappa:
#        iw = 0
        sys.stdout.flush()
        phialphab = 0 # b x alpha matrix
        for alpha in dims:
            temp2 = 0

            for b in range(atomspercell):
                temp = 1j*np.zeros(len(time))
                for iN_T in range(N_T):
                    nxyz = cmap[iN_T]
                    aid = krnl.amap[nxyz[0],nxyz[1],nxyz[2],b]
                    if aid < 0: continue
                    kdotr = ( k[0]*r[nxyz[0],nxyz[1],nxyz[2],0] +
                              k[1]*r[nxyz[0],nxyz[1],nxyz[2],1] +
                              k[2]*r[nxyz[0],nxyz[1],nxyz[2],2] )
                    temp += mtyp[int(krnl.atype[aid])-1]*krnl.getAtomData(nxyz,alpha,b)*np.exp(1j*kdotr)

                temp = np.fft.fft(temp)
                temp  = np.real(np.abs(temp)**2)
                temp2 += temp

            phialphab += temp2

        phi[ik] = phialphab[:nw]
        if ik%steps == 0:
            print '\b.',
            sys.stdout.flush()
        ik += 1
    print '\b]  Done!'
    phi *= 1.0/4/np.pi/N_T/tau
    write2file(phi,outname)

def write2file(phi,fname):
    fout = open(fname,'w')
    for iw in range(len(phi[0])-1,-1,-1):
        line = ""
        for ik in range(len(phi)):
            line += str(phi[ik,iw]) + " "
        fout.write(line + "\n")
    print "...done writing to file "+ fname
    fout.close()

def buildmap(n):
    """
    sequence of numbers to loop over unit cells in the data
    """
    comb = []
    for i1 in range(n[0]):
        for i2 in range(n[1]):
            for i3 in range(n[2]):
                comb.append([i1,i2,i3])
    return comb

if __name__=="__main__":
    main()


