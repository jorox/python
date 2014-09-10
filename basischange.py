import numpy as np
import sys

def main():
    """
    changes the unit cell of a certain crystal
    takes a special file format
    """
    
    fname = sys.argv[1]
    fin = open(fname)
    a123 = []
    batms = []
##### Read in old basis and vectors
    for line in fin:
        if line[0] == "#": continue
        
        line = line.split()
        line = [ float(x.strip()) for x in line[:3] ]
        
        if len(a123) == 3: batms.append(line); continue
        a123.append(line)
    
    fname = sys.argv[2]
    fin = open(fname)
    b123 = []
    for line in fin:
        if line[0] == "#": continue
        
        line = line.split()
        line = [ float(x.strip()) for x in line[:3] ]
        
        b123.append(line)
        if len(b123) == 3: break
        
    print "... lattice vectors \n          old                     new "
    for i in range(3):
        print ("     %1.4f | %1.4f | %1.4f     %1.4f | %1.4f | %1.4f " % 
               (a123[0][i], a123[1][i], a123[2][i], b123[0][i], b123[1][i], b123[2][i]) )
    
    print "... basis atoms = "
    for i in range(len(batms)):
        print "      %1.4f %1.4f %1.4f" % (batms[i][0], batms[i][1], batms[i][2])
        

##### Read in new basis that you want to switch to
##### Take any point q_A = (q1,q2,q3) then q_E = q1*a1_E + q2*a2_E + q3*a3_E = (x, y ,z)
#####     Hence, we can say that q_B = (p1,p2,p3) = q1*a1_B + q2*a2_B + q3*a3_B

#####     Writing in matrix form we can say that [ a1_E | a2_E | a3_E ]*q_A = q_E
#####     Apply the same logic to vector ai we see  [b1_E | b2_E | b3_E ]*ai_B = ai_E
#####     Hence, --> ////ai_B = cbm*ae_E|\\\\
            
    a123 = [np.array(x) for x in a123] #old basis
    b123 = [np.array(x) for x in b123] #new basis
    
    B = np.transpose(b123)
    
    invB = np.linalg.inv(B)
    a123_B = [np.dot(x,invB) for x in a123]
    A_B = np.transpose(a123_B)  #representation of old vectors in new space (colum  wise)
    print " ... representation of old vectors in the new basis = "
    for i in range(3):
        print "      %1.7f | %1.7f | %1.7f" % (A_B[0][i], A_B[1][i], A_B[2][i])
    
##### Build 5 unit cells all around
    comb = []  #array containing unit cell coordinates
    for i1 in range(-2,2):
        for i2 in range(-2,2):
            for i3 in range(-2,2):
                comb.append([i1,i2,i3])
                
    nuc = len(comb)
    b2atms = []     #new basis atoms 
    b2map = []      #new basis map

    for uc in comb:
        for i in range(len(batms)):
            tmp = [ uc[0]+batms[i][0], uc[1]+batms[i][1], uc[2]+batms[i][2] ] # add all basis atoms in each unit cell
            prcs = 4   # significat figures for rounding
            tmp = np.array(tmp)
            tmp = np.dot(A_B,tmp) # matrix multiplication
            tmp = np.round(tmp,prcs) 
            eps = 0 #needed for round off error
            if  -eps<=tmp[0]<1+eps and -eps<=tmp[1]<1+eps and -eps<=tmp[2]<1+eps: # if in first unit cell
                b2atms.append(tmp.tolist())
                b2map.append( [uc[0],uc[1],uc[2],i] )  
    
    print "--> New basis has " + str(len(b2atms)) + " atoms in fractional coordinates:"
    for i in range(len(b2atms)):
        print ( "     %1.4f %1.4f %1.4f <-- %1.0f %1.0f %1.0f|%1.0f" % 
        (b2atms[i][0], b2atms[i][1], b2atms[i][2], b2map[i][0], b2map[i][1], b2map[i][2], b2map[i][3]) )    
    

if __name__ == "__main__":
    main()
        
        
        
        
        
        