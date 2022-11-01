import statistics
import math
import numpy as np



def getPconfs(f):


    conf_file = open(f, 'r')
    pconfs = conf_file.readlines()
    conf_file.close()
    
    pconfs = [x.strip() for x in pconfs]

    setup = pconfs[:3]
    pconfs = pconfs[3:]
    pconfs = [[float(y) for y in x.split(" ")] for x in pconfs]
    
    numP = len(pconfs)
    
    r = [] # positions of center of mass
    b = [] # base vectors
    n = [] # base normal vectors
    v = [] # velocities
    L = [] # angular velocities


    for i in range(numP):
        p = pconfs[i]
        r.append(np.array(p[:3]))
        b.append(np.array(p[3:6]))
        n.append(np.array(p[6:9]))
        v.append(np.array(p[9:12]))
        L.append(np.array(p[12:]))
        
    return r, b, n, v, L
    
    
def getTopInfo(f):
    
    # read top info
    top_file = open(f,'r')
    top = top_file.readlines()
    top_file.close()
    
    return top
    
def dist(r1,r2):
    return math.sqrt((r1[0]-r2[0])**2 + (r1[1]-r2[1])**2 + (r1[2]-r2[2])**2)
    
def getBlocs(r, b):
    bLocs = []
    for i in range(numP):
        bLocs.append(np.add(r[i],b[i]))
    return bLocs
    
def findH(bLocs, top):
    h_bonds = []
    for i in range(len(top)-1):
        strandID = top[i+1][0]
        center = bLocs[i]
        closest = 99999
        closestIndex = i
        for z in range(len(top)-1):
            if strandID != top[z+1][0]:
                dist2center = dist(center,bLocs[z])
                if dist2center < closest:
                    closest = dist2center
                    closestIndex = z
        h_bonds.append((i,closestIndex))
    return h_bonds
    
    
def getHbonds(confFile, topFile):
    r, b, n, v, L = getPconfs(confFile)
    top = getTopInfo(topFile)
    bLocs = getBlocs(r,b)
    return findH(bLocs,top)
    
    
    