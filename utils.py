import statistics
import math
import numpy as np
import pandas as pd


# get clean conf fiel data
def getPconfs(f):

    ''' test mkdocs'''

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
    
    
    confs_dict = {}
    
    confs_dict['r'] = r
    confs_dict['b'] = b
    confs_dict['n'] = n
    confs_dict['v'] = v
    confs_dict['L'] = L
    
    confs = pd.DataFrame.from_dict(confs_dict, index = pd.RangeIndex(numP))
    
    return confs



 #returns a clean toplogy info data   
def getTopInfo(f):
    top_file = open(f,'r')
    top = top_file.readlines()
    top_file.close()
    
    top = [t.strip() for t in top]
    top = [t.split() for t in top]
    return top
    
# simple 3D distance function
def dist(r1,r2):
    return math.sqrt((r1[0]-r2[0])**2 + (r1[1]-r2[1])**2 + (r1[2]-r2[2])**2)
   
#gets the location of every base for every nucleotide   
def getBlocs(r, b):
    bLocs = []
    numP = len(r)
    for i in range(numP):
        bLocs.append(np.add(r[i],b[i]))
    return bLocs

#Wil break on nicks and go to neighbor.  P
#takes in the location of each base and the topology info of the model
def findH(bLocs, top):
    h_bonds = []
    
    # Iterate through all particles to attempt to find their H bond neighbor
    for i in range(1,len(top)):
        strandID = top[i][0]
        center = bLocs[i-1]
        closest = 69.0
        closestIndex = -1
        for z in range(1,len(top)):
            if strandID != top[z][0]:
                dist2center = dist(center,bLocs[z-1])
                if dist2center < closest and dist2center > 0.500:
                    closest = dist2center
                    closestIndex = z-1
        if closest > 0.85:
            closestIndex = -1
        h_bonds.append((i-1,closestIndex))
    return h_bonds
    

# wrapper function for convenient file input
def getHbonds(confFile, topFile):
    r, b, n, v, L = getPconfs(confFile)
    top = getTopInfo(topFile)
    bLocs = getBlocs(r,b)
    return findH(bLocs,top)
    
    
    