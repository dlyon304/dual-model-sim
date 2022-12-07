import utils
import pandas as pd
import math
import numpy as np



class Model:
    
    #numP: number of particles in the model
    
    #pconfs_header: simulation parameters for model
    #pconfs: pandas dataframe of all properties of the particles
    
    #top_header: number of strands, number of particles
    #top_df: pandas dataframe for all particles interactions with eachother
    
    #bLocs: the location of each base of each particle
    #h_bonds: the expected hydrogen bonding partner of each particles, -1 if single stranded
    
    
    
    def __init__(self, conf_file, top_file):
        
        # Retrieve the particle configurations and make pandas dataframe containing particles attributes
        #************************************
        temp_file = open(conf_file, 'r')
        confs = temp_file.readlines()
        temp_file.close()
        
        confs = [x.strip() for x in confs]
        self.pconfs_header = confs[:3]
        confs = [[float(y) for y in x.split(" ")] for x in confs]
        confs = confs[3:]
        
        self.numP = len(confs)
        
        r = [] # positions of center of mass
        b = [] # base vectors
        n = [] # base normal vectors
        v = [] # velocities
        L = [] # angular velocities
        
        for i in range(self.numP):
            p = confs[i]
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
        
        self.pconfs = pd.DataFrame(confs_dict, index = pd.RangeIndex(self.numP))
        
        
        #Retrieve top info and create pandas dataframe to store particle interactions
        #***************************
        temp_file = open(top_file, 'r')
        top = temp_file.readlines()
        temp_file.close()
        
        top = [t.strip() for t in top]
        top = [t.split() for t in top]
        self.top_header = top[0]
        top = top[1:]
        top = [[int(line[0]),line[1],int(line[2]),int(line[3])] for line in top]
        
        bl = []
        for i in range(self.numP):
            bl.append(np.add(self.pconfs['r'].iloc[i]), self.pconfs['b'].iloc[i])
        self.bLocs = bl
        
        
        #Find all hydrogen bonds, -1 index if single stranded
        hb = []
        for i in range(self.numP-1):
            strandID = top[i][0]
            center = self.bLocs[i]
            closest = 69.0
            closestIndex = -1
            for z in range(self.numP):
                if strandID != top[z][0]:
                    dist2center = math.sqrt((center[0]-self.bLocs[z][0])**2 + (center[1]-self.bLocs[z][1])**2 + (center[2]-self.bLocs[2][z])**2)
                    if dist2center < closest and dist2center > 0.500:
                        closest = dist2center
                        closestIndex = z
            if closest > 0.85:
                closestIndex = -1
            hb.append((i-1,closestIndex))
        self.h_bonds = hb
        
        df = df = pd.DataFrame(top, 
                      columns = ['strand', 'base', 'downStream', 'upStream'],
                      index = range(self.numP)
                      )
        df['h_bond'] = [b[1] for b in self.h_bonds]
        
        self.top_df = df
        
        
        
    





