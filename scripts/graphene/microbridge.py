#!/usr/bin/python
import numpy as np

# Method to make geometric calculations for semiconductor microbridges
#
# A temperature dependent measurement of the carrier velocity vs. electric field 
# characteristic for as-grown and H-intercalated epitaxial graphene on SiC
#
# Journal of Applied Physics 113, 193708 (2013); https://doi.org/10.1063/1.4807162

class microbridge_geometry: 

    def __init__(self, DIM): 
        self.l = 28.0
        self.w = 20.9
        self.n = 500
        self.DIM = DIM
        self.bridge(DIM)
        self.form_factor()

    # Build dimensions for the bridge
    def bridge(self, d):
        x,y   = np.linspace(0.0, self.l, self.n), []
        DELTA = x[1]-x[0] 
    
        for i, _ in enumerate(x[0:self.n/2]): 
            if _ <= 3.80:  y.append( self.w/2 )
            else:  
                tmp = y[i-1] - DELTA
                if tmp > d/2.0: y.append(tmp)
                else: y.append(d/2.0)
      
        self.X = x
        self.Y = y+y[::-1]

    def form_factor(self):
        DELTA = self.X[1]-self.X[0] 
        self.W = [2.0*_ for _ in self.Y]
        self.F = 0.0
        for _ in self.W: self.F += DELTA/_
        
    def transform_iv(self, VD, ID, n=1.0e13, RC=200, W0=20.9):
        e,s1,s2 = 1.602e-19, 10, 1e-4 
        eD = [s1*(V-I*RC/W0)/(self.F*self.DIM) for V,I in zip(VD,ID)] 
        vD = [I/(2*e*n*self.DIM*s2) for I in ID]
        return eD, vD
