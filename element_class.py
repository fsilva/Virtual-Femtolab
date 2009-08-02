
#Master Element Class, return unaltered electric field/spectra
#                      calculate ABCD propagation

from numpy import *

import pulseBeam

class Element:
    def __init__(self,length):
        self.length = length
        self.name = 'Dummy'
        self.n = 1
        pass
        
    def is_discrete(self):
        #discrete elements can only return calculations after the whole length
        return False 
        
    def propagate(self,z,input_beam,lambdaZero):
        return input_beam

        
    def calc_pulseBeam(self,z,input_beam):
        #calculate pulseBeam after propagation through z 
        if(z < 0 or z > self.length):
            print 'Error with z - outside of length'
            raise
        if(self.is_discrete() and z < self.length):
            return input_beam
        
        return input_beam
