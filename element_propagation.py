
from numpy import *

import element_class
import pulseBeam


# name, GVD,TOD,QOD
materials = [ ['Glass - Fused Silica',1.45332,36160e-30,27495e-45,0],
              ['Glass - BK7',1.51078,44650e-30,32099e-45,0],
              ['Air - 1 bar',1.0,20.4e-30,10.896e-45,0],
              ['Vacuum',1.0,0,0,0]]

class Element_Propagation(element_class.Element):
    def __init__(self,length,material,lambdaZero):
        self.length = length
        self.material = material[:]
        self.lambdaZero = lambdaZero
        self.name = material[0]
        self.n = material[1]
        
    def is_discrete(self):
        #discrete elements can only return calculations after the whole length
        return False 
        
    def calc_pulseBeam(self,z,input_beam):
        #calculate pulseBeam after propagation through z 
        if(z < 0 or (abs(z-self.length) > 1e-10 and z > self.length)): #i.e. z > self.length
            print 'Error with z - outside of length - z=',z
            raise Exception, 'Error with z - outside of length - z=%e'%z
        if(self.is_discrete() and z < self.length):
            return input_beam
        
        #1) material propagation is in the spectral domain
        input_beam.field_to_spectrum()
        
       
        #2) apply propagation
        input_beam.apply_dispersion(z*self.material[2],z*self.material[3],z*self.material[4])
                
        #3) recalculate time domain
        input_beam.spectrum_to_field()

        