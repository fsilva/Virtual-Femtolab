
from numpy import *

import pulseBeam

class Element_ThinLens:
    def __init__(self,f,lambdaZero):
        self.f = f
        self.name = 'Thin Lens'
        self.n = 0 #0 means this element calculates the refraction
        self.length = 0
        self.lambdaZero = lambdaZero

        
    def is_discrete(self):
        #discrete elements can only return calculations after the whole length
        return True
        
    def calc_pulseBeam(self,z,input_beam):
        #calculate pulseBeam after propagation through z 
        if(z < 0 or z > self.length):
            print 'Error with z - outside of length'
            raise
        if(self.is_discrete() and z < self.length):
            return input_beam
            
        input_beam.beam_apply_thinlens(self.f)
            
        return input_beam
        
        
    def open_edit_dialog(self,refresh_callback):
        self.refresh_callback = refresh_callback
        
        #TODO
        

        
