
from numpy import *

import element_class
import pulseBeam


class Element_Aperture(element_class.Element):
    def __init__(self,new_spot):
        self.length = 0
        self.new_spot = new_spot
        self.name = 'Aperture' 
        self.n = 0

    def is_discrete(self):
        #discrete elements can only return calculations after the whole length
        return True
        
    def calc_pulseBeam(self,z,input_beam):
        #calculate pulseBeam after propagation through z 
        if(z < 0 or z > self.length):
            print 'Error with z - outside of length - z=',z
            raise Exception, 'Error with z - outside of length - z=%e'%z
        if(self.is_discrete() and z < self.length):
            return input_beam
        
        input_beam.beam_apply_fake_aperture(self.new_spot)
        
    def open_edit_dialog(self,refresh_callback):
        self.refresh_callback = refresh_callback
        
        #TODO
