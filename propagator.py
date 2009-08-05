
import pulseBeam
#import initialpulse_dialog
#import addelement_dialog

import element_propagation
import element_thinlens
import element_chirpedmirror
import element_filter

class Propagator:
    def __init__(self,NT,deltaT,lambdaZero):
        self.initialPulseBeam = pulseBeam.pulseBeam(NT,deltaT,lambdaZero)
        self.pulseBeam = pulseBeam.pulseBeam(NT,deltaT,lambdaZero)
        
        self.lambdaZero = lambdaZero
        
        self.elements = []
        self.current_z = 0
    
    def change_computational_window(self,NT,deltaT):
        self.initialPulseBeam.change_window(NT,deltaT)
        del self.pulseBeam
        self.pulseBeam = pulseBeam.pulseBeam(NT,deltaT,self.lambdaZero)
        self.change_z(self.current_z)
        
    def get_initialPulseBeam(self):
        return self.initialPulseBeam
        
    def get_pulseBeam(self):
        return self.pulseBeam
        
    def example_pulseBeam(self):
        self.initialPulseBeam.initialize_pulse(10e-15,0,0,0,1e-4,1e10,0,1e-3,1000)
        
        self.initialPulseBeam.calculate_autoco()
        self.initialPulseBeam.calculate_FROG()
        
        self.pulseBeam.copy(self.initialPulseBeam)
        
        

    
    def example_elements(self):
        self.elements.append(element_thinlens.Element_ThinLens(1e-3,self.lambdaZero))
        self.elements.append(element_propagation.Element_Propagation(0.001,element_propagation.materials[1],self.lambdaZero))
        self.elements.append(element_thinlens.Element_ThinLens(1,self.lambdaZero))
        self.elements.append(element_propagation.Element_Propagation(0.001,element_propagation.materials[1],800e-9))
        self.elements.append(element_thinlens.Element_ThinLens(1,800e-9))
        self.elements.append(element_propagation.Element_Propagation(0.001,element_propagation.materials[0],800e-9))
        #self.elements.append(element_propagation.Element_Propagation(0.002,element_propagation.materials[1],800e-9))
        #self.elements.append(element_propagation.Element_Propagation(0.001,element_propagation.materials[1],800e-9))
        #self.elements.append(element_propagation.Element_Propagation(0.002,element_propagation.materials[1],800e-9))
        #self.elements.append(element_propagation.Element_Propagation(0.005,element_propagation.materials[0],800e-9))
        #self.elements.append(element_chirpedmirror.Element_ChirpedMirror(element_chirpedmirror.mirrors[0],10,800e-9))
        #self.elements.append(element_propagation.Element_Propagation(0.005,element_propagation.materials[0],800e-9))
        
    def get_elements(self):
        return self.elements
        
    def get_materials(self):
        return element_propagation.materials
        
        
    def get_max_z(self):
        z0 = 0
        for element in self.elements:
            z0 += element.length
        return z0
        
    def change_z(self,z):
        max_z = self.get_max_z()
        if(z > max_z):
            z = max_z
            return
        
        self.pulseBeam.copy(self.initialPulseBeam)
        self.current_z = z
    
        z0 = 0
        i = 0
        lastn = 0
        while(z0 < z and i < len(self.elements)):
            if(z0+self.elements[i].length < z):
                #refraction at the first material interface
                if(self.elements[i].n != 0): #if n = 0, the element doesn't have refraction (e.g. a thin lens)
                    if(lastn != 0 and self.elements[i].n != 0 and self.elements[i].n != lastn):
                        self.pulseBeam.beam_apply_refraction(lastn,self.elements[i].n)    
                        lastn = self.elements[i].n
                
                self.elements[i].calc_pulseBeam(self.elements[i].length,self.pulseBeam)
                self.pulseBeam.beam_apply_propagation(self.elements[i].length,self.elements[i].n) #should this move into calc_pulseBeam?

                z0 += self.elements[i].length
                i += 1
            else:
               self.elements[i].calc_pulseBeam(z-z0,self.pulseBeam)
               self.pulseBeam.beam_apply_propagation(z-z0,self.elements[i].n)
               z0 = z
    
        self.pulseBeam.calculate_autoco()
        self.pulseBeam.calculate_FROG()

        
    def get_spots(self):
        spots = []
        
        spot      = self.initialPulseBeam.get_beam_spot()
        curvature = self.initialPulseBeam.get_beam_curvature()

        lastn = 0 #for refraction calc
        for element in self.elements:
            spots.append(spot)
            element_name = str(element.__class__)
            if(element_name == 'element_propagation.Element_Propagation'):
                
                #refraction at the first material interface
                if(element.n != 0 and lastn != 0 and element.n != 0 and elements.n != lastn):  #if n = 0, the element doesn't have refraction (e.g. a thin lens)
                    spot_tmp,curvature_tmp = self.initialPulseBeam.beam_apply_refraction(lastn,self.elements[i].n)    
                    lastn = self.elements[i].n
                    new_spot,new_curvature = self.initialPulseBeam.beam_calc_propagation(element.n, element.length,spot_tmp,curvature_tmp)
                else:        
                    new_spot,new_curvature = self.initialPulseBeam.beam_calc_propagation(element.n, element.length,spot,curvature)
                    
            elif(element_name == 'element_thinlens.Element_ThinLens'):
                
                new_spot,new_curvature = self.initialPulseBeam.beam_calc_thinlens(element.f,spot,curvature)                
            else:
                pass
            if(new_curvature > 0 and curvature < 0):
                spots.append(-new_spot)
            else:
                spots.append(new_spot)
            
            spot = new_spot
            curvature =  new_curvature
              
        return spots
        
        
    def add_element(self,element_type,position):
        if(element_type == 'Material Propagation'):
           element = element_propagation.Element_Propagation(0.001,element_propagation.materials[0],800e-9)
        elif(element_type == 'Thin Lens'):
           element = element_thinlens.Element_ThinLens(0.1,800e-9)
        elif(element_type == 'Chirped Mirror'):
           element = element_chirpedmirror.Element_ChirpedMirror(element_chirpedmirror.mirrors[0],22,800e-9)
        elif(element_type == 'Bandstop Spectral Filter'):
           element = element_filter.Element_Filter(780e-9,820e-9)
        else:
            return
            
        self.elements.insert(position,element)
        
#        element.open_edit_dialog(self.refresh_callback)
        
#        self.refresh_callback() 
        
    def remove_element(self,selected):
        if(selected >= len(self.elements)):
            print 'remove_element: invalid selected data'
            return 
        
        del self.elements[selected] 
        
    def delete_frogs_for_output(self):
        #delete frog data in pulsebeams so that each output file is much smaller
        self.initialPulseBeam.delete_frog()
        self.pulseBeam.delete_frog() 
        
    def recreate_frogs(self):
        self.initialPulseBeam.recreate_frog()
        self.pulseBeam.recreate_frog()
            
