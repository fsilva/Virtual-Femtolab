
from numpy import *
from numpy.fft import *

class pulseBeam:

    def __init__(self,NT,deltaT,lambdaZero):
        self.NT         = NT
        self.deltaT     = deltaT
        self.lambdaZero = lambdaZero
        self.freqZero   = 3e8/self.lambdaZero
        #self.t                     = zeros((self.NT),dtype=complex)#linspace(0,self.deltaT,self.NT)-self.deltaT/2.
        #self.frequencies           = zeros((self.NT/2),dtype=complex)
        self.wavelengths           = zeros((self.NT),dtype=complex)
        self.ElectricField         = zeros((self.NT),dtype=complex)
        self.FFT                   = zeros((self.NT),dtype=complex)
        self.Spectrum              = zeros((self.NT),dtype=complex)
        self.BeamProfile_spot      = 0
        self.BeamProfile_curvature = 0
        self.AutoCo                =  zeros((self.NT))
        self.AutoCoFFT             =  zeros((self.NT/2))
        self.IntensiometricAutoCo  =  zeros((self.NT))
        self.FROG                  = zeros((self.NT/2,self.NT))
        #fill the frequencies and wavelengths arrays and t
        #self.frequencies           = arange(0,self.NT/self.deltaT/2,1/self.deltaT)
        #self.frequencies             = zeros((self.NT/2))
        #self.frequencies[:self.NT/4] = arange(-self.NT/self.deltaT/4.,0,1/self.deltaT)
        #self.frequencies[self.NT/4:] = arange(0,self.NT/self.deltaT/4.,1/self.deltaT)
        self.frequencies           = 2*pi*fftfreq(self.NT,self.deltaT/self.NT)
        self.wavelengths           = 3e8/self.frequencies   
        self.t                     = arange(-self.deltaT/2,self.deltaT/2,self.deltaT/self.NT) 
        self.energy                = 0
        self.peak_power            = 0
        self.rate                  = 0
        self.amplitudeE            = 0
        
    def copy(self,pulseBeam):
        self.ElectricField[:]     = pulseBeam.ElectricField[:]
        self.FFT[:]               = pulseBeam.FFT[:]
        self.Spectrum[:]          = pulseBeam.Spectrum[:]
        self.BeamProfile_spot       = pulseBeam.BeamProfile_spot
        self.BeamProfile_curvature  = pulseBeam.BeamProfile_curvature
        self.AutoCo[:]              = pulseBeam.AutoCo[:]
        self.AutoCoFFT[:]           = pulseBeam.AutoCoFFT[:]
        self.IntensiometricAutoCo[:]= pulseBeam.IntensiometricAutoCo[:]
        self.FROG[:,:]              = pulseBeam.FROG[:,:]
        
        self.FROGxmin   = pulseBeam.FROGxmin
        self.FROGxmax   = pulseBeam.FROGxmax
        self.FROGdeltax = pulseBeam.FROGdeltax
        self.FROGymin   = pulseBeam.FROGymin
        self.FROGymax   = pulseBeam.FROGymax
        self.FROGdeltay = pulseBeam.FROGdeltay
        
        self.BeamProfile_spot      = pulseBeam.BeamProfile_spot 
        self.BeamProfile_curvature = pulseBeam.BeamProfile_curvature
        
        self.amplitudeE = pulseBeam.amplitudeE
        
        self.energy = pulseBeam.energy
        self.peak_power = pulseBeam.peak_power
        self.rate = pulseBeam.rate

    def update_spectrum(self):
        self.Spectrum[:] = roll(self.FFT,self.NT)**2
        
    def field_to_spectrum(self):
        self.FFT[:] = fft(self.ElectricField) #fft(ifftshift(self.ElectricField[:])) 
        self.update_spectrum()
        #self.Spectrum[:self.NT/4] = self.FFT[3*self.NT/4:]**2
        #self.Spectrum[self.NT/4:] = self.FFT[:self.NT/4]**2
    
    def spectrum_to_field(self): 
        #self.FFT[:]=zeros((self.NT))
        #self.FFT[3*self.NT/4:] = sqrt(self.Spectrum[:self.NT/4])
        #self.FFT[:self.NT/4]   = sqrt(self.Spectrum[self.NT/4:])
        self.ElectricField[:] = ifft(self.FFT[:])
        
        
    def initialize_pulse(self, timeFWHM, GVD,TOD,QOD,FOD, spot, curvature, peak_power, energy, rate):
        self.BeamProfile_spot = spot
        self.BeamProfile_curvature = curvature
        self.energy = energy
        self.peak_power = peak_power
        self.rate = rate
        
        t_sigma = timeFWHM/2.35
        
        #phase = self.freqZero*2*3.14*self.t
        #self.ElectricField[:]   = (cos(phase)+1j*sin(phase))*exp(-self.t**2/t_sigma**2/4)
        self.ElectricField[:]   = exp(-self.t**2/t_sigma**2/4)
        
        self.rescale_field_to_energy_or_peak_power()
        
        self.field_to_spectrum()
        self.apply_dispersion(GVD,TOD,QOD,FOD)
        self.spectrum_to_field()
        
        
        
    
    def initialize_spectrum(self, deltaFreq, GVD, TOD, QOD, FOD, spot, curvature, peak_power, energy, rate):   
        self.BeamProfile_spot = spot
        self.BeamProfile_curvature = curvature 
        self.energy = energy
        self.peak_power = peak_power
        self.rate = rate
        
        freq = self.frequencies#-self.freqZero
        self.Spectrum[:] = exp(-freq**2/2/deltaFreq**2)
                        
        self.apply_dispersion(GVD,TOD,QOD,FOD) 
        self.spectrum_to_field()
        
        self.rescale_field_to_energy_or_peak_power()
            
        
        
        
        
    def initialize_spectrum_loaded(self, loader, phaseloader, GVD, TOD, QOD,FOD, spot,curvature, peak_power, energy, rate):  
        self.BeamProfile_spot = spot
        self.BeamProfile_curvature = curvature 
        self.energy = energy
        self.peak_power = peak_power
        self.rate = rate

        freq = self.frequencies
        lambdas = 3e8/freq
        spec = loader.get_data_for_wavelengths(lambdas)
        spec /= max(spec)

        if(phaseloader is None):
            phase = zeros((len(self.Spectrum)))
        else:
            phase = phaseloader.get_data_for_wavelengths(lambdas)
        
        self.Spectrum[:] = spec*exp(1j*phase/2.)
            
        self.apply_dispersion(GVD,TOD,QOD,FOD) 
        self.spectrum_to_field()
        
        self.rescale_field_to_energy_or_peak_power()
        
        
    def rescale_field_to_energy_or_peak_power(self):
        self.ElectricField[:] /= max(abs(self.ElectricField[:]))
        
        if(self.energy > 0):
            self.amplitudeE = sqrt(self.energy/(sum(abs(self.ElectricField[:])**2)*self.deltaT/self.NT)) 
        else:
            self.amplitudeE = sqrt(self.peak_power) 
 
        self.ElectricField[:]   *= self.amplitudeE
        
        
        
    def apply_dispersion(self,GVD,TOD,QOD,FOD):
        freq = self.frequencies#-self.freqZero)
        delta = exp(1j*(freq**2*GVD/2.+ \
                           freq**3*TOD/2./3.+ \
                           freq**4*QOD/2./3./4.+ \
                           freq**5*FOD/2./3./4./5.))
        #self.Spectrum *= delta
        self.FFT *= delta 
        self.update_spectrum()


    def apply_spectral_filter(self,lambda1,lambda2,mult_factor):
        if(lambda1 > lambda2):
            tmp = lambda1
            lambda1 = lambda2
            lambda2 = tmp
        for i in xrange(self.NT/2):
            wavelength = 3e8/self.frequencies[i]
            if(lambda1 < wavelength and wavelength < lambda2):
                self.Spectrum[i] *= mult_factor
                self.Spectrum[i] *= mult_factor
            
    def beam_apply_propagation(self,length,n):
        A = 1
        B = length*n
        C = 0 
        D = 1
        
        self.beam_apply_ABCD(A,B,C,D)
        
    def beam_apply_thinlens(self,f):
        A = 1
        B = 0
        C = -1/f 
        D = 1
        
        self.beam_apply_ABCD(A,B,C,D)
        
    def beam_apply_refraction(self,n1,n2):
        A = 1
        B = 0
        C = 0
        D = n1/n2
        
        self.beam_apply_ABCD(A,B,C,D)
        
    def beam_apply_ABCD(self,A,B,C,D):
        R = self.BeamProfile_curvature
        w = self.BeamProfile_spot
        
        invQ1 = 1/R-1j*self.lambdaZero/pi/(w**2)
        
        invQ2 = (C+D*invQ1)/(A+B*invQ1)
        
        R = 1/real(invQ2)
        w = sqrt(-1/imag(invQ2)/pi*self.lambdaZero)
        
        self.BeamProfile_curvature = R
        self.BeamProfile_spot = w
        
        
    def calculate_autoco(self):
    #autoco calculation from Joao Silva
        print 'autoco calculation wrong.'
        E = self.ElectricField[:]
        I = abs(E)**2
        
        self.AutoCo[:] = 2.0*sum(I**2)

        #Re[[4I(t)E(t)]E*(t-tau)]
        self.AutoCo += 4.0*real(correlate(I*E,E.conj(),"full"))[self.NT/2:self.NT*3/2]

        #4Re[[E(t)][I(t-tau)E*(t-tau)]]
        self.AutoCo += 4.0*real(correlate(E,I*E.conj(),"full"))[self.NT/2:self.NT*3/2]

        #2Re[E^2(t)E*^2(t-tau)]
        self.AutoCo += 2.0*real(correlate(E**2,E.conj()**2,"full"))[self.NT/2:self.NT*3/2]

        #4I(t)I(t-tau)
        self.IntensiometricAutoCo[:] = 4.0*correlate(I,I,"full")[self.NT/2:self.NT*3/2]
        self.AutoCo += self.IntensiometricAutoCo
        
    
#negative delay
        #for i in xrange(0,self.NT/2+1,1):
        #    total = sum(abs((real(self.ElectricField[:self.NT-i,0])+real(self.ElectricField[i:,0]))**2)**2) + \
        #            sum(abs((real(self.ElectricField[:self.NT-i,0]))**2)**2)
        #    #total = sum(real(self.ElectricField[:self.NT-i]**3*conjugate(self.ElectricField[i:])**3))
        #    self.AutoCo[self.NT/2-i] = total
        #    total = sum(abs(((self.ElectricField[:self.NT-i,0])*(self.ElectricField[i:,0]))**2))
        #    self.IntensiometricAutoCo[self.NT/2-i] = total
#positive delay
        #for i in xrange(0,self.NT/2,1):
        #    total = sum(abs((real(self.ElectricField[i:,0])+real(self.ElectricField[:self.NT-i,0]))**2)**2) + \
        #            sum(abs((real(self.ElectricField[i:,0]))**2)**2)
        #    #total = sum(real(self.ElectricField[i:]**3*conjugate(self.ElectricField[:self.NT-i])**3))
        #    self.AutoCo[i+self.NT/2] = total
        #    total = sum(abs(((self.ElectricField[i:,0])*(self.ElectricField[:self.NT-i,0]))**2))
        #    self.IntensiometricAutoCo[self.NT/2+i] = total
            
     
        self.AutoCo -= min(abs(self.AutoCo))
        self.AutoCo /= max(abs(self.AutoCo))
        self.AutoCo *= 8
        self.IntensiometricAutoCo -= min(abs(self.IntensiometricAutoCo))
        self.IntensiometricAutoCo /= max(abs(self.IntensiometricAutoCo))
        self.IntensiometricAutoCo *= 2
        self.IntensiometricAutoCo += 1

        self.AutoCoFFT[:] = fft(self.AutoCo)[:self.NT/2]
        self.AutoCoFFT /= max(self.AutoCoFFT)
        



    def calculate_FROG(self):
        print 'SHGFROG calc wrong'
        self.FROGxmin   = self.t[0]
        self.FROGxmax   = self.t[-1]
        self.FROGdeltax = (self.t[-1]-self.t[0])/self.NT
        self.FROGymin   = self.frequencies[0]
        self.FROGymax   = self.frequencies[-1]
        self.FROGdeltay = (self.frequencies[-1]-self.frequencies[0])/self.NT

#negative delay
        for i in xrange(self.NT/2):
            self.field     =  zeros((self.NT))
            self.field[i:] = (self.ElectricField[:self.NT-i]*self.ElectricField[i:])**2
            self.field_fft = fft(self.field)       
            self.field_fft[:self.NT/2] = zeros((self.NT/2))   #no negative freq
            self.FROG[:,self.NT/2-i] = abs(self.field_fft[self.NT/2:])**2#fftshift(abs(self.field_fft[self.NT/2:])**2)
            
#positive delay
        for i in xrange(self.NT/2):
            self.field     = zeros((self.NT))
            self.field[i:] = (self.ElectricField[i:]*self.ElectricField[:self.NT-i])**2
            self.field_fft = fft(self.field)  
            self.field_fft[:self.NT/2] = zeros((self.NT/2))      #no negative freq
            self.FROG[:,self.NT/2+i] = abs(self.field_fft[self.NT/2:])**2#fftshift(abs(self.field_fft[self.NT/2:])**2)
            
        self.FROG /= ma.max(abs(self.FROG))
        
        
    def calc_fwhm(self): # return two approximations to FWHM
    
        #calculate intensity envelope
        envelope = abs(self.ElectricField[:])**2
        
        # calculate FWHM by gaussian fit
        X = arange(len(envelope))
        x = sum(X*envelope)/sum(envelope)
        width = sqrt(abs(sum((X-x)**2*envelope)/sum(envelope)))
        fwhm1 = 2.35*width*self.deltaT/self.NT
        
        #calculate FWHM by finding the half maximum level intersections
        max_envelope = max(envelope)
        i = 0
        #find first intersection
        while(envelope[i] < max_envelope/2. and i < len(envelope)):
           i+=1
        #interpolate
        i_ = i+(max_envelope/2.-envelope[i-1])/(envelope[i]-envelope[i-1])
           
        #find second intersection
        j = len(envelope)-1
        while(envelope[j] < max_envelope/2. and j > 0):
           j -= 1
        #interpolate
        j_ = j-(max_envelope/2.-envelope[j])/(envelope[j-1]-envelope[j])
        
        fwhm2 = (j_-i_)*self.deltaT/self.NT
        
        return fwhm1,fwhm2
        
    def calc_energy(self):
        return sum(abs(self.ElectricField[:]**2))*self.deltaT/self.NT
        
    def calc_peak_power(self):
        return max(abs(self.ElectricField[:]))**2
            
    def calc_peak_intensity(self):
        return max(abs(self.ElectricField[:]))**2/(self.BeamProfile_spot**2*pi)
        
    def calc_CW_power(self):
        return self.calc_energy()*self.rate
        
    def get_t(self):
        return self.t
    
    def get_frequencies(self):
        return roll(self.frequencies,self.NT/2)
    
    def get_temporal_intensity(self):
        return abs(self.ElectricField)**2
        
    def get_temporal_phase(self):
        phase = unwrap(angle(self.ElectricField)) #TODO: should we multiply by 2, because of intensity? (**2)
        return phase - phase[self.NT/2]
    
    def get_spectral_intensity(self):
        return abs(roll(self.Spectrum,self.NT/2))**2

    def get_spectral_phase(self):
        phase = unwrap(angle(roll(self.Spectrum,self.NT/2))) #TODO: should we multiply by 2, because of intensity? (**2)
        return phase - phase[self.NT/2]
        
    def get_autoco(self):
        return self.AutoCo
        
    def get_SHGFROG(self):
        return self.FROG
        
    def phase_blank(self,t,data_array,phase_array,threshold):
        absarray = abs(data_array)
        absarray /= max(absarray) 
        i = 0
        while(i < self.NT and absarray[i] < threshold):
            i += 1
        
        if(i == self.NT):
            return #no information here
        
        j = self.NT-1
        
        while(j >= 0 and absarray[j] < threshold):
            j -= 1
        
        if(j == -1):
            return #no information here
            
        return t[i:j], phase_array[i:j]
