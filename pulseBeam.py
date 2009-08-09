
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
        self.InterferometricAutoCo =  zeros((self.NT))
        self.InterferometricAutoCoEnvelope =  zeros((self.NT))
        self.AutoCoFFT             =  zeros((self.NT/2))
        self.IntensiometricAutoCo  =  zeros((self.NT))
        self.FROG                  = zeros((self.NT,self.NT))
        self.frequencies           = 2*pi*fftfreq(self.NT,self.deltaT/self.NT)
        self.wavelengths           = 3e8/(self.frequencies/2/pi+self.freqZero) 
        self.t                     = arange(-self.deltaT/2,self.deltaT/2,self.deltaT/self.NT) 
        self.energy                = 0
        self.peak_power            = 0
        self.rate                  = 0
        self.amplitudeE            = 0
        
        self.initialGVD            = 0
        self.initialTOD            = 0
        self.initialFOD            = 0
        self.initialTemporalFWHM = 0
        self.initialSpectralFWHM = 0
        self.initialPulseShapeLoader = None
        self.initialSpectralShapeLoader = None
        
    def copy(self,pulseBeam):
        self.ElectricField[:]     = pulseBeam.ElectricField[:]
        self.FFT[:]               = pulseBeam.FFT[:]
        self.Spectrum[:]          = pulseBeam.Spectrum[:]
        self.BeamProfile_spot       = pulseBeam.BeamProfile_spot
        self.BeamProfile_curvature  = pulseBeam.BeamProfile_curvature
        self.InterferometricAutoCoEnvelope[:]   = pulseBeam.InterferometricAutoCoEnvelope[:]
        self.InterferometricAutoCo[:]   = pulseBeam.InterferometricAutoCo[:]
#        self.AutoCoFFT[:]           = pulseBeam.AutoCoFFT[:]
        self.IntensiometricAutoCo[:]= pulseBeam.IntensiometricAutoCo[:]
        self.FROG[:,:]              = pulseBeam.FROG[:,:]
        
        self.FROGxmin   = pulseBeam.FROGxmin
        self.FROGxmax   = pulseBeam.FROGxmax
        self.FROGymin   = pulseBeam.FROGymin
        self.FROGymax   = pulseBeam.FROGymax
        
        self.BeamProfile_spot      = pulseBeam.BeamProfile_spot 
        self.BeamProfile_curvature = pulseBeam.BeamProfile_curvature
        
        self.amplitudeE = pulseBeam.amplitudeE
        
        self.initialGVD            = pulseBeam.initialGVD 
        self.initialTOD            = pulseBeam.initialTOD
        self.initialFOD            = pulseBeam.initialFOD 
        
        self.initialTemporalFWHM = pulseBeam.initialTemporalFWHM
        self.initialSpectralFWHM = pulseBeam.initialSpectralFWHM
        
        self.initialPulseShapeLoader = pulseBeam.initialPulseShapeLoader
        self.initialSpectralShapeLoader = pulseBeam.initialSpectralShapeLoader
        
        self.energy = pulseBeam.energy
        self.peak_power = pulseBeam.peak_power
        self.rate = pulseBeam.rate

    def update_spectrum(self):
        self.Spectrum[:] = self.FFT**2
        
    def field_to_spectrum(self):
        self.FFT[:] = fft(self.ElectricField) #fft(ifftshift(self.ElectricField[:])) 
        self.update_spectrum()
        #self.Spectrum[:self.NT/4] = self.FFT[3*self.NT/4:]**2
        #self.Spectrum[self.NT/4:] = self.FFT[:self.NT/4]**2
    
    def spectrum_to_field(self): 
        #self.FFT[:]=zeros((self.NT))
        #self.FFT[3*self.NT/4:] = sqrt(self.Spectrum[:self.NT/4])
        #self.FFT[:self.NT/4]   = sqrt(self.Spectrum[self.NT/4:])
        self.ElectricField[:] = ifft(self.FFT)
        
    def change_window(self, NT, deltaT):
        if(self.initialTemporalFWHM != 0):
            #save data
            tmp_timeFWHM = self.initialTemporalFWHM
            tmp_GVD = self.initialGVD
            tmp_TOD = self.initialTOD
            tmp_FOD = self.initialFOD
            tmp_spot = self.BeamProfile_spot
            tmp_curvature = self.BeamProfile_curvature
            tmp_peak_power = self.peak_power
            tmp_energy = self.energy
            tmp_rate = self.rate
            self.__init__(NT,deltaT,self.lambdaZero)
            self.initialize_pulse(tmp_timeFWHM , tmp_GVD , tmp_TOD , tmp_FOD , tmp_spot, tmp_curvature, tmp_peak_power , tmp_energy , tmp_rate)
        else:
            tmp_spectralFWHM = self.initialSpectralFWHM
            tmp_GVD = self.initialGVD
            tmp_TOD = self.initialTOD
            tmp_FOD = self.initialFOD
            tmp_spot = self.BeamProfile_spot
            tmp_curvature = self.BeamProfile_curvature
            tmp_peak_power = self.peak_power
            tmp_energy = self.energy
            tmp_rate = self.rate
            self.__init__(NT,deltaT,self.lambdaZero)
            self.initialize_spectrum(tmp_spectralFWHM , tmp_GVD , tmp_TOD , tmp_FOD , tmp_spot, tmp_curvature, tmp_peak_power , tmp_energy , tmp_rate)
    
    def change_central_wavelength(self, lambdaZero):
        if(self.initialTemporalFWHM != 0):
            #save data
            tmp_timeFWHM = self.initialTemporalFWHM
            tmp_GVD = self.initialGVD
            tmp_TOD = self.initialTOD
            tmp_FOD = self.initialFOD
            tmp_spot = self.BeamProfile_spot
            tmp_curvature = self.BeamProfile_curvature
            tmp_peak_power = self.peak_power
            tmp_energy = self.energy
            tmp_rate = self.rate
            self.__init__(self.NT,self.deltaT,lambdaZero)
            self.initialize_pulse(tmp_timeFWHM , tmp_GVD , tmp_TOD , tmp_FOD , tmp_spot, tmp_curvature, tmp_peak_power , tmp_energy , tmp_rate)
        else:
            tmp_spectralFWHM = self.initialSpectralFWHM
            tmp_GVD = self.initialGVD
            tmp_TOD = self.initialTOD
            tmp_FOD = self.initialFOD
            tmp_spot = self.BeamProfile_spot
            tmp_curvature = self.BeamProfile_curvature
            tmp_peak_power = self.peak_power
            tmp_energy = self.energy
            tmp_rate = self.rate
            self.__init__(self.NT,self.deltaT,lambdaZero)
            self.initialize_spectrum(tmp_spectralFWHM , tmp_GVD , tmp_TOD , tmp_FOD , tmp_spot, tmp_curvature, tmp_peak_power , tmp_energy , tmp_rate)
                
        
        
    def initialize_pulse(self, timeFWHM, GVD,TOD,FOD, spot, curvature, peak_power, energy, rate):
        self.BeamProfile_spot = spot
        self.BeamProfile_curvature = curvature
        self.energy = energy
        self.peak_power = peak_power
        self.rate = rate
        
        self.initialGVD = GVD
        self.initialTOD = TOD 
        self.initialFOD = FOD
        
        self.initialTemporalFWHM = timeFWHM
        self.initialSpectralFWHM = 0
        
        t_sigma = timeFWHM/2.35
        
        #phase = self.freqZero*2*3.14*self.t
        #self.ElectricField[:]   = (cos(phase)+1j*sin(phase))*exp(-self.t**2/t_sigma**2/4)
        self.ElectricField[:]   = exp(-self.t**2/t_sigma**2/4)
        
        self.rescale_field_to_energy_or_peak_power()
        
        self.field_to_spectrum()
        self.apply_dispersion(GVD,TOD,FOD)
        self.spectrum_to_field()
        



    def initialize_load_envelope(self, loader, GVD, TOD, FOD, spot, curvature, peak_power, energy, rate):   
        self.BeamProfile_spot = spot
        self.BeamProfile_curvature = curvature 
        self.energy = energy
        self.peak_power = peak_power
        self.rate = rate

        self.initialGVD = GVD
        self.initialTOD = TOD 
        self.initialFOD = FOD

        self.initialTemporalFWHM = 0
        self.initialSpectralFWHM = 0

        
        loader.detect_time_scale()
        
        self.ElectricField[:]   = loader.get_data_for_x_series(self.t)
        
        print self.ElectricField
        
        self.rescale_field_to_energy_or_peak_power()
        
        self.field_to_spectrum()
        self.apply_dispersion(GVD,TOD,FOD)
        self.spectrum_to_field()
        
        
        
        
    
    def initialize_spectrum(self, spectralfwhm, GVD, TOD, FOD, spot, curvature, peak_power, energy, rate):   
        self.BeamProfile_spot = spot
        self.BeamProfile_curvature = curvature 
        self.energy = energy
        self.peak_power = peak_power
        self.rate = rate
        
        self.initialGVD = GVD
        self.initialTOD = TOD 
        self.initialFOD = FOD
        
        self.initialTemporalFWHM = 0
        self.initialSpectralFWHM = spectralfwhm
        
        deltaFreq = 3e8/self.lambdaZero**2*spectralfwhm/2.35*2
        
        freq = self.frequencies/2/pi
        self.FFT[:] = exp(-freq**2/2/deltaFreq**2)
        
        self.apply_dispersion(GVD,TOD,FOD) 
        self.spectrum_to_field()

        #field is shifted due to fft, shift it back
        self.ElectricField = fftshift(self.ElectricField)
        
        self.rescale_field_to_energy_or_peak_power()
            
        
        
        
        
    def initialize_load_spectrum(self, loader, GVD, TOD, FOD, spot, curvature, peak_power, energy, rate):   
        self.BeamProfile_spot = spot
        self.BeamProfile_curvature = curvature 
        self.energy = energy
        self.peak_power = peak_power
        self.rate = rate

        self.initialGVD = GVD
        self.initialTOD = TOD 
        self.initialFOD = FOD

        self.initialTemporalFWHM = 0
        self.initialSpectralFWHM = 0

        wavelengths = self.wavelengths[:]

    
        loader.detect_wavelength_scale()
        
        self.FFT[:] = sqrt(loader.get_data_for_x_series(wavelengths[:]))
        
        self.apply_dispersion(GVD,TOD,FOD) 
        self.spectrum_to_field()

        #field is shifted due to fft, shift it back
        self.ElectricField = fftshift(self.ElectricField)

        self.rescale_field_to_energy_or_peak_power()
    
        
        
        
        
        
    def rescale_field_to_energy_or_peak_power(self):
        self.ElectricField[:] /= max(abs(self.ElectricField[:]))
        
        if(self.energy > 0):
            self.amplitudeE = sqrt(self.energy/(sum(abs(self.ElectricField[:])**2)*self.deltaT/self.NT)) 
        else:
            self.amplitudeE = sqrt(self.peak_power) 
 
        self.ElectricField[:]   *= self.amplitudeE
        
        
        
    def apply_dispersion(self,GVD,TOD,FOD):
        freq = self.frequencies#-self.freqZero)
        delta = exp(1j*(freq**2*GVD/2.+ \
                           freq**3*TOD/2./3.+ \
                           freq**4*FOD/2./3./4.))
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
        w, R = self.beam_calc_propagation(length,n)

        self.BeamProfile_curvature = R
        self.BeamProfile_spot = w
    
    def beam_calc_propagation(self,length,n,w=-1,R=0):
        A = 1
        B = length*n
        C = 0 
        D = 1
        
        if(w <= 0):
            return self.beam_calc_ABCD(self.BeamProfile_spot,self.BeamProfile_curvature,A,B,C,D)
        else:
            return self.beam_calc_ABCD(w,R,A,B,C,D)
        
    def beam_apply_thinlens(self,f):
        w, R = self.beam_calc_thinlens(f)

        self.BeamProfile_curvature = R
        self.BeamProfile_spot = w
        
    def beam_calc_thinlens(self,f,w=-1,R=0):
        A = 1
        B = 0
        C = -1/f 
        D = 1
        
        if(w <= 0):
            return self.beam_calc_ABCD(self.BeamProfile_spot,self.BeamProfile_curvature,A,B,C,D)
        else:
            return self.beam_calc_ABCD(w,R,A,B,C,D)
        
    def beam_apply_refraction(self,n1,n2):
        w, R = self.beam_calc_refraction(n1,n2)

        self.BeamProfile_curvature = R
        self.BeamProfile_spot = w
        
    def beam_calc_refraction(self,n1,n2,w=-1,R=0):
        A = 1
        B = 0
        C = 0
        D = n1/n2
        
        if(w <= 0):
            return self.beam_calc_ABCD(self.BeamProfile_spot,self.BeamProfile_curvature,A,B,C,D)
        else:
            return self.beam_calc_ABCD(w,R,A,B,C,D)
                
    def beam_calc_ABCD(self,w,R,A,B,C,D):
        invQ1 = 1/R-1j*self.lambdaZero/pi/(w**2)
        invQ2 = (C+D*invQ1)/(A+B*invQ1)
        
        R = 1/real(invQ2)
        w = sqrt(-1/imag(invQ2)/pi*self.lambdaZero)
        
        return w,R
        
        
        
    def calculate_autoco_helper(self,E):
        #interferometric autoco Calculation by Joao Silva
        I = abs(E)**2
        
        AutoCo = 2.0*sum(I**2)

        #Re[[4I(t)E(t)]E*(t-tau)]
        AutoCo += 4.0*real(correlate(I*E,E.conj(),"full"))[self.NT/2:self.NT*3/2]

        #4Re[[E(t)][I(t-tau)E*(t-tau)]]
        AutoCo += 4.0*real(correlate(E,I*E.conj(),"full"))[self.NT/2:self.NT*3/2]

        #2Re[E^2(t)E*^2(t-tau)]
        AutoCo += 2.0*real(correlate(E**2,E.conj()**2,"full"))[self.NT/2:self.NT*3/2]

        #4I(t)I(t-tau)
        AutoCo += 4.0*correlate(I,I,"full")[self.NT/2:self.NT*3/2]
        
        return AutoCo
        
        
    def calculate_autoco(self):
    #calculates 2 things:
    #   -interferometric autoco (if sampling rate is sufficient)
    #   -intensiometric autoco

    
    #interferometric autoco calculation from Joao Silva
        E = self.get_real_electric_field()
        I = abs(self.ElectricField)**2
        
        if(not E is None):
            self.InterferometricAutoCo[:] = self.calculate_autoco_helper(E)
        else:
            self.InterferometricAutoCo[:] = zeros((self.NT))

        #4I(t)I(t-tau)
        self.IntensiometricAutoCo[:] = 4.0*correlate(I,I,"full")[self.NT/2:self.NT*3/2]
    
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
            
     
        self.InterferometricAutoCo -= min(abs(self.InterferometricAutoCo))
        self.InterferometricAutoCo /= max(abs(self.InterferometricAutoCo))
        self.InterferometricAutoCo *= 8
        self.IntensiometricAutoCo -= min(abs(self.IntensiometricAutoCo))
        self.IntensiometricAutoCo /= max(abs(self.IntensiometricAutoCo))
        self.IntensiometricAutoCo *= 2
        self.IntensiometricAutoCo += 1
        
        

        #self.AutoCoFFT[:] = fft(self.AutoCo)[:self.NT/2]
        #self.AutoCoFFT /= max(self.AutoCoFFT)
        



    def calculate_FROG(self):
        
        self.FROGxmin   = self.t[0]
        self.FROGxmax   = self.t[-1]
        self.FROGymin   = self.get_frequencies()[0]
        self.FROGymax   = self.get_frequencies()[-1]
        
        ElectricField = real(self.ElectricField)
        # complex field
        
#negative delay
        for i in xrange(self.NT/2):
            self.field     =  zeros((self.NT))
            self.field[i:] = (ElectricField[:self.NT-i]*ElectricField[i:])**2
            self.field_fft = fft(self.field)       
            self.FROG[:,self.NT/2-i] = fftshift(abs(self.field_fft)**2)
            
#positive delay
        for i in xrange(self.NT/2):
            self.field     = zeros((self.NT))
            self.field[i:] = (ElectricField[i:]*ElectricField[:self.NT-i])**2
            self.field_fft = fft(self.field)  
            self.FROG[:,self.NT/2+i] = fftshift(abs(self.field_fft)**2)
            
        self.FROG /= ma.max(abs(self.FROG))
        
        
        
    def fwhm_helper(self,envelope,divisor):
        max_envelope = max(envelope)
        i = 0
        
        
        #find first intersection
        while(envelope[i] < max_envelope/divisor and i < len(envelope)):
           i+=1
        #interpolate
    #    i_ = i+(envelope[i-1]-max_envelope/2.)/(envelope[i]-envelope[i-1])
        i_ = i+(max_envelope/divisor-envelope[i])/(envelope[i+1]-envelope[i])
           
           
        #find second intersection
        j = len(envelope)-1
        while(envelope[j] < max_envelope/divisor and j > 0):
           j -= 1
        #interpolate
        #j_ = j+(envelope[j]-max_envelope/2.)/(envelope[j]-envelope[j+1])
        j_ = j+(max_envelope/divisor-envelope[j])/(envelope[j+1]-envelope[j])
        
        if(j_-i_ < 0):
            return 0 #strangely, this happens with very odd spectral shapes 
        return j_-i_
        
    def get_temporal_fwhm(self): # return two approximations to FWHM
        #calculate intensity envelope
        envelope = abs(self.ElectricField[:])**2
        
        # calculate FWHM by gaussian fit
        X = arange(len(envelope))
        x = sum(X*envelope)/sum(envelope)
        width = sqrt(abs(sum((X-x)**2*envelope)/sum(envelope)))
        fwhm1 = 2.35*width*self.deltaT/self.NT
        
        #calculate FWHM by iteractively checking for intersections at max/2
        fwhm2 = self.fwhm_helper(envelope,2)*self.deltaT/self.NT
        
        return fwhm1,fwhm2
        
    def get_spectral_fwhm(self): 
        envelope = self.get_spectral_intensity()
        
        fwhm = self.fwhm_helper(envelope,2)/self.deltaT 
        
        fwhm = 3e8/self.freqZero**2*fwhm
        
        return fwhm
        
    def get_spectral_width(self): 
        envelope = self.get_spectral_intensity()

        width = self.fwhm_helper(envelope,2.718*2.718)/self.deltaT 

        width *= 3e8/self.freqZero**2

        return width
        
        
    def calc_energy(self):
        return sum(abs(self.ElectricField[:]**2))*self.deltaT/self.NT
        
    def calc_peak_power(self):
        return max(abs(self.ElectricField[:]))**2
            
    def calc_peak_intensity(self):
        return max(abs(self.ElectricField[:]))**2/(self.BeamProfile_spot**2*pi)
        
    def calc_CW_power(self):
        return self.calc_energy()*self.rate
        
    def get_beam_spot(self):
        return self.BeamProfile_spot
    
    def get_beam_curvature(self):
        return self.BeamProfile_curvature
        
    def get_rep_rate(self):
        return self.rate
        
    def get_t(self):
        return self.t
    
    def get_frequencies(self):
        return roll(self.frequencies,self.NT/2)/2/pi
    
    def get_wavelengths(self):
        return roll(self.wavelengths,self.NT/2)
    
    def get_temporal_envelope(self): #rename?
        return abs(self.ElectricField)
        
    def get_temporal_phase(self):
        phase = unwrap(angle(self.ElectricField)) 
        return phase - phase[self.NT/2]
    
    def get_real_electric_field(self):
        if(self.deltaT/self.NT > 1./self.freqZero/4.):
            return None #no point calculating the real electric field. it would be undersampled.
        else:
            return real(self.ElectricField*exp(1j*(2*pi*self.freqZero*self.t+self.get_temporal_phase())))
    
    def get_spectral_intensity(self):
        return abs(roll(self.Spectrum,self.NT/2))**2
    
    def get_spectral_phase(self):
        phase = unwrap(angle(roll(self.Spectrum,self.NT/2))) 
        return phase - phase[self.NT/2]
        
    def get_spectral_intensity_and_phase_vs_wavelength(self,wavelength_limit):
        intensity = self.get_spectral_intensity()
        phase = self.get_spectral_phase()
        wavelengths = roll(self.wavelengths,self.NT/2)
        i = 0
        while(wavelengths[i] <= 0 or wavelengths[i] > wavelength_limit and i < self.NT): 
            i += 1 

        return intensity[i:],phase[i:] ,wavelengths[i:]       
        
    def get_interferometric_autoco(self):
        return self.InterferometricAutoCo

    def get_intensiometric_autoco(self):
        return self.IntensiometricAutoCo
        
    def get_SHGFROG(self):
        frog_limits = self.FROGxmin,self.FROGxmax,self.FROGymin,self.FROGymax
        return self.FROG,frog_limits
    
    def delete_frog(self):
        del self.FROG
        
    def recreate_frog(self):
        self.FROG = zeros((self.NT,self.NT))
        
        
    def phase_blank(self,t,data_array,phase_array,threshold): 
        length = len(t) 
        absarray = abs(data_array)
        absarray /= max(absarray) 
        i = 0
        while(i < length and absarray[i] < threshold):
            i += 1
        
        if(i == length):
            return #no information here
        
        j = length-1
        
        while(j >= 0 and absarray[j] < threshold):
            j -= 1
        
        if(j == -1):
            return #no information here
            
        return t[i:j], phase_array[i:j]
