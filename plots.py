#!/usr/bin/env python
# -*- coding: us-ascii -*-

import wx
    
import matplotlib as mpl
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as Toolbar
from numpy import array,vstack,hstack,transpose



class FourPlots(wx.Panel):
    def __init__(self, parent, id = -1, dpi = None, **kwargs):
        wx.Panel.__init__(self, parent, id=id, **kwargs)
        self.figure = mpl.figure.Figure(dpi=75, figsize=(4,3),facecolor='#dfdfdf')
        self.figure.subplots_adjust(hspace = 0.4,wspace=0.4)
        self.canvas = Canvas(self, -1, self.figure)
        self.toolbar = Toolbar(self.canvas)
        self.toolbar.Realize()
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas,1,wx.EXPAND)
        sizer.Add(self.toolbar, 0 , wx.LEFT | wx.EXPAND)
        self.SetSizer(sizer)
        
        self.create_plots()
    
    def create_plots(self):
        # create the four plots
        self.plot1 = self.figure.add_subplot(221)
        self.axes1 = self.figure.gca()
        self.axes1.set_title('Temporal Profile')#,fontsize=12)
        
        self.plot2 = self.figure.add_subplot(223)
        self.axes2 = self.figure.gca()
        self.axes2.set_title('Spectral Profile')#,fontsize=12)
        
        self.plot3 = self.figure.add_subplot(222)
        self.axes3 = self.figure.gca()
        self.axes3.set_title('Autocorrelations')#,fontsize=12)
        
        self.plot4 = self.figure.add_subplot(224)
        self.axes4 = self.figure.gca()
        self.axes4.set_title('SHG FROG')#,fontsize=12)

        self.plot_init = False

    def redraw(self, t_envelope, envelope, electric_field, t_phase, temporal_phase, wavelength, spectrum, wavelength_phase, spectral_phase, t_autoco, inter_autoco, inten_autoco, frog, frog_limits):
        if(self.plot_init == False):
            #temporal profile
            self.plot1_line1, = self.plot1.plot(t_envelope, envelope, 'r')
            if(not electric_field is None):
                self.plot1_line2, = self.plot1.plot(t_envelope, electric_field, 'b')
            self.plot1_twinx = self.plot1.twinx()
            self.plot1_line3, = self.plot1_twinx.plot(t_phase, temporal_phase, 'k')
            
            self.plot1.set_xlabel('Time(s)',fontsize='small')
            self.plot1.set_ylabel('Electric Field(sqrt(Ws^-1))',fontsize='small')
            self.plot1_twinx.set_ylabel('Phase (rad))',fontsize='small')
            
            self.plot1.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
            
            if(len(temporal_phase) > 0):
                max_phase = max(temporal_phase)
                min_phase = min(temporal_phase)
                delta = max_phase-min_phase
                if(delta < 0.1):
                    self.plot1_twinx.set_ylim((-0.09,0.01))
            
            #spectral profile
            self.plot2_line1, = self.plot2.plot(wavelength, spectrum, 'g')
            self.plot2_twinx = self.plot2.twinx()
            self.plot2_line2, = self.plot2_twinx.plot(wavelength_phase, spectral_phase, 'k')
            
            self.plot2.set_xlabel('Wavelength (m)',fontsize='small')
            self.plot2.set_ylabel('Intensity (Wm^-1))',fontsize='small')
            self.plot2_twinx.set_ylabel('Phase (rad))',fontsize='small')
            
            self.plot2.ticklabel_format(style='sci', scilimits=(0,0), axis='x')
            self.plot2_twinx.ticklabel_format(style='sci', scilimits=(0,0), axis='x')
            
            #autoco
            self.plot3_line1, = self.plot3.plot(t_autoco, inter_autoco, 'b')
            self.plot3_line2, = self.plot3.plot(t_autoco, inten_autoco, 'r')            
            self.plot3.set_xlabel('Delay(s)',fontsize='small')
#            self.plot3_twinx = self.plot2.twinx()  #TODO: add autocoFFT
#            self.plot3_line2, = self.plot2_twinx.plot(autoco_fft, 'k')

            #SHG FROG
            self.plot4_imshow = self.plot4.imshow(frog, interpolation='bilinear',extent=frog_limits,aspect="auto")
            self.plot4.set_xlabel('Delay (s)',fontsize='small')
            self.plot4.set_ylabel('Envelope Frequency(rad s^-1)',fontsize='small')
            #TODO: change bilinear to better, but also fast, interpolation
            #extent = pulsebeam.FROGxmin+k*pulsebeam.FROGdeltax, pulsebeam.FROGxmin+l*pulsebeam.FROGdeltax, \
            #             pulsebeam.FROGymin+i*pulsebeam.FROGdeltay, pulsebeam.FROGymin+j*pulsebeam.FROGdeltay
            
            self.plot_init = True
        else:
            #temporal profile
            #self.plot1_line1.set_ydata(envelope)
            self.plot1_line1.set_data(t_envelope,envelope)
            self.plot1.set_xlim((min(t_envelope),max(t_envelope)))
            if(not electric_field is None):
#                self.plot1_line2.set_ydata(electric_field)
                self.plot1_line2.set_data(t_envelope,electric_field)
                self.plot1.set_ylim((min(electric_field),max(envelope)))
            else:
                self.plot1.set_ylim((min(envelope),max(envelope)))
                
            self.plot1_line3.set_data(t_phase,temporal_phase)
            max_phase = max(temporal_phase)
            min_phase = min(temporal_phase)
            delta = max_phase-min_phase
            if(delta < 0.1):
                self.plot1_twinx.set_ylim((-0.09,0.01))
            else:
                self.plot1_twinx.set_ylim((min_phase-delta/10.,max_phase+delta/10.))
            
            #spectral profile
            #self.plot2_line1.set_ydata(spectrum)
            self.plot2_line1.set_data(wavelength,spectrum)
            self.plot2.set_ylim((min(spectrum),max(spectrum)))
            self.plot2_line2.set_data(wavelength_phase,spectral_phase)
            
            if(len(spectral_phase) != 0):
                max_phase = max(spectral_phase)
                min_phase = min(spectral_phase)
                delta = max_phase-min_phase
                self.plot2_twinx.set_ylim((min_phase-delta/10.,max_phase+delta/10.))
            
            #Autoco
            self.plot3.set_xlim((min(t_autoco),max(t_autoco)))
            self.plot3_line1.set_data(t_autoco,inter_autoco)
            self.plot3_line2.set_data(t_autoco,inten_autoco)
            #TODO: add autocofft
            #SHG FROG
            self.plot4_imshow.set_data(frog)
            
        self.figure.canvas.draw()
        
    def reset(self):
        self.plot_init = False
        
        self.figure.clf()
        
        self.create_plots()
        
