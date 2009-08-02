#!/usr/bin/env python
# -*- coding: us-ascii -*-
# generated by wxGlade 0.6.3 on Sun Aug  2 17:33:21 2009

import wx

# begin wxGlade: extracode
# end wxGlade



class EditInitialPulse(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: EditInitialPulse.__init__
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.notebook_1 = wx.Notebook(self, -1, style=0)
        self.notebook_1_pane_4 = wx.Panel(self.notebook_1, -1)
        self.notebook_1_pane_3 = wx.Panel(self.notebook_1, -1)
        self.notebook_1_pane_2 = wx.Panel(self.notebook_1, -1)
        self.notebook_1_pane_1 = wx.Panel(self.notebook_1, -1)
        self.radio_btn_1 = wx.RadioButton(self.notebook_1_pane_1, -1, "Temporal FWHM", style=wx.RB_GROUP)
        self.temporalFWHM = wx.TextCtrl(self.notebook_1_pane_1, -1, "6.6", style=wx.TE_RIGHT)
        self.label_1 = wx.StaticText(self.notebook_1_pane_1, -1, "fs")
        self.radio_btn_2 = wx.RadioButton(self.notebook_1_pane_1, -1, "Spectral FWHM")
        self.spectralFWHM = wx.TextCtrl(self.notebook_1_pane_1, -1, "200", style=wx.TE_RIGHT)
        self.label_3 = wx.StaticText(self.notebook_1_pane_1, -1, "nm")
        self.radio_btn_3 = wx.RadioButton(self.notebook_1_pane_1, -1, "Load Electric Field")
        self.radio_btn_4 = wx.RadioButton(self.notebook_1_pane_1, -1, "Load Spectrum")
        self.label_13 = wx.StaticText(self.notebook_1_pane_2, -1, "GVD")
        self.gvd = wx.TextCtrl(self.notebook_1_pane_2, -1, "0.0", style=wx.TE_RIGHT)
        self.label_14 = wx.StaticText(self.notebook_1_pane_2, -1, "fs ^2", style=wx.ALIGN_CENTRE)
        self.label_15 = wx.StaticText(self.notebook_1_pane_2, -1, "TOD")
        self.tod = wx.TextCtrl(self.notebook_1_pane_2, -1, "0.0", style=wx.TE_RIGHT)
        self.label_16 = wx.StaticText(self.notebook_1_pane_2, -1, "fs ^3")
        self.label_17 = wx.StaticText(self.notebook_1_pane_2, -1, "FOD")
        self.fod = wx.TextCtrl(self.notebook_1_pane_2, -1, "0.0", style=wx.TE_RIGHT)
        self.label_18 = wx.StaticText(self.notebook_1_pane_2, -1, "fs ^4")
        self.label_4 = wx.StaticText(self.notebook_1_pane_3, -1, "Spot Size")
        self.spotsize = wx.TextCtrl(self.notebook_1_pane_3, -1, "10", style=wx.TE_RIGHT)
        self.label_6 = wx.StaticText(self.notebook_1_pane_3, -1, "mm")
        self.label_5 = wx.StaticText(self.notebook_1_pane_3, -1, "Curvature")
        self.curvature = wx.TextCtrl(self.notebook_1_pane_3, -1, "1000", style=wx.TE_RIGHT)
        self.label_7 = wx.StaticText(self.notebook_1_pane_3, -1, "m")
        self.label_8 = wx.StaticText(self.notebook_1_pane_3, -1, "Waist at\n(vacuum)")
        self.waist_distance = wx.StaticText(self.notebook_1_pane_3, -1, "")
        self.label_10 = wx.StaticText(self.notebook_1_pane_3, -1, "m")
        self.label_11 = wx.StaticText(self.notebook_1_pane_3, -1, "Spot at waist")
        self.spot_waist = wx.StaticText(self.notebook_1_pane_3, -1, "")
        self.label_12 = wx.StaticText(self.notebook_1_pane_3, -1, "mm")
        self.label_2 = wx.StaticText(self.notebook_1_pane_4, -1, "Rate")
        self.rate = wx.TextCtrl(self.notebook_1_pane_4, -1, "1", style=wx.TE_RIGHT)
        self.label_1_copy = wx.StaticText(self.notebook_1_pane_4, -1, "KHz")
        self.radio_btn_2_copy = wx.RadioButton(self.notebook_1_pane_4, -1, "Pulse Energy")
        self.pulseenergy = wx.TextCtrl(self.notebook_1_pane_4, -1, "3", style=wx.TE_RIGHT)
        self.label_3_copy = wx.StaticText(self.notebook_1_pane_4, -1, "uJ")
        self.radio_btn_3_copy = wx.RadioButton(self.notebook_1_pane_4, -1, "Peak Power")
        self.peakpower = wx.TextCtrl(self.notebook_1_pane_4, -1, "1", style=wx.TE_RIGHT)
        self.label_9 = wx.StaticText(self.notebook_1_pane_4, -1, "MW")
        self.button_1 = wx.Button(self, -1, "Refresh")
        self.button_2 = wx.Button(self, -1, "Done")

        self.__set_properties()
        self.__do_layout()

        self.Bind(wx.EVT_RADIOBUTTON, self.load_electric_field, self.radio_btn_3)
        self.Bind(wx.EVT_RADIOBUTTON, self.load_spectrum, self.radio_btn_4)
        self.Bind(wx.EVT_BUTTON, self.refresh_click, self.button_1)
        self.Bind(wx.EVT_BUTTON, self.done_click, self.button_2)
        # end wxGlade

    def __set_properties(self):
        # begin wxGlade: EditInitialPulse.__set_properties
        self.SetTitle("Edit Initial Pulse")
        self.radio_btn_1.SetValue(1)
        self.radio_btn_2_copy.SetValue(1)
        self.notebook_1.SetMinSize((573, 164))
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: EditInitialPulse.__do_layout
        sizer_1 = wx.BoxSizer(wx.VERTICAL)
        sizer_2 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_6 = wx.BoxSizer(wx.HORIZONTAL)
        grid_sizer_1_copy = wx.FlexGridSizer(4, 3, 5, 5)
        sizer_5 = wx.BoxSizer(wx.HORIZONTAL)
        grid_sizer_2 = wx.FlexGridSizer(4, 3, 5, 5)
        sizer_3 = wx.BoxSizer(wx.HORIZONTAL)
        grid_sizer_3 = wx.FlexGridSizer(3, 3, 5, 5)
        sizer_4 = wx.BoxSizer(wx.HORIZONTAL)
        grid_sizer_1 = wx.FlexGridSizer(4, 3, 5, 5)
        sizer_4.Add((20, 20), 1, 0, 0)
        grid_sizer_1.Add(self.radio_btn_1, 0, 0, 0)
        grid_sizer_1.Add(self.temporalFWHM, 0, 0, 0)
        grid_sizer_1.Add(self.label_1, 0, 0, 0)
        grid_sizer_1.Add(self.radio_btn_2, 0, 0, 0)
        grid_sizer_1.Add(self.spectralFWHM, 0, 0, 0)
        grid_sizer_1.Add(self.label_3, 0, 0, 0)
        grid_sizer_1.Add(self.radio_btn_3, 0, 0, 0)
        grid_sizer_1.Add((20, 20), 0, 0, 0)
        grid_sizer_1.Add((20, 20), 0, 0, 0)
        grid_sizer_1.Add(self.radio_btn_4, 0, 0, 0)
        grid_sizer_1.Add((20, 20), 0, 0, 0)
        grid_sizer_1.Add((20, 20), 0, 0, 0)
        sizer_4.Add(grid_sizer_1, 2, wx.EXPAND, 0)
        sizer_4.Add((20, 20), 1, 0, 0)
        self.notebook_1_pane_1.SetSizer(sizer_4)
        sizer_3.Add((20, 20), 1, 0, 0)
        grid_sizer_3.Add(self.label_13, 0, 0, 0)
        grid_sizer_3.Add(self.gvd, 0, 0, 0)
        grid_sizer_3.Add(self.label_14, 0, 0, 0)
        grid_sizer_3.Add(self.label_15, 0, 0, 0)
        grid_sizer_3.Add(self.tod, 0, 0, 0)
        grid_sizer_3.Add(self.label_16, 0, 0, 0)
        grid_sizer_3.Add(self.label_17, 0, 0, 0)
        grid_sizer_3.Add(self.fod, 0, 0, 0)
        grid_sizer_3.Add(self.label_18, 0, 0, 0)
        sizer_3.Add(grid_sizer_3, 1, wx.EXPAND, 0)
        sizer_3.Add((20, 20), 1, 0, 0)
        self.notebook_1_pane_2.SetSizer(sizer_3)
        sizer_5.Add((20, 20), 1, 0, 0)
        grid_sizer_2.Add(self.label_4, 0, 0, 0)
        grid_sizer_2.Add(self.spotsize, 0, 0, 0)
        grid_sizer_2.Add(self.label_6, 0, 0, 0)
        grid_sizer_2.Add(self.label_5, 0, 0, 0)
        grid_sizer_2.Add(self.curvature, 0, 0, 0)
        grid_sizer_2.Add(self.label_7, 0, 0, 0)
        grid_sizer_2.Add(self.label_8, 0, 0, 0)
        grid_sizer_2.Add(self.waist_distance, 0, 0, 0)
        grid_sizer_2.Add(self.label_10, 0, 0, 0)
        grid_sizer_2.Add(self.label_11, 0, 0, 0)
        grid_sizer_2.Add(self.spot_waist, 0, 0, 0)
        grid_sizer_2.Add(self.label_12, 0, 0, 0)
        sizer_5.Add(grid_sizer_2, 1, wx.EXPAND, 0)
        sizer_5.Add((20, 20), 1, 0, 0)
        self.notebook_1_pane_3.SetSizer(sizer_5)
        sizer_6.Add((20, 20), 1, 0, 0)
        grid_sizer_1_copy.Add(self.label_2, 0, 0, 0)
        grid_sizer_1_copy.Add(self.rate, 0, 0, 0)
        grid_sizer_1_copy.Add(self.label_1_copy, 0, 0, 0)
        grid_sizer_1_copy.Add(self.radio_btn_2_copy, 0, 0, 0)
        grid_sizer_1_copy.Add(self.pulseenergy, 0, 0, 0)
        grid_sizer_1_copy.Add(self.label_3_copy, 0, 0, 0)
        grid_sizer_1_copy.Add(self.radio_btn_3_copy, 0, 0, 0)
        grid_sizer_1_copy.Add(self.peakpower, 0, 0, 0)
        grid_sizer_1_copy.Add(self.label_9, 0, 0, 0)
        sizer_6.Add(grid_sizer_1_copy, 1, wx.EXPAND, 0)
        sizer_6.Add((20, 20), 1, 0, 0)
        self.notebook_1_pane_4.SetSizer(sizer_6)
        self.notebook_1.AddPage(self.notebook_1_pane_1, "Temporal/Spectral Shape")
        self.notebook_1.AddPage(self.notebook_1_pane_2, "Spectral Phase")
        self.notebook_1.AddPage(self.notebook_1_pane_3, "Beam Parameters")
        self.notebook_1.AddPage(self.notebook_1_pane_4, "Energy Parameters")
        sizer_1.Add(self.notebook_1, 1, wx.LEFT|wx.RIGHT|wx.EXPAND|wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 5)
        sizer_2.Add((20, 20), 1, 0, 0)
        sizer_2.Add(self.button_1, 0, wx.ALL, 5)
        sizer_2.Add((20, 20), 0, 0, 0)
        sizer_2.Add(self.button_2, 0, wx.ALL, 5)
        sizer_2.Add((20, 20), 0, 0, 0)
        sizer_1.Add(sizer_2, 0, wx.EXPAND, 0)
        self.SetSizer(sizer_1)
        sizer_1.Fit(self)
        self.Layout()
        # end wxGlade
        
    def set_info(self,initialPulseBeam,refresh_function):
        self.initialPulseBeam = initialPulseBeam
        self.refresh_function = refresh_function
        
        #load values into GUI
        
        ipb = initialPulseBeam
        
        #shape
        if(ipb.initialTemporalFWHM != 0):
            self.radio_btn_1.SetValue(1)
        else:
            self.radio_btn_2.SetValue(1)
        self.temporalFWHM.SetValue(str(ipb.initialTemporalFWHM*1e15))
        self.spectralFWHM.SetValue(str(ipb.initialSpectralFWHM*1e9))
        
        #phase
        self.gvd.SetValue(str(ipb.initialGVD*1e30))
        self.tod.SetValue(str(ipb.initialTOD*1e45))
        self.fod.SetValue(str(ipb.initialFOD*1e60))
        
        #beam spot and curvature
        self.spotsize.SetValue(str(ipb.BeamProfile_spot*1000))
        self.curvature.SetValue(str(ipb.BeamProfile_curvature))
        self.calculate_waist()
        
        #energy/etc
        self.rate.SetValue(str(ipb.rate*0.001))
        self.peakpower.SetValue(str(ipb.peak_power*1e-6))
        self.pulseenergy.SetValue(str(ipb.energy*1e6))
        if(ipb.energy == 0):
            self.radio_btn_3_copy.SetValue(1)
        else:
            self.radio_btn_2_copy.SetValue(1)
            
    def calculate_waist(self):
        print 'TODO: calculate waist'

            

    def load_electric_field(self, event): # wxGlade: EditInitialPulse.<event_handler>
        print "Event handler `load_electric_field' not implemented!"
        event.Skip()

    def load_spectrum(self, event): # wxGlade: EditInitialPulse.<event_handler>
        print "Event handler `load_spectrum' not implemented!"
        event.Skip()

    def refresh_click(self, event): # wxGlade: EditInitialPulse.<event_handler>
        
        #energy 
        try:
            self.initialPulseBeam.rate = float(self.rate.GetValue())*1000
        except:
            self.rate.SetValue(str(self.initialPulseBeam.rate))    
        if(self.radio_btn_3_copy.GetValue()==1): #peak power
            try:
                self.initialPulseBeam.peak_power = float(self.peakpower.GetValue())*1000
                self.initialPulseBeam.energy = 0
            except:
                self.peakpower.SetValue(str(self.initialPulseBeam.peak_power*1e-6))
        else: #energy
            try:
                self.initialPulseBeam.energy = float(self.pulseenergy.GetValue())*1e-6
                self.initialPulseBeam.peak_power = 0
            except:
                self.pulseenergy.SetValue(str(self.initialPulseBeam.energy*1e-6))
                
        #beam spot and curvature
        try:
            self.initialPulseBeam.BeamProfile_spot = float(self.spotsize.GetValue())*1e-3
        except: 
            self.spotsize.SetValue(str(self.initialPulseBeam.BeamProfile_spot*1e3))
        try:
            self.initialPulseBeam.BeamProfile_curvature = float(self.curvature.GetValue())
        except: 
            self.curvature.SetValue(str(self.initialPulseBeam.BeamProfile_curvature))
            
        #phase
        try:
            self.initialPulseBeam.initialGVD = float(self.gvd.GetValue())*1e-30
        except:
            self.gvd.SetValue(str(self.initialPulseBeam.initialGVD*1e-30))
        try:
            self.initialPulseBeam.initialTOD = float(self.tod.GetValue())*1e-45
        except:
            self.tod.SetValue(str(self.initialPulseBeam.initialTOD*1e-30))
        try:
            self.initialPulseBeam.initialFOD = float(self.fod.GetValue())*1e-60
        except:
            self.fod.SetValue(str(self.initialPulseBeam.initialFOD*1e-30))    
            
        #shape
        if(self.radio_btn_1.GetValue()==1):
            #temporal fwhm
            try:
                self.initialPulseBeam.initialTemporalFWHM = float(self.temporalFWHM.GetValue())*1e-15
                self.initialPulseBeam.initialize_pulse(self.initialPulseBeam.initialTemporalFWHM,
                                                       self.initialPulseBeam.initialGVD,
                                                       self.initialPulseBeam.initialTOD,
                                                       self.initialPulseBeam.initialFOD,
                                                       self.initialPulseBeam.BeamProfile_spot,
                                                       self.initialPulseBeam.BeamProfile_curvature,
                                                       self.initialPulseBeam.peak_power,
                                                       self.initialPulseBeam.energy,
                                                       self.initialPulseBeam.rate)
            except:
                self.temporalFWHM.SetValue(str(self.initialPulseBeam.initialTemporalFWHM*1e15))
        elif(self.radio_btn_2.GetValue()==1):
            try:
                print 'yeah'
                self.initialPulseBeam.initialSpectralFWHM = float(self.spectralFWHM.GetValue())*1e-9
                self.initialPulseBeam.initialize_spectrum(self.initialPulseBeam.initialSpectralFWHM,
                                                       self.initialPulseBeam.initialGVD,
                                                       self.initialPulseBeam.initialTOD,
                                                       self.initialPulseBeam.initialFOD,
                                                       self.initialPulseBeam.BeamProfile_spot,
                                                       self.initialPulseBeam.BeamProfile_curvature,
                                                       self.initialPulseBeam.peak_power,
                                                       self.initialPulseBeam.energy,
                                                       self.initialPulseBeam.rate)
            except:
                self.spectralFWHM.SetValue(str(self.initialPulseBeam.initialSpectralFWHM*1e9)) 
        else:
            pass
            
            
               
        
        self.refresh_function()
        if(not event is None):
            event.Skip()

    def done_click(self, event): # wxGlade: EditInitialPulse.<event_handler>
        self.refresh_click(None)
        self.Close()
        event.Skip()

# end of class EditInitialPulse
