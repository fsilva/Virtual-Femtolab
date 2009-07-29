#!/usr/bin/env python
# -*- coding: us-ascii -*-
# generated by wxGlade 0.6.3 on Mon Jul 27 14:34:59 2009

import wxversion
wxversion.ensureMinimal('2.8')
import wx
import wx.grid
    
import matplotlib as mpl
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as Toolbar

import propagator

class FourPlots(wx.Panel):
    def __init__(self, parent, id = -1, dpi = None, **kwargs):
        wx.Panel.__init__(self, parent, id=id, **kwargs)
        self.figure = mpl.figure.Figure(dpi=75, figsize=(4,3),facecolor='#dfdfdf')
        self.figure.subplots_adjust(hspace = 0.33,wspace=0.33)
        self.canvas = Canvas(self, -1, self.figure)
        self.toolbar = Toolbar(self.canvas)
        self.toolbar.Realize()
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas,1,wx.EXPAND)
        sizer.Add(self.toolbar, 0 , wx.LEFT | wx.EXPAND)
        self.SetSizer(sizer)

        # create the four plots
        self.plot1 = self.figure.add_subplot(221)
        self.axes1 = self.figure.gca()
        self.axes1.set_title('Temporal Profile')#,fontsize=12)
        
        self.plot2 = self.figure.add_subplot(223)
        self.axes2 = self.figure.gca()
        self.axes2.set_title('Spectral Profile')#,fontsize=12)
        
        self.plot3 = self.figure.add_subplot(222)
        self.axes3 = self.figure.gca()
        self.axes3.set_title('Autocorrelation')#,fontsize=12)
        
        self.plot4 = self.figure.add_subplot(224)
        self.axes4 = self.figure.gca()
        self.axes4.set_title('SHG FROG')#,fontsize=12)

        self.plot_init = False

    def redraw(self, t, intensity, tphase, freq, spectrum, sphase, autoco, autoco_fft, frog):
        if(self.plot_init == False):
            #temporal profile
            self.plot1_line1, = self.plot1.plot(t, intensity, 'r')
            self.plot1_twinx = self.plot1.twinx()
            self.plot1_line2, = self.plot1_twinx.plot(t, tphase, 'k')
            
            #spectral profile
            self.plot2_line1, = self.plot2.plot(freq, spectrum, 'g')
            self.plot2_twinx = self.plot2.twinx()
            self.plot2_line2, = self.plot2_twinx.plot(freq, sphase, 'k')
            
            #autoco
            self.plot3_line1, = self.plot3.plot(t, autoco, 'b')
#            self.plot3_twinx = self.plot2.twinx()  #TODO: add autocoFFT
#            self.plot3_line2, = self.plot2_twinx.plot(autoco_fft, 'k')

            #SHG FROG
            extent = 0,1,0,1
            self.plot4_imshow = self.plot4.imshow(frog, interpolation='bilinear',extent=extent,aspect="auto")
            #TODO: change bilinear to better, but also fast, interpolation
            #extent = pulsebeam.FROGxmin+k*pulsebeam.FROGdeltax, pulsebeam.FROGxmin+l*pulsebeam.FROGdeltax, \
            #             pulsebeam.FROGymin+i*pulsebeam.FROGdeltay, pulsebeam.FROGymin+j*pulsebeam.FROGdeltay
            
            self.plot_init = True
        else:
            #temporal profile
            self.plot1_line1.set_ydata(intensity)
            self.plot1_line2.set_ydata(tphase)
            #spectral profile
            self.plot2_line1.set_ydata(spectrum)
            self.plot2_line2.set_ydata(sphase)
            #Autoco
            self.plot3_line1.set_ydata(autoco)
            #TODO: add autocofft
            #SHG FROG
            self.plot4_imshow.set_data(frog)
            
        self.figure.canvas.draw()
        



# begin wxGlade: extracode
# end wxGlade



class VFFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: VFFrame.__init__
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        
        # Menu Bar
        self.MainFrame_menubar = wx.MenuBar()
        wxglade_tmp_menu = wx.Menu()
        wxglade_tmp_menu.Append(wx.NewId(), "Open", "", wx.ITEM_NORMAL)
        wxglade_tmp_menu.Append(wx.NewId(), "Save", "", wx.ITEM_NORMAL)
        wxglade_tmp_menu.AppendSeparator()
        wxglade_tmp_menu.Append(wx.NewId(), "Exit", "", wx.ITEM_NORMAL)
        self.MainFrame_menubar.Append(wxglade_tmp_menu, "File")
        wxglade_tmp_menu = wx.Menu()
        self.MainFrame_menubar.Append(wxglade_tmp_menu, "Preferences")
        wxglade_tmp_menu = wx.Menu()
        wxglade_tmp_menu.Append(wx.NewId(), "Export Plots", "", wx.ITEM_NORMAL)
        wxglade_tmp_menu.Append(wx.NewId(), "Export Animation", "", wx.ITEM_NORMAL)
        self.MainFrame_menubar.Append(wxglade_tmp_menu, "Export")
        wxglade_tmp_menu = wx.Menu()
        self.MainFrame_menubar.Append(wxglade_tmp_menu, "About")
        self.SetMenuBar(self.MainFrame_menubar)
        # Menu Bar end

        self.AddButton = wx.Button(self, wx.ID_ADD, "")
        self.EditButton = wx.Button(self, wx.ID_PROPERTIES, "")
        self.RemoveButton = wx.Button(self, wx.ID_REMOVE, "")
        self.VFData = wx.grid.Grid(self, -1, size=(1, 1))
        self.SchematicPanel = wx.Panel(self, -1)
        self.DistanceSlider = wx.Slider(self, -1, 0, 0, 1000)

        # Add Matplotlib widgets to window
        self.plot = FourPlots(self)

        self.__set_properties()

        #Set Parameter Values in Grid
        self.refresh_grid_information()

        self.__do_layout()

        self.Bind(wx.EVT_SIZE, self.resize)
        self.Bind(wx.EVT_MENU, self.menu_open_click, id=-1)
        self.Bind(wx.EVT_MENU, self.menu_save_click, id=-1)
        self.Bind(wx.EVT_MENU, self.menu_exit_click, id=-1)
        self.Bind(wx.EVT_MENU, self.menu_exportplots_click, id=-1)
        self.Bind(wx.EVT_MENU, self.menu_animation_click, id=-1)
        self.Bind(wx.EVT_BUTTON, self.addbutton_click, self.AddButton)
        self.Bind(wx.EVT_BUTTON, self.editbutton_click, self.EditButton)
        self.Bind(wx.EVT_BUTTON, self.removebutton_click, self.RemoveButton)
        self.Bind(wx.EVT_COMMAND_SCROLL, self.distanceslider_change, self.DistanceSlider)
        # end wxGlade
        
        self.init_calculations()
        self.refresh_interface()

    def __set_properties(self):
        # begin wxGlade: VFFrame.__set_properties
        self.SetTitle("Virtual Femtolab")
        self.SetSize((900,675))
        self.VFData.CreateGrid(10, 3)
        self.VFData.SetRowLabelSize(0)
        self.VFData.SetColLabelSize(0)
        self.VFData.EnableEditing(0)
        self.VFData.SetColLabelValue(0, "Name")
        self.VFData.SetColLabelValue(1, "Value")
        self.VFData.SetColLabelValue(2, "Units")
        self.SchematicPanel.SetMinSize((800, 200))
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: VFFrame.__do_layout
        sizer_1 = wx.BoxSizer(wx.VERTICAL)
        sizer_2 = wx.BoxSizer(wx.VERTICAL)
        sizer_9 = wx.BoxSizer(wx.VERTICAL)
        sizer_6 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_7 = wx.BoxSizer(wx.VERTICAL)
        sizer_8 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_7.Add(self.plot, 1, wx.EXPAND, 0)
        sizer_8.Add(self.AddButton, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        sizer_8.Add(self.EditButton, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        sizer_8.Add(self.RemoveButton, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        sizer_7.Add(sizer_8, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        sizer_6.Add(sizer_7, 3, wx.EXPAND, 0)
        sizer_6.Add(self.VFData, 1, wx.EXPAND|wx.ALIGN_RIGHT, 0)
        sizer_2.Add(sizer_6, 3, wx.EXPAND, 0)
        sizer_9.Add(self.SchematicPanel, 1, wx.EXPAND, 0)
        sizer_9.Add(self.DistanceSlider, 0, wx.EXPAND, 0)
        sizer_2.Add(sizer_9, 1, wx.EXPAND, 0)
        sizer_1.Add(sizer_2, 1, wx.EXPAND, 0)
        self.SetSizer(sizer_1)
        self.Layout()
        # end wxGlade

    def menu_open_click(self, event): # wxGlade: VFFrame.<event_handler>
        print "Event handler `menu_open_click' not implemented!"
        event.Skip()

    def menu_save_click(self, event): # wxGlade: VFFrame.<event_handler>
        print "Event handler `menu_save_click' not implemented!"
        event.Skip()

    def menu_exit_click(self, event): # wxGlade: VFFrame.<event_handler>
        print "Event handler `menu_exit_click' not implemented!"
        event.Skip()

    def menu_exportplots_click(self, event): # wxGlade: VFFrame.<event_handler>
        print "Event handler `menu_exportplots_click' not implemented!"
        event.Skip()

    def menu_animation_click(self, event): # wxGlade: VFFrame.<event_handler>
        print "Event handler `menu_animation_click' not implemented!"
        event.Skip()

    def addbutton_click(self, event): # wxGlade: VFFrame.<event_handler>
        print "Event handler `addbutton_click' not implemented!"
        event.Skip()

    def editbutton_click(self, event): # wxGlade: VFFrame.<event_handler>
        print "Event handler `editbutton_click' not implemented!"
        event.Skip()

    def removebutton_click(self, event): # wxGlade: VFFrame.<event_handler>
        print "Event handler `removebutton_click' not implemented!"
        event.Skip()

    def distanceslider_change(self, event): # wxGlade: VFFrame.<event_handler>
        distance = self.DistanceSlider.GetValue()/1000.*self.propagator.get_max_z()
        self.change_distance(distance)

        event.Skip()
        
    def init_calculations(self):
        self.propagator = propagator.Propagator(128,64.e-15,8.e-7)
        self.propagator.example_pulseBeam()

    def change_distance(self, distance):
        # Recalculate everything
        self.propagator.change_z(distance)
    
        # Refresh interface/redraw
        self.refresh_interface()
    
    def refresh_interface(self):
        pulseBeam = self.propagator.get_pulseBeam()
        t = pulseBeam.get_t()
        intensity = pulseBeam.get_temporal_intensity()
        tphase = pulseBeam.get_temporal_phase()
        freq = pulseBeam.get_frequencies()
        spectrum = pulseBeam.get_spectral_intensity()
        sphase = pulseBeam.get_spectral_phase()
        autoco = pulseBeam.get_autoco()
        frog = pulseBeam.get_SHGFROG()
        self.plot.redraw(t,intensity,tphase,freq,spectrum,sphase,autoco,0,frog)
        
    def refresh_grid_information(self):
        pass
        
    def resize(self,event):
        print 'resize not implemented'
        #TODO: resize properties table
        event.Skip()
        


# end of class VFFrame

def main():
    app = wx.PySimpleApp(0)
    wx.InitAllImageHandlers()
    MainFrame = VFFrame(None, -1, "")
    app.SetTopWindow(MainFrame)
    MainFrame.Show()
    app.MainLoop()
