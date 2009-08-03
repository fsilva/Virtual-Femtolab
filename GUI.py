#!/usr/bin/env python
# -*- coding: us-ascii -*-
# generated by wxGlade 0.6.3 on Mon Jul 27 14:34:59 2009

import wxversion
wxversion.ensureMinimal('2.8')
import wx
import wx.grid
import wx.lib.scrolledpanel
    
import matplotlib as mpl
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as Toolbar
from numpy import array

import propagator
import draw_schematic

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

    def redraw(self, t_envelope, envelope, electric_field, t_phase, temporal_phase, freq, spectrum, freq_phase, spectral_phase, inter_autoco, inten_autoco, autoco_fft, frog, frog_limits):
        if(self.plot_init == False):
            #temporal profile
            self.plot1_line1, = self.plot1.plot(t_envelope, envelope, 'r')
            if(not electric_field is None):
                self.plot1_line2, = self.plot1.plot(t_envelope, electric_field, 'b')
            self.plot1_twinx = self.plot1.twinx()
            self.plot1_line3, = self.plot1_twinx.plot(t_phase, temporal_phase, 'k')
            
            if(len(temporal_phase) > 0):
                max_phase = max(temporal_phase)
                min_phase = min(temporal_phase)
                delta = max_phase-min_phase
                if(delta < 0.1):
                    self.plot1_twinx.set_ylim((-0.09,0.01))
            
            #spectral profile
            self.plot2_line1, = self.plot2.plot(freq, spectrum, 'g')
            self.plot2_twinx = self.plot2.twinx()
            self.plot2_line2, = self.plot2_twinx.plot(freq_phase, spectral_phase, 'k')
            
            #autoco
            self.plot3_line1, = self.plot3.plot(t_envelope, inter_autoco, 'b')
            self.plot3_line2, = self.plot3.plot(t_envelope, inten_autoco, 'r')            
#            self.plot3_twinx = self.plot2.twinx()  #TODO: add autocoFFT
#            self.plot3_line2, = self.plot2_twinx.plot(autoco_fft, 'k')

            #SHG FROG
            self.plot4_imshow = self.plot4.imshow(frog, interpolation='bilinear',extent=frog_limits,aspect="auto")
            #TODO: change bilinear to better, but also fast, interpolation
            #extent = pulsebeam.FROGxmin+k*pulsebeam.FROGdeltax, pulsebeam.FROGxmin+l*pulsebeam.FROGdeltax, \
            #             pulsebeam.FROGymin+i*pulsebeam.FROGdeltay, pulsebeam.FROGymin+j*pulsebeam.FROGdeltay
            
            self.plot_init = True
        else:
            #temporal profile
            self.plot1_line1.set_ydata(envelope)
            if(not electric_field is None):
                self.plot1_line2.set_ydata(electric_field)
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
            self.plot2_line1.set_ydata(spectrum)
            self.plot2.set_ylim((min(spectrum),max(spectrum)))
            self.plot2_line2.set_data(freq_phase,spectral_phase)
            
            if(len(spectral_phase) != 0):
                max_phase = max(spectral_phase)
                min_phase = min(spectral_phase)
                delta = max_phase-min_phase
                self.plot2_twinx.set_ylim((min_phase-delta/10.,max_phase+delta/10.))
            
            #Autoco
            self.plot3_line1.set_ydata(inter_autoco)
            self.plot3_line2.set_ydata(inten_autoco)
            #TODO: add autocofft
            #SHG FROG
            self.plot4_imshow.set_data(frog)
            
        self.figure.canvas.draw()
        
    def reset(self):
        self.plot_init = False
        
        self.figure.clf()
        
        self.create_plots()
        



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
        self.menu_exit = wx.NewId()
        wxglade_tmp_menu.Append(self.menu_exit, "Exit", "", wx.ITEM_NORMAL)
        self.MainFrame_menubar.Append(wxglade_tmp_menu, "File")
        wxglade_tmp_menu = wx.Menu()
        self.menu_compwindow_id = wx.NewId()
        wxglade_tmp_menu.Append(self.menu_compwindow_id, "Computational Window", "", wx.ITEM_NORMAL)
        wxglade_tmp_menu.Append(wx.NewId(), "Preferences", "", wx.ITEM_NORMAL)
        self.MainFrame_menubar.Append(wxglade_tmp_menu, "Options")
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
        self.VFData = wx.grid.Grid(self)#, -1,  size=(255,400))
        self.SchematicPanel = wx.Panel(self, -1)
        self.DistanceSlider = wx.Slider(self, -1, 0, 0, 1000)

        # Add Matplotlib widgets to window
        self.plot = FourPlots(self)

        self.__set_properties()

        #Set Parameter Values in Grid
        self.init_grid_information()

        self.__do_layout()
        
        self.Bind(wx.EVT_MENU, self.menu_computational_window_click, id=self.menu_compwindow_id)

        #self.Bind(wx.EVT_MENU, self.menu_open_click, id=-1)
        #self.Bind(wx.EVT_MENU, self.menu_save_click, id=-1)
        self.Bind(wx.EVT_MENU, self.menu_exit_click, id=self.menu_exit)
        #self.Bind(wx.EVT_MENU, self.menu_exportplots_click, id=-1)
        #self.Bind(wx.EVT_MENU, self.menu_animation_click, id=-1)
        self.Bind(wx.EVT_BUTTON, self.addbutton_click, self.AddButton)
        self.Bind(wx.EVT_BUTTON, self.editbutton_click, self.EditButton)
        self.Bind(wx.EVT_BUTTON, self.removebutton_click, self.RemoveButton)
        self.Bind(wx.EVT_COMMAND_SCROLL, self.distanceslider_change, self.DistanceSlider)
        
        #self.SchematicPanel.Bind(wx.EVT_PAINT, self.repaint_schematic) 
        self.Bind(wx.EVT_PAINT, self.paint_event) 
        self.SchematicPanel.Bind(wx.EVT_LEFT_UP, self.click_schematic) 
    
        # end wxGlade
        
        self.init_calculations()
        self.refresh_interface()
        self.refresh_grid_information()
        self.repaint_schematic()

    def __set_properties(self):
        # begin wxGlade: VFFrame.__set_properties
        self.SetTitle("Virtual Femtolab")
        self.SetSize((900,675))
        self.VFData.CreateGrid(13, 3)
        self.VFData.SetRowLabelSize(0)
        self.VFData.SetColLabelSize(0)
        self.VFData.EnableEditing(0)
        self.VFData.SetColLabelValue(0, "Name")
        self.VFData.SetColLabelValue(1, "Value")
        self.VFData.SetColLabelValue(2, "Units")
        self.VFData.SetSize((255, 150)) #150 is wrong, but it works anyway. 255 is hardcoded - TODO:fix
        self.SchematicPanel.SetSize((890, 130))
        size = self.SchematicPanel.GetClientSize()
        self.buffer = wx.EmptyBitmap(size.width, size.height)
        self.dc = wx.BufferedDC(None, self.buffer )
        #self.dc.SetBackground(wx.Brush('gray'))
        #self.dc.Clear()
        
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: VFFrame.__do_layout
        self.sizer_1 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_2 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_9 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_6 = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer_7 = wx.BoxSizer(wx.VERTICAL)
        self.sizer_8 = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer_7.Add(self.plot, 1, wx.EXPAND, 0)
        self.sizer_8.Add(self.AddButton, 0, wx.ALL|wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 5)
        self.sizer_8.Add(self.EditButton, 0, wx.ALL|wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 5)
        self.sizer_8.Add(self.RemoveButton, 0, wx.ALL|wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 5)
        self.sizer_7.Add(self.sizer_8,0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.sizer_6.Add(self.sizer_7, 3, wx.EXPAND, 0)
        self.sizer_6.Add(self.VFData,0,wx.EXPAND,0)
        self.sizer_2.Add(self.sizer_6, 1, wx.EXPAND, 0)
        self.sizer_9.Add(self.SchematicPanel, 5, wx.EXPAND, 0)
        self.sizer_9.Add(self.DistanceSlider, 1, wx.EXPAND, 0)
        self.sizer_2.Add(self.sizer_9, 0, wx.EXPAND, 0)
        self.sizer_1.Add(self.sizer_2, 1, wx.EXPAND, 0)
        self.SetSizer(self.sizer_1)
        self.Layout()
        # end wxGlade

    def menu_open_click(self, event): # wxGlade: VFFrame.<event_handler>
        print "Event handler `menu_open_click' not implemented!"
        event.Skip()

    def menu_save_click(self, event): # wxGlade: VFFrame.<event_handler>
        print "Event handler `menu_save_click' not implemented!"
        event.Skip()

    def menu_exit_click(self, event): # wxGlade: VFFrame.<event_handler>
        self.Close()

    def menu_exportplots_click(self, event): # wxGlade: VFFrame.<event_handler>
        print "Event handler `menu_exportplots_click' not implemented!"
        event.Skip()

    def menu_animation_click(self, event): # wxGlade: VFFrame.<event_handler>
        print "Event handler `menu_animation_click' not implemented!"
        event.Skip()
        
    def menu_computational_window_click(self, event):
        import edit_computationalwindow
        dialog = edit_computationalwindow.EditComputationalWindow(self)
        dialog.set_info(self.NT,self.deltaT,self.change_computational_window,self.refresh_everything)
        dialog.Show()
        event.Skip()

    def addbutton_click(self, event): # wxGlade: VFFrame.<event_handler>
        event.Skip()
        # Open add dialog
        import add_dialog
        dialog = add_dialog.AddDialog(self)
        dialog.set_info(self.selected, self.propagator.add_element,self.refresh_everything) 
        dialog.Show()

    def editbutton_click(self, event): # wxGlade: VFFrame.<event_handler>
        if(self.selected == 0):
            import edit_initialpulse
            dialog = edit_initialpulse.EditInitialPulse(self)
            dialog.set_info(self.propagator.get_initialPulseBeam(),self.refresh_everything)
            dialog.Show()
            return
    
        name = str(self.propagator.get_elements()[self.selected-1].__class__)

        if(name == 'element_propagation.Element_Propagation'):
            import edit_materialpropagation
            dialog = edit_materialpropagation.EditMaterialPropagation(self)
            dialog.set_info(self.propagator.get_elements()[self.selected-1],self.propagator.get_materials(),self.refresh_everything)
            dialog.Show()
        elif(name == 'element_thinlens.Element_ThinLens'):
            import edit_thinlens
            dialog = edit_thinlens.EditThinLens(self)
            dialog.set_info(self.propagator.get_elements()[self.selected-1],self.refresh_everything)
            dialog.Show()
            
        event.Skip()

    def removebutton_click(self, event): # wxGlade: VFFrame.<event_handler>
        if(self.selected == 0 or self.selected > len(self.propagator.get_elements())):
            event.Skip()
            return
            
        element = self.propagator.get_elements()[self.selected-1]
        dialog = wx.MessageDialog(None, 'Are you sure you want to remove %s ?'%element.name, 'Confirmation', 
                    wx.YES_NO | wx.NO_DEFAULT | wx.ICON_QUESTION)
        if(dialog.ShowModal() == wx.ID_YES):
            self.propagator.remove_element(self.selected-1)
            self.selected = 0
            self.refresh_everything()
        
        event.Skip()
        
        

    def distanceslider_change(self, event): # wxGlade: VFFrame.<event_handler>
        distance = self.DistanceSlider.GetValue()/1000.*self.propagator.get_max_z()
        if(abs(self.distance-distance) > 1e-15):
            self.change_distance(distance)

        event.Skip()
        
        
        
    def init_calculations(self):
        #TODO: move this into a separate config file
        self.NT = 1024
        self.deltaT = 64.e-15
        self.lambdaZero = 8.0e-7
        self.propagator = propagator.Propagator(self.NT,self.deltaT,self.lambdaZero)
        self.propagator.example_pulseBeam()
        self.propagator.example_elements()
        
        self.distance = 0
        self.selected = 0
        
        

    def change_distance(self, distance):
        # Recalculate everything
        
        #keep within bounds
        if(distance < 0):
            distance = 0
        max_distance = self.propagator.get_max_z()
        if(distance > max_distance):
            distance = max_distance
            
        #recalculate
        self.distance = distance 
        self.propagator.change_z(distance)
    
        # Refresh interface/redraw
        self.refresh_interface()
        self.refresh_grid_information()
        self.repaint_schematic()
        
        
        
    def refresh_everything(self):
        self.change_distance(self.distance)
        
        
    
    def refresh_interface(self):
        pulseBeam = self.propagator.get_pulseBeam()
        
        t = pulseBeam.get_t()
        envelope = pulseBeam.get_temporal_envelope()
        electric_field = pulseBeam.get_real_electric_field()
        
        t_phase = t[:]
        temporal_phase = pulseBeam.get_temporal_phase()
        t_phase,temporal_phase = pulseBeam.phase_blank(t_phase,envelope,temporal_phase,1e-2)

        
        freq = pulseBeam.get_frequencies()
        spectrum = pulseBeam.get_spectral_intensity()
        
        freq_phase = freq[:]
        spectral_phase = pulseBeam.get_spectral_phase()
        freq_phase,spectral_phase = pulseBeam.phase_blank(freq_phase,spectrum,spectral_phase,1e-2)
        
        inter_autoco     = pulseBeam.get_interferometric_autoco()
        inten_autoco     = pulseBeam.get_intensiometric_autoco()
        frog, frog_limits = pulseBeam.get_SHGFROG()
        self.plot.redraw(t,envelope,electric_field,t_phase,temporal_phase,freq,spectrum,freq_phase,spectral_phase,inter_autoco,inten_autoco,0,frog,frog_limits)
        
        
        
    def init_grid_information(self):
        self.VFData.SetDefaultCellFont(wx.Font(12, wx.FONTFAMILY_SWISS, wx.NORMAL, wx.FONTWEIGHT_NORMAL))
        self.VFData.SetCellValue(0,0,'Num.Points')
        self.VFData.SetCellValue(1,0,'Delta T')
        self.VFData.SetCellValue(1,2,'fs')
        self.VFData.SetCellValue(2,0,'FWHM (iterative)')
        self.VFData.SetCellValue(2,2,'fs')
        self.VFData.SetCellValue(3,0,'FWHM (gaussian fit)')
        self.VFData.SetCellValue(3,2,'fs')        
        self.VFData.SetCellValue(4,0,'Spectral FWHM')
        self.VFData.SetCellValue(4,2,'nm')
        self.VFData.SetCellValue(5,0,'Beam Spot')
        self.VFData.SetCellValue(5,2,'mm')
        self.VFData.SetCellValue(6,0,'Beam Curvature')
        self.VFData.SetCellValue(6,2,'m')
        self.VFData.SetCellValue(7,0,'Peak Power')
        self.VFData.SetCellValue(7,2,'W')
        self.VFData.SetCellValue(8,0,'Peak Intensity')
        self.VFData.SetCellValue(8,2,'Wm^-2')
        self.VFData.SetCellValue(9,0,'Pulse Energy')
        self.VFData.SetCellValue(9,2,'J')
        self.VFData.SetCellValue(10,0,'Rep. Rate')
        self.VFData.SetCellValue(10,2,'MHz')
        self.VFData.SetCellValue(11,0,'CW Power')
        self.VFData.SetCellValue(11,2,'W')
        self.VFData.SetCellValue(12,0,'z')
        self.VFData.SetCellValue(12,2,'m')

        
    def change_computational_window(self, NT, deltaT):
        self.NT = NT
        self.deltaT = deltaT
        
        self.propagator.change_computational_window(NT, deltaT) 
        self.plot.reset()
        
        
        
    def refresh_grid_information(self):
        pulseBeam = self.propagator.get_pulseBeam()
        self.VFData.SetCellValue(0,1,'%d'%pulseBeam.NT)
        self.VFData.SetCellValue(1,1,'%3.1f'%(pulseBeam.deltaT*1e15))
        fwhm1,fwhm2 = pulseBeam.calculate_fwhm()
        self.VFData.SetCellValue(2,1,'%3.3f'%(fwhm1*1e15))
        self.VFData.SetCellValue(3,1,'%3.3f'%(fwhm2*1e15))
        self.VFData.SetCellValue(4,1,'%3.3f'%(pulseBeam.get_spectral_fwhm()*1e9))
        spot = pulseBeam.get_beam_spot()
        if(spot < 1e-3):
            self.VFData.SetCellValue(5,1,'%3.3e'%(spot*1e6))
            self.VFData.SetCellValue(5,2,'um')
        else:
            self.VFData.SetCellValue(5,1,'%3.3e'%(spot*1e3))
            self.VFData.SetCellValue(5,2,'mm')
        self.VFData.SetCellValue(6,1,'%3.3e'%(pulseBeam.get_beam_curvature()))
        self.VFData.SetCellValue(7,1,'%3.3e'%(pulseBeam.calc_peak_power()))
        self.VFData.SetCellValue(8,1,'%3.3e'%(pulseBeam.calc_peak_intensity()))
        self.VFData.SetCellValue(9,1,'%3.3e'%(pulseBeam.calc_energy()))
        rep_rate = pulseBeam.get_rep_rate()
        if(rep_rate < 1e6):
            self.VFData.SetCellValue(10,1,'%3.3e'%(rep_rate*1e-3))
            self.VFData.SetCellValue(10,2,'KHz')
        else:
            self.VFData.SetCellValue(10,1,'%3.3e'%(rep_rate*1e-6))
            self.VFData.SetCellValue(10,2,'MHz')
        self.VFData.SetCellValue(11,1,'%3.3e'%(pulseBeam.calc_CW_power()))
        self.VFData.SetCellValue(12,1,'%3.3e'%(self.distance))

        self.VFData.Fit()
    
    def click_schematic(self,event):
        
        width = 120

        x = event.GetPosition()[0]
        
        elements = self.propagator.get_elements()
        
        #calculate how much to scroll
        total_width = (width+5)*(len(elements)+1)
        available_width = self.SchematicPanel.GetSizeTuple()[0]-10

        if(available_width > total_width):
            x -= 5
        else:
            x -= 5-self.distance/self.propagator.get_max_z()*(total_width-available_width)
            
        i = int(x/(width+5))
        
        print 'TODO: click_schematic() - 120 needs to become global'

        if(i <= len(self.propagator.get_elements())):
            self.selected = i
            
        self.repaint_schematic()
        
    def paint_event(self,event):
        dc = wx.BufferedPaintDC(self.SchematicPanel, self.buffer)
    
        
    def repaint_schematic(self):
        dc = wx.BufferedPaintDC(self.SchematicPanel, self.buffer)
        #dc.SetBackground(wx.Brush(self.GetBackgroundColour()))
        dc.SetBackground(wx.Brush('#efefef'))
        dc.BeginDrawing()
        dc.Clear()

        width = 120
        height = 100   
        
        elements = self.propagator.get_elements()
        
        #calculate how much to scroll
        total_width = (width+5)*(len(elements)+1)
        available_width = self.SchematicPanel.GetSizeTuple()[0]-10

        if(available_width > total_width):
            x = 5
        else:
            x = 5-self.distance/self.propagator.get_max_z()*(total_width-available_width)
        
        text = '6.6 fs' #TODO: fix
        if(self.selected == 0):
            selected = True
        else:
            selected = False
            
        draw_schematic.draw_initial_pulse(dc,x,width,height,text,selected) 
        x += width+5
        
        
        spots = self.propagator.get_spots()
        spots /= max(abs(array(spots)))
        
        z = 0
        
        for i in xrange(len(elements)):
            spot_in  = spots[2*i+0]
            spot_out = spots[2*i+1]
            if(self.selected == i+1):
                selected = True
            else:
                selected = False

            draw_schematic.draw_element(dc,x,width,height,elements[i],str(elements[i].__class__),selected,spot_in,spot_out)
            
            
            #should we draw the position line in this element?
            if(z <= self.distance and z+elements[i].length > self.distance):
                #calc position to draw the line
                line_x = x+(self.distance-z)/elements[i].length*width
                draw_schematic.draw_line(dc,line_x,height)
            
            z += elements[i].length
            x += width+5

        #event.Skip()
        dc.EndDrawing()
        



# end of class VFFrame

def main():
    app = wx.PySimpleApp(0)
    wx.InitAllImageHandlers()
    MainFrame = VFFrame(None, -1, "")
    app.SetTopWindow(MainFrame)
    MainFrame.Show()
    app.MainLoop()
