#!/usr/bin/env python
# -*- coding: us-ascii -*-
# generated by wxGlade 0.6.3 on Sat Aug  1 22:40:27 2009

import wx

# begin wxGlade: extracode
# end wxGlade



class AddDialog(wx.Dialog):
    def __init__(self, *args, **kwds):
        # begin wxGlade: AddDialog.__init__
        kwds["style"] = wx.DEFAULT_DIALOG_STYLE
        wx.Dialog.__init__(self, *args, **kwds)
        self.element = wx.StaticText(self, -1, "Element Type:", style=wx.ALIGN_CENTRE)
        self.list_box_1 = wx.ListBox(self, -1, choices=["Material Propagation", "Thin Lens"], style=wx.LB_SINGLE)
        self.button_1 = wx.Button(self, wx.ID_ADD, "")

        self.__set_properties()
        self.__do_layout()

        self.Bind(wx.EVT_BUTTON, self.add_clicked, self.button_1)
        # end wxGlade
        
        
    def __set_properties(self):
        # begin wxGlade: AddDialog.__set_properties
        self.SetTitle("Add New Element")
        self.list_box_1.SetMinSize((234, 139))
        self.list_box_1.SetSelection(0)
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: AddDialog.__do_layout
        sizer_1 = wx.BoxSizer(wx.VERTICAL)
        sizer_2 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_2.Add(self.element, 2, wx.ALL|wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_2.Add(self.list_box_1, 3, wx.EXPAND, 0)
        sizer_1.Add(sizer_2, 1, wx.EXPAND, 0)
        sizer_1.Add(self.button_1, 0, wx.ALL|wx.ALIGN_RIGHT|wx.SHAPED, 5)
        self.SetSizer(sizer_1)
        sizer_1.Fit(self)
        self.Layout()
        # end wxGlade
        
    def set_info(self, position, add_function, refresh_function):
        self.position = position
        self.add_function = add_function
        self.refresh_function = refresh_function

    def add_clicked(self, event): # wxGlade: AddDialog.<event_handler>
        event.Skip()
        
        x = self.list_box_1.GetSelection()
        
        if(x == 0):
            self.add_function('Material Propagation',self.position)
        elif(x == 1):
            self.add_function('Thin Lens',self.position)
        else:
            #do nothing
            pass 
        self.Close()
        
        self.refresh_function()

# end of class AddDialog