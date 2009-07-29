
import wx
import element_propagation
import element_thinlens

from numpy import cos, exp

def draw_initial_pulse(dc,x,width,height,text):
    #draw box
    dc.SetPen(wx.Pen('brown', 2))
    dc.SetBrush(wx.Brush('lightgray'))
    dc.DrawRoundedRectangle(x,5,width,height,2)

    #draw beam propagation
    
    #draw initial pulse

    lastx = x0 = x+5
    lasty = y0 = height*0.5
    lines1 = []
    lines2 = []
    N = 50
    f1,f2 = 1,10
    
    for i in xrange(N):
        x = x0 + width*0.5*i/N
        y = y0 - height*0.25*(cos((i-N/2.)*f1)*exp(-(i-N/2.)**2/f2**2))
        
        lines1.append((lastx,lasty,x,y))
        lastx = x
        lasty = y
        
    lastx = x0
    lasty = y0    
    for i in xrange(N):
            x = x0 + width*0.5*i/N
            y = y0 - height*0.25*exp(-(i-N/2.)**2/f2**2)

            lines2.append((lastx,lasty,x,y))
            lastx = x
            lasty = y
         
    dc.SetPen(wx.Pen('blue', 1))     
    dc.DrawLineList(lines1)
    dc.SetPen(wx.Pen('red', 1))
    dc.DrawLineList(lines2)
    
    #draw text
    dc.DrawLabel(text,wx.Rect(x,height*0.5,width,height*0.5),wx.ALIGN_CENTER)
    print text


def draw_element(dc,x,width,height,element,type):
    pass


