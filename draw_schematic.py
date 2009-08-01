
import wx
import element_propagation
import element_thinlens

from numpy import cos, exp


def helper_draw_box(dc,x,width,height,selected):
    #draw box
    if(selected == True):
        dc.SetPen(wx.Pen('black', 2))
    else:
        dc.SetPen(wx.Pen('gray', 2))
    dc.SetBrush(wx.Brush('gray',style=wx.TRANSPARENT))
    dc.DrawRoundedRectangle(x,5,width,height,2)
    
    
def helper_draw_text(dc,x,width,height,text_top,text_bottom):
    #draw text
    dc.DrawLabel(text_top,wx.Rect(x,height*0.1,width,height*0.1),wx.ALIGN_CENTER)
    dc.DrawLabel(text_bottom,wx.Rect(x,height*0.85,width,height*0.15),wx.ALIGN_CENTER)

def helper_draw_beam(dc,x,width,height,size_in,size_out):
    dc.SetPen(wx.Pen('red', 1))
    dc.SetBrush(wx.Brush('#ff4f4f'))
    points = [wx.Point(x,height*(0.5+size_in*0.2)),wx.Point(x,height*(0.5-size_in*0.2)), \
              wx.Point(x+width,height*(0.5-size_out*0.2)),wx.Point(x+width,height*(0.5+size_out*0.2))]
    dc.DrawPolygon(points)




def draw_initial_pulse(dc,x,width,height,text,selected):
    helper_draw_box(dc,x,width,height,selected)
    
    #draw initial pulse

    lastx = x0 = x+10
    lasty = y0 = height*0.55
    lines1 = []
    lines2 = []
    N = 50
    f1,f2 = 1,8.
    
    for i in xrange(N):
        px = x0 + (width-20)*i/N
        py = y0 - height*0.25*(cos((i-N/2.)*f1)*exp(-(i-N/2.)**2/f2**2))
        
        lines1.append((lastx,lasty,px,py))
        lastx = px
        lasty = py
        
    lastx = x0
    lasty = y0    
    for i in xrange(N):
            px = x0 + (width-20)*i/N
            py = y0 - height*0.25*exp(-(i-N/2.)**2/f2**2)

            lines2.append((lastx,lasty,px,py))
            lastx = px
            lasty = py
         
    dc.SetPen(wx.Pen('blue', 1))     
    dc.DrawLineList(lines1)
    dc.SetPen(wx.Pen('red', 1))
    dc.DrawLineList(lines2)
    
    helper_draw_text(dc,x,width,height,'Initial Pulse',text)




def draw_element(dc,x,width,height,element,type,selected,beam_size_in,beam_size_out):
    if(type == 'element_propagation.Element_Propagation'):
        draw_propagation(dc,x,width,height,element,selected,beam_size_in,beam_size_out)
    elif(type == 'element_thinlens.Element_ThinLens'):
        draw_thinlens(dc,x,width,height,element,selected,beam_size_in,beam_size_out)
    
    
def draw_propagation(dc,x,width,height,element,selected,beam_size_in,beam_size_out):

    helper_draw_box(dc,x,width,height,selected)

    #draw beam propagation
    helper_draw_beam(dc,x,width,height,beam_size_in,beam_size_out)

    #draw text
    if(element.length >= 1):
        text = 'L= %1.1e m'%element.length
    elif(element.length >= 1e-2):
        text = 'L= %1.1e cm'%(element.length*100)
    elif(element.length >= 1e-4):
        text = 'L= %1.1e mm'%(element.length*1000)
    else:
        text = 'L= %1.1e um'%(element.length*1e6)
    helper_draw_text(dc,x,width,height,element.name,text) 



def draw_thinlens(dc,x,width,height,element,selected,beam_size_in,beam_size_out):

    helper_draw_box(dc,x,width,height,selected)

    #draw beam propagation
    helper_draw_beam(dc,x,width,height,beam_size_in,beam_size_out)
    
    #draw thin lens
    x0 = x+width*0.5
    y0 = height*0.3
    y1 = height*0.7
    l = 5
    
    if(element.f > 0):
        lines = [(x0,y0,x0,y1),(x0,y0,x0-l,y0+l),(x0,y0,x0+l,y0+l),(x0,y1,x0+l,y1-l),(x0,y1,x0-l,y1-l)]
    else:
        lines = [(x0,y0,x0,y1),(x0,y0,x0-l,y0-l),(x0,y0,x0+l,y0-l),(x0,y1,x0+l,y1+l),(x0,y1,x0-l,y1+l)]
    dc.SetPen(wx.Pen('black', 1))     
    dc.DrawLineList(lines)

    #draw text
    helper_draw_text(dc,x,width,height,element.name,'f= %1.1e m'%element.f)


