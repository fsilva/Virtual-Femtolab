from glob import glob
from py2deb import Py2deb


version="0.8"
changelog='See repository at github for details'

p=Py2deb("virtual-femtolab")
p.author="Francisco Silva"
p.mail="fsilvaportugal@gmail.com"
p.description="""Application to calculate propagation effects in ultrafast pulsed laser beams. Uses Python + wxWidgets + matplotlib."""
p.url = "http://github.com/fsilva/Virtual-Femtolab"
p.depends="python-wxgtk2.8, python-numpy, python-matplotlib, python"
p.license="gpl"
p.section="science"
p.arch="all"

#p["/usr/share/applications"]=["data/virtualfemtolab.desktop|virtualfemtolab.desktop"]
p["/usr/share/virtualfemtolab"]=[i for i in glob("*.py")]
p["/usr/bin"]=["virtualfemtolab",]
#p["/usr/share/doc/fricorder"]=["README","COPYING",]

p.generate(version)

