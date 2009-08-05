
import numpy
import csv

def saveCSVtable(filename,header,table):
    f = open(filename,'w')
    f.write(header+'\n')
    obj = csv.writer(f)
    obj.writerows(table)