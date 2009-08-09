
import numpy
import csv

def saveCSVtable(filename,header,table):
    f = open(filename,'w')
    f.write(header+'\n')
    obj = csv.writer(f)
    obj.writerows(table)
    

#loads a csv file
# the separator can be ; , tab space


class csv_loader:

    def __init__(self,filename):
        
        csvfile = open(filename,"rU")
        try:
            sample = csvfile.read(1024)
            dialect = csv.Sniffer().sniff(sample)
        except:
            print 'Error - Could not load ',filename,' - no delimiter found'
            raise csv.Error
        has_header = csv.Sniffer().has_header(sample)
        csvfile.seek(0)
        

        reader = csv.reader(csvfile, dialect)
        
        if(has_header):
            reader.next()
    
        dataX = []
        dataY = []
        for line in reader:
            dataX.append(float(line[0]))
            dataY.append(float(line[1]))

        self.x = numpy.array(dataX)
        self.y = abs(numpy.array(dataY)) #TODO: fix


    def get_data_for_x_series(self,desiredX): #lambdas must be increasing
        #print desiredX
        #numpy.set_printoptions(threshold=1e50)
        
        #print self.x
        #print self.y
        return numpy.interp(desiredX,self.x,self.y)
        
    def detect_wavelength_scale(self): #detect if in nm or m and make sure y is positive
        if(self.x[0] > 1 and self.x[1] > 1):
            self.x *= 1e-9
        self.y = abs(self.y)
            
    def detect_time_scale(self): #detect if in fs or s
        if(self.x[0] > 1 and self.x[1] > 1):
            self.x *= 1e-15

