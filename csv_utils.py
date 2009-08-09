
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

    def __init__(self,filename,num_columns):
        
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
        for i in xrange(num_columns):
            dataY.append([])
        
        try:
            for line in reader:
                dataX.append(float(line[0]))
                for i in xrange(num_columns):
                    dataY[i].append(float(line[1+i]))
        except:
            raise csv.Error

        self.x = numpy.array(dataX)
        self.y = numpy.array(dataY) 

    def get_data_for_x_series(self,desiredX,column): #lambdas must be increasing
        if(len(self.y) <= column):
            return numpy.zeros(len(desiredX)) #in case we are being asked for an inexistent column, return 0 (probably the phase)
        return numpy.interp(desiredX,self.x,self.y[column])
        
    def detect_wavelength_scale(self): #detect if in nm or m and make sure y is positive
        if(self.x[0] > 1 and self.x[1] > 1):
            self.x *= 1e-9
        self.y = abs(self.y)
            
    def detect_time_scale(self): #detect if in fs or s
        if(self.x[0] > 1 and self.x[1] > 1):
            self.x *= 1e-15

