#!/usr/bin/env python
"""ReadNDVI - For a given month, reads a NDVI csv"""

"""For a prescribed month, reads a CSV of NDVI information, and outputs a 1D list"""

from datetime import datetime as dt
from datetime import timedelta as td
import math
import calendar
import csv
import numpy as np

class d:
    """A class to hold one set of 1D data"""
    
    def __init__(self,val,name,description=""):
        self.unit = ""
        self.name = name
        self.val = val
        self.meta = {}
        if description == "":
            self.description = name
        else:
            self.description = description
    
    def __str__(self):
        return(self.name)
    
    #def __repr__(self): #not sure this is good code
    #    return self.val    
        
    def add_meta(self,key,value):
        self.meta[key]= value
           
class d_all:
    """A class to hold a load of d objects"""
    
    def __init__(self,lat,lon,time=None):
        """Set the core data lists when setting up"""
        self.lat  = lat
        self.lon  = lon
        self.data = {}
        self.meta = {}
        #the time option is optional 
        #(will be used for indiv but not binned  
        self.time = time
    
    def valid(self):
        """Checks all lists are same length"""
        #first check lat, lon, time
        x  = len(self.lat)        
        if x != len(self.lon):
            return False
        if self.time != None:
            if x != len(self.time):
                return False
            
        #now go through every val in the data dictionary
        for key in self.data:
            if x != len(self.data[key].val):
                return False
        
        #if we get here, things are OK
        return True
    
    def list_all_datasets(self):
        list_of_datasets = []
        for key in self.data:
            list_of_datasets.append(self.data[key].description)
        return(list_of_datasets)
        
    def filter_all(self,filter_list):
        """For filer_list the same length as the data, keeps only datapoints where filter_list == True"""
        
        filter_as_np = np.array(filter_list)
        
        print "Filtering"
        print "Start: %i" %len(self.lat)
        all_indexes = range(0,len(self.lat))
        
        #core sets
        self.lat = list(np.array(self.lat)[filter_as_np])
        self.lon = list(np.array(self.lon)[filter_as_np])
        if self.time != None:
            self.time = list(np.array(self.time)[filter_as_np])
        
        #data sets
        for key in self.data:
            #print key    
            self.data[key].val = list(np.array(self.data[key].val)[filter_as_np])
        
        print "End: %i" %len(self.lat)
            
    def add_data(self,new_data):
        #adds a new d object, or multiple d objects
        if type(new_data) == list:
            for this_new_data in new_data:
                self.data[this_new_data.name] = this_new_data
        else:
            #if it's just one
            self.data[new_data.name] = new_data       

    def make_menu(self,num_only=False):
        """Create menu options of all the datasets in data"""
        
        menu_text = []
        answers_dict = {}
        iii = 0
        for key in self.data:
            if isinstance(self.data[key].val[0],Number) or not(num_only):
                menu_text.append([str(iii),self.data[key].description])
                answers_dict[str(iii)] = key
                iii += 1
        
        return(menu_text,answers_dict)

def add_months_dt(sourcedate,months):
    month = sourcedate.month - 1 + months
    year = int(sourcedate.year + month / 12 )
    month = month % 12 + 1
    day = min(sourcedate.day,calendar.monthrange(year,month)[1])
    return dt(year,month,day,sourcedate.hour,sourcedate.minute,sourcedate.second)


def readLANDCOVER(landcoverfolder,xres,yres):
    #opens and reads LANDCOVER.xyz, a space-separated, column-arranged file
    
    filename = '%s/2011LANDCOVER.xyz' %landcoverfolder
    print filename
    f = open(filename,'r')
    
    lats=[]
    lons=[]
    LC=[]
    for line in f:
        line = line.strip()
        columns = line.split() #use whitespace as delimiter
        #remove floating point errors in 
        lons.append(int(float(columns[0])/xres) * xres)
        lats.append(int(float(columns[1])/yres) * yres)
        LC.append(int(columns[2]))
    f.close()
    return(lats,lons,LC)


def readNDVI(NDVIfolder,year,month,xres,yres,compass,datatype="NDVI"):
    
    if datatype=="NDVI":
        #data from http://neo.sci.gsfc.nasa.gov/view.php?datasetId=MOD13A2_M_N
        #work out the file name to be read
        filename = '%s/MOD13A2_M_NDVI_%i-%02d-01_rgb_3600x1800.CSV' % (NDVIfolder,year,month)
    elif datatype=="LAI":
        #work out the file name to be read
        filename = '%s/MOD15A2_M_LAI_%i-%02d-01_rgb_3600x1800.CSV' % (NDVIfolder,year,month)
    elif datatype == "Fire count":
        #work out the file name to be read
        filename = '%s/MOD14A1_M_FIRE_%i-%02d-01_rgb_3600x1800.CSV' % (NDVIfolder,year,month)
    elif datatype == "Land classification":
        (lats,lons,NDVI) = readLANDCOVER(NDVIfolder,xres,yres) #this one is different
    elif datatype == "Population density":
        filename = '%s/SEDAC_POP_2000-01-01_rgb_3600x1800.CSV' %NDVIfolder
        
    if datatype != "Land classification":  
                
        print filename
        
        f = open(filename,'r')
        
        numberofcols = int(360./xres)
        
        latcounter = 90. - 0.5*yres
        lats = []
        lons = []
        NDVI = []
        l = 0
        for line in f:
            l+=1
            line = line.strip() #get rid of the newline character
            columns = line.split(",") #split by columns, comma separated
            for i in range(0,numberofcols):
                try:
                    #"{0:.3f}".format(x) keeps things at 3dp, avoids proiferation of floating point errors.
                    lats.append(float("{0:.3f}".format(latcounter)))
                    lons.append(float("{0:.3f}".format(-180.+(xres*(i+0.5)))))
                    NDVI.append(float(columns[i]))
                except IndexError:
                    print "Index error at Line number %i, column number %i" %(l,i)
                    sys.exit() 
            latcounter -= yres
    # End if datatype != "Land classification"
        
    north = compass.n
    south = compass.s
    east  = compass.e
    west  = compass.w
    
    cut_lats = [lats[i] for i in range(0,len(NDVI)) if (south < lats[i] < north) and (west < lons[i] < east)]
    cut_lons = [lons[i] for i in range(0,len(NDVI)) if (south < lats[i] < north) and (west < lons[i] < east)] 
    cut_NDVI = [NDVI[i] for i in range(0,len(NDVI)) if (south < lats[i] < north) and (west < lons[i] < east)]
    
    lats = cut_lats
    lons = cut_lons
    NDVI = cut_NDVI
    del cut_lats,cut_lons,cut_NDVI
    
    NDVI_sea0 = []
    for i in range(0,len(NDVI)):
        if NDVI[i] == 99999.:
            NDVI_sea0.append( 0. )
        elif datatype == "Fire count" and NDVI[i] == 0.1: #for fires, 0.1 means 0.
            NDVI_sea0.append( 0. )
        else:
            NDVI_sea0.append( NDVI[i] )
    
    return(lats,lons,NDVI_sea0)
    
def CSV_months(startdate,enddate,NDVIfolder,xres,yres,compass,CSV_name):
    """Returns NDVIs for the months between two dates, plus time & location info"""
    #should be generic to all CSV form data now.
        
    lats_all = []
    lons_all = []
    NDVI_all = []
    year_all = []
    month_all= []
    
    if CSV_name != "Land classification":  
        monthcounter = startdate
        while monthcounter <= enddate:
            year  = monthcounter.year
            month = monthcounter.month
            print "loading %s data for %i-%02d" %(CSV_name,year,month)
            try:
                (lats_thismonth,lons_thismonth,NDVI_thismonth) = readNDVI(NDVIfolder,year,month,xres,yres,compass,datatype=CSV_name)
            except IOError: #if file doesn't exist
                print "No file for %i-%i" %(year,month)
                monthcounter = add_months_dt(monthcounter,1)
                continue
                
            lenofmonth = len(NDVI_thismonth)
            year_thismonth = [year]*lenofmonth
            month_thismonth = [month]*lenofmonth
            
            lats_all.extend(lats_thismonth)
            lons_all.extend(lons_thismonth)
            NDVI_all.extend(NDVI_thismonth)
            year_all.extend(year_thismonth)
            month_all.extend(month_thismonth)
            del lats_thismonth,lons_thismonth,NDVI_thismonth,year_thismonth,month_thismonth
            #print year_all
            print len(NDVI_all)
            monthcounter = add_months_dt(monthcounter,1)
    else: #if it IS land class
        (lats_all,lons_all,NDVI_all) = readNDVI(NDVIfolder,2011,1,xres,yres,compass,datatype=CSV_name) #2011 and 1 are not used
        year_all = None
        month_all = None
        
    return(lats_all,lons_all,NDVI_all,year_all,month_all)
        
def readNDVI_v2(NDVIfolder,year,month,xres,yres,compass,datatype="NDVI"):
    if datatype == "NDVI":
        filename = '%s/MOD13A2_M_NDVI_%i-%02d-01_rgb_3600x1800.CSV' % (NDVIfolder,year,month)
    elif datatype == "Fire count":
        filename = '%s/MOD14A1_M_FIRE_%i-%02d-01_rgb_3600x1800.CSV' % (NDVIfolder,year,month)
    f = open(filename,"rb")
    data = list(csv.reader(f, delimiter=','))
    f.close()
    print "Full Shape.."
    print np.array(data).shape
    
    north = compass.n
    south = compass.s
    east  = compass.e
    west  = compass.w
    
    #north = 89.95
    #south = -89.95
    #east  = 179.95
    #west  = -179.95
    
    
    north_j = int((90 - north) / yres)
    south_j = int((90 - south) / yres)+1
    east_i  = int((180 + east) / xres)
    west_i  = int((180 + west) / xres)+1
    
    print "north_j = %i" %north_j
    print "south_j = %i" %south_j
    print "east_i = %i" %east_i
    print "west_i = %i" %west_i
                
    
    #slice
    data = [row[west_i:east_i] for row in data[north_j:south_j]]
    
    print "Cut Shape.."
    print np.array(data).shape
    
    return(data)

def NDVI_months_v2(startdate,enddate,NDVIfolder,xres,yres,compass,datatype="NDVI"):

    dateclock = startdate
    m = 0
    NDVI_3D = []
    while dateclock < enddate:
        year  = dateclock.year
        month = dateclock.month
        print "loading %s data for %i-%02d" %(datatype,year,month)
        NDVI_3D.append(readNDVI_v2(NDVIfolder,year,month,xres,yres,compass,datatype=datatype))
        dateclock = add_months_dt(dateclock,1)    
    return(NDVI_3D)
    
def associate_NDVI_v2(ida,NDVI_3D,xres,yres,compass,startdate,datatype="NDVI"):

    north=compass.n
    south=compass.s
    east=compass.e
    west=compass.w
    
    #north = 89.95
    #south = -89.95
    #east  = 179.95
    #west  = -179.95
    
    print "Associating with individual observations"
    
    print "3D Shape.."
    print np.array(NDVI_3D).shape  
    
    NDVI = []
    
    max_i = int((east-west)/xres)
    max_j = int((north-south)/yres)
    
    if datatype == "Fire screen":
        NDVI_3D = np.array(NDVI_3D)
        filter_range = float(raw_input("Minimum acceptable distance from fire in degrees?-->"))
    
    for k in range(0,len(ida.lat)):

        i =int((ida.lon[k]-west)/xres)
        j =int((north-ida.lat[k])/yres)
        m = 12*(ida.time[k].year - startdate.year) + (ida.time[k].month - startdate.month)
        
        if datatype != "Fire screen":
        
            if k % 1000 == 0:
                print "%i of %i" %(k,len(ida.lat))
                print "lat=%g, lon=%g" %(ida.lat[k],ida.lon[k])
                print "j  =%i, i  =%i" %(j,i)
                print "%s = %g" %(datatype,float(NDVI_3D[m][j][i]))
            try:
                NDVI.append(float(NDVI_3D[m][j][i]))
            except IndexError:
                print "Out of range, assigning zero"
                print "lat=%g, lon=%g" %(ida.lat[k],ida.lon[k])
                print "m=%i, j=%i, i=%i" %(m,j,i)            
                NDVI.append(0.)
                print "================"
        
        else:
            j_low = j-int(filter_range/yres)
            j_high = j+int(filter_range/yres)
            i_low = i-int(filter_range/xres)
            i_high = i+int(filter_range/xres)
            #print(j_low)
            #type(j_low)
            #subset = np.array([float(i) for i in list(np.flatten(np.array(NDVI_3D[m][j_low:j_high,i_low:i_high])))])
            subset = np.array(NDVI_3D[m][j_low:j_high,i_low:i_high]).astype(float)
            #print subset
            
            #if k % 1000 == 0:
            #    print "%i of %i" %(k,len(ida.lat))
            #    print "lat=%g, lon=%g" %(ida.lat[k],ida.lon[k])
            #    print "j  =%i, i  =%i" %(j,i)
            #    print "%s = %g" %(datatype,float(np.max(subset)))
            try:
                NDVI.append(float(np.max(subset)))
            except IndexError:
                print "Out of range, assigning 99999"
                print "lat=%g, lon=%g" %(ida.lat[k],ida.lon[k])
                print "m=%i, j=%i, i=%i" %(m,j,i)            
                NDVI.append(99999.)
                print "================"
            
    
    #seas from 99999 -> 0
    NDVI_sea0 = []
    if datatype == "NDVI":
        for i in range(0,len(NDVI)):
            if NDVI[i] == 99999.:
                NDVI_sea0.append( 0. )
            else:
                NDVI_sea0.append( float(NDVI[i]) ) #also, make sure float
        NDVI_d = d(NDVI_sea0,"NDVI",description="Normalised Diffusive Vegetation Index")
    elif datatype == "Fire count":
        for i in range(0,len(NDVI)):
            if NDVI[i] == 99999.:
                NDVI_sea0.append( 0. )
            elif NDVI[i] == 0.1: #for fires, 0.1 = 0
                NDVI_sea0.append( 0. )                
            else:
                NDVI_sea0.append( float(NDVI[i]) ) #also, make sure float    
        NDVI_d = d(NDVI_sea0,"Fire count",description="MODIS fire count (fire pixels/1000km2/day)")
    elif datatype == "Fire screen":
        for i in range(0,len(NDVI)):
            if NDVI[i] == 99999.:
                NDVI_sea0.append( 1. )
            elif NDVI[i] == 0.1: #NO FIRES -- for fires, 0.1 = 0
                NDVI_sea0.append( 0. )                
            else: #anything else
                NDVI_sea0.append( 1. ) #also, make sure float    
        NDVI_d = d(NDVI_sea0,"Fire screen",description="Fire in month in %f degree range?"%filter_range)    
    
    ida.add_data(NDVI_d)
