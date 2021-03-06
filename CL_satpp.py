#!/usr/bin/env python
"""CL_satpp.py -- A command-line interface for interogating the sat_pp dataset"""

from datetime import datetime as dt
from datetime import timedelta as td
from datetime import date as date
import sys
try:
    from matplotlib.patches import Polygon
except ImportError:
    print "Need to be in a (virtual) environment with several modules pip installed to work!"
    sys.exit()
    
import matplotlib.pyplot as plt
import matplotlib.colorbar as colorbar
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np
import math
import os
from geopy.distance import great_circle
import scipy.stats as scistats
from numbers import Number
from os import listdir
from os.path import isfile, join
from bpch import bpch
from dateutil.relativedelta import relativedelta
import copy
import pickle
import cPickle
from WYIBF import *
from readNDVI import *
import brewer2mpl
from netCDF4 import Dataset

def lat_str(y):
    """Returns string for latitude"""
    deg = u"\u00b0"
    if y >= 0: #north
        return "%.2f%sN" %(abs(y),deg)
    else: #south
        return "%.2f%sS" %(abs(y),deg)

def lon_str(x):
    """Returns string for latitude"""
    deg = u"\u00b0"
    if x >= 0: #east
        return "%.2f%sE" %(abs(x),deg)
    else: #west
        return "%.2f%sW" %(abs(x),deg)        

class box:
    """A definition of a box"""
    
    def __init__(self,n,s,e,w):
        self.n = n
        self.s = s
        self.e = e
        self.w = w
           
    def __str__(self):
        return lat_str(self.n)+" - "+lat_str(self.s)+" latitude; "+lon_str(self.e)+" - "+lon_str(self.w)+" longitude"
        
    def valid(self):
        """Checks north is north of south and east is east of west"""
        if self.n>self.s and self.e>self.w:
            return True
        else:
            print "Not a valid spatial extent"
            return False
    
    def within(self,lat_test,lon_test):
        """Checks if lat lon pair are in the box"""
        if type(lat_test) == list:
            #if we have lists, get a corresponding list output.
            output = []
            for i in range(0,len(lat_test)):
                if self.s <= lat_test[i] <= self.n and \
                       self.w <= lon_test[i] <= self.e:
                    #print "list item true"
                    output.append(True)
                else:
                    #print "list item false"
                    output.append(False)
            return(output)
        else: #if scalar
            if self.s <= lat_test <= self.n and \
                   self.w <= lon_test <= self.e:
                return(True)
            else:
                return(False)   

class circle:
    """An object describing a circlular area"""
    
    def __init__(self,lat,lon,radius):
        self.lat = lat
        self.lon = lon
        self.radius = radius
        
    def __str__(self):
        return "Circle centred at "+lat_str(self.lat)+", "+lon_str(self.lon)+", with radius "+str(self.radius)+"km."

    def within(self,lat_test,lon_test):
       """Checks in lat lon pair are within the circle"""
       if type(lat_test) == list:
            #if we have lists, get a corresponding list output.
            output = []
            for i in range(0,len(lat_test)):
                if great_circle((lat_test[i],lon_test[i]),
                                (self.lat,   self.lon   )).km \
                                <= self.radius :
                    #print "list item true"
                    output.append(True)
                else:
                    #print "list item false"
                    output.append(False)
            return(output)
       else: #if scalar
           if great_circle((lat_test,   lon_test   ),
                           (self.lat,   self.lon   )).km \
                           <= self.radius :
               return(True)
           else:
               return(False)    
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

def build_ida(ULN,lat,lon,time,geos_VC,sat_VC,sat_DVC,AMF):
    """Creates an d_all object containing ind objects of basic data"""
    ida = d_all(lat,lon,time=time)
    ida.data['ULN'] = d(ULN,'ULN','Unique Line Number')
    ida.data['geos_VC'] = d(geos_VC,'geos_VC','GEOS Chem modelled vertical column')
    ida.data['geos_VC'].unit = 'molec/cm2'
    ida.data['sat_VC'] = d(sat_VC,'sat_VC','Observed satellite column')
    ida.data['sat_VC'].unit = 'molec/cm2'
    ida.data['sat_DVC'] = d(sat_DVC,'sat_DVC','Uncertainty in observed satellite column')
    ida.data['sat_DVC'].unit = 'molec/cm2'
    ida.data['AMF'] = d(AMF,'AMF','Air mass factor')
    return(ida)


# do we want a build_bda def here? 
   
class geos_data():
    """A class for objects of data imported from geos chem output"""
    
    def __init__(self):
       self.data = [] #the data
       self.time = [] #time of each block of data
       self.time_bounds = [] #upper and lower bounds of time-relevancy
       self.lat = [] #the latitudes
       self.lon = [] #the longitudes
       
       
    def add_data(self, new_data, new_time, new_time_bounds=None):
        self.data.append(new_data)
        self.time.append(new_time)
        self.time_bounds.append(new_time_bounds)
    
    def set_name(self, name): #to name the dataset (name in geos output)
        self.name = name
        
    def set_human_name(self, human_name): #to name the dataset (human name)
        self.human_name = human_name
        
    def set_unit(self, unit): #to set a unit for the dataset
        self.unit = unit
        
    def set_lat_lon(self,lat,lon): #to set latitudes and longitudes
        self.lat = lat
        self.lat_spacing = abs(lat[2]-lat[1])
        self.lon = lon
        self.lon_spacing = abs(lon[2]-lon[1])
    


def clearscreen():
    """Clears the screen. Checks the OS to deliver the correct command"""
    os.system('cls' if os.name == 'nt' else 'clear')


def clear_screen():
    """Runs clearscreen()"""
    clearscreen()

def basic_menu(title,menu_options,quit_option=True):
    """Loads a basic menu to screen and returns the output"""
    clearscreen() #clears the screen
    menu_letters = [] #an initially empty list for the user choice options
    print title
    for i in range(0,len(menu_options)):
        #add option to menu_letters
        menu_letters.append( menu_options[i][0].upper() )
        print "[" + menu_options[i][0] + "] " + menu_options[i][1]
    if quit_option: #add an option to leave the menu (go up a level)
        menu_letters.append( "Z" )
        print "[Z] Exit this menu"
    choice = ""
    while choice not in menu_letters : #only allow choices in menu_letters
        choice = raw_input("-->").upper()
    return(choice)
    

def initial():
    choice = ""
    while not(choice.upper() in ["Y","N"]):
        os.system('clear')
        print "Welcome to the sat_pp visualisation user interface"
        print "Load default options and datasets? Y/N"
        choice = raw_input("-->")
    return(choice.upper())

#Simple menus replaced with calls to basic menu.

#def top_level_menu(current_pickle,current_emfiles):
#    choice = ""
#    while not(choice.upper() in ["1","2","G","P","R","E","X","Z"]):
#        os.system('clear')
#        print "TOP LEVEL MENU"
#        print "Current pickles loaded: %s" %current_pickle
#        print "Current emissions files are: %sXXX%s" %(current_emfiles[0],current_emfiles[1])
#        print "[1] Use binned data"
#        print "[2] Use individual observations data"
#        print "[G] Geographically select data to use"
#        print "[P] Change pickle"
#        print "[R] Reload pickle"
#        print "[E] Change current emissions files"   
#        print "[X] Inspect and change global options"
#        print "[Z] Quit"
#        choice = raw_input("-->")
#    return(choice.upper())
#
#def use_binned_data_menu():
#    choice = ""
#    while not(choice.upper() in ["1","2","Z"]):
#        os.system('clear')
#        print "USING BINNED DATA"
#        print "[1] Plot gridded data on map"
#        print "[2] Compute basic statistics"
#        print "[Z] Return to top level menu" 
#        choice = raw_input("-->")
#    return(choice.upper())
#
#def use_indiv_data_menu():
#    choice = ""
#    while choice.upper() not in ["1","2","3","4","Z"]:
#        os.system('clear')
#        print "USING INDIVIDUAL DATA"
#        print "[1] Simple dots-on-map"
#        print "[2] Compute basic statistics"
#        print "[3] Compute timeseries statistics"
#        print "[4] Compare two datasets"
#        print "[Z] Return to previous menu"
#        choice = raw_input("-->")
#    return(choice.upper())        
#    
#def binned_map_menu():
#    choice = ""
#    while not(choice.upper() in ["1","2","3","4","5","6","Z"]):
#        os.system('clear')
#        print "Plotting map of binned data"
#        print "Loaded pickle and global options will be used"
#        print "Print which variable?"
#        print "[1] Mean satellite observations"  
#        print "[2] Calculated uncertainty in satellite obsservations" 
#        print "[3] Mean model values" 
#        print "[4] Standard deviation of model values" 
#        print "[5] Mean NDVI value" 
#        print "[6] Mean deviation of model from observations (sigmas)" 
#        print "[Z] Return to previous menu"
#        choice = raw_input("-->")
#    return(choice.upper())
#
#def dots_on_map_menu():
#    choice = ""
#    while not(choice.upper() in ["1","2","3","Z"]):
#        os.system('clear')
#        print "Plotting map of binned data"
#        print "Loaded pickle and global options will be used"
#        print "Print which variable?"
#        print "[1] Mean satellite observations"  
#        print "[2] Calculated uncertainty in satellite obsservations" 
#        print "[3] Mean model values" 
#        print "[Z] Return to previous menu"
#        choice = raw_input("-->")
#    return(choice.upper())

def basic_statistics_menu(data_type,quickinput=["0","0"]):
    """A menu for selecting a basic statistical operator for a given dataset"""
    #data_type is binned or indiv    
    while True: #this loop is broken by return
        
        dataset_choice = ""
        if data_type == "binned":
            while dataset_choice not in ["1","2","3","4","5","6","Z"]:
                os.system('clear')
                print "BASIC STATISTICS MENU"
                print "Select dataset to run statisics on:"
                print "[1] Mean satellite observations"  
                print "[2] Calculated uncertainty in satellite obsservations" 
                print "[3] Mean model values" 
                print "[4] Standard deviation of model values" 
                print "[5] Mean NDVI value" 
                print "[6] Mean deviation of model from observations (sigmas)"
                print "[Z] Return to previous menu"
                if quickinput[0] == "0":
                    dataset_choice = raw_input("-->")
                else:
                    dataset_choice = quickinput[0]
                if dataset_choice == "Z": #if returning to previous menu
                    return (["Z","Z"])
        elif data_type == "indiv":
            while dataset_choice not in ["1","2","3","Z"]:
                os.system('clear')
                print "BASIC STATISTICS MENU"
                print "Select dataset to run statisics on:"
                print "[1] Satellite observations"  
                print "[2] Uncertainty in satellite obsservations" 
                print "[3] Model values" 
                print "[Z] Return to previous menu"
                if quickinput[0] == "0":
                    dataset_choice = raw_input("-->")
                else:
                    dataset_choice = quickinput[0]
                if dataset_choice == "Z": #if returning to previous menu
                    return (["Z","Z"])
                    
        stat_choice = ""
        while stat_choice not in ["1","2","3","Z"]:
            os.system('clear')                
            print "BASIC STATISTICS MENU"
            print "Select statistic to process:"
            print "[1] mean"
            print "[2] standard deviation"
            print "[3] count"
            print "[Z] return to previous menu"
            if quickinput[1] == "0":
                stat_choice = raw_input("-->")
            else:
                stat_choice = quickinput[1]
        if stat_choice in ["1","2","3"]:
            return([dataset_choice,stat_choice])
        else:
            pass #Z will return to previous menu       

def timeseries_statistics_menu(quickinput=["0","0"]):
    """A menu for selecting statistical operator and dataset to compute over timeseries"""
    #data_type is binned or indiv    
    while True: #this loop is broken by return

        dataset_choice = ""
        while dataset_choice not in ["1","2","3","Z"]:
            os.system('clear')
            print "TIMESERIES STATISTICS MENU"
            print "Select dataset to run statisics on:"
            print "[1] Satellite observations"
            print "[2] Uncertainty in satellite obsservations"
            print "[3] Model values"
            print "[Z] Return to previous menu"
            if quickinput[0] == "0":
                dataset_choice = raw_input("-->")
            else:
                dataset_choice = quickinput[0]
            if dataset_choice == "Z": #if returning to previous menu
                return (["Z","Z",0])

        stat_choice = ""
        while stat_choice not in ["1","2","3","Z"]:
            os.system('clear')
            print "TIMESERIES STATISTICS MENU"
            print "Select statistic to process:"
            print "[1] mean"
            print "[2] standard deviation"
            print "[3] count"
            print "[Z] return to previous menu"
            if quickinput[1] == "0":
                stat_choice = raw_input("-->")
            else:
                stat_choice = quickinput[1]
            if stat_choice != "Z":
                print "Lastly, number of days per time step:"
                timestep = float(raw_input("-->"))
                return(dataset_choice,stat_choice,timestep)

def frange(x, y, jump=1.0):
  while x < y:
    yield x
    return
    x += jump
  

def calc_statistic(dataset,stat_choice):
    if stat_choice == "1": #mean
        return np.nanmean(dataset)
    elif stat_choice == "2": #standard deviation
        return np.nanstd(dataset)
    elif stat_choice == "3": #count
        return np.count_nonzero(~np.isnan(dataset))
    elif stat_choice == "4": #percentiles
        return np.nanpercentile(dataset,range(0,101))
        
    
def map_preplot_menu(dataset_name,title="",vmin=0.,vmax=3.e16,unit="molec/cm2",do_labels=False):
    choice = "xx"    
    if title == "":
        title = dataset_name
    save_filename = title + ".png"
    while not(choice.upper() in ["P","S","Z"]): #menu-leaving options       
        os.system('clear')
        print "Preparing to map data"
        print "Dataset to plot  = %s" %dataset_name
        print "Title for figure = %s" %title
        print "Colourbar minimum = %g" %vmin
        print "Colourbar maximum = %g" %vmax
        print "Units text        = %s" %unit
        if do_labels:
            print "Label lat+lon     = True"
        else:
            print "Label lat+lon     = False"
        print "OPTIONS:"
        print "[T] Change title"
        print "[C] Change colourbar min/max"
        print "[U] Change units text"
        print "[L] Switch lat+lon label option"
        print "[P] Plot figure on screen"
        print "[S] Save figure."
        print "[Z] Return to previous menu"
        choice = raw_input("-->")
        if choice.upper() == "T":
            title = change_var(title,"Figure title")
        elif choice.upper() == "C":
            vmin = change_var(vmin,"Colourbar minimum")
            vmax = change_var(vmax,"Colourbar maximum")
        elif choice.upper() == "U":
            unit = change_var(unit,"Units text")
        
        if choice.upper() == "S":
            save_filename = title + ".png"
            save_filename = change_var(save_filename,"Image filename to save (include extension)")
        if choice.upper() == "L":
            do_labels = not(do_labels)
    return(choice.upper(),title,vmin,vmax,unit,save_filename,do_labels)
    
def geo_select_menu(type_geo_select,geo_selection):

    #type_geo_select can be either:
    # latlon
    choice = ""
    while choice not in ["1","2","3","4","Z"]    :
        os.system('clear')
        print "GEOGRAPHIC SELECTION OPTION"
        if type_geo_select == "latlon":
            print "Currently, geographic selection done by rectangular area,"
            print "%s to %s latitude, %s to %s longitude" %(
                   lat_str(geo_selection.s),lat_str(geo_selection.n),lon_str(geo_selection.w),lon_str(geo_selection.e)
                   )
        elif type_geo_select == "circle":
            print "Currently, geographic selection done by circular area,"
            print "Centre point: %s , %s , Radius %g km" %(
                   lat_str(geo_selection.lat),lon_str(geo_selection.lon),geo_selection.radius)
        elif type_geo_select == "country":
            print "Currenty, geographic selection done by country"            
            #if type(selection) == str: #if it's one country, it should be a 1-member list
            #    geo_selection = [geo_selection]
            print "Included countries are:"
            for a in geo_selection:
                print a
        elif type_geo_select == "state":
            print "Currently, geographic selection done by state"
            #if type(geo_selection) == str: #if it's one state, it should be a 1-member list
            #    geo_selection = [geo_selection]
            print "Included states are:"
            for a in geo_selection:
                print a
        print "OPTIONS"
        print "[1] Define new rectangular selection area"
        print "[2] Define new circular selection area"
        print "[3] Define new selection area by country/countries"
        print "[4] Define new selection area by state(s)"
        print "[Z] Return to previous menu"
        choice = raw_input("-->").upper()
        if choice == "1":
            valid = False
            while valid == False: #use the in-built validation method for this class
                type_geo_select = "rectangle"
                new_north = float(raw_input("Enter lat for N boundary (N +ve, S -ve) -->"))
                new_south = float(raw_input("Enter lat for S boundary (N +ve, S -ve) -->"))
                new_east  = float(raw_input("Enter lon for E boundary (E +ve, W -ve) -->"))
                new_west  = float(raw_input("Enter lon for W boundary (E +ve, W -ve) -->"))
                geo_selection = box(new_north,new_south,new_east,new_west)
                valid = geo_selection.valid()
                
        elif choice == "2":
            type_geo_select = "circle"
            new_cen_lat = float(raw_input("Enter lat for circle centre (N +ve, S -ve) -->"))
            new_cen_lon = float(raw_input("Enter lon for circle centre (E +ve, W -ve) -->"))
            new_radius  = float(raw_input("Enter radius of circle (km) -->"))
            geo_selection = circle(new_cen_lat,new_cen_lon,new_radius)
        elif choice == "3":
            type_geo_select = "country"
            geo_selection = []
            print "Enter the countries to be included. An empty entry will finalise the list"
            while True:
                add_country= raw_input("-->")
                if add_country == "":
                    break
                else:
                    geo_selection.append(add_country)
        elif choice == "4":
            type_geo_select = "state"
            geo_selection = []
            print "Enter the states to be included. An empty entry will finalise the list"
            while True:
                add_state= raw_input("-->")
                if add_state == "":
                    break
                else:
                    geo_selection.append(add_state)
        elif choice == "Z":
            print "No change made"
          
    return(type_geo_select,geo_selection)
                                    
                                 
def global_options(map_box,xdim,ydim,startdate,enddate,winter_flag):
    leave = False
    winter_flag = False
    while leave == False :
        os.system('clear')
        print "GLOBAL OPTIONS"
        print "Spatial viewing boundaries for maps:"
        print "NORTH = "+lat_str(map_box.n)
        print "SOUTH = "+lat_str(map_box.s)
        print "EAST  = "+lon_str(map_box.e)
        print "WEST  = "+lon_str(map_box.w)
        print "Longitudinal spacing of gridded data = "+str(xdim)
        print "Latitudinal  spacing of gridded data = "+str(ydim)
        print "Start time = " + str(startdate.year) +"-"+str(startdate.month) + "-"+str(startdate.day)
        print "End  time = " + str(enddate.year) +"-"+str(enddate.month) + "-"+str(enddate.day)
        print "[B] Change boundary options"
        print "[S] Change spacing options"  
        print "[T] Change time frame"
        print "[W] Toggle special winter months in this time-frame option (currently %r)"%winter_flag
        print "[Z] Return to top menu"
        choice = raw_input("-->").upper()
        if choice == "B" or choice == "b":
            new_north = float(raw_input("Enter new value for NORTH:"))
            new_south = float(raw_input("Enter new value for SOUTH:"))
            new_east  = float(raw_input("Enter new value for EAST :"))
            new_west  = float(raw_input("Enter new value for WEST :"))
            map_box = box(new_north,new_south,new_east,new_west)
        elif choice == "S" or choice == "s":
            xdim = change_var(xdim,"Longitudinal spacing of gridded data")
            ydim = change_var(ydim,"Latitudinal  spacing of gridded data")
        elif choice == "T" or choice == "t":
            new_s_year = change_var(startdate.year,"Start time year")
            new_s_month= change_var(startdate.month,"Start time month")            
            new_s_day  = change_var(startdate.day,"Start time day")
            startdate  = dt(new_s_year,new_s_month,new_s_day,0,0,0)
            new_e_year = change_var(  enddate.year,"End time year")
            new_e_month= change_var(  enddate.month,"End time month")
            new_e_day  = change_var(  enddate.day,"End time day")
            enddate    = dt(new_e_year,new_e_month,new_e_day,0,0,0)
        elif choice == "W":
            winter_flag = not winter_flag
        elif choice == "Z" or choice == "z":
            leave = True
        
    return(map_box,xdim,ydim,startdate,enddate,winter_flag)

def change_pickle(current_pickle):
    #os.system('clear')
    print "Currently, pickles named %s* are loaded." %current_pickle
    print "To change this, enter a new path and pickle name (excluding suffixes)."
    print "Otherwise, enter Z to keep current option and return totop menu"
    choice = raw_input("-->")
    if choice == "Z" or choice == "z":
        changed = False
    else:
        current_pickle = choice
        changed = True
    return(changed,current_pickle)
    
def change_emfiles(current_emfiles):
    os.system('clear')
    print "Currently, pickles from folder %s are loaded."%current_emfiles
    print "To change this, enter a new folder."
    print "Otherwise, enter Z to keep current option and return to top menu"
    choice = raw_input("-->")
    if choice == "Z" or choice == "z":
        changed = False
    else:
        current_emfiles = choice
        changed = True
    return(changed,current_emfiles)
    
def change_var(var,var_name="[unspecified]"):
    former_type = type(var)
    var = former_type(raw_input("Enter new value for "+var_name+"-->"))
    return(var)    

#def load_new_pickles_all(current_pickle,verbose=True,NDVI=True,indv=True,binn=True):
#    if NDVI:
#        NDVI_data = load_new_pickles_NDVI(current_pickle,verbose=verbose)
#    if indv:
#        indv_data = load_new_pickles_indv(current_pickle,verbose=verbose)
#    if binn:
#        binn_data = load_new_pickles_binn(current_pickle,verbose=verbose)
#    
#    if NDVI:
#        if indv:
#            if binn:
#                return(NDVI_data,indv_data,binn_data) #NIB
#            else:
#                return(NDVI_data,indv_data)           #NIx
#        else:
#            if binn:
#                return(NDVI_data,binn_data)           #NxB
#            else:
#                return(NDVI_data)                     #Nxx
#    else:
#        if indv:
#            if binn:
#                return(indv_data,binn_data)           #xIB
#            else:
#                return(indv_data)                     #xIx
#        else:
#            if binn:
#                return(binn_data)                     #xxB
#            else:
#                return()                              #xxx                    

    
def load_new_pickles_NDVI(current_pickle,verbose=True):
    #NDVIs
    if verbose:
        print "Loading pickled NDVI data from %s" %(current_pickle + "_NDVI.p")
    (NDVI_lat,NDVI_lon,NDVI,NDVI_year,NDVI_month) = pickle.load( open(current_pickle + "_NDVI.p","rb") )
    #(NDVI_lat,NDVI_lon,NDVI,NDVI_year,NDVI_month) = stripallnans(NDVI_lat,NDVI_lon,NDVI,NDVI_year,NDVI_month)
    if verbose:
        print "Pickled NDVI data now loaded as %i-long lists" %len(NDVI)
        print("NDVI latitudes        : NDVI_lat  ")
        print("NDVI longitudes       : NDVI_lon  ")
        print("NDVIs                 : NDVI      ") 
        print("NDVI years            : NDVI_year ")
        print("NDVI months           : NDVI_month")
        print("-----------------------------")
    
    NDVI_data = (NDVI_lat,NDVI_lon,NDVI,NDVI_year,NDVI_month)
    return(NDVI_data)
    
def load_new_pickles_da(current_pickle,suffix):    
    #Individual observations
    print("Loading pickled individual observations")
    da = cPickle.load( open(current_pickle + suffix,"rb") )
    print("Loaded the following datasets")
    for text in da.list_all_datasets():
        print text
    
    return(da)

def stripallnans(*a):
    """For multiple 1D lists of floats of the same size, indentifies locations of nan entries, and cuts these indexes from all lists"""   
    #a is a tuple of lists   
    n = 0
    #first, lets check all lists are the same length. Do this by cycling through and checking each is the same length as the previous one.
    for a_entry in a:
       if n == 0 :
          checklen = len(a_entry)
          n += 1
       else :
          new_checklen = len(a_entry)
          print new_checklen
          n += 1
          if checklen != new_checklen :
             raise ValueError("In stripall nans, len of entry #%i not equal to len of entry #%i" %(n,n-1))
          else :
             checklen = new_checklen   
    #Now that's done, loop through indexes
    truefalse = []
    
    for i in range(0,checklen) :
        flag = 0
        for a_entry in a:
            if type(a_entry[i]) == str:
                pass
            elif type(a_entry[i]) == float :               
                if a_entry[i] == float("nan") or np.isnan(a_entry[i]):
                    flag += 1
            elif type(a_entry[i]) == np.float64 :
                if np.isnan(a_entry[i]):
                    flag += 1
            elif a_entry[i] == np.nan:
                flag += 1                    
        if flag == 0 :
            truefalse.append(True)
        else:
            truefalse.append(False)            
    b = []
    bb = 0
    for a_entry in a:
        b.append([a_entry[i] for i in range(0,checklen) if truefalse[i] ])
        print len(b)
        print len(b[bb])
        bb+=1     
    return(tuple(b))



def draw_screen_poly( rec_lats, rec_lons, m, c, color_index,colmap ):
    #rainbow = plt.get_cmap('rainbow')
    rainbow = colmap.mpl_colormap
    x, y = m( rec_lons, rec_lats )
    xy = zip(x,y)
    this_col=color_index(c,clip=True)
    #poly = Polygon( xy, facecolor=[min(1.,0.+this_col*2),1.-2*abs(0.5-this_col),min(1.,2.-this_col*2)], edgecolor='none', alpha=1.0 )
    poly = Polygon( xy, facecolor=rainbow(this_col), edgecolor='none', alpha=1.0 )
    plt.gca().add_patch(poly)
    
def prepare_map(map_box,boundary=0.1,labels=False):
    """Creates a basemap object ready for plotting"""
    m = Basemap(projection='merc',
            llcrnrlat=map_box.s-boundary,urcrnrlat=map_box.n+boundary,
            llcrnrlon=map_box.w-boundary,urcrnrlon=map_box.e+boundary,
            lat_ts=(map_box.n+map_box.s)/2.,resolution='h')
    
    #get appropriate grid line spacing for this map
    min_dimension = min(map_box.n-map_box.s,map_box.e-map_box.w)
    if min_dimension < 1:
        gl_spacing = 0.1
    elif min_dimension < 5 :
        gl_spacing = 0.5
    elif min_dimension < 10:
        gl_spacing = 1.0
    elif min_dimension < 50:
        gl_spacing = 5.0
    else:
        gl_spacing = 10.
       
    #Draw lines
    m.drawcoastlines()
    m.drawcountries(linewidth=0.25)
    if labels:
        m.drawparallels(np.arange(-90.,90.,gl_spacing),labels=[True,False,False,False])
        m.drawmeridians(np.arange(-180.,181.,gl_spacing),labels=[False,False,False,True]) 
    else:
        m.drawparallels(np.arange(-90.,90.,gl_spacing),labels=[False,False,False,False])
        m.drawmeridians(np.arange(-180.,181.,gl_spacing),labels=[False,False,False,False]) 
    return m
    
def free_colorbar(vmin,vmax,label="no label selected",colmap=None,logplot=False):
    #if coltype == "bwr":
    #    cols_for_bar = [[0.,0.,1.],[1.,1.,1.],[1.,0.,0.]]
    #elif coltype == "lg-dg":
    #    cols_for_bar = [[0.,0.2,0.],[0.8,0.2,0.8]]
    #else:
    #    rainbow = plt.get_cmap(coltype)
    #    cols_for_bar = []
    #    for i in range(0,256):
    #        cols_for_bar.append(list(rainbow(i)[0:3]))
            
    cols_for_bar = colmap.mpl_colors        
    #can add in more cols_for_bar options here.
    cm1 = colors.LinearSegmentedColormap.from_list("MyCmapName",cols_for_bar)
    if logplot:
        cnorm = colors.LogNorm(vmin=vmin,vmax=vmax)
    else:
        cnorm = colors.Normalize(vmin=vmin,vmax=vmax)
    cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
    cpick.set_array([])
    plt.colorbar(cpick,label=label)

def plot_grid_from_list(lat,lon,var,xdim,ydim,map_box,title="Unnamed plot",vmin=0.0e16,vmax=3.0e16,lab="HCHO column molec/cm2",
                        save=False,save_filename="nofilenamespecified",
                        colmap=None,do_plot=True,labels=False,logplot=False):
    #Define basemap
    m = prepare_map(map_box,labels=labels)
    if logplot:
        color_index = colors.LogNorm(vmin=vmin,vmax=vmax)
    else:
        color_index = colors.Normalize(vmin,vmax)
    #(lat,lon,var) = stripallnans(lat,lon,var) #cut out nan values (stops crashes)    
    for i in range(0,len(lat)):
        this_w = lon[i]-xdim/2
        this_e = lon[i]+xdim/2
        this_s = lat[i]-ydim/2
        this_n = lat[i]+ydim/2
        this_c = var[i]
        if math.isnan(this_c) or math.isnan(this_w) or math.isnan(this_e) or math.isnan(this_s) or math.isnan(this_n):
            continue #don't plot if there's no data
        draw_screen_poly( [this_s,this_n,this_n,this_s],[this_w,this_w,this_e,this_e], m, this_c, color_index, colmap )        
    plt.title(title)
    free_colorbar(vmin,vmax,label=lab,colmap=colmap,logplot=logplot)
    fig = plt.gcf()
    if save:    
        fig.savefig(save_filename, dpi=600)
    if do_plot:
        plt.show()
    else:
        plt.close(fig)

def plot_dots_on_map(lat,lon,var,map_box,vmin=0.,vmax=3.0E16,title="Unnamed plot",lab="HCHO column molec/cm2",
                     save=False,save_filename="nofilenamespecified",
                     colmap=None):
    """Draws dots on a map at lat/lon positions, color-coded based on values"""
    m = prepare_map(map_box)
    x, y = m(lon,lat)
    m.scatter(x,y,18,marker='.',edgecolors='none',c=var, vmin=vmin, vmax=vmax, cmap=colmap.mpl_colormap)
    plt.title(title)
    free_colorbar(vmin,vmax,label=lab,colmap=colmap)
    fig = plt.gcf()
    if save:    
        fig.savefig(save_filename, dpi=600)
    plt.show()

    
def geo_select_regional(region_list,text_to_match,*datasets):
    """Reduces datasets down to only where the associated region matches a given text string"""
    output = []
    for dataset in datasets :  
        selected = [dataset[i] for i in range(0,len(dataset)) if region_list[i] in text_to_match]
        output.append(selected)        
    return(tuple(output))
    
def geo_select_regional_i(region_list,text_to_match,da):
    """Reduces datasets down to only where the associated region matches a given text string"""
    selector = []
    for region in region_list :  
        selector.append(region in text_to_match)
    
    da.filter_all(selector)      
    return(da)    

def da_select_by_match(field_to_match,val_to_match,da):
    """Reduces data_all object down to where a dataset matches a value"""
    #useful for filtering by region
    
    if type(val_to_match) != list: #if it's a single value, make 1-member list
        list_to_match = [val_to_match]
    else:
        list_to_match = val_to_match
        
    len_da = len(da.lat)
    
    selector = [ da.data[field_to_match].val[i] in list_to_match
                 for i in range(0,len_da) ]
    da.filter_all(selector)
    return(da)

#def in_box(lat,lon,north,south,east,west):
#    if lat < south :
#        return False
#    if lat > north :
#        return False
#    if lon < west :
#        return False
#    if lon > east :
#        return False
#    return True
    
def geo_select_rectangle(lat,lon,geo_selection,*datasets):
    """Reduces datasets down to only those points within designated bounds"""
    #[northbound,southbound,eastbound,westbound] = geo_selection
    output = []
    for dataset in datasets :
        selected = [dataset [i] for i in range(0,len(dataset)) if geo_selection.within(lat[i],lon[i]) ]
        output.append(selected)
    return(tuple(output))
    
def geo_select_rectangle_i(geo_selection,ida):
    """Reduces datasets down to only those points within designated bounds"""

    selector = geo_selection.within(ida.lat,ida.lon)
    ida.filter_all(selector)
    return(ida)    

def da_select_by_shape(select_shape,da):
    """Reduces data_all objects down to only those points within designated bounds"""
    selector = select_shape.within(da.lat,da.lon)
    da.filter_all(selector)
    #return(da)
    
def geo_select_circle(lat,lon,geo_selection,*datasets):
    """Reduces datasets down to only those points within a specified radius of a specified point"""
    [centre_lat,centre_lon,radius] = geo_selection
    output = []
    distances=[]
    for i in range(0,len(lat)):
       distances.append(great_circle((lat[i],lon[i]),(centre_lat,centre_lon)).km)
    
    for dataset in datasets :
        selected = [dataset[i] for i in range(0,len(dataset)) if distances[i] <= radius ]
        output.append(selected)
    return(tuple(output))
      
    
def geo_select_circle_i(geo_selection,ida):
    """Reduces datasets down to only those points within a specified radius of a specified point"""
    [centre_lat,centre_lon,radius] = geo_selection
    output = []
    selector = []
    for i in range(0,len(ida.lat)):
       distance = great_circle( (ida.lat[i],ida.lon[i]),
                                (centre_lat,centre_lon)
                              ).km
       selector.append(distance <= radius)
    
    ida.filter_all(selector)    
    return(ida)
    
def time_select(startdate,enddate,ida,winter_flag=False):
    """Reduces ida down to only points within two dates"""
    selector = []
    if winter_flag:
        for t in ida.time:
            if t.month not in [1,2,3,12]:
                selector.append(False) #automatic filter out if not JFM or D
            else:
                selector.append(startdate <= t <= enddate)
    else:
        for t in ida.time:
            selector.append(startdate <= t <= enddate)
    ida.filter_all(selector)
    
    return(ida)
    
def time_select_list(startdate,enddate,time,data):
    """Reduces data down to only points where associated time within two dates"""
    data_filtered = []
    for i in range(0,len(data)):
        if startdate <= time[i] <= enddate:
            data_filtered.append(data[i])
    
    return(data_filtered)
    
def time_cycle(ida,data_key,stat_choice,
               step,month_flag=False,
               plot=False):
    
    time = ida.time
    data = ida.data[data_key].val           
    alpha = 1.0
       
    #sort the data into time bins
    earliest_dt = min(time)
    last_dt     = max(time)
    
    if month_flag: #if we're doing calendar months        
        clocklow = dt(earliest_dt.year,earliest_dt.month,1,0,0,0)
        
        step = relativedelta(months=step)
    else: #normal steps
        clocklow = dt(earliest_dt.year,earliest_dt.month,earliest_dt.day,0,0,0)
        step = relativedelta(days=step)
        
    stat_collection = []
    time_collection = []
    time_bounds_collection = []    

    print "Computing..."
    
    while clocklow <= last_dt:
        clockhigh = clocklow + step
        #print clocklow
        #print clockhigh
        #print step
        #print "-----"
        this_dataset = time_select_list(clocklow,clockhigh,time,data)
        this_stat = calc_statistic(this_dataset,stat_choice)
        central_time = clocklow + 0.5*step
        stat_collection.append(this_stat)
        time_collection.append(central_time)
        time_bounds_collection.append([clocklow,clockhigh])
        #now continue the cycle
        clocklow = clocklow + step
        
    if plot:
        plot_title = ida.data[data_key].description
        x_label = "date"
        if stat_choice == "1":
            y_label = "mean value / " + ida.data[data_key].unit
        elif stat_choice ==  "2":
            y_label = "standard deviation of values / " + ida.data[data_key].unit
        elif stat_choice == "3":
            y_label = "number of values"

        
        #plot these:
        time_scatter(time_collection,stat_collection,
                     title=plot_title,
                     x_label=x_label,y_label=y_label,alpha=1.0)
                     
        print time_collection
        print stat_collection
        _ = raw_input("Press enter to continue")
    
    return(stat_collection,time_collection,time_bounds_collection)

def time_scatter(time,y,yerr=[],title="UNNAMED PLOT",y_label="",x_label="",alpha=1.0):
    
    if yerr == []:
        plt.plot(time, y, alpha=alpha, color='g',marker="o")
    else:
        plt.errorbar(time, y, yerr=yerr, alpha=alpha, color='g',marker="o")

    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(y_label, fontsize=16)
    plt.grid(b=True, which='major', color='0.65')

    plt.title(title)
    plt.show()

def get_mode(inlist):
    return(max(set(list(inlist)), key=list(inlist).count))

def get_meanerr(inlist):
    return(np.sqrt(np.divide(1.,np.sum(np.divide(1.,inlist)))))
       
def binner(lat,lon,vals,map_box,stat="mean",xdim=0.3125,ydim=0.25,do_extras=False):
    """Divides the region into a grid, and computes a statistic for the values falling within it"""

    num_xbins = int((map_box.e-map_box.w)/xdim) + 1
    num_ybins = int((map_box.n-map_box.s)/ydim) + 1
    
    #use scistats.binned_statistic_2d with the np function chosen.
    if stat == "mean" : #mean of values in each box
        (binned_stat,xedges,yedges,binnumber) = \
            scistats.binned_statistic_2d(lon,lat,vals,np.nanmean,bins=[num_xbins,num_ybins],
            range=[[map_box.w-xdim/2,map_box.w+(num_xbins-1)*xdim+xdim/2],[map_box.s-ydim/2,map_box.s+(num_ybins-1)*ydim+ydim/2]],
            expand_binnumbers=True)
    elif stat == "stdev" : #standard deviation of values in each box
        (binned_stat,xedges,yedges,binnumber) = \
            scistats.binned_statistic_2d(lon,lat,vals,np.nanstd,bins=[num_xbins,num_ybins],
            range=[[map_box.w-xdim/2,map_box.w+(num_xbins-1)*xdim+xdim/2],[map_box.s-ydim/2,map_box.s+(num_ybins-1)*ydim+ydim/2]],
            expand_binnumbers=True)
    elif stat == "sum" : #sum of all values in each box
        (binned_stat,xedges,yedges,binnumber) = \
            scistats.binned_statistic_2d(lon,lat,vals,np.nansum,bins=[num_xbins,num_ybins],
            range=[[map_box.w-xdim/2,map_box.w+(num_xbins-1)*xdim+xdim/2],[map_box.s-ydim/2,map_box.s+(num_ybins-1)*ydim+ydim/2]],
            expand_binnumbers=True)
    elif stat == "count" : #count of number of values in each box
        (binned_stat,xedges,yedges,binnumber) = \
            scistats.binned_statistic_2d(lon,lat,vals,notnancount,bins=[num_xbins,num_ybins],
            range=[[map_box.w-xdim/2,map_box.w+(num_xbins-1)*xdim+xdim/2],[map_box.s-ydim/2,map_box.s+(num_ybins-1)*ydim+ydim/2]],
            expand_binnumbers=True)
    elif stat == "quadsum" : #sum of squares of all values
        qvals = list(np.multiply(vals,vals))
        (binned_stat,xedges,yedges,binnumber) = \
            scistats.binned_statistic_2d(lon,lat,qvals,np.nansum,bins=[num_xbins,num_ybins],
            range=[[map_box.w-xdim/2,map_box.w+(num_xbins-1)*xdim+xdim/2],[map_box.s-ydim/2,map_box.s+(num_ybins-1)*ydim+ydim/2]],
            expand_binnumbers=True)
    elif stat == "mode": #most common value in each cell
        (binned_stat,xedges,yedges,binnumber) = \
            scistats.binned_statistic_2d(lon,lat,vals,get_mode,bins=[num_xbins,num_ybins],
            range=[[map_box.w-xdim/2,map_box.w+(num_xbins-1)*xdim+xdim/2],[map_box.s-ydim/2,map_box.s+(num_ybins-1)*ydim+ydim/2]],
            expand_binnumbers=True)
    elif stat == "combined measurement error": #error in mean
        #(binned_stat,xedges,yedges,binnumber) = \
        #    scistats.binned_statistic_2d(lon,lat,vals,get_meanerr,bins=[num_xbins,num_ybins],
        #    range=[[map_box.w-xdim/2,map_box.w+(num_xbins-1)*xdim+xdim/2],[map_box.s-ydim/2,map_box.s+(num_ybins-1)*ydim+ydim/2]],
        #    expand_binnumbers=True)            
        qvals = list(np.multiply(vals,vals))
        (binned_stat_q,xedges,yedges,binnumber) = \
            scistats.binned_statistic_2d(lon,lat,qvals,np.nansum,bins=[num_xbins,num_ybins],
            range=[[map_box.w-xdim/2,map_box.w+(num_xbins-1)*xdim+xdim/2],[map_box.s-ydim/2,map_box.s+(num_ybins-1)*ydim+ydim/2]],
            expand_binnumbers=True)
        (binned_stat_c,xedges,yedges,binnumber) = \
            scistats.binned_statistic_2d(lon,lat,vals,notnancount,bins=[num_xbins,num_ybins],
            range=[[map_box.w-xdim/2,map_box.w+(num_xbins-1)*xdim+xdim/2],[map_box.s-ydim/2,map_box.s+(num_ybins-1)*ydim+ydim/2]],
            expand_binnumbers=True)
        binned_stat = np.divide(np.sqrt(binned_stat_q),binned_stat_c)
    
    #can also work out lats, lons and index map of bins.
    #As this function is run several times, and these data would be the same
    #on each iteration, this is made optional.
    #Return of xedges,yedges also dependant on this option.
    if do_extras:
            
        numberofdata = len(binnumber[0])
        #index map is a list the same length as the individual data input
        #it refers to which bin that data point has been assigned.
        index_map = []
        for i in range(0,numberofdata):
            index_map.append(binnumber[0][i]*num_ybins + binnumber[1][i])
        
        #binned_stat is a 2D data array. The below makes this a 1D list
        #Also creates 1D list of bins' lats and lons
        lat_list = []
        lon_list = []
        out_list = []
        for i in range(0,len(xedges)-1):
           this_w = xedges[i]
           this_e = xedges[i+1]
           for j in range(0,len(yedges)-1):
               this_s = yedges[j]
               this_n = yedges[j+1]
               lat_list.append((this_s+this_n)/2.)
               lon_list.append((this_e+this_w)/2.)
               out_list.append(binned_stat[i][j])
        return(lat_list,lon_list,out_list,xedges,yedges,index_map)
    else:
        out_list = []
        for i in range(0,len(xedges)-1):
            for j in range(0,len(yedges)-1):
                out_list.append(binned_stat[i][j])
        return(out_list)
    
def create_binned_set(ida,map_box,xdim,ydim):
    """Given the individual data, returns a full set of binned data"""
    
    #pull variables out of ida
    lat=ida.lat
    lon=ida.lon
    
    #get lat_ & lon_ binned
    (lat_binned,lon_binned,null_1,xedges,yedges,index_map) = \
        binner(lat,lon,lat,map_box,stat="mean",xdim=xdim,ydim=ydim,
               do_extras=True)   
    
    bda = d_all(lat_binned,lon_binned)            
    bda.meta["Data type"] = "Binned data"
    bda.meta["Area"] = map_box
    bda.meta["Binning dimensions"] = [xdim,ydim]
         
    try:
        sat_VC=ida.data['sat_VC'].val
        sat_DVC=ida.data['sat_DVC'].val
        geos_VC=ida.data['geos_VC'].val

        #Mean satellite observation 
        sat_VC_mean_binned = \
            binner(lat,lon,sat_VC,map_box,stat="mean",xdim=xdim,ydim=ydim)
        #Standard deviation in satellite observation
        sat_VC_stdev_binned = \
            binner(lat,lon,sat_VC,map_box,stat="stdev",xdim=xdim,ydim=ydim)
        #Count of number of observations
        sat_VC_count_binned = \
            binner(lat,lon,sat_VC,map_box,stat="count",xdim=xdim,ydim=ydim)
        #Mean error in satellite observation
        sat_DVC_mean_binned = \
            binner(lat,lon,sat_DVC,map_box,stat="mean",xdim=xdim,ydim=ydim)    
        #Mean model observation
        geos_VC_mean_binned = \
            binner(lat,lon,geos_VC,map_box,stat="mean",xdim=xdim,ydim=ydim)
        #Standard deviation in model values
        geos_VC_stdev_binned = \
            binner(lat,lon,geos_VC,map_box,stat="stdev",xdim=xdim,ydim=ydim)

        #push these lists into bda
        vals         = [sat_VC_mean_binned,
                        sat_VC_stdev_binned,
                        sat_VC_count_binned,
                        sat_DVC_mean_binned,
                        geos_VC_mean_binned,
                        geos_VC_stdev_binned]
        names        = ["sat_VC_mean",
                        "sat_VC_stdev",
                        "sat_VC_count",
                        "sat_DVC_mean",
                        "geos_VC_mean",
                        "geos_VC_stdev"]
        descriptions = ["Mean satellite vertical column",
                        "Standard deviation of satellite vertical column",
                        "Number of observations in bin",
                        "Mean error in satellite vertical column",
                        "Mean modelled vertical column",
                        "Standard deviation of modelled vertical column"]
        units        = ["molec/cm2",
                        "molec/cm2",
                        ""          ,
                        "molec/cm2",
                        "molec/cm2",
                        "molec/cm2"]
        
        for i in range(0,len(names)):
            bda.data[names[i]] = d(vals[i],names[i],description=descriptions[i])
            bda.data[names[i]].unit = units[i]
    
    except: #error means that sat_VC doesn't exist, will happen for BRUG datasets.
        print "Not standard form data, use 'bin additional datasets' option"
        _ = raw_input("Press enter to continue-->")
    
    return(bda)

def load_regiondata(filename,states=True):
    """Reads a csv file which reports the country (and, optionally, the state) on a grid"""
    cs = open(filename,"r")
    #Empty lists to populate
    cs_lat = []
    cs_lon = []
    cs_country = []
    if states:
        cs_state = []   
    for line in cs:
        line = line.strip()
        columns = line.split(",")
        cs_lat.append(float(columns[0]))
        cs_lon.append(float(columns[1]))
        cs_country.append( columns[2] )
        if states:
            cs_state.append( columns[3]    )    
    cs.close() #we can close the file here
    if states:
        return(cs_lat,cs_lon,cs_country,cs_state)
    else:
        return(cs_lat,cs_lon,cs_country)

def find_nearest(array,value):
    """A function that returns the nearest value in a sorted array""" 
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

def notnancount(data):
    """Returns count of values of data that are not NaN""" 
    return(np.count_nonzero(~np.isnan(data)))


def region_matcher(cs_lat,cs_lon,cs_country,cs_state,lat,lon,states=True):
    """For each lat,lon pair, assigns a country/state based on the output of load_regiondata"""
    
    num_cs = len(cs_lat)
    
    grid_lat = sorted(list(set(cs_lat)))
    grid_lon = sorted(list(set(cs_lon)))
    
    #Now, for each entry in lat_list etc., get the corresponding index in cs
    #and populate lists country_out_list and state_out   
    num_grid = len(lat)       
    country_out = []
    if states:
        state_out = []
    for line in range(0,num_grid):
        nearest_lat = find_nearest(grid_lat,lat[line])
        nearest_lon = find_nearest(grid_lon,lon[line])
        index = [i for i in range(0,num_cs) if ((nearest_lat == cs_lat[i]) and (nearest_lon == cs_lon[i]))]
        #print index
        country_out.append(cs_country[index[0]])
        if states:
            state_out.append(cs_state[index[0]])
    if states:
        return(country_out,state_out)
    else:
        return(country_out)
        
def region_matcher_fast(cs_lat,cs_lon,cs_country,cs_state,lat,lon,states=True):

    """A fast version of region_matcher that only works for country_state.csv"""
    
    num_grid = len(lat)
    
    #np arrays
    (cs_lat,cs_lon,lat,lon) = (np.array(cs_lat),np.array(cs_lon),np.array(lat),np.array(lon))
    
    #get i and j indexes
    
    j = np.around(np.divide((lat-2.),0.25))
    j = np.maximum(j,np.zeros_like(j))
    j = np.minimum(j,np.zeros_like(j)+144)
    
    i = np.around(np.divide((lon-65.),0.3125))
    i = np.maximum(i,np.zeros_like(i))
    i = np.minimum(i,np.zeros_like(i)+112)
    ind = np.add(113*j,i)
    
    print len(j)
    print len(ind)
    print np.min(ind)
    print np.max(ind)
    print len(cs_country)
    
    print "got indexes"
    country_out = []
    if states:
        state_out = []
    
    for line in range(0,num_grid):
        if line % 1000 == 0:
            print "Assigning %i of %i" %(line,num_grid)
        try:    
            country_out.append(cs_country[int(ind[line])])
            if states:
                state_out.append(cs_state[int(ind[line])])
        except IndexError:
            print "line = %i" %line
            print "ind[line] = %g" %ind[line]
            raise ValueError
    if states:
        return(country_out,state_out)
    else:
        return(country_out)        

def add_slash(string):
    """Adds a trailing / character if string does not currently end with one"""
    if string.endswith("/"):
        return(string)
    else:
        return(string+"/")

def save_pickle(da,save_path=None,suffix="_x.p"):
    """A routine to save a data_all object"""
    
    pickle_name = save_path+suffix  
    save_location = pickle_name
    
    cPickle.dump(da,open(save_location,"wb"))
    print "Pickle saved to %s" %save_location                 
    
def select_a_var(var_names,vartuple,purpose="",numbers_only=True):
    """menu to select one of many variables supplied in a tuple"""
    menu_options = []
    menu_counter = 0
    for var in vartuple:
        #don't include non-numerical datasets if numbers_only is True
        if isinstance(var[0],Number) or not numbers_only: 
            menu_options.append([str(menu_counter),var_names[menu_counter]])
        menu_counter += 1

    if purpose == "":
        #if purpose isn't chosen, have a generic title for this menu
        menu_title = "Choose data:"
    else:
        menu_title = "Choose data for %s:" %purpose
    var_choice_no = int(basic_menu(menu_title,menu_options,quit_option=False))
    return(vartuple[var_choice_no],var_names[var_choice_no])

def two_var_comparison(ida):
    """Compare two variables and plot a scatter"""
    while True: #
        tvc1_menu_title = "DATASET COMPARISON:\n"\
                          "Select first (x-axis) dataset:"
        [tvc1_menu_text,tcv1_answers_dict] = ida.make_menu()
        tvc1_menu_choice = basic_menu(tvc1_menu_title,
                                      tvc1_menu_text,
                                      quit_option=True)
        if tvc1_menu_choice == "Z":
            break
        x_key =      tcv1_answers_dict[tvc1_menu_choice]    
        x_var =      ida.data[x_key].val
        x_var_name = ida.data[x_key].description+" ("+ida.data[x_key].unit+")"
            
            
        while True:
            tvc2_menu_title = "DATASET COMPARISON:\n"\
                              "Select second (y-axis) dataset:"
            [tvc2_menu_text,tcv2_answers_dict] = ida.make_menu()
            tvc2_menu_choice = basic_menu(tvc2_menu_title,
                                          tvc2_menu_text,
                                          quit_option=True)
            if tvc2_menu_choice == "Z":
                break
            
            y_key =      tcv2_answers_dict[tvc2_menu_choice]   
            y_var =      ida.data[y_key].val    
            y_var_name = ida.data[y_key].description+" ("+ida.data[y_key].unit+")"
            
            while True:
                tvc3_menu_title = "Choose type of comparison:"
                tvc3_menu_text = [["1","Scatter plot"],
                                  ["2","Simple correlation statistics"],
                                  ["3","Error bar plot"],
                                  ["4","Advanced correlation statistics"],
                                  ["5","write out datasets as CSV"]]
                type_of_comparison = basic_menu(tvc3_menu_title,
                                                tvc3_menu_text,
                                                quit_option=True)
                if type_of_comparison == "Z":
                    break
                                                    
                #if 3 or 4, will need to define errors on one or both datasets.
                if type_of_comparison in ["3","4"]:
                    x_var_errs = get_errors(ida,x_key)
                    if x_var_errs == "Z":
                        continue
                    
                    y_var_errs = get_errors(ida,y_key)
                    if y_var_errs == "Z":
                        continue
                        
                #Need to make sure there are no nans in the data going in.
                if type_of_comparison in ["1","2","5"]:
                    (nonan_x_var,nonan_y_var) = \
                            stripallnans(x_var,y_var)
                else:                        
                    if x_var_errs != None:
                        if y_var_errs != None:
                            (nonan_x_var,nonan_y_var,nonan_x_var_errs,nonan_y_var_errs) = \
                            stripallnans(x_var,y_var,x_var_errs,y_var_errs)
                        else:
                            (nonan_x_var,nonan_y_var,nonan_x_var_errs) = \
                            stripallnans(x_var,y_var,x_var_errs)
                    else:
                        if y_var_errs != None:
                            (nonan_x_var,nonan_y_var,nonan_y_var_errs) = \
                            stripallnans(x_var,y_var,y_var_errs)
                        else:
                            (nonan_x_var,nonan_y_var) = \
                            stripallnans(x_var,y_var)                                                    
                                                                
                if type_of_comparison == "1":
                    goplot = False
                    axes_lims = [np.nanmin(nonan_x_var),np.nanmax(nonan_x_var),
                                 np.nanmin(nonan_y_var),np.nanmax(nonan_y_var)]
                    alpha = "auto" 
                    while goplot == False:
                        print "x-axis min: %g \nx-axis max: %g\ny-axis min: %g \ny-axis max: %g\nalpha: %s"\
                              %(axes_lims[0],axes_lims[1],axes_lims[2],axes_lims[3],str(alpha))
                        print "Type C to change, or press ENTER to plot"
                        opt = raw_input("-->").upper()
                        if opt == "C":
                            axes_lims[0] = input("New value for x axis minimum-->")
                            axes_lims[1] = input("New value for x axis maximum-->")
                            axes_lims[2] = input("New value for y axis minimum-->")
                            axes_lims[3] = input("New value for y axis maximum-->")
                            while True:
                                try:
                                    alpha = input("Alpha value ('auto' for automatic) -->")
                                    break
                                except NameError:
                                    print "Not valid input. Put 'single quotes' around text input"
                                except SyntaxError:
                                    print "Not valid input. Put 'single quotes' around text input"
                        elif opt == "":
                            goplot = True
                        
                    if type(alpha) == str:
                        alpha = None    
                    error_bar_scatter(nonan_x_var,nonan_y_var,
                                      x_error=None,y_error=None,
                                      title="Scatter plot",
                                      x_label=x_var_name,y_label=y_var_name,
                                      alpha=alpha,
                                      x_min=axes_lims[0],x_max=axes_lims[1],
                                      y_min=axes_lims[2],y_max=axes_lims[3],do_best_fit=True)
                elif type_of_comparison == "2":
                    print "Linear regression:"
                    slope, intercept, r_value, p_value, std_err = scistats.linregress(nonan_x_var,nonan_y_var)
                    print "Best-fit line: y = %+.2gx%+.2g " %(slope,intercept)
                    r_squared = r_value * r_value
                    print "R-squared    : %.4g " %r_squared
                    print "P-value      : %.4g " %p_value
                    print "Standard err : %g "   %std_err        
                    pearsonr = scistats.pearsonr(nonan_x_var,nonan_y_var)
                    print "Pearson's correlation coefficient:\n %f, 2-tailed p-value: %f" %(pearsonr[0],pearsonr[1])
                    polyfit = np.polyfit(nonan_x_var, nonan_y_var, 1)
                    inv_opt = raw_input("Enter to continue or I to use these as inputs to an inversion -->").upper()
                    if inv_opt == "I":
                        ida = second_stage(slope,std_err,intercept,ida)
                        _ = raw_input("Inversion done. Press enter to continue -->")
    
                elif type_of_comparison == "3":
                    goplot = False
                    axes_lims = [np.nanmin(nonan_x_var),np.nanmax(nonan_x_var),
                                 np.nanmin(nonan_y_var),np.nanmax(nonan_y_var)]
                    while goplot == False:
                        print "x-axis min: %g \nx-axis min: %g\ny-axis min: %g \ny-axis min: %g"\
                              %(axes_lims[0],axes_lims[1],axes_lims[2],axes_lims[3])
                        print "Type C to change, or press ENTER to plot"
                        opt = raw_input("-->").upper()
                        if opt == "C":
                           axes_lims[0] = input("New value for x axis minimum-->")
                           axes_lims[1] = input("New value for x axis maximum-->")
                           axes_lims[2] = input("New value for y axis minimum-->")
                           axes_lims[3] = input("New value for y axis maximum-->")
                        elif opt == "":
                            goplot = True
                    error_bar_scatter(nonan_x_var,nonan_y_var,
                          x_error=nonan_x_var_errs,y_error=nonan_y_var_errs,
                          title="Error bar plot",
                          x_label=x_var_name,y_label=y_var_name,
                          x_min=axes_lims[0],x_max=axes_lims[1],
                          y_min=axes_lims[2],y_max=axes_lims[3])
                elif type_of_comparison == "4":
                    #WYIBF analysis
                    (m,b,merr,berr,gof) = WYIBF(
                        nonan_x_var     , nonan_y_var,
                        nonan_x_var_errs, nonan_y_var_errs,
                        err_type="stdev",iterate=True)
                    print "Best-fit line: y = (%+.2g +/- %+.2g)x + (%+.2g +/- %+.2g)" %(m,merr,b,berr)
                    print "Goodness of fit: %g" %gof 
                    _ = raw_input("Press enter to continue-->")
                elif type_of_comparison == "5":
                    #write as csv
                    write_filename = raw_input("Enter path to save to, including extension -->")
                    f = open(write_filename,"wb")
                    for i in range(len(nonan_x_var)):
                        f.write(str(nonan_x_var[i]))
                        f.write(",")
                        f.write(str(nonan_y_var[i]))
                        f.write(",")
                        f.write("\n")
                    f.close()

def get_errors(ida,var_key):
    """Selects errors based on user options"""
    
    var_name = ida.data[var_key].description
    var = ida.data[var_key].val
       
    e_choice = basic_menu("How is error in %s determined?" %var_name,[
                        ["1","No error"],
                        ["2","Fixed value for all data"],
                        ["3","Fixed fraction for all data"],
                        ["4","Use existing dataset"]],
                        quit_option=True)
    if e_choice == "Z":
        var_errs = "Z"
    elif e_choice == "1":
        var_errs = list(np.zeros_like(var))
    elif e_choice == "2":                                  
        var_errs = list(np.zeros_like(var))
        fixed_error_value = float(raw_input("Enter fixed error value\n"
                                        "-->"))
        for i in range(0,len(var_errs)):
            var_errs[i] = fixed_error_value

    elif e_choice == "3":
        var_errs = list(np.zeros_like(var))
        fraction_error_value =float(raw_input("Enter fractional error value\n"
                                              "-->"))
        for i in range(0,len(var_errs)):
            var_errs[i] = var[i]*fraction_error_value
        
        
    elif e_choice == "4":
        e_d_menu_title = "DATASET COMPARISON:\n"\
                      "Select second (y-axis) dataset:"
        [e_d_menu_text,e_d_answers_dict] = ida.make_menu()
        e_d_menu_choice = basic_menu(e_d_menu_title,
                                  e_d_menu_text,
                                  quit_option=False)
        
        var_errs = ida.data[e_d_answers_dict[e_d_menu_choice]].val
     
    return(var_errs)    

        
def error_bar_scatter(x_var,y_var,
                      x_error=None,y_error=None,
                      x_min=None,x_max=None,y_min=None,y_max=None,
                      x_label="",y_label="",
                      alpha=None,
                      title="",
                      do_best_fit=False,
                      do_best_fit_equation=False,
                      show=True):
    """Draws a scatter plot of two data sets. Best fit line and error bars optional"""
    
    #if alpha (transparency) undefined, estimate a good alpha based on the dataset size
    if alpha == None:
        alpha = float(len(x_var))**-0.6
        
    #if plotting mins and maxes are not defined, use mins and maxes of data,
    if x_min == None:
        x_min = np.nanmin(x_var)
    if x_max == None:
        x_max = np.nanmax(x_var)
    if y_min == None:
        y_min = np.nanmin(y_var)
    if y_max == None:
        y_max = np.nanmin(y_var)
        
    #consider all 0 error to be no errors defined    
    if x_error == None and y_error == None: #no error bars
        plt.plot(x_var, y_var,'o',
                    alpha=alpha)
    elif x_error != None and y_error == None: #bars for x, none for y
        plt.errorbar(x_var, y_var, xerr=x_error,fmt='o',
                     alpha=alpha)
    elif x_error == None and y_error != None: #none for x, bars for y
        plt.errorbar(x_var, y_var, yerr=y_error,fmt='o',
                     alpha=alpha)
    elif x_error != None and y_error != None: #bars for x, bars for y
        plt.errorbar(x_var, y_var, xerr=x_error, yerr=y_error,fmt='o',
                     alpha=alpha)   
    
    if do_best_fit or do_best_fit_equation:
    
        slope, intercept, r_value, p_value, std_err = scistats.linregress(x_var, y_var)
        
        line_xs = [x_min, x_max]
        line_ys = [x_min*slope + intercept,x_max*slope + intercept]
        if do_best_fit:
            plt.plot(line_xs,line_ys,'-r')
        if do_best_fit_equation:
            plt.text((x_max-x_min)*0.1 + x_min,(y_max-y_min)*0.9 + y_min,"fit: "+str("%.2f" %slope)+"x "+str("%+.2e" %intercept))
    
    #axis labels       
    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(y_label, fontsize=16)
    #axis limits
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    #grid
    plt.grid(b=True, which='major', color='0.65')
    #title
    plt.title(title)
    
    if show:
        plt.show()



def load_geosfile(emfiles_folder,file_tag="trac_avg",columns=False):
    """Loads up emissions from one or several files"""
    
    export_dict = {} #eventually we'll populate this dictionary and return it
      
    #make sure provided folder string ends with a slash
    if not emfiles_folder.endswith("/"):
        emfiles_folder = emfiles_folder + "/"
        
    print "Will search for files in %s"%emfiles_folder    
        
    files_in_folder = [f for f in listdir(emfiles_folder) if isfile(join(emfiles_folder, f))]
    #print files_in_folder
    #filter this list down to just trac_avg files.
    files_in_folder = [f for f in files_in_folder if f.startswith(file_tag)]
    
    #add path to filenames
    for i in range(0,len(files_in_folder)):
        files_in_folder[i] = emfiles_folder + files_in_folder[i]
    
    option = ""
    while option != "G":
        #clear_screen()
        print "Indentifed the following geos_chem output files:"
        number_options = []
        for i in range(0,len(files_in_folder)):
            print "[%i] %s" %(i,files_in_folder[i])
            number_options.append(i)
        print "To remove a file from this list, type X then its number (i.e. X0)"
        print "To process only one file from this list, type Y then its number (i.e. Y0)"
        print "If you are happy with this list, type G"
        option = raw_input("-->").upper()
        if option == "G":
            pass
        elif option.startswith("X"):
            #excising a file from the list
            option_int = int(option[1:])
            del files_in_folder[option_int]
        elif option.startswith("Y"):
            option_int = int(option[1:])
            files_in_folder = [files_in_folder[option_int]]
            option = "G" #function will proceed
    
    #Quite often these bpch files cannot be opened. What we'll do is try to load each
    #and see what works
    good_files = []
    for bpch_filename in files_in_folder:
        print "Attempting to open %s" %bpch_filename
        
        try:
            this_bpch = bpch(bpch_filename)
            good_files.append(bpch_filename)
        except:
            print("Unexpected error:", sys.exc_info()[0])
            print "Error opening %s. This file will be ignored" %bpch_filename
    
    print "Found %i valid bpch files" %len(good_files)
    
    raw_input("Press enter to proceed-->")
    
    if good_files == []:
        #if there are no valid bpch files
        return(export_dict) #leave the function with empty dictionary.

    #Now, inspect the first file in the list to get dimensions and variables.
    
    f = bpch(good_files[0])    
    variables_list = list(f.variables)
    
    num_variables = len(variables_list)
    using_variables_list = []
    done = False
    while done == False:
        clearscreen()
        print "There are %i different variables in the GEOS Chem output" %num_variables
        print "Write the name of a variable to add it to the list of variables to read"
        print "To see the entire list of variables, write 'V'"
        print "To see the different groups of variables, write 'G'"
        print "To see a list of all variables within a group (and, optionally, add them all), write the group name followed by a *"
        print "You have currently chosen the following variables:"
        print using_variables_list
        print "To remove an item from this list, type @ followed by its name"
        print "To clear this list, type 'C'"
        print "Latitude, longitude and time variables will be used automatically"
        print "Once this list is complete, type 'Y' to proceed"
        option = raw_input("-->")
        if   option.upper() == "V":
            for variable in variables_list:
                print variable
            raw_input("Press enter to continue-->")
        elif option.upper() == "G":
            groups = list(f.groups)
            for group in groups:
                print group
            raw_input("Press enter to continue-->")
        elif option.endswith("*"):
            group_variables_list = [variables_list[i] for i in range(0,num_variables) if variables_list[i].startswith(option[:(len(option)-1)])]
            if group_variables_list == []:
                print "This is not a valid group"
            else:
                for variable in group_variables_list:
                    print variable
            add_all = raw_input("To add all of these variables, type A, otherwise, press enter to continue-->").upper()
            if add_all == "A":
                for variable in group_variables_list:
                     if not(option in using_variables_list): #don't allow duplication 
                        using_variables_list.append(variable)
        elif option.upper() == "Y":
            done = True
        elif option.upper() == "C":
            using_variables_list = []
        elif option.startswith("@"):
            try:
                using_variables_list.remove(option[1:])
            except ValueError:
                print option[1:] + " is not in the list"
                raw_input("Press enter to continue-->")
        else: #if adding a variable
            if option in variables_list:
                if not(option in using_variables_list): #don't allow duplication 
                    using_variables_list.append(option)
            else:
                print "%s is not a valid variable" %option
                raw_input("Press enter to continue-->")
    
    del f
    
    #now to process this list
    
    #if the user hasn't chosen any variables
    if using_variables_list == []:
        print "No variables chosen"
        raw_input("Press enter to continue-->")
        return(export_dict)
    
    #initialise the dataset holders and assign name
    geos_datasets = []
    i = 0
    for variable in using_variables_list:
        geos_datasets.append(geos_data())
        geos_datasets[i].name = variable #set machine name of variables
        #prompt the user for "human name"
        print "Please enter a short 'human readable' name for the variable %s, (press enter to just use '%s'" %(variable,variable)
        a = raw_input("-->")
        if a == "":
            geos_datasets[i].human_name = variable
        else:
            geos_datasets[i].human_name = a
        del a
        i += 1
        
    del i
    
    if columns:
        #print "For 3D data, compute the column (C) or take a level (enter number)?"
        while True:
            c_or_lev = raw_input("For 3D data, compute the column (C) or take a level (enter number)? -->").upper()
            if c_or_lev == "C":
                break
            else:
                try:
                    _ = int(c_or_lev)
                    break
                except ValueError:
                    print "not a valid option!"
                    pass
                                  
    else:
        c_or_lev = "NA"
        
    
    #Right, now let's get onto the serious stuff
    for this_file in good_files: #access each file in turn
        print "Accessing %s" %this_file
        f = bpch(this_file)
        tau_time = list(f.variables['time'])
        lat = list(f.variables['latitude'])
        lon = list(f.variables['longitude'])
        num_times = len(tau_time)
        
        #convert times to datetime
        time = []
        for i in range(0,len(tau_time)):
            time.append(tau_to_datetime(tau_time[i]))
        
        print "There are %i different times recorded in this file" %num_times
        for geos_dataset in geos_datasets: 
            this_name = geos_dataset.name
            print "Reading variable %s" %this_name
            this_data = f.variables[this_name]
            
            #set units and lat/lon
            #we might end up setting these multiple times but
            #that's not a problem big deal
            if columns:
                geos_dataset.unit = "molec.cm-2"
            else:
                geos_dataset.unit = this_data.units
            tau_time_bounds = f.variables['time_bounds']
            dt_time_bounds = [ [tau_to_datetime(a),tau_to_datetime(b)] for [a,b] in tau_time_bounds]
            geos_dataset.set_lat_lon(lat,lon)
            
            if len(this_data[0]) != 1: #if there is vertical information
                if columns and c_or_lev == "C":
                    print "Calculating column of this variable"                   
                    this_data_vv = np.array(this_data) * 1e-9 #convert to v/v
                    box_height = np.array(f.variables['BXHGHT-$_BXHEIGHT']) * 100 #box height in cm
                    air_den =    np.array(f.variables['TIME-SER_AIRDEN'])  # air density molec.cm-3
                    air_amount = np.multiply(box_height,air_den) #air per grid box molec.cm-2
                    this_data_col = np.sum(np.multiply(this_data_vv,air_amount),axis=1) #column molec.cm-3                    
                else:    
                    print "For this 3D variable, only one layer information will be read"
            else:
                print "Variable is 2D"
            for t in range(0,num_times):
                if columns and c_or_lev == "C":
                    geos_dataset.add_data(this_data_col[t],time[t],new_time_bounds=dt_time_bounds[t])
                elif columns:
                    geos_dataset.add_data(this_data[t][int(c_or_lev)],time[t],new_time_bounds=dt_time_bounds[t])
                else:                
                    geos_dataset.add_data(this_data[t][0],time[t],new_time_bounds=dt_time_bounds[t])
    
    print "Reading geos_chem variables done"
    #OK. let's turn this into a dictionary
    
    print "Dataset created with the following variables:" 
    for geos_dataset in geos_datasets:
        print geos_dataset.human_name
        export_dict[geos_dataset.human_name] = geos_dataset
        
    return(export_dict)

def tau_to_datetime(tau_time):
    """Converts time from the TAU used in GEOS_Chem files to datetime"""
    
    #TAU is count of hours since 1985-01-01 00:00:00
    datum = dt(1985,01,01,0,0,0)
    hours_to_add = td(hours=tau_time)
    return(datum+hours_to_add)

def straighten_geos(geos_chem_var,time_option):
    """Converts GEOS Chem variable into 1D lists"""
    lat = geos_chem_var.lat
    lenlat = len(lat)
    lat1D = []
    
    lon = geos_chem_var.lon
    lenlon = len(lon)
    lon1D = []
    
    data = geos_chem_var.data[time_option]
    data1D = []
    for j in range(0,lenlat):
        for i in range(0,lenlon):
            data1D.append(data[j,i])
            lat1D.append(lat[j])
            lon1D.append(lon[i])
    
    return(lat1D,lon1D,data1D)

            
def plot_geos_chem(geos_chem_var_dict,colmap):
    """Plots GEOS Chem data using plt_grid_from_list"""
    clearscreen()
    valid_input = False
    while valid_input == False:
        
        print "Enter human name of GEOS Chem var to plot"
        print "The following datasets are saved:"
        for key in geos_chem_var_dict.keys():
            print key
        name_choice = raw_input("-->") 
        try:
            geos_chem_var = geos_chem_var_dict[name_choice]
            valid_input = True
        except(KeyError): #if not a valid key option
            print "%s is not a currently saved GEOS Chem dataset" %name_choice
    
    clearscreen()
    print "Preparing to plot GEOS Chem variable %s" %geos_chem_var.human_name
    
    #if there is more than one valid time, ask the user to select one
    if len(geos_chem_var.time) > 1:
        print "Select a time to plot"
        i = 0
        for i in range(0,len(geos_chem_var.time)):
            print "[%i] %s" %(i,str(geos_chem_var.time[i]))
        time_choice = int(raw_input("-->"))
    else:
        time_choice = 0
    
       
    (lat1D,lon1D,data1D) = straighten_geos(geos_chem_var,time_choice)
    
    title_for_plot = geos_chem_var.human_name + " at " + str(geos_chem_var.time[time_choice])
    
    plot_grid_from_list(lat1D,lon1D,data1D,
                        geos_chem_var.lon_spacing,geos_chem_var.lat_spacing,
                        box(max(lat1D),min(lat1D),max(lon1D),min(lon1D)),
                        title=geos_chem_var.human_name,
                        vmin=min(data1D),vmax=max(data1D),
                        lab=geos_chem_var.unit,
                        save=False,save_filename="nofilenamespecified",
                        colmap=colmap)

def stats_from_da(da):
    """calculate simple stats from a data_all object"""
      
    while True: #allow for break to exit
        s1_menu_title = "Which dataset do you wish to calculate statistics for?"
        [s1_menu_text,s1_answers_dict] = da.make_menu()
        s1_menu_choice = basic_menu(s1_menu_title,
                                    s1_menu_text,
                                    quit_option=True)
        if s1_menu_choice == "Z":
            break
        
        data_key = s1_answers_dict[s1_menu_choice]
        
        while True: #allow for break to exit   
            s2_menu_title = "Which statistical operator do you want to compute?"
            s2_menu_text=[["1","mean"],
                           ["2","standard deviation"],
                           ["3","count"],
                           ["4","percentiles"]]
            s2_menu_choice= basic_menu(s2_menu_title,
                                       s2_menu_text,
                                       quit_option=True)
            if s2_menu_choice == "Z":
                break
            
            #a cumbersome but robust way of setting stat_choice_text
            for pair in s2_menu_text:
                if pair[0] == s2_menu_choice:
                    stat_choice_text = pair[1]
            
            stat = calc_statistic(da.data[data_key].val,s2_menu_choice)
            
            print "Dataset: %s" %da.data[data_key].description
            print "Statistic: %s" %stat_choice_text
            if s2_menu_choice == "4":
                for i in range(0,101):
                    print "%i percentile: %g" %(i,stat[i])
            else:
                print "Result: %g" %stat
            
            _ = raw_input("Press enter to continue-->")    
    
    
def match_index(bounds,item_to_match):
    """Finds the first index where item_to_match is within bounds"""
    #bounds -> list of pairs [lower,higher] of values
    for i in range(0,len(bounds)):
        if bounds[i][0] <= item_to_match < bounds[i][1]:
            return(i)
    #if no match
    return(None)


def associate_ind_geos(ida,geos_dict):
    """associates geos_chem data with points in ida"""
    
    ida_len = len(ida.lat)
    
    #will assume that all data in geos_dict has the same
    #lat/lon/time arrangement
    
    geos_keys = geos_dict.keys()
    
    print "Getting indicies"
    geos_lat = geos_dict[geos_keys[0]].lat    
    geos_lat_bounds = []
    for i in range(0,len(geos_lat)):
        geos_lat_bounds.append([geos_lat[i]-0.5*geos_dict[geos_keys[0]].lat_spacing,
                                geos_lat[i]+0.5*geos_dict[geos_keys[0]].lat_spacing]
                              )
    
    geos_lon = geos_dict[geos_keys[0]].lon    
    geos_lon_bounds = []
    for j in range(0,len(geos_lon)):
        geos_lon_bounds.append([geos_lon[j]-0.5*geos_dict[geos_keys[0]].lon_spacing,
                                geos_lon[j]+0.5*geos_dict[geos_keys[0]].lon_spacing]
                              )
                    
    
    #time bounds should already be in the file
    geos_time = geos_dict[geos_keys[0]].time
    geos_time_bounds = geos_dict[geos_keys[0]].time_bounds
    
    #Now, for every point in ida, find the lat,lon and time indexes
    #return None if there's no match
    match_indexes = []
    success_match = 0 #keep a counter of sucessful matches 
    for p in range(0,ida_len):
        
        match_indexes.append(
            [match_index(geos_time_bounds,ida.time[p]),
             match_index(geos_lat_bounds ,ida.lat[p] ),
             match_index(geos_lon_bounds ,ida.lon[p] )]
             )
        if not(None in match_indexes[p]):
            success_match += 1 
    
    
    print "Matched %i out of %i datapoints" %(success_match,ida_len)
    
    #now add these data to ida         
    for geos_key in geos_keys:
        print "Associating "+geos_key
        this_val = []
        for p in range(0,ida_len):
            if None in match_indexes[p]: #if there's no matching point
                this_val.append(np.nan)
            else:
                this_val.append(geos_dict[geos_key].data[ match_indexes[p][0] ][ match_indexes[p][1] ][ match_indexes[p][2] ])
                
        new_d = d(this_val,
                  geos_dict[geos_key].name,
                  geos_dict[geos_key].human_name)
        new_d.unit = geos_dict[geos_key].unit
        new_d.add_meta('Source','GEOS Chem output associated with core dataset')
        ida.data[geos_key] = copy.deepcopy(new_d) #add to ida
        del this_val
        del new_d 
    
    
    _=raw_input("Press enter to continue-->")    
    return(ida)
    
def bin_extra(ida,bda):
    """Allows the user to create a binned dataset additional to the ones automatically created"""
    #useful for associated sets.
    while True: #allow for break to exit
        be1_menu_title = "Which dataset do you wish to bin?"
        [be1_menu_text,be1_answers_dict] = ida.make_menu()
        be1_menu_choice = basic_menu(be1_menu_title,
                                    be1_menu_text,
                                    quit_option=True)
        if be1_menu_choice == "Z":
            break    
            
        data_key = be1_answers_dict[be1_menu_choice]
        while True: #allow for break to exit   
            be2_menu_title = "Which statistical operator do you want to compute for each bin?"
            be2_menu_text=[["1","mean"],
                           ["2","stdev"],
                           ["3","count"],
                           ["4","combined measurement error"]]
            be2_menu_choice= basic_menu(be2_menu_title,
                                       be2_menu_text,
                                       quit_option=True)
            if be2_menu_choice == "Z":
                break
            
            #a cumbersome but robust way of setting stat_choice_text
            for pair in be2_menu_text:
                if pair[0] == be2_menu_choice:
                    stat_choice_text = pair[1]
           
            (lat_list,lon_list,this_binned,null1,null2,null3) = \
               binner(ida.lat,ida.lon,
                      ida.data[data_key].val,
                      bda.meta["Area"],
                      stat=stat_choice_text,
                      xdim=bda.meta["Binning dimensions"][0],
                      ydim=bda.meta["Binning dimensions"][1],
                      do_extras=True)
            
            del null1,null2,null3
               
            this_name = ida.data[data_key].name
            this_description = ida.data[data_key].description
            
            this_binned = grid_match(bda,this_binned,lat_list,lon_list)
            if this_binned==None:
                _ = raw_input("Matching failed. Press enter to continue -->")
                return(bda)
            
            
            bda.data[this_name] = d(this_binned,this_name,description=this_description+": "+stat_choice_text)
            #units are same as in ida, unless this is a count (then it's unitless)
            if stat_choice_text in ["mean","standard deviation"]:
                bda.data[this_name].unit = ida.data[data_key].unit
            else:
                bda.data[this_name].unit = ""
            
            return(bda)           

def grid_match(bda,new_data,new_lat,new_lon):
    #if new_data is created by gridding onto the grid of bda, there will be issues if bda does not cover a complete rectangular area.
    #so, this routine goes point-by-point in bda, and finds the associated point in new_data
    
    matched_data = []
    j = 0
    new_data_len = len(new_data)
    for i in range(len(bda.lat)):
        (lat_to_match,lon_to_match) = (bda.lat[i],bda.lon[i])
        cycle_round=0
        while True:
            if bda.lat[i] == new_lat[j] and bda.lon[i] == new_lon[j]:
                matched_data.append(new_data[j])
                cycle_round=0
                break
            else:
                j+=1
                if j >= new_data_len:
                    j=0
                    cycle_round+=1
                    if cycle_round>=2:
                        print 'Cannot match grid point for lat=%g, lon=%g' %(bda.lat[i],bda.lon[i])
                        return(None)
    #return new data as list
    return(matched_data)
    
                    
        

def criteria_filtering(ida):
    """Filter based on values"""
    cf1_menu_title = "FILTER:\n"\
                      "Filter based on which dataset?:"
    [cf1_menu_text,cf1_answers_dict] = ida.make_menu()
    cf1_menu_choice = basic_menu(cf1_menu_title,
                                cf1_menu_text,
                                quit_option=True)
    if cf1_menu_choice == "Z":
        return(ida)
    cf1_key =     cf1_answers_dict[cf1_menu_choice]    
    cf_var =      ida.data[cf1_key].val
    cf_var_name = ida.data[cf1_key].description
    
    cf2_menu_choice = basic_menu(
                      "FILERING SETTING:",
                      [["1","Keep values of %s which MATCH a specified value"%cf_var_name],
                       ["2","Keep values of %s which DO NOT MATCH a specified value"%cf_var_name],
                       ["3","Keep values of %s which are LESS THAN a specified value"%cf_var_name],
                       ["4","Keep values of %s which are MORE THAN a specified value"%cf_var_name]]
                       )
    if cf1_menu_choice.upper() == "Z":
        return(ida)
    
    cf3_menu_choice = input("Enter the critical value-->")
    pre_f_n = len(ida.lat)
    print "%i values prior to filtering."%pre_f_n
    
    filter_TF = []
    if cf2_menu_choice == "1":     
        for i in range(0,len(ida.lat)):
            if ida.data[cf1_key].val[i] == cf3_menu_choice:
                filter_TF.append(True)
            else:
                filter_TF.append(False)
    elif cf2_menu_choice == "2":     
        for i in range(0,len(ida.lat)):
            if ida.data[cf1_key].val[i] != cf3_menu_choice:
                filter_TF.append(True)
            else:
                filter_TF.append(False)
    elif cf2_menu_choice == "3":     
        for i in range(0,len(ida.lat)):
            if ida.data[cf1_key].val[i] < cf3_menu_choice:
                filter_TF.append(True)
            else:
                filter_TF.append(False)
    elif cf2_menu_choice == "4":     
        for i in range(0,len(ida.lat)):
            if ida.data[cf1_key].val[i] > cf3_menu_choice:
                filter_TF.append(True)
            else:
                filter_TF.append(False)
        
    ida.filter_all(filter_TF)
    
    post_f_n = len(ida.lat)
    print "%i values post filtering."%post_f_n
    _=raw_input("Press enter to continue-->")
    return(ida)

def betwixt(a,b,s):
    """Gets the substring between two substrings in s"""    
    return(s[s.index(a)+len(a):s.index(b)])

def only_months(ida):
    chosen_month = input("Enter month (as number) to highlight-->")
    keep = [ida.time[i].month == chosen_month for i in range(0,len(ida.lat)) ]
    ida.filter_all(keep)
    return(ida)

def associate_CSV(CSV_dir,map_box,bda,earliest,latest,CSV_name,CSV_desc,new_units=""):
    
    #check time bounds
    #earliest = min(ida.time)
    #latest   = max(ida.time)
    
    
    #read NDVI data
    (lat_CSV,lon_CSV,data_CSV,year_CSV,month_CSV) = CSV_months(earliest,latest,CSV_dir,0.1,0.1,map_box,CSV_name)
    
    if CSV_name != "Land classification":
        CSV_dt = []
        for i in range(0,len(data_CSV)):
            CSV_dt.append(dt(year_CSV[i],month_CSV[i],1,0,0,0))
    
        CSV_da = d_all(lat_CSV,lon_CSV,time=CSV_dt)
    else:
        CSV_da = d_all(lat_CSV,lon_CSV)
        
    CSV_da.add_data(d(data_CSV,CSV_name,description=CSV_desc)) 
    if CSV_name=="Land classification":
        (lat_list,lon_list,this_binned,null1,null2,null3) = \
           binner(CSV_da.lat,CSV_da.lon,
                  CSV_da.data[CSV_name].val,
                  bda.meta["Area"],
                  stat="mode",
                  xdim=bda.meta["Binning dimensions"][0],
                  ydim=bda.meta["Binning dimensions"][1],
                  do_extras=True)
    else:
        (lat_list,lon_list,this_binned,null1,null2,null3) = \
               binner(CSV_da.lat,CSV_da.lon,
                      CSV_da.data[CSV_name].val,
                      bda.meta["Area"],
                      stat="mean",
                      xdim=bda.meta["Binning dimensions"][0],
                      ydim=bda.meta["Binning dimensions"][1],
                      do_extras=True)
    
    this_binned = grid_match(bda,this_binned,lat_list,lon_list)
    if this_binned==None:
        _ = raw_input("Matching failed. Press enter to continue -->")
        return(bda)   
       
    this_name = CSV_da.data[CSV_name].name
    this_description = CSV_da.data[CSV_name].description
    
    if this_name in bda.data: #if we aready have a dataset with this name, delete it.
        del bda.data[this_name]
    bda.data[this_name] = d(this_binned,this_name,description=this_description+": mean")
    bda.data[this_name].unit = new_units
               
    return(bda)

def write_as_csv(xda):
    csv_out = raw_input("Enter full path of location for CSV file-->")
    wf = open(csv_out,'w')
    lat_to_write = np.array(xda.lat)
    lon_to_write = np.array(xda.lon)
    
    do_time = False
    #if len(xda.time) != 0:
    #    do_time = True
    #    time_to_write = np.array(xda.time)
    #else:
    #    do_time = False
        
    selectvar_menu_title = "Choose variable:"
    [selectvar_menu_text,selectvar_answers_dict] = xda.make_menu()
    selectvar_menu_choice = basic_menu(selectvar_menu_title,
                                  selectvar_menu_text,
                                  quit_option=True)
    if selectvar_menu_choice == "Z":
        return
    
    the_key =      selectvar_answers_dict[selectvar_menu_choice]   
    the_var =      np.array(xda.data[the_key].val)   
    for i in range(len(xda.lat)):
        if do_time:
            this_time_to_write = time_to_write[i]
            wf.write(str(this_time_to_write.year))
            wf.write(",")
            wf.write(str(this_time_to_write.month))
            wf.write(",")
            wf.write(str(this_time_to_write.day))
            wf.write(",")
            wf.write(str(this_time_to_write.hour))
            wf.write(",")
            wf.write(str(this_time_to_write.minute))
            wf.write(",")
            wf.write(str(this_time_to_write.second))
            wf.write(",")
        wf.write(str(lat_to_write[i]))
        wf.write(",")
        wf.write(str(lon_to_write[i]))
        wf.write(",")
        wf.write(str(the_var[i]))
        wf.write("\n")
    wf.close()
    return
    
def binned_to_indv(ida,bda):
    """For each individual data point, get a set of data from the the closest binned point"""
    
    bda_lats = np.array(bda.lat)
    bda_lons = np.array(bda.lon)
    
    print "Identifying datasets"
    assoc_dict = {}
    for key in bda.data.keys():
        print key
        assoc_dict["fb_"+key] = []
      
    print "Computing indexes"
    li = len(ida.lat)
    for a in range(li):
        #print "%i of %i" %(a,li)
        #get index of smallest distance from lat/lom
        b = np.argmin(np.add(np.square(bda_lats - ida.lat[a]),np.square(bda_lons - ida.lon[a])))
        
        for key in bda.data.keys():
             assoc_dict["fb_"+key].append(bda.data[key].val[b])
             
    print "Adding to store"
    for key in bda.data.keys():
        ida.add_data(d(assoc_dict["fb_"+key],"fb_"+key,description=bda.data[key].description+" (from binned)"))
        
    return(ida)
        
def second_stage(gradient,err,intercept,ida):
    
    HCHO_obs = ida.data['sat_VC'].val
    HCHO_errs = ida.data['sat_DVC'].val
    
    isop_inverted = np.divide(np.subtract(HCHO_obs,intercept),gradient)
    isop_inverted_low = np.divide(np.subtract(HCHO_obs,intercept),gradient-err)
    isop_inverted_high = np.divide(np.subtract(HCHO_obs,intercept),gradient+err)
    
    #isop_inverted_errs = np.multiply(isop_inverted,
    #                       np.divide(HCHO_errs,np.absolute(HCHO_obs)))
                           
    ida.add_data(d(isop_inverted,"isop_inv",description="Isoprene emissions derived from inversion"))
    ida.data["isop_inv"].units="C.cm-2.s-1"
    ida.add_data(d(isop_inverted_low,"isop_inv_low",description="Isoprene emissions derived from inversion - low estimate"))
    ida.data["isop_inv_low"].units="C.cm-2.s-1"
    ida.add_data(d(isop_inverted_high,"isop_inv_high",description="Isoprene emissions derived from inversion - high estimate"))
    ida.data["isop_inv_high"].units="C.cm-2.s-1"
    
    #ida.add_data(d(isop_inverted,"isop_inv_errs",description="Error in isoprene emissions derived from inversion"))
    #ida.data["isop_inv_errs"].units="C.cm-2.s-1"
   
    return(ida)        
    
def daily_fire_filter(ida,fire_file,start_date,end_date,map_box,xdim,ydim):
    
    do_a_save = raw_input("If you want to automatically save the output from this, enter the full filepath here. Otherwise, enter for normal operations-->")
    
    #filter_range = float(raw_input("Minimum acceptable distance from fire in degrees?-->"))
    #filter_t_range = float(raw_input("Minimum acceptable time from fire in days?-->"))
    filter_t_range = 1
        
    fire_lats = []
    fire_lons = []
    fire_dates = []
    
    confidence_threshold = 0.
    print "Considering all detections with a confidence of %g or greater"%confidence_threshold
    
    #Read in fire file
    f = open(fire_file,"rb")
    data =csv.reader(f, delimiter=',')
    
    l = 0
    for line in data: 
        if l == 0: #skip first line (headers)
            l+=1
                        
        elif float(line[9]) < confidence_threshold:
            l+=1
            
        else:
            fire_lats.append(float(line[0])) #lat
            fire_lons.append(float(line[1])) #lon
            fire_dates.append(dt.strptime(line[5],"%Y-%m-%d")) #dates
            l+=1
        if l%10000 == 0:
            print "line %i"%l
            
    f.close()
    
    num_fires = len(fire_lats)
    
    print "%i fires in file above confidence limit" %num_fires
    
    #fire_lats = np.array(fire_lats)
    #fire_lons = np.array(fire_lons)
    #fire_dates = np.array(fire_dates)
    
    #now get a grid for every date from the beginning to end of the date range (plus buffer)
    
    list_of_dates = []
    
    list_of_grids = []
        
    date_clock = start_date
    
    print "Binning"
    while date_clock <= end_date:
        list_of_dates.append(date_clock)
        print date_clock
        
        accepted_dates = [date_clock + td(days=x) for x in range(-filter_t_range,filter_t_range+1)]
        fire_indexes = [i for i in range(num_fires) if fire_dates[i] in accepted_dates]
        
        these_lats = [fire_lats[i] - ydim for i in fire_indexes] + [fire_lats[i] - ydim for i in fire_indexes] + [fire_lats[i] - ydim for i in fire_indexes] + \
                     [fire_lats[i]        for i in fire_indexes] + [fire_lats[i]        for i in fire_indexes] + [fire_lats[i]        for i in fire_indexes] + \
                     [fire_lats[i] + ydim for i in fire_indexes] + [fire_lats[i] + ydim for i in fire_indexes] + [fire_lats[i] + ydim for i in fire_indexes]     
        these_lons = [fire_lons[i] - xdim for i in fire_indexes] + [fire_lons[i]        for i in fire_indexes] + [fire_lons[i] + xdim for i in fire_indexes] + \
                     [fire_lons[i] - xdim for i in fire_indexes] + [fire_lons[i]        for i in fire_indexes] + [fire_lons[i] + xdim for i in fire_indexes] + \
                     [fire_lons[i] - xdim for i in fire_indexes] + [fire_lons[i]        for i in fire_indexes] + [fire_lons[i] + xdim for i in fire_indexes]
        dummy      = [1.           for i in fire_indexes] * 9
        
        list_of_grids.append(binner(these_lats,these_lons,dummy,map_box,stat="count",xdim=xdim,ydim=ydim,do_extras=False))
        
        if date_clock == start_date:
            (bin_lats,bin_lons,dummy2,xedges,yedges,index_map) = (binner(these_lats,these_lons,dummy,map_box,stat="count",xdim=xdim,ydim=ydim,do_extras=True))
               
        date_clock += td(days=1)
    
    
    print "Now matching up observations"
    li = len(ida.lat)
    fire_count = []
    for a in range(li):
        if a%100 == 0:
            print "data point %i of %i"%(a,li)
            
        #print "%i of %i" %(a,li)
        #get index of smallest distance from lat/lom
        b  = np.argmin(np.add(np.square(bin_lats - (ida.lat[a]       )),np.square(bin_lons - (ida.lon[a]       ))))
        
        t = (ida.time[a]-start_date).days
        
        fire_count.append(list_of_grids[t][b])
        
    ida.add_data(d(fire_count,"Fire screen",description="Fire within range 1 grid box and 1 day?"))
    
    if do_a_save != "":
       print "Saving ida as %s" %do_a_save
       cPickle.dump(ida,open(do_a_save,"wb"))
       print "Pickle saved to %s" %do_a_save
       
       
def seasonal(perm_ida,map_box,xdim,ydim,colmap,save_path,
             do_options=["obs_map","mod_map","inversion"]):
    print "Doing seasonal processing on data"
    
    #if we are doing inversions, ask user for X and Y variables
    if "inversion" in do_options:
         #X       
         invx_menu_title = "For inversion, which is the independent (X axis) varible?"
         [invx_menu_text,invx_answers_dict] = perm_ida.make_menu()
         invx_menu_choice = basic_menu(invx_menu_title,
                                        invx_menu_text,
                                        quit_option=False)
         invx_key = invx_answers_dict[invx_menu_choice]
         
         #Y
         invy_menu_title = "For inversion, which is the dependent (Y axis) varible?"
         [invy_menu_text,invy_answers_dict] = perm_ida.make_menu()
         invy_menu_choice = basic_menu(invy_menu_title,
                                        invy_menu_text,
                                        quit_option=False)
         invy_key = invx_answers_dict[invy_menu_choice]
    
    for season in ["all","winter","pre-monsoon","monsoon","post-monsoon"]:
        print "========"
        #define months in season as season_months
        if season == "all":
            print "No seasonal selection"
            season_months = [1,2,3,4,5,6,7,8,9,10,11,12]
        elif season == "winter":
            print "Winter season"
            season_months = [1,2]
        elif season == "pre-monsoon":
            print "Pre-monsoon season"
            season_months = [3,4,5]
        elif season == "monsoon":
            print "Monsoon season"
            season_months = [6,7,8,9]
        elif season == "post-monsoon":
            print "Post-monsoon season"
            season_months = [10,11,12]
        
        #make a deep copy of ida
        print "making a copy..."
        this_ida = copy.deepcopy(perm_ida)
        
        #now filter it
        print "filtering for season..."
        filter_list = [this_ida.time[i].month in season_months for i in range(len(this_ida.time)) ]
        this_ida.filter_all(filter_list)
        
        #are we doing either mapping options?
        if "obs_map" in do_options or "mod_map" in do_options:
            #if so, standard binning.
            print "Binning..."
            this_bda = create_binned_set(this_ida,map_box,xdim,ydim)
            
            #make obs_map
            if "obs_map" in do_options:
                save_filename = "%s/obscol_%s.png" %(save_path,season)
                print "Creating %s" %save_filename
                plot_grid_from_list(this_bda.lat,this_bda.lon,
                                    this_bda.data["sat_VC_mean"].val,
                                    this_bda.meta["Binning dimensions"][0],
                                    this_bda.meta["Binning dimensions"][1],
                                    map_box,title="Observed HCHO vertical column, season = %s"%season,
                                    vmin=0,vmax=2e16,lab="molec.cm-2",
                                    save=True,
                                    save_filename=save_filename,
                                    colmap=colmap,do_plot=False)
            
            #make mod_map
            if "mod_map" in do_options:
                save_filename = "%s/modcol_%s.png" %(save_path,season)
                print "Creating %s" %save_filename
                plot_grid_from_list(this_bda.lat,this_bda.lon,
                                    this_bda.data["geos_VC_mean"].val,
                                    this_bda.meta["Binning dimensions"][0],
                                    this_bda.meta["Binning dimensions"][1],
                                    map_box,title="Modelled HCHO vertical column, season = %s"%season,
                                    vmin=0,vmax=2e16,lab="molec.cm-2",
                                    save=True,
                                    save_filename=save_filename,
                                    colmap=colmap,do_plot=False)
            
        
        #are we doing inversions?
        if "inversion" in do_options:
            print "Inversion"
            
            (nonan_x_var,nonan_y_var) = \
                            stripallnans(this_ida.data[invx_key].val,this_ida.data[invy_key].val)
            
            print "Linear regression:"
            slope, intercept, r_value, p_value, std_err = scistats.linregress(nonan_x_var,nonan_y_var)
            print "Best-fit line: y = %+.2gx%+.2g " %(slope,intercept)
            r_squared = r_value * r_value
            print "R-squared    : %.4g " %r_squared
            print "P-value      : %.4g " %p_value
            print "Standard err : %g "   %std_err        
            pearsonr = scistats.pearsonr(nonan_x_var,nonan_y_var)
            print "Pearson's correlation coefficient:\n %f, 2-tailed p-value: %f" %(pearsonr[0],pearsonr[1])
            polyfit = np.polyfit(nonan_x_var, nonan_y_var, 1)
            
            this_ida = second_stage(slope,std_err,intercept,this_ida)
            
            print "Prior emission rate: %g" %np.nanmean(this_ida.data[invx_key].val)
            
            frac_err = np.nanmean(this_ida.data["sat_DVC"].val) / np.nanmean(this_ida.data["sat_VC"].val)
            count = np.count_nonzero(~np.isnan(this_ida.data["isop_inv"].val))
            sqrtcount = math.sqrt(count)
            
            pm = np.nanmean(this_ida.data["isop_inv"].val)*(frac_err/sqrtcount)
            
            #print "Post emission rate, low: %g" %np.nanmean(this_ida.data["isop_inv_low"].val)
            print "Post emission rate: %g" %np.nanmean(this_ida.data["isop_inv"].val)
            print "plus/minus --> %g" %pm
            print pm
            #print "Post emission rate, low: %g" %np.nanmean(this_ida.data["isop_inv_high"].val)            
        del this_ida
            
def simple_maths(xda):
    """Simple mathematical operators on one or two fields"""
    
    #should work for both ida and bda objects
    
    #operation
    sm_menu_title = "Choose type of operation (A=first, B=second)"
    sm_menu_text=[["1","A + B"],
                  ["2","A - B"],
                   ["3","A * B"],
                   ["4","A / B"],
                   ["5","A + number"],
                   ["6","A * number"]]
    sm_menu_choice= basic_menu(sm_menu_title,
                                sm_menu_text,
                                quit_option=True)
    if sm_menu_choice == "Z":
        return xda
    
    #First dataset
    d1_menu_title = "Select dataset A:"
    [d1_menu_text,d1_answers_dict] = xda.make_menu()
    d1_menu_choice = basic_menu(d1_menu_title,
                                  d1_menu_text,
                                  quit_option=False)
    d1_key =      d1_answers_dict[d1_menu_choice]    
    d1_var =      xda.data[d1_key].val
    d1_unit = xda.data[d1_key].unit
    
    if sm_menu_choice in ["1","2","3","4"]:
        #Second dataset
        d2_menu_title = "Select dataset B:"
        [d2_menu_text,d2_answers_dict] = xda.make_menu()
        d2_menu_choice = basic_menu(d2_menu_title,
                                      d2_menu_text,
                                      quit_option=False)
        d2_key =      d2_answers_dict[d2_menu_choice]    
        d2_var =      xda.data[d2_key].val
        d2_unit = xda.data[d2_key].unit
    elif sm_menu_choice == ["5"]:
        flat_num = input("Enter number to ADD to each datapoint")
    else:
        flat_num = input("Enter number to MULTIPLY each datapoint by")
        
    #Perform operation
    if sm_menu_choice == "1":
        new_var = np.add(d1_var,d2_var)
        new_var_sname = d1_key+"_PLUS_"+d2_key
    elif sm_menu_choice == "2":
        new_var = np.subtract(d1_var,d2_var)
        new_var_sname = d1_key+"_MINUS_"+d2_key
    elif sm_menu_choice == "3":
        new_var = np.multiply(d1_var,d2_var)
        new_var_sname = d1_key+"_TIMES_"+d2_key
    elif sm_menu_choice == "4":
        new_var = np.divide(d1_var,d2_var)
        new_var_sname = d1_key+"_DIV_"+d2_key
    elif sm_menu_choice == "5":
        new_var = np.add(d1_var,flat_num)
        new_var_sname = d1_key+"_PLUS_"+str(flat_num)
    elif sm_menu_choice == "6":
        new_var = np.multiply(d1_var,flat_num)
        new_var_sname = d1_key+"_TIMES_"+str(flat_num)
    
    print "Operation successful"
    
    #Get new name and units
    new_var_lname = ""
    while new_var_lname == "": #keep asking the question if a blank name is specified
        new_var_lname = raw_input("Type a name for this new field --> ")
    
    if sm_menu_choice in ["1","2","6"]:
        new_var_unit = d1_unit
    else:
        new_var_unit = raw_input("Type units for this new field --> ")
    
    #create object and add to container variable
    
    new_var_sname = "smo_"+new_var_lname #hack to avoid overwriting names
    
    d_obj = d(new_var,new_var_sname,description=new_var_lname)
    d_obj.unit = new_var_unit
    xda.add_data(d_obj)
    
    return(xda)
        
def enter_date(prompt='date in YYYY-MM-DD format ->'):
    date_entry = raw_input(prompt)
    year, month, day = map(int, date_entry.split('-'))
    return(date(year, month, day))

#WORK IN PROGRESS    
#def write_netCDF(bda,map_box,xdim,ydim):
#    
#   #file to save to
#   save_file = raw_input("Enter full path for netCDF to be saved to (.nc extension recommended) ->")
#   #open file in write mode
#   dataset = Dataset(save_file,'w',format='NETCDF4_CLASSIC')
#   
#   #work out the lats and lons to have in the file - choose the mapping boundaries
#   num_lats = int((map_box.n-map_box.s)/ydim + 1)
#   num_lons = int((map_box.e-map_box.w)/xdim + 1)
#   
#   lat = dataset.createDimension('lat', num_lats)
#   lon = dataset.createDimension('lon', num_lons)
   
def bda_add_frombda(bda_old,bda_new):
    bda_old_len = len(bda_old.lat)
    bda_new_len = len(bda_new.lat)
    
    d_menu_title = "Which dataset do you wish to bring in? Quit to return to main menu"
    [d_menu_text,d_answers_dict] = bda_new.make_menu()
    d_menu_choice = basic_menu(d_menu_title,
                            d_menu_text,
                            quit_option=True)
    if d_menu_choice == "Z":
        return(bda_old)
    else:
        d_key =      d_answers_dict[d_menu_choice] 
    
    if d_key not in bda_new.data:
        print "%s not a variable in new dataset" %d_key
        _ = raw_input("Press enter to continue -->")
        return(bda_old)
    
    new_vals = np.array([])
    for i in range(bda_old_len):
    
        if (i+1)%100 == 0: #ticker
            print "%i of %i" %(i+1,bda_old_len)
            
        lat_to_find = bda_old.lat[i]
        lon_to_find = bda_old.lon[i]
        new_index = [j for j in range(bda_new_len) if bda_new.lat[j] == lat_to_find and bda_new.lon[j] == lon_to_find]
        if new_index == []:
            new_val = np.nan
            print "nan at point %i" %i
        elif len(new_index) > 1:
            print "Multiple lat/lon matches, this shouldn't happen!"
            new_val = np.nan
        else:
            new_val = bda_new.data[d_key].val[new_index[0]]
        new_vals = np.append(new_vals,new_val)
    
    print new_vals
    
    bda_old.add_data(d(list(new_vals),d_key+"_aux",description=bda_new.data[d_key].description+"_aux"))
        
    bda_old.data[d_key+"_aux"].unit = bda_new.data[d_key].unit 

           
    print "%s added" %d_key+"_aux"   
    return(bda_old)
