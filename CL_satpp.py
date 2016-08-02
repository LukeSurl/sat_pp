#!/usr/bin/env python
"""CL_satpp.py -- A command-line interface for interogating the sat_pp dataset"""

from datetime import datetime as dt
from datetime import timedelta as td
import pickle
from matplotlib.patches import Polygon
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
import sys
from bpch import bpch

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


def calc_statistic(dataset,stat_choice):
    if stat_choice == "1": #mean
        return np.nanmean(dataset)
    elif stat_choice == "2": #standard deviation
        return np.nanstd(dataset)
    elif stat_choice == "3": #count
        return np.count_nonzero(~np.isnan(dataset))
        
    
def map_preplot_menu(dataset_name,title="",vmin=0.,vmax=3.e16,units="molec.cm-2"):
    choice = "xx"    
    if title == "":
        title = dataset_name
    save_filename = title + ".png"
    while not(choice.upper() in ["P","S","Z"]): #menu-leaving options
        while not(choice.upper() in ["T","C","U","P","S","Z"]):
            os.system('clear')
            print "Preparing to map data"
            print "Dataset to plot  = " + dataset_name
            print "Title for figure = " + title
            print "Colourbar minimum = " + str(vmin)
            print "Colourbar maximum = " + str(vmax)
            print "Units text        = " + units
            print "OPTIONS:"
            print "[T] Change title"
            print "[C] Change colourbar min/max"
            print "[U] Change units text"
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
                units = change_var(units,"Units text")
            
            if choice.upper() == "S":
                save_filename = title + ".png"
                save_filename = change_var(save_filename,"Image filename to save (include extension)")
    return(choice.upper(),title,vmin,vmax,units,save_filename)
    
def geo_select_menu(type_geo_select,geo_selection):

    #type_geo_select can be either:
    # latlon
    choice = ""
    while choice not in ["1","2","3","4","Z"]    :
        os.system('clear')
        print "GEOGRAPHIC SELECTION OPTION"
        if type_geo_select == "latlon":
            print "Currently, geographic selection done by rectangular area,"
            print "%f to %f latitude, %f to %f longitude" %(geo_selection[1],geo_selection[0],geo_selection[3],geo_selection[2])
        elif type_geo_select == "circle":
            print "Currently, geographic selection done by circular area,"
            print "Centre point: %f N, %f E, Radius %f km" %(geo_selection[0],geo_selection[1],geo_selection[2])
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
            type_geo_select = "latlon"
            new_north = float(raw_input("Enter lat for N boundary (N +ve, S -ve) -->"))
            new_south = float(raw_input("Enter lat for S boundary (N +ve, S -ve) -->"))
            new_west  = float(raw_input("Enter lon for W boundary (E +ve, W -ve) -->"))
            new_east  = float(raw_input("Enter lon for E boundart (E +ve, W -ve) -->"))
            geo_selection = [new_north,new_south,new_west,new_east]
        elif choice == "2":
            type_geo_select = "circle"
            new_cen_lat = float(raw_input("Enter lat for circle centre (N +ve, S -ve) -->"))
            new_cen_lon = float(raw_input("Enter lon for circle centre (E +ve, W -ve) -->"))
            new_radius  = float(raw_input("Enter radius of circle (km) -->"))
            geo_selection = [new_cen_lat,new_cen_lon,new_radius]
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
                                    
                                 
def global_options(north,south,east,west,xdim,ydim,startdate,enddate):
    leave = False
    while leave == False :
        os.system('clear')
        print "GLOBAL OPTIONS"
        print "Current boundaries for maps etc.:"
        print "NORTH = "+str(north)
        print "SOUTH = "+str(south)
        print "EAST  = "+str(east)
        print "WEST  = "+str(west)
        print "Longitudinal spacing of gridded data = "+str(xdim)
        print "Latitudinal  spacing of gridded data = "+str(ydim)
        print "Start time = " + str(startdate.year) +"-"+str(startdate.month) + "-"+str(startdate.day)
        print "End  time = " + str(enddate.year) +"-"+str(enddate.month) + "-"+str(enddate.day)
        print "[B] Change boundary options"
        print "[S] Change spacing options"  
        print "[T] Change time frame"
        print "[Z] Return to top menu"
        choice = raw_input("-->")
        if choice == "B" or choice == "b":
            north = float(raw_input("Enter new value for NORTH:"))
            south = float(raw_input("Enter new value for SOUTH:"))
            east  = float(raw_input("Enter new value for EAST :"))
            west  = float(raw_input("Enter new value for WEST :"))
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
        elif choice == "Z" or choice == "z":
            leave = True
        
    return(north,south,east,west,xdim,ydim,startdate,enddate)

def change_pickle(current_pickle):
    os.system('clear')
    print "Currently, pickles with prefix %s are loaded." %current_pickle
    print "To change this, enter a new prefix."
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

def load_new_pickles_all(current_pickle,verbose=True,NDVI=True,indv=True,binn=True):
    if NDVI:
        NDVI_data = load_new_pickles_NDVI(current_pickle,verbose=verbose)
    if indv:
        indv_data = load_new_pickles_indv(current_pickle,verbose=verbose)
    if binn:
        binn_data = load_new_pickles_binn(current_pickle,verbose=verbose)
    
    if NDVI:
        if indv:
            if binn:
                return(NDVI_data,indv_data,binn_data) #NIB
            else:
                return(NDVI_data,indv_data)           #NIx
        else:
            if binn:
                return(NDVI_data,binn_data)           #NxB
            else:
                return(NDVI_data)                     #Nxx
    else:
        if indv:
            if binn:
                return(indv_data,binn_data)           #xIB
            else:
                return(indv_data)                     #xIx
        else:
            if binn:
                return(binn_data)                     #xxB
            else:
                return()                              #xxx                    

    
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
    
def load_new_pickles_indv(current_pickle,verbose=True):    
    #Individual observations
    if verbose:
        print("Loading pickled individual observations")
    (lat,lon,sat_VC,sat_DVC,geos_VC,time,country,state,index_map) = pickle.load( open(current_pickle + "_2.p","rb") )
    #(lat,lon,sat_VC,sat_DVC,geos_VC,time,country,state,index_map) = stripallnans(lat,lon,sat_VC,sat_DVC,geos_VC,time,country,state,index_map)
    if verbose:
        print "Pickled individual points now loaded as %i-long lists" %len(sat_VC)
        print("Latitudes             : lat       ")
        print("Longitudes            : lon       ")
        print("Timestamps            : time      ") 
        print("Observed columns      : sat_VC    ")
        print("Error in obs.         : sat_DVC   ")
        print("Modelled columns      : geos_VC   ")
        print("Country assignments   : country   ")
        print("State assignments     : state     ") 
        print("-----------------------------")

    indv_data = (lat,lon,sat_VC,sat_DVC,geos_VC,time,country,state,index_map)
    return(indv_data)

def load_new_pickles_binn(current_pickle,verbose=True):     
    #Binned data
    if verbose:
        print("Loading pickled binned data (full date range)")      
        
    (lat_binned,lon_binned,
     sat_VC_mean_binned,sat_VC_stdev_binned,
     sat_DVC_mean_binned,
     geos_VC_mean_binned,geos_VC_stdev_binned,
     NDVI_binned,
     country_binned,state_binned,
     index_map2) = pickle.load( open(current_pickle + "_binned.p","rb") )
    
    if verbose:
        print "Pickled binned data now loaded as %i-long lists" %len(lat_binned)
        print("Binned latitudes     : lat_binned         ")
        print("Binned_longitudes    : lon_binned         ")
        print("Binned obs., mean    : sat_VC_mean_binned ") 
        print("Binned obs., stdev   : sat_VC_stdev_binned")
        print("Binned obs. av. error: sat_DVC_mean_binned")
        print("Binned modelled mean : geos_VC_mean_binned")
        print("Binned modelled stdev: geos_VC_stdev_binned")
        print("Country strings      : country_binned ")
        print("State strings        : state_binned    ")
        print("-----------------------------")
    binn_data = (lat_binned,lon_binned,
     sat_VC_mean_binned,sat_VC_stdev_binned,
     sat_DVC_mean_binned,
     geos_VC_mean_binned,geos_VC_stdev_binned,
     country_binned,state_binned)
    
    print("All data loaded. Have fun with it!")
    return(binn_data)

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
                if a_entry[i] == float("nan"):
                    flag += 1
            elif type(a_entry[i]) == np.float64 :
                if np.isnan(a_entry[i]):
                    flag += 1                    
        if flag == 0 :
            truefalse.append(True)
        else:
            truefalse.append(False)            
    b = []
    for a_entry in a:
        b.append([a_entry[i] for i in range(0,checklen) if truefalse[i] ])     
    return(tuple(b))



def draw_screen_poly( rec_lats, rec_lons, m, c, color_index ):
    x, y = m( rec_lons, rec_lats )
    xy = zip(x,y)
    this_col=color_index(c,clip=True)
    poly = Polygon( xy, facecolor=[min(1.,0.+this_col*2),1.-2*abs(0.5-this_col),min(1.,2.-this_col*2)], edgecolor='none', alpha=1.0 )
    plt.gca().add_patch(poly)
    
def prepare_map(north,south,east,west):
    
    m = Basemap(projection='merc',
            llcrnrlat=south-0.1,urcrnrlat=north+0.1,
            llcrnrlon=west -0.1,urcrnrlon=east +0.1,
            lat_ts=(north+south)/2.,resolution='i')
    
    #get appropriate grid line spacing for this map
    min_dimension = min(north-south,east-west)
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
    m.drawparallels(np.arange(-90.,90.,gl_spacing))
    m.drawmeridians(np.arange(-180.,181.,gl_spacing)) 
    return m
    
def free_colorbar(vmin,vmax,label="no label selected",coltype="bwr"):
    if coltype == "bwr":
        cols_for_bar = [[0.,0.,1.],[1.,1.,1.],[1.,0.,0.]]
    elif coltype == "lg-dg":
        cols_for_bar = [[0.,0.2,0.],[0.8,0.2,0.8]]
    #can add in more cols_for_bar options here.
    cm1 = colors.LinearSegmentedColormap.from_list("MyCmapName",cols_for_bar)
    cnorm = colors.Normalize(vmin=vmin,vmax=vmax)
    cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
    cpick.set_array([])
    plt.colorbar(cpick,label=label)

def plot_grid_from_list(lat,lon,var,xdim,ydim,north,south,east,west,title="Unnamed plot",vmin=0.0e16,vmax=3.0e16,lab="HCHO column molec.cm-3",save=False,save_filename="nofilenamespecified"):
    #Define basemap
    m = prepare_map(north,south,east,west)
    color_index = colors.Normalize(vmin,vmax)
    (lat,lon,var) = stripallnans(lat,lon,var) #cut out nan values (stops crashes)    
    for i in range(0,len(lat)):
        this_w = lon[i]-xdim/2
        this_e = lon[i]+xdim/2
        this_s = lat[i]-ydim/2
        this_n = lat[i]+ydim/2
        this_c = var[i]
        if math.isnan(this_c):
            continue #don't plot if there's no data
        draw_screen_poly( [this_s,this_n,this_n,this_s],[this_w,this_w,this_e,this_e], m, this_c, color_index )        
    plt.title(title)
    free_colorbar(vmin,vmax,label=lab)
    fig = plt.gcf()
    if save:    
        fig.savefig(save_filename)
    plt.show()

def plot_dots_on_map(lat,lon,var,north,south,east,west,vmin=0.,vmax=3.0E16,title="Unnamed plot",lab="HCHO column molec.cm-3",save=False,save_filename="nofilenamespecified"):
    """Draws dots on a map at lat/lon positions, color-coded based on values"""
    m = prepare_map(north,south,east,west)
    x, y = m(lon,lat)
    m.scatter(x,y,18,marker='.',edgecolors='none',c=var, vmin=vmin, vmax=vmax)
    plt.title(title)
    free_colorbar(vmin,vmax,label=lab)
    fig = plt.gcf()
    if save:    
        fig.savefig(save_filename)
    plt.show()

    
def geo_select_regional(region_list,text_to_match,*datasets):
    """Reduces datasets down to only where the associated region matches a given text string"""
    output = []
    for dataset in datasets :  
        selected = [dataset[i] for i in range(0,len(dataset)) if region_list[i] in text_to_match]
        output.append(selected)        
    return(tuple(output))

def in_box(lat,lon,north,south,east,west):
    if lat < south :
        return False
    if lat > north :
        return False
    if lon < west :
        return False
    if lon > east :
        return False
    return True
    
def geo_select_rectangle(lat,lon,geo_selection,*datasets):
    """Reduces datasets down to only those points within designated bounds"""
    [northbound,southbound,eastbound,westbound] = geo_selection
    output = []
    for dataset in datasets :
        selected = [dataset [i] for i in range(0,len(dataset)) if in_box(lat[i],lon[i],northbound,southbound,eastbound,westbound) ]
        output.append(selected)
    return(tuple(output))

    
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
    
def time_select(startdate,enddate,time,*datasets):
    """Reduces datasets down to only where the associated region matches a given text string"""
    output = []
    for dataset in datasets :  
        selected = [dataset[i] for i in range(0,len(dataset)) if startdate <= dt.fromtimestamp(time[i]) <= enddate ]
        output.append(selected)        
    return(tuple(output))
    
def time_cycle(startdate,enddate,step_days,time,dataset,statistic):
    alpha = float(raw_input("alpha value -->"))
    time_scatter(time,dataset,alpha=alpha)
    
    clocklow = startdate
    stat_collection = []
    time_collection = []
    while clocklow <= enddate:
        clockhigh = clocklow + td(days=step_days)
        this_dataset = time_select(clocklow,clockhigh,time,dataset)
        this_stat = calc_statistic(this_dataset,statistic)
        central_time = clocklow + td(days=step_days*0.5)
        stat_collection.append(this_stat)
        time_collection.append(central_time)
        clocklow = clockhigh
    
    return(stat_collection,time_collection)

def time_scatter(time,y,yerr=[],title="UNNAMED PLOT",y_label="",x_label="",alpha=0.2):
    if yerr == []:
        plt.scatter(time, y, alpha=alpha, color='g')
    else:
        plt.errorbar(time, y, yerr=yerr, alpha=alpha, color='g')

    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(y_label, fontsize=16)
    plt.grid(b=True, which='major', color='0.65')

    plt.title(title)
    plt.show()

def save_pickle_indiv(name,ULN,sat_VC,sat_DVC,geos_VC,lat,lon,time):
        print("Pickling individual data")
        pikname = name + "_1.p"
        pickle.dump( (ULN,sat_VC,sat_DVC,geos_VC,lat,lon,time), open(pikname,"wb") )
        
def binner(lat,lon,vals,north,south,east,west,stat="mean",xdim=0.3125,ydim=0.25,do_extras=False):
    """Divides the region into a grid, and computes a statistic for the values falling within it"""

    num_xbins = int((east-west)/xdim) + 1
    num_ybins = int((north-south)/ydim) + 1
    
    #use scistats.binned_statistic_2d with the np function chosen.
    if stat == "mean" : #mean of values in each box
        (binned_stat,xedges,yedges,binnumber) = \
            scistats.binned_statistic_2d(lon,lat,vals,np.nanmean,bins=[num_xbins,num_ybins],
            range=[[west-xdim/2,east+xdim/2],[south-ydim/2,north+ydim/2]],
            expand_binnumbers=True)
    elif stat == "stdev" : #standard deviation of values in each box
        (binned_stat,xedges,yedges,binnumber) = \
            scistats.binned_statistic_2d(lon,lat,vals,np.nanstd,bins=[num_xbins,num_ybins],
            range=[[west-xdim/2,east+xdim/2],[south-ydim/2,north+ydim/2]],
            expand_binnumbers=True)
    elif stat == "sum" : #sum of all values in each box
        (binned_stat,xedges,yedges,binnumber) = \
            scistats.binned_statistic_2d(lon,lat,vals,np.nansum,bins=[num_xbins,num_ybins],
            range=[[west-xdim/2,east+xdim/2],[south-ydim/2,north+ydim/2]],
            expand_binnumbers=True)
    elif stat == "count" : #count of number of values in each box
        (binned_stat,xedges,yedges,binnumber) = \
            scistats.binned_statistic_2d(lon,lat,vals,notnancount,bins=[num_xbins,num_ybins],
            range=[[west-xdim/2,east+xdim/2],[south-ydim/2,north+ydim/2]],
            expand_binnumbers=True)
    elif stat == "quadsum" : #sum of squares of all values
        qvals = list(np.multiply(vals,vals))
        (binned_stat,xedges,yedges,binnumber) = \
            scistats.binned_statistic_2d(lon,lat,qvals,np.nansum,bins=[num_xbins,num_ybins],
            range=[[west-xdim/2,east+xdim/2],[south-ydim/2,north+ydim/2]],
            expand_binnumbers=True)
    
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
    
def create_binned_set(lat,lon,geos_VC,sat_VC,sat_DVC,north,south,east,west,xdim,ydim):
    """Given the individual data, returns a full set of binned data"""

    #Mean satellite observation 
    (lat_binned,lon_binned,sat_VC_mean_binned,xedges,yedges,index_map) = \
        binner(lat,lon,sat_VC,north,south,east,west,stat="mean",xdim=xdim,ydim=ydim,
               do_extras=True)
    #Standard deviation in satellite observation
    sat_VC_stdev_binned = \
        binner(lat,lon,sat_VC,north,south,east,west,stat="stdev",xdim=xdim,ydim=ydim)
    #Count of number of observations
    sat_VC_count_binned = \
        binner(lat,lon,sat_VC,north,south,east,west,stat="count",xdim=xdim,ydim=ydim)
    #Mean error in satellite observation
    sat_DVC_mean_binned = \
        binner(lat,lon,sat_DVC,north,south,east,west,stat="mean",xdim=xdim,ydim=ydim)    
    #Mean model observation
    geos_VC_mean_binned = \
        binner(lat,lon,geos_VC,north,south,east,west,stat="mean",xdim=xdim,ydim=ydim)
    #Standard deviation in model values
    geos_VC_stdev_binned = \
        binner(lat,lon,geos_VC,north,south,east,west,stat="stdev",xdim=xdim,ydim=ydim)

    return(lat_binned,lon_binned,xedges,yedges,index_map,
           sat_VC_mean_binned,sat_VC_stdev_binned,sat_VC_count_binned,
           sat_DVC_mean_binned,
           geos_VC_mean_binned,geos_VC_stdev_binned)

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
        
def save_pickles(lat,lon,sat_VC,sat_DVC,geos_VC,time,
                 country,state,index_map,
                 lat_binned,lon_binned,
                 sat_VC_mean_binned,sat_VC_stdev_binned,
                 sat_DVC_mean_binned,
                 geos_VC_mean_binned,geos_VC_stdev_binned,
                 NDVI_binned,
                 country_binned,state_binned,index_map2):
    print "Enter path for pickles to saved to (blank entry for script folder):"
    pickle_save_path = raw_input("-->")
    
    #add in a final slash if the user has missed it
    #if it's blank it defaults to the script directory
    if pickle_save_path != "": 
        if not pickle_save_path.endswith("/"):
            pickle_save_path = pickle_save_path + "/"
            
    print "Enter prefix for pickle file names:"
    pickle_prefix = raw_input("-->")
    
    pickle.dump((lat,lon,sat_VC,sat_DVC,geos_VC,time,
                 country,state,index_map),
                 open(pickle_save_path + pickle_prefix + "_2.p","wb"))
    #"None" is space for NDVI data             
    pickle.dump((lat_binned,lon_binned,
                 sat_VC_mean_binned,sat_VC_stdev_binned,
                 sat_DVC_mean_binned,
                 geos_VC_mean_binned,geos_VC_stdev_binned,
                 NDVI_binned,
                 country_binned,state_binned,
                 index_map2),
                 open(pickle_save_path + pickle_prefix + "_binned.p","wb"))
    current_pickle = pickle_save_path + pickle_prefix
    return(current_pickle)
    
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

def two_var_comparison(var_names,vartuple):
    (xvar,xvar_name) = select_a_var(var_names,vartuple,
                                    "independant (x-axis) variable",numbers_only=True)
    (yvar,yvar_name) = select_a_var(var_names,vartuple,
                                    "dependant (y-axis) variable",numbers_only=True)
    type_of_comparison = ""
    while type_of_comparison != "Z":
        type_of_comparison = basic_menu("Choose type of comparison:",[
                                        ["1","Scatter plot"],
                                        ["2","Simple correlation statistics"],
                                        ["3","Error bar plot"],
                                        ["4","Advanced correlation statistics"]],
                                        quit_option=True)
                                        
        #if 3 or 4, will need to define errors on one or both datasets.
        if type_of_comparison in ["3","4"]:
            [error_x_option,error_y_option] = ["",""]
            while "Z" not in [error_x_option,error_y_option]:
                
                #for x
                error_x_option = basic_menu("How is error in %s determined?" %xvar_name,[
                                            ["1","No error"],
                                            ["2","Fixed value for all data"],
                                            ["3","Fixed fraction for all data"],
                                            ["4","Use existing dataset"]],
                                            quit_option=True)
                if error_x_option == "Z":
                    type_of_comparison = "backtomenu" #go back to type_of_comparison menu
                    continue
                elif error_x_option == "1":
                    xvar_errs = None
                elif error_x_option == "2":                                  
                    xvar_errs = list(np.zeros_like(xvar))
                    fixed_error_value = float(raw_input("Enter fixed error value\n"
                                                        "-->"))
                    for i in range(0,len(xvar_errs)):
                        xvar_errs[i] = fixed_error_value
                    del fixed_error_value
                elif error_x_option == "3":
                    xvar_errs = list(np.zeros_like(xvar))
                    fraction_error_value =float(raw_input("Enter fractional error value\n"
                                                          "-->"))
                    for i in range(0,len(xvar_errs)):
                        xvar_errs[i] = xvar[i]*fraction_error_value
                    del fractional_error_value
                elif error_x_option == "4":
                    (_,xvar_errs) = select_a_var(var_names,vartuple,
                                                 "errors in independant (x-axis) variable",
                                                 numbers_only=True)   
                
                #for y
                error_y_option = basic_menu("How is error in %s determined?" %yvar_name,[
                                            ["1","No error"],
                                            ["2","Fixed value for all data"],
                                            ["3","Fixed fraction for all data"],
                                            ["4","Use existing dataset"]],
                                            quit_option=True)
                if error_y_option == "Z":
                    type_of_comparison = "backtomenu" #go back to type_of_comparison menu
                    continue
                elif error_y_option == "1":
                    yvar_errs = None
                elif error_y_option == "2":                                  
                    yvar_errs = list(np.zeros_like(yvar))
                    fiyed_error_value = float(raw_input("Enter fixed error value\n"
                                                        "-->"))
                    for i in range(0,len(yvar_errs)):
                        yvar_errs[i] = fixed_error_value
                    del fixed_error_value
                elif error_y_option == "3":
                    yvar_errs = list(np.zeros_like(yvar))
                    fraction_error_value =float(raw_input("Enter fractional error value\n"
                                                          "-->"))
                    for i in range(0,len(yvar_errs)):
                        yvar_errs[i] = yvar[i]*fraction_error_value
                    del fractional_error_value
                elif error_y_option == "4":
                    (_,yvar_errs) = select_a_var(var_names,vartuple,
                                                 "errors in independant (y-axis) variable",
                                                 numbers_only=True)
                #if we've got this far we can leave the while loop
                break                 
                                        
        if type_of_comparison == "1":
            error_bar_scatter(xvar,yvar,
                              x_error=None,y_error=None,
                              title="Scatter plot",
                              x_label=xvar_name,y_label=yvar_name)
        elif type_of_comparison == "2":
            print "Linear regression:"
            slope, intercept, r_value, p_value, std_err = scistats.linregress(xvar,yvar)
            print "Best-fit line: y = %+.2gx%+.2g " %(slope,intercept)
            r_squared = r_value * r_value
            print "R-squared    : %.4g " %r_squared
            print "P-value      : %.4g " %p_value
            print "Standard err : %g "   %std_err        
            pearsonr = scistats.pearsonr(xvar,yvar)
            print "Pearson's correlation coefficient:\n %f, 2-tailed p-value: %f" %(pearsonr[0],pearsonr[1])
            polyfit = np.polyfit(xvar, yvar, 1)
        elif type_of_comparison == "3":
            error_bar_scatter(xvar,yvar,
                  x_error=xvar_errs,y_error=yvar_errs,
                  title="Error bar plot",
                  x_label=xvar_name,y_label=yvar_name)
        elif type_of_comparison == "4":
            #WYIBF analysis
            pass
        _ = raw_input("Press enter to continue -->")
        
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
        alpha = float(len(x_var))^0.2
        
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
        plt.scatter(x_var, y_var,
                    alpha=alpha, fmt="o", color='g')
    elif x_error != None and y_error == None: #bars for x, none for y
        plt.errorbar(x_data, y_data, xerr=x_error,
                     alpha=alpha, fmt="o", color='g')
    elif x_error == None and y_error != None: #none for x, bars for y
        plt.errorbar(x_data, y_data, yerr=y_error,
                     alpha=alpha, fmt="o", color='g')
    elif x_error != None and y_error != None: #bars for x, bars for y
        plt.errorbar(x_data, y_data, xerr=x_error, yerr=y_error,
                     alpha=alpha, fmt="o", color='g')   
    
    if do_best_fit or do_best_fit_equation:
        par = np.polyfit(x_data, y_data, 1, full=True)
        slope=par[0][0]
        intercept=par[0][1]
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


class geos_data():
    """A class for objects of data imported from geos chem output"""
    
    def __init__(self):
       self.data = [] #the data
       self.time = [] #time of each block of data
       self.lat = [] #the latitudes
       self.lon = [] #the longitudes
       
       
    def add_data(self, new_data, new_time):
        self.data.append(new_data)
        self.time.append(new_time)
    
    def set_name(self, name): #to name the dataset (name in geos output)
        self.name = name
        
    def set_human_name(self, human_name): #to name the dataset (human name)
        self.human_name = human_name
        
    def set_unit(self, unit): #to set a unit for the dataset
        self.unit = unit
        
    def set_lat_lon(self,lat,lon): #to set latitudes and longitudes
        self.lat = lat
        self.lon = lon

def load_emfile(emfiles_folder,file_tag="trac_avg"):
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
        print "To see a list of all variables within a group, write the group name followed by a *"
        print "You have currently chosen the following variables:"
        print using_variables_list
        print "To clear this list, type 'C'"
        print "Latitude, longitude and time variables will be used automatically"
        print "Once this list is complete, type 'Y' to proceed"
        option = raw_input("-->").upper()
        if   option == "V":
            for variable in variables_list:
                print variable
            raw_input("Press enter to continue-->")
        elif option == "G":
            groups = list(f.groups)
            for group in groups:
                print group
            raw_input("Press enter to continue-->")
        elif option.endswith("*"):
            group_variables_list = [variables_list[i] for i in range(0,num_variables) if variables_list[i].startswith(option[:1])]
            if group_variables_list == []:
                print "This is not a valid group"
            else:
                for variable in group_variables_list:
                    print variable
            raw_input("Press enter to continue-->")
        elif option == "Y":
            done = True
        elif option == "C":
            using_variables_list = []
        else: #if adding a variable
            if option in variables_list:
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
        print "Please enter a short 'human readable' name for the variable %s" %variable
        geos_datasets[i].human_name = raw_input("-->")
        i += 1
    del i
        
    
    #Right, now let's get onto the serious stuff
    for this_file in good_files: #access each file in turn
        print "Accessing %s" %this_file
        f = bpch(this_file)
        time = list(f.variables['time'])
        lat = list(f.variables['latitude'])
        lon = list(f.variables['longitude'])
        num_times = len(time)
        print "There are %i different times recorded in this file" %num_times
        for geos_dataset in geos_datasets: 
            this_name = geos_dataset.name
            print "Reading variable %s" %this_name
            this_data = f.variables[this_name]
            
            #set units and lat/lon
            #we might end up setting these multiple times but
            #that's not a problem big deal
            geos_dataset.unit = this_data.units
            geos_dataset.set_lat_lon(lat,lon)
            
            if len(this_data[0]) != 1: #if there is vertical information
                print "For this 3D variable, only lowest layer information will be read"
            else:
                print "Variable is 2D"
            for t in range(0,num_times):                
                geos_dataset.add_data(this_data[t][0],time[t])
    
    print "Reading geos_chem variables done"
    #OK. let's turn this into a dictionary
    
    print "Dataset created with the following variables:" 
    for geos_dataset in geos_datasets:
        print geos_dataset.human_name
        export_dict[geos_dataset.human_name] = geos_dataset
        
    return(export_dict)
            

    
    
    
