#!/usr/bin/env python
"""CL_satpp.py -- A command-line interface for interogating the sat_pp dataset"""

from datetime import datetime as dt
from datetime import timedelta
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

def clearscreen():
    """Clears the screen. Checks the OS to deliver the correct command"""
    os.system('cls' if os.name == 'nt' else 'clear')


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

def top_level_menu(current_pickle,current_emfiles):
    choice = ""
    while not(choice.upper() in ["1","2","G","P","R","E","X","Z"]):
        os.system('clear')
        print "TOP LEVEL MENU"
        print "Current pickles loaded: %s" %current_pickle
        print "Current emissions files are: %sXXX%s" %(current_emfiles[0],current_emfiles[1])
        print "[1] Use binned data"
        print "[2] Use individual observations data"
        print "[G] Geographically select data to use"
        print "[P] Change pickle"
        print "[R] Reload pickle"
        print "[E] Change current emissions files"   
        print "[X] Inspect and change global options"
        print "[Z] Quit"
        choice = raw_input("-->")
    return(choice.upper())

def use_binned_data_menu():
    choice = ""
    while not(choice.upper() in ["1","2","Z"]):
        os.system('clear')
        print "USING BINNED DATA"
        print "[1] Plot gridded data on map"
        print "[2] Compute basic statistics"
        print "[Z] Return to top level menu" 
        choice = raw_input("-->")
    return(choice.upper())

def use_indiv_data_menu():
    choice = ""
    while choice.upper() not in ["1","2","3","Z"]:
        os.system('clear')
        print "USING INDIVIDUAL DATA"
        print "[1] Simple dots-on-map"
        print "[2] Compute basic statistics"
        print "[3] Compute timeseries statistics"
        print "[Z] Return to previous menu"
        choice = raw_input("-->")
    return(choice.upper())        
    
def binned_map_menu():
    choice = ""
    while not(choice.upper() in ["1","2","3","4","5","6","Z"]):
        os.system('clear')
        print "Plotting map of binned data"
        print "Loaded pickle and global options will be used"
        print "Print which variable?"
        print "[1] Mean satellite observations"  
        print "[2] Calculated uncertainty in satellite obsservations" 
        print "[3] Mean model values" 
        print "[4] Standard deviation of model values" 
        print "[5] Mean NDVI value" 
        print "[6] Mean deviation of model from observations (sigmas)" 
        print "[Z] Return to previous menu"
        choice = raw_input("-->")
    return(choice.upper())

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
        

    
def dots_on_map_menu():
    choice = ""
    while not(choice.upper() in ["1","2","3","Z"]):
        os.system('clear')
        print "Plotting map of binned data"
        print "Loaded pickle and global options will be used"
        print "Print which variable?"
        print "[1] Mean satellite observations"  
        print "[2] Calculated uncertainty in satellite obsservations" 
        print "[3] Mean model values" 
        print "[Z] Return to previous menu"
        choice = raw_input("-->")
    return(choice.upper())
    
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
            print "[P] or [ ] Plot figure on screen"
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
            new_e_year = change_var(  enddate.year,"Start time year")
            new_e_month= change_var(  enddate.month,"Start time month")
            new_e_day  = change_var(  enddate.day,"Start time day")
            enddate    = dt(new_e_year,new_e_month,new_e_day,0,0,0)
        elif choice == "Z" or choice == "z":
            leave = True
        
    return(north,south,east,west,xdim,ydim,startdate,enddate)

def change_pickle(current_pickle):
    os.system('clear')
    print "Currently, pickles with prefix %s are loaded."
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
    print "Currently, pickles with prefix %s are loaded."
    print "To change this, enter a new prefix."
    print "Otherwise, enter Z to keep current option and return to top menu"
    choice = raw_input("-->")
    if choice == "Z" or choice == "z":
        changed = False
    else:
        current_pickle = choice
        changed = True
    return(changed,current_emfiles)
    
def change_var(var,var_name="[unspecified]"):
    former_type = type(var)
    var = former_type(raw_input("Enter new value for "+var_name+"-->"))
    return(var)    

def load_new_pickles_all(current_pickle,verbose=True):
    NDVI_data = load_new_pickles_NDVI(current_pickle,verbose=verbose)
    indv_data = load_new_pickles_indv(current_pickle,verbose=verbose)
    binn_data = load_new_pickles_binn(current_pickle,verbose=verbose)
    
    return(NDVI_data,indv_data,binn_data)

    
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
    (lat_binned,lon_binned,sat_mean_binned,geos_mean_binned,sat_uncer_binned,geos_stdev_binned,dev_mean_binned,NDVI_mean_binned,countries_binned,states_binned,index_map) = pickle.load( open(current_pickle + "_binned.p","rb") )
    #(lat_binned,lon_binned,sat_mean_binned,geos_mean_binned,sat_uncer_binned,geos_stdev_binned,dev_mean_binned,NDVI_mean_binned,countries_binned,states_binned) = stripallnans(lat_binned,lon_binned,sat_mean_binned,geos_mean_binned,sat_uncer_binned,geos_stdev_binned,dev_mean_binned,NDVI_mean_binned,countries_binned,states_binned)
    
    if verbose:
        print "Pickled binned data now loaded as %i-long lists" %len(lat_binned)
        print("Binned latitudes     : lat_binned       ")
        print("Binned_longitudes    : lon_binned       ")
        print("Binned obs., mean    : sat_mean_binned  ") 
        print("Binned obs., uncer   : sat_uncer_binned ")
        print("Binned modelled mean : geos_mean_binned ")
        print("Binned modelled stdev: geos_stdev_binned")
        print("Binned NDVI mean     : NDVI_mean_binned ")
        print("Binned deviation mean: dev_mean_binned  ")
        print("Country strings      : countries_binned ")
        print("State strings        : states_binned    ")
        print("-----------------------------")
    binn_data = (lat_binned,lon_binned,sat_mean_binned,geos_mean_binned,sat_uncer_binned,geos_stdev_binned,dev_mean_binned,NDVI_mean_binned,countries_binned,states_binned)
    
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
    m = Basemap(projection='merc',llcrnrlat=south+1,urcrnrlat=north-1,\
            llcrnrlon=west+1,urcrnrlon=east-1,lat_ts=(north+south)/2.,resolution='i')
    
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
    [northbound,southbound,westbound,eastbound] = geo_selection
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
        clockhigh = clocklow + timedelta(days=step_days)
        this_dataset = time_select(clocklow,clockhigh,time,dataset)
        this_stat = calc_statistic(this_dataset,statistic)
        central_time = clocklow + timedelta(days=step_days*0.5)
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
