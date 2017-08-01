#!/usr/bin/env python

from CL_satpp import *
from CL_corrections import *
from CL_oversampler import *
import os
from datetime import datetime as dt
import copy
import pickle
import brewer2mpl
import pick_colmap

#You will need to make this file executable
#chmod +x CL.py
#Excute this from the command line:
#./CL.py

#Here begins the program.

yesno = [["Y","Yes"],["N","No"]] #helps with simple questions

#"default" variable nameds for data 
indiv_varnames = ["latitude",
                  "longitude",
                  "Observed vertical column",
                  "Error in observed vertical column",
                  "Modelled vertical column",
                  "Time",
                  "Country",
                  "State"]
binned_varnames= ["latitude",
                  "longitude",
                  "Observed vertical column, mean in bin",
                  "Observed vertical column, standard deviation in bin",
                  "Error in observed vertical column, mean in bin",
                  "Modelled vertical column, mean in bin",
                  "Modelled vertical column, standard deviation in bin",
                  "Country",
                  "State"]

#We'll look for settings.p in the current working directory
#The idea of this file is to save you having to set the same variables
#repeatedly from one session to another.

#load settings (if settings.p file exists)
if os.path.isfile("settings.p"):
    (country_statefile,
    map_box,
    xdim,ydim,
    startdate,enddate,
    type_geo_select,
    geo_selection,
    current_pickle,current_geosfolder
    ) = cPickle.load(open("settings.p","rb"))
else: #otherise, some default options
    country_statefile="country_state.csv"    
    map_box = box(38., #north
                  2.,  #south
                  100.,#east
                  65.) #west
    xdim = 0.3125
    ydim = 0.25
    startdate = dt(2014,01,01, 0, 0, 0)
    enddate   = dt(2014,12,31,23,59,59)
    type_geo_select = "latlon"
    geo_selection = box(38., #north
                        2.,  #south
                        100.,#east
                        65.) #west                              
    current_pickle = "NONE SELECTED"
    current_geosfolder = "/group_workspaces/jasmin/geoschem/local_users/lsurl/runs/geosfp_025x03125_tropchem_in_2014/" 

#default colormap
colmap = brewer2mpl.get_map("RdYlBu","Diverging", 6, reverse=True)

#==Loading data==
#The main data sets for this utility are saved as pickles:
#a) All the individual data (data for each sateliite pixel)
#Such files are saved with the suffix "_indiv.p", though we can also 
#load data in the correct format with any filename (such as the
# output from the pre-processor)
# This data is loaded to the variable "ida"
#b) The binned data (binned onto a grid)
#This file has the suffix "_binned.p"
#The data is loaded into the variable "bda"

#This script asks the user if they want to load the last set of data used
#as recorded in settings.p

if current_pickle == "NONE SELECTED": #if settings.p does not exist
    initial_choice = "2"
else:
    initial_choice = basic_menu(
                    "INITIAL MENU",
                    [
                     ["1","Load pickles: %s" %current_pickle],
                     ["2","Proceed without loading"]
                    ],
                    quit_option=False)
 

if initial_choice == "1":
    print "Loading default options"       
    
    #ida is the main data holder for all individual data
    #it is of type ind_data_all and holds ind type objects
    try:
        ida_p = load_new_pickles_da(current_pickle,"_2.p") #legacy
    except IOError: #if "_2.p" file does not exist
        try:
            ida_p = load_new_pickles_da(current_pickle,"_indiv.p")
        except IOError:
            print "No ida file? Will try to proceed anyway"
    ida = copy.deepcopy(ida_p) #ida_p is kept static while ida gets modified.
    
    #bda is the main data holder for all binned data
    #it is of type bin_data_all and holds bnd type objects 
    try:   
        bda_p = load_new_pickles_da(current_pickle,"_binned.p") #bda_p is kept static while bda gets modified.
        bda = copy.deepcopy(bda_p)
    except IOError:
        print "%s_binned.p does not exist. If you want binned data, you'll need to do the binning" %current_pickle
        _ = raw_input("Press enter to continue")
    

#Main menu loop.
top_level_menu_choice = ""
while top_level_menu_choice != "Z": #loop unless ordered to quit

    #every time we get here, quickly re-save settings
    cPickle.dump(
    (country_statefile,
    map_box,
    xdim,ydim,
    startdate,enddate,
    type_geo_select,
    geo_selection,
    current_pickle,current_geosfolder
    ),
    open("settings.p","wb"))

    #The menu
    
    top_level_menu_title = "TOP LEVEL MENU\n" \
               "Current pickles loaded: %s\n" \
               "Current GEOS chem directory: %s" \
               %(current_pickle,current_geosfolder)
    top_level_menu_text = [
            ["1","Use binned data"],
            ["2","Use individual observations data"],
            ["G","Geographically select data to use"],
            ["F","Filter indvidual data based on values"],
            ["M","Filter down to a defined month"],
            ["S","Save current dataset as new pickle"],
            ["P","Change pickle"],
            ["O","Create binned from individual data"],
            ["C","Assign country/state to data"],
            ["CF","Assign country/state to data FAST"],                        
            ["B","Bin additional datasets"],
            ["BI","Add binned data to individual data"],
            ["N","Associate CSV-form data"],
            ["Ni","Associate NDVI data (individual)"],
            ["Li","Associate LAI data (individual)"],            
            ["Fi","Associate Fire data (individual)"],
            ["FM","Associate Fire mask"],
            ["R","Reload pickle"],
            ["E","Change GEOS Chem directory"],
            ["A","Load GEOS Chem trac_avg data"],
            ["T","Load GOES Chem ND51 data"],
            ["K","Plot GOES Chem data"],
            ["COL","Change colour bar choice"],   
            ["X","Inspect and change global options"],
            ["Z","Quit"]
        ]
    top_level_menu_choice = basic_menu(top_level_menu_title,
                                       top_level_menu_text,
                                       quit_option=False)
    
    if top_level_menu_choice == "P": #change pickle
        found_something = False
        while found_something == False:
            (changed,current_pickle) = change_pickle(current_pickle)
            if changed == False:
                break
            #reloading is time-consuming so only do if change made
            #delete current
            if 'ida' in locals():
                del ida
            if 'bda' in locals():
                del bda
                             
            if os.path.isfile(current_pickle+"_2.p"): #try with _2.p suffix
                print "Found individual data: " + current_pickle+"_2.p"
                ida_p = load_new_pickles_da(current_pickle,"_2.p")
                ida = copy.deepcopy(ida_p) #ida_p is kept static while ida gets modified.
                found_something = True
            elif os.path.isfile(current_pickle+"_indiv.p"): #try with _indiv.p suffix
                print "Found individual data: " + current_pickle+"_indiv.p"
                ida_p = load_new_pickles_da(current_pickle,"_indiv.p")
                ida = copy.deepcopy(ida_p) #ida_p is kept static while ida gets modified.
                found_something = True
            elif os.path.isfile(current_pickle+".p"): #try with no suffix
                print "Found individual data: " + current_pickle+".p"
                ida_p = load_new_pickles_da(current_pickle,".p")
                ida = copy.deepcopy(ida_p) #ida_p is kept static while ida gets modified.
                found_something = True
            else:
                print "No individual data found!"
            
            if os.path.isfile(current_pickle+"_binned.p"): #look for binned
                print "Found binned data: "+current_pickle+"_binned.p"
                bda_p = load_new_pickles_da(current_pickle,"_binned.p")
                bda = bda_p #bda_p is kept static while bda gets modified
                found_something = True       
            else:
                print "No binned data detected!"
            if found_something == False:
                print "ALERT! No data was found!"
                current_pickle = "NONE"
                    
        
    elif top_level_menu_choice == "E": #change emissions
        (changed,current_geosfolder) = change_emfiles(current_geosfolder)
        #currently no loading of the emisisons files are done here
    
    elif top_level_menu_choice == "A": #load geos chem data
        geos_dict = load_geosfile(current_geosfolder,file_tag="trac_avg")
        option = basic_menu("Done. Do you wish to associate"
                            " these data with the observations?",
                            [["Y","yes"],["N","no"]],quit_option=False)
        if option == "Y":
            associate_ind_geos(ida,geos_dict)
    
    elif top_level_menu_choice == "T": #load psuedo-satellite data
        geos_dict = load_geosfile(current_geosfolder,file_tag="ts_",columns=True)
        option = basic_menu("Done. Do you wish to associate"
                            " these data with the observations?",
                            [["Y","yes"],["N","no"]],quit_option=False)
        if option == "Y":
            ida = associate_ind_geos(ida,geos_dict)
            
    elif top_level_menu_choice == "B": #Bin additional datasets
        bda = bin_extra(ida,bda)
        
    elif top_level_menu_choice == "BI": #Bin additional datasets
        ida = binned_to_indv(ida,bda)  
    
    elif top_level_menu_choice == "N": #associate NDVI
        #NDVI_dir = raw_input("Enter directory containing NDVI data-->")
        
        option=basic_menu("Which type of gridded data to associate?",
                          [["N","NDVI"],["L","LAI"],["F","Fire count"],["LC","Land classification"],
                           ["P","Population density"]                                                      ],
                          quit_option=True)
        if option=="N":
            CSV_dir = "/group_workspaces/cems2/nceo_generic/nceo_ed/NDVI"
            CSV_name = "NDVI"
            CSV_desc = "Normalised Diffusive Vegitation Index (NDVI)"
        elif option =="L":
            CSV_dir = "/group_workspaces/cems2/nceo_generic/nceo_ed/LAI"
            CSV_name = "LAI"
            CSV_desc = "Leaf Area Index (LAI)"
        elif option =="F":
            CSV_dir = "/group_workspaces/cems2/nceo_generic/nceo_ed/fire"
            CSV_name = "Fire count"
            CSV_desc= "MODIS fire count (fire pixels/1000km2/day)"
        elif option =="LC":
            CSV_dir = "/group_workspaces/cems2/nceo_generic/nceo_ed/LC"
            CSV_name = "Land classification"
            CSV_desc= "Land cover classification"
        elif option=="P":
            CSV_dir = "/group_workspaces/cems2/nceo_generic/nceo_ed/population"
            CSV_name = "Population density"
            CSV_desc= "Population per km2"
                         
        elif option=="Z":
            continue #quit menu    
        bda = associate_CSV(CSV_dir,map_box,bda,startdate,enddate,CSV_name,CSV_desc)
    
    elif top_level_menu_choice == "NI": #associate NDVI ida
        NDVI_dir = "/group_workspaces/cems2/nceo_generic/nceo_ed/NDVI"
        NDVI_3D = NDVI_months_v2(startdate,enddate,NDVI_dir,0.1,0.1,map_box)
        associate_NDVI_v2(ida,NDVI_3D,0.1,0.1,map_box,startdate)
        
    elif top_level_menu_choice == "LI": #associate LAI ida
        LAI_dir = "/group_workspaces/cems2/nceo_generic/nceo_ed/LAI"
        LAI_3D = NDVI_months_v2(startdate,enddate,LAI_dir,0.1,0.1,map_box,datatype="LAI")
        associate_NDVI_v2(ida,LAI_3D,0.1,0.1,map_box,startdate,datatype="LAI")        
    elif top_level_menu_choice == "FI": #associate fires ida
        FIRE_dir = "/group_workspaces/cems2/nceo_generic/nceo_ed/fire"
        FIRE_3D = NDVI_months_v2(startdate,enddate,FIRE_dir,0.1,0.1,map_box,datatype="Fire count")
        associate_NDVI_v2(ida,FIRE_3D,0.1,0.1,map_box,startdate,datatype="Fire count")
    
    elif top_level_menu_choice == "FM": #associate fire mask 
        FIRE_file = "/group_workspaces/cems2/nceo_generic/nceo_ed/fire_daily/fire_archive_M6_7826.csv"
        daily_fire_filter(ida,FIRE_file,startdate,enddate,map_box,xdim,ydim)
                
    elif top_level_menu_choice == "K": #plot geos chem data
        plot_geos_chem(geos_dict)
    
    elif top_level_menu_choice == "X": #change options
        if 'winter_flag' not in locals():
            winter_flag = False
        datesbefore = (startdate,enddate)
        (map_box,
         xdim,ydim,
         startdate,enddate,winter_flag) = global_options(
                             map_box,
                             xdim,ydim,
                             startdate,enddate,winter_flag)
        ida = time_select(startdate,enddate,ida,winter_flag=winter_flag)
            
    elif top_level_menu_choice == "R": #reload pickle
        #re-unpack data sets to clear any previous geographic selection
        print "Reloading..."
        ida = copy.deepcopy(ida_p) #ida_p is kept static while ida gets modified.
        try:
            bda = copy.deepcopy(bda_p) #bda_p is kept static while bda gets modified.
        except NameError: #bda doesn't exist
            pass
        #dev_mean_binned = list(np.subtract(geos_VC_mean_binned,sat_VC_mean_binned))
        
        
    elif top_level_menu_choice == "G": #geographic selection option
        before_geo = (type_geo_select,geo_selection)
        (type_geo_select,geo_selection) = geo_select_menu(type_geo_select,geo_selection)
        
        print "Updating geographic selection..."
        if type_geo_select in ["country","state"]:
            try: #will fail if ida doesn't have country/state data
                ida = da_select_by_match(type_geo_select,geo_selection,ida)
            except:
                print "Cannot process individual data, country/state not assigned"
            if 'bda' in locals(): #only if bda exits
                try: #will fail if ida doesn't have country/state data
                    bda = da_select_by_match(type_geo_select,geo_selection,bda)
                except:
                    print "Cannot process individual data, country/state not assigned"    

        elif type_geo_select in ["rectangle","circle"]:
            da_select_by_shape(geo_selection,ida)
            if 'bda' in locals(): #only if bda exits
                da_select_by_shape(geo_selection,bda)                   
 
    
    elif top_level_menu_choice == "F": #values filtering
        ida = criteria_filtering(ida)
    
    elif top_level_menu_choice == "M": #month filtering
        ida = only_months(ida)
    
    elif top_level_menu_choice == "O": #create binned data
        "Binning data to %fx%f grid cells" %(ydim,xdim)   
        bda = create_binned_set(ida,map_box,xdim,ydim)
    
    elif top_level_menu_choice in ["C","CF"]: #Assign country/state to newly-loaded pickles
                
        #Country and state matching
        print "Matching countries and states"
        print "Countries and state assignments will be loaded from %s" %country_statefile
        print "Press enter to use this file or enter new location"
        new_option = raw_input("-->")
        if new_option != "":
            country_statefile == new_option
        del new_option
        print "Reading %s" %country_statefile
        (csv_lat,csv_lon,csv_country,csv_state) = \
            load_regiondata(country_statefile,states=True)
        do_ind_countrystate = "x"
        while do_ind_countrystate not in ["Y","y","N","n"]:
            do_ind_countrystate = raw_input("Assign to individual data points? Y/N -->").upper()
        if do_ind_countrystate in ["Y","y"]:
            print "Assigning country+state to indiv. data points"
            if top_level_menu_choice == "C":
                (country,state) =\
                    region_matcher(csv_lat,csv_lon,csv_country,csv_state,
                                   ida.lat,ida.lon)
            elif top_level_menu_choice == "CF":
                (country,state) =\
                    region_matcher_fast(csv_lat,csv_lon,csv_country,csv_state,
                                   ida.lat,ida.lon)
                               
            #add these to ida properly
            ida.data['country'] = d(country,'country','Country observation falls within')
            del country
            ida.data['state'  ] = d(state  ,'state'  ,'State observation falls within')
            del state
            
        if 'bda' in locals():                   
            print "Assigning country+state to binned data points"
            if top_level_menu_choice == "C":
                (country_binned,state_binned) =\
                    region_matcher(csv_lat,csv_lon,csv_country,csv_state,
                                   bda.lat,bda.lon)
            elif top_level_menu_choice == "CF":
                (country_binned,state_binned) =\
                    region_matcher_fast(csv_lat,csv_lon,csv_country,csv_state,
                                   bda.lat,bda.lon)
            bda.data["country"] = d(country_binned,"country","Country")
            del country_binned
            bda.data["state"  ] = d(state_binned  ,"state"  ,"State"  )
            del state_binned

                  
    elif top_level_menu_choice == "S":
        print "Enter path for pickles to saved to (blank entry for script folder):"
        save_path = raw_input("-->")

        #add in a final slash if the user has missed it
        #if it's blank it defaults to the script directory
        if save_path != "": 
            if not save_path.endswith("/"):
                save_path =save_path + "/"
        
        print "Enter name for pickle file. Suffixes will be appended automatically:"
        prefix = raw_input("-->")
                                           
        save_pickle(ida,save_path=save_path,prefix=prefix,suffix="_2.p")
        ida_p = copy.deepcopy(ida) #update the static store
        try:
            save_pickle(bda,save_path=save_path,prefix=prefix,suffix="_binned.p")
            bda_p = copy.deepcopy(bda) #update the static store
        except NameError: #bda doesn't exist
            pass
        current_pickle = save_path+prefix
    
    elif top_level_menu_choice == "COL": #change color bar
        colmap = pick_colmap.pick_colmap(preview=False)   
        
    elif top_level_menu_choice == "1": #Binned data
        #ubd = use binned data
        ubd_menu_choice = ""
        while ubd_menu_choice != "Z":

            ubd_menu_title = "USING BINNED DATA"
            ubd_menu_text  = [
                ["1","Plot gridded data on map"],
                ["2","Compute basic statistics"],
                ["3","Compare datasets"],
                ["4","Write out as CSV"]
                ]
            
            ubd_menu_choice = basic_menu(ubd_menu_title,
                                         ubd_menu_text,
                                         quit_option="True")
                   
            if ubd_menu_choice == "1": #Binned map menu (bmm)
                while True:                    
                    os.system('clear')
                    
                    bmm_title = "Plotting map of binned data\n"\
                                "Loaded pickle and global options will be used\n"\
                                "Plot which variable?"
                    [bmm_menu_text,bmm_answers_dict] = bda.make_menu()                        
                    bmm_menu_choice = basic_menu(bmm_title,
                                            bmm_menu_text,
                                            quit_option=True)
                    if bmm_menu_choice != "Z": #unless we're quitting
                        data_key = bmm_answers_dict[bmm_menu_choice]                     
                        
                        (map_preplot_menu_choice,title,vmin,vmax,unit,save_filename) =\
                            map_preplot_menu(bda.data[data_key].description,
                            vmin=np.nanmin(bda.data[data_key].val),
                            vmax=np.nanmax(bda.data[data_key].val),
                            unit=bda.data[data_key].unit)
                    else: #go up a menu
                        break
                    
                    if map_preplot_menu_choice == "Z":
                        break
                                                                                                                     
                    #At this point binned_map_preplot_menu_choice will be P (plot) or S (save)                    
                    do_save = map_preplot_menu_choice.upper() == "S" #saving the figure?
                    
                    plot_grid_from_list(bda.lat,bda.lon,bda.data[data_key].val,
                                        bda.meta["Binning dimensions"][0],bda.meta["Binning dimensions"][1],
                                        map_box,title=title,vmin=vmin,vmax=vmax,lab=unit,
                                        save=do_save,save_filename=save_filename,colmap=colmap)
                    break
            
            #bda-ification done to this point
                    
            elif ubd_menu_choice == "2": #binned statisics
                stats_from_da(bda)
                
            elif ubd_menu_choice == "3": #compare_datasets 
                two_var_comparison(bda)
                
            elif ubd_menu_choice == "4": #write as CSV
                write_as_csv(bda)
                
            elif ubd_menu_choice == "5": #simple maths
                bda = simple_maths(bda)
                     
    elif top_level_menu_choice == "2": #Indiv data
        #uid = use individual data
        uid_menu_choice = ""
        while uid_menu_choice != "Z":
            uid_menu_title = "USING INDIVIDUAL DATA"
            uid_menu_text  = [
                ["1","Simple dots-on-map"],
                ["2","Compute basic statistics"],
                ["3","Timeseries"],
                ["4","Compare two datasets"],
                ["5","Oversample"],
                ["6","Overlay"],
                ["7","SEASONAL maps"],
                ["8","SEASONAL inversions"],
                ["9","Simple mathetatical operations"]
                ]
            uid_menu_choice = basic_menu(uid_menu_title,
                                         uid_menu_text,
                                         quit_option=True)
            if uid_menu_choice == "1": #dots on map
                #dom = dots on map
                dom_menu_choice = ""
                while dom_menu_choice != "Z":
                    dom_menu_title = "Plotting map of data\n"\
                                     "Loaded pickle and global options will be used\n"\
                                     "Print which variable?"
                    [dom_menu_text,dom_answers_dict] = ida.make_menu()

                    dom_menu_choice = basic_menu(dom_menu_title,
                                                 dom_menu_text,
                                                 quit_option=True)
                    if dom_menu_choice != "Z": #unless we're quitting
                        data_key = dom_answers_dict[dom_menu_choice]
                        (map_preplot_menu_choice,title,vmin,vmax,unit,save_filename) = \
                            map_preplot_menu(ida.data[data_key].description,unit=ida.data[data_key].unit)                  
                        if   map_preplot_menu_choice.upper() == "S": #saving the figure
                            plot_dots_on_map(ida.lat,ida.lon,ida.data[data_key].val,map_box,vmin=vmin,vmax=vmax,title=title,lab=unit,save=True,save_filename=save_filename)
                        elif map_preplot_menu_choice.upper() == "P": #plotting the figure
                            plot_dots_on_map(ida.lat,ida.lon,ida.data[data_key].val,map_box,vmin=vmin,vmax=vmax,title=title,lab=unit,save=False,colmap=colmap)
                        else: #no plot
                            pass     
                                       
            elif uid_menu_choice == "2": #statisics
                stats_from_da(ida)                
            
            elif uid_menu_choice == "3": #timeseries statistics
            
                while True: #allow for break to exit
                    ts1_menu_title = "Which dataset do you wish to get a timeseries for?"
                    [ts1_menu_text,ts1_answers_dict] = ida.make_menu()
                    ts1_menu_choice = basic_menu(ts1_menu_title,
                                                 ts1_menu_text,
                                                 quit_option=True)
                    if ts1_menu_choice == "Z":
                        break
                    
                    data_key = ts1_answers_dict[ts1_menu_choice]
                      
                    while True: #allow for break to exit   
                        ts2_menu_title = "Which statistical operator do you want to compute for each time period?"
                        ts2_menu_text=[["1","mean"],
                                       ["2","standard deviation"],
                                       ["3","count"]]
                        ts2_menu_choice= basic_menu(ts2_menu_title,
                                                    ts2_menu_text,
                                                    quit_option=True)
                        if ts2_menu_choice == "Z":
                            break
                        
                        #a cumbersome but robust way of setting stat_choice_text
                        for pair in ts2_menu_text:
                            if pair[0] == ts2_menu_choice:
                                stat_choice_text = pair[1]
                                stat_choice = ts2_menu_choice
                        
                        step = 0 #assign 0 before assignment
                        month_flag = False #False unless set to True       
                        while True: #allow for break to exit
                            print "How long should the time series be for each individual point?"
                            print "Enter a number to choose that many days"
                            print "Alternatively, enter M then a number (i.e. M1) to set timestep to that many calendar months"
                            print "[Z] Return to previous menu"
                            ts3_menu_choice = raw_input("-->").upper()
                            if ts3_menu_choice == "Z":
                                break
                            elif ts3_menu_choice.startswith("M"):
                                try:
                                    step = int( ts3_menu_choice[1:] )
                                    month_flag = True
                                except ValueError:
                                    print "Invalid selection."
                                    continue
                            else:
                                try:
                                    step = int(ts3_menu_choice)
                                except ValueError:
                                    print "Invalid selection."
                                    continue
                            
                            if step > 0: #if we've assigned a valid option
                                clearscreen()
                                # visual readout
                                print "Dataset:   " + ida.data[data_key].description
                                print "Statistic: " + stat_choice_text
                                if month_flag:
                                    print "Timesteps: "+ str(step) + " calendar months"
                                else:
                                    print "Timesteps: " + str(step) + " days"
                                                                        
                                this_ts = \
                                    time_cycle(ida,data_key,stat_choice,
                                    step,month_flag=month_flag,plot=True)  
                                
                            break
                        break
                   
            elif uid_menu_choice == "4": #compare datasets
                two_var_comparison(ida)
            
            elif uid_menu_choice == "5": #oversampler
                #print len(ida.lat)
                #raw_input("halt-->")
                
                os_data = oversample(ida)
                save_location = raw_input("Dump this data as a pickle? Enter full path or press enter to skip this-->")
                if save_location != "":
                    cPickle.dump(os_data,open(save_location,"wb"))
                #print len(ida.lat)
                #raw_input("halt-->")
           
            elif uid_menu_choice == "6": #overlayer
                ov_data = overlay(ida)
                save_location = "/home/users/lsurl/CL/OL/allindia_2011-2015_DJF_OL.p"
                if save_location != "":
                    cPickle.dump(ov_data,open(save_location,"wb"))
            
            elif uid_menu_choice == "7": #SEASONAL MAPS
                save_location = "/home/users/lsurl/CL/imgs"
                save_location_o = raw_input("Enter path to save images to, or press enter for %s --> "%save_location)
                if save_location_o != "":
                    save_location = save_location_o
                seasonal(ida,map_box,xdim,ydim,colmap,save_location,
                         do_options=["obs_map","mod_map"])
            elif uid_menu_choice == "8": #SEASONAL INVERSIONS
                save_location = "/home/users/lsurl/CL/imgs" #not used 
                seasonal(ida,map_box,xdim,ydim,colmap,save_location,
                         do_options=["inversion"])
            
            elif uid_menu_choice == "9": #Simple maths
                ida = simple_maths(ida)
        
    elif top_level_menu_choice == "Z": #Quit
        reallyquit = ""
        while reallyquit not in ["Y","N"]:
            reallyquit = raw_input("Are you sure you want to quit? Y/N").upper()
            if reallyquit == "Y":
                sys.exit()
            else:
                top_level_menu_choice = ""
                       
        
       
