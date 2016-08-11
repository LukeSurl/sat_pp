#!/usr/bin/env python

from CL_satpp import *
from CL_corrections import *
from CL_oversampler import *
import os
from datetime import datetime as dt
import copy

#"shorthand" variables
yesno = [["Y","Yes"],["N","No"]] 

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


#initialisation
initial_choice = basic_menu(
                "INITIAL MENU",
                [
                 ["1","Load default data"],
                 ["2","Proceed without loading"]
                ],
                quit_option=False)

#"Default settings"
filtered_startdate = dt(1970,1,1)
filtered_enddate   = dt(2099,12,31)
filtered_folder    = \
    "/group_workspaces/jasmin/geoschem/local_users/lsurl/sat_pp/filtered/"
corrections_folder= \
    "/group_workspaces/jasmin/geoschem/local_users/lsurl/sat_pp/corrections/"
corr_type="basic"
country_statefile="/home/users/lsurl/CL/country_state.csv"    
#south = 2.
#north = 38.
#west = 65.
#east = 100.
#compass = [north,south,east,west]

map_box = box(38., #north
              2.,  #south
              100.,#east
              65.) #west

xdim = 0.3125
ydim = 0.25
startdate = dt(2014,03,01, 0, 0, 0)
enddate   = dt(2014,12,31,23,59,59)
type_geo_select = "latlon"
geo_selection = box(38., #north
                    2.,  #south
                    100.,#east
                    65.) #west
                           
current_pickle = "NONE SELECTED"
current_geosfolder = "/group_workspaces/jasmin/geoschem/local_users/lsurl/runs/geosfp_025x03125_tropchem_in_2014/"   

if initial_choice == "1":
    print "Loading default options"       
    current_pickle =   '/home/users/lsurl/CL/pickles/2014-11-B'
    current_geosfolder = '/group_workspaces/jasmin/geoschem/local_users/lsurl/runs/geosfp_025x03125_tropchem_in_2014'
    
    #ida is the main data holder for all individual data
    #it is of type ind_data_all and holds ind type objects
    ida_p = load_new_pickles_da(current_pickle,"_2.p")
    ida = copy.deepcopy(ida_p) #ida_p is kept static while ida gets modified.
    
    #bda is the main data holder for all binned data
    #it is of type bin_data_all and holds bnd type objects    
    bda_p = load_new_pickles_da(current_pickle,"_binned.p") #bda_p is kept static while bda gets modified.
    bda = copy.deepcopy(bda_p) 
   
    #(NDVI_lat,NDVI_lon,NDVI,NDVI_year,NDVI_month) = NDVI_data_all
    #(lat,lon,
    #    sat_VC,sat_DVC,
    #    geos_VC,time,
    #    country,state,
    #    index_map) = indv_data_all
    #(lat_binned,lon_binned,
    #    sat_VC_mean_binned,sat_VC_stdev_binned,
    #    sat_DVC_mean_binned,
    #    geos_VC_mean_binned,geos_VC_stdev_binned,
    #    country_binned,state_binned) = binn_data_all
    #dev_mean_binned = list(np.subtract(geos_VC_mean_binned,sat_VC_mean_binned))
    

#main loop
top_level_menu_choice = ""
while top_level_menu_choice != "Z": #loop unless ordered to quit

    top_level_menu_title = "TOP LEVEL MENU\n" \
               "Current pickles loaded: %s\n" \
               "Current GEOS chem directory: %s" \
               %(current_pickle,current_geosfolder)
    top_level_menu_text = [
            ["1","Use binned data"],
            ["2","Use individual observations data"],
            ["G","Geographically select data to use"],
            ["C","Generate new dataset from filtered data and pacific corrections"],
            ["S","Save current dataset as new pickle"],
            ["P","Change pickle"],
            ["B","Bin additional datasets"],
            ["R","Reload pickle"],
            ["E","Change GEOS Chem directory"],
            ["A","Load GEOS Chem trac_avg data"],
            ["T","Load GOES Chem ND51 data"],
            ["K","Plot GOES Chem data"],   
            ["X","Inspect and change global options"]
        ]
    top_level_menu_choice = basic_menu(top_level_menu_title,
                                       top_level_menu_text)
    
    if top_level_menu_choice == "P": #change pickle
        (changed,current_pickle) = change_pickle(current_pickle)
        if changed:
            #reloading is time-consuming so only do if change made
            ida_p = load_new_pickles_da(current_pickle,"_2.p")
            ida = ida_p #ida_p is kept static while ida gets modified.
            bda_p = load_new_pickles_da(current_pickle,"_binned.p")
            bda = bda_p #bda_p is kept static while bda gets modified
        
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
        geos_dict = load_geosfile(current_geosfolder,file_tag="ts_")
        option = basic_menu("Done. Do you wish to associate"
                            " these data with the observations?",
                            [["Y","yes"],["N","no"]],quit_option=False)
        if option == "Y":
            ida = associate_ind_geos(ida,geos_dict)
            
    elif top_level_menu_choice == "B": #Bin additional datasets
        bda = bin_extra(ida,bda)
                
    elif top_level_menu_choice == "K": #plot geos chem data
        plot_geos_chem(geos_dict)
    
    elif top_level_menu_choice == "X": #change options
        datesbefore = (startdate,enddate)
        (map_box,
         xdim,ydim,
         startdate,enddate) = global_options(
                             map_box,
                             xdim,ydim,
                             startdate,enddate)
        #if time selection has changed, reselect data.
        if datesbefore != (startdate,enddate): 
            print "Start and/or end time has changed. Selecting data..."                
            ida = time_select(startdate,enddate,ida)
            
    elif top_level_menu_choice == "R": #reload pickle
        #re-unpack data sets to clear any previous geographic selection
        print "Reloading..."
        ida = copy.deepcopy(ida_p) #ida_p is kept static while ida gets modified.
        bda = copy.deepcopy(bda_p) #bda_p is kept static while bda gets modified.
        #dev_mean_binned = list(np.subtract(geos_VC_mean_binned,sat_VC_mean_binned))
        
        
    elif top_level_menu_choice == "G": #geographic selection option
        before_geo = (type_geo_select,geo_selection)
        (type_geo_select,geo_selection) = geo_select_menu(type_geo_select,geo_selection)
        if before_geo != (type_geo_select,geo_selection): #if geo selection has changed
            print "Updating geographic selection..."
            if type_geo_select in ["country","state"]:
                ida = da_select_by_match(type_geo_select,geo_selection,ida)
                bda = da_select_by_match(type_geo_select,geo_selection,bda)
            #    dev_mean_binned = list(np.subtract(geos_VC_mean_binned,sat_VC_mean_binned))

            elif type_geo_select in ["rectangle","circle"]:
                ida = da_select_by_shape(geo_selection,ida)
                bda = da_select_by_shape(geo_selection,bda)                   

            #    dev_mean_binned = list(np.subtract(geos_VC_mean_binned,sat_VC_mean_binned))
    
    elif top_level_menu_choice == "C": #load and correct a new filtered dataset (lcnfd)
        lcnfd_menu_choice = "" 
        while lcnfd_menu_choice != "Z": 
        
            lcnfd_menu_title = "CURRENT SETTINGS:\n" \
                    "Loading filtered data from: %s\n" \
                    "Loading pacific correction from: %s\n" \
                    "Pacific correction type: %s\n" \
                    "Date range for data (set wide to get all data):\n"\
                    "Begin: %s\n"\
                    "End: %s"\
                    %(filtered_folder,corrections_folder,corr_type,
                    str(filtered_startdate),str(filtered_enddate))
            lcnfd_menu_text = [
                    ["1","Change filtered data folder"],
                    ["2","Change pacific correction folder"],
                    ["3","Change pacific correction type"],
                    ["4","Change dates"],
                    ["5","Load using these settings"]
                ]
            lcnfd_menu_choice = basic_menu(lcnfd_menu_title,
                                          lcnfd_menu_text)
            
            if lcnfd_menu_choice == "1": #change filtered data folder
                filtered_folder = \
                    change_var(filtered_folder,"Filtered data folder")
            elif lcnfd_menu_choice == "2": #change pacific corrections folder
                corrections_folder = \
                    change_var(corrections_folder,"Pacific corrections folder")
            elif lcnfd_menu_choice == "3": #switch correction type
                if corr_type == "basic":
                    corr_type = "advanced"
                elif corr_type == "advanced":
                    corr_type = "basic"

            elif lcnfd_menu_choice == "4": #change dates
                new_s_year = change_var(filtered_startdate.year,"Start time year")
                new_s_month= change_var(filtered_startdate.month,"Start time month")
                new_s_day  = change_var(filtered_startdate.day,"Start time day")
                filtered_startdate  = dt(new_s_year,new_s_month,new_s_day,0,0,0)
                new_e_year = change_var(  filtered_enddate.year,"End time year")
                new_e_month= change_var(  filtered_enddate.month,"End time month")
                new_e_day  = change_var(  filtered_enddate.day,"End time day")
                filtered_enddate    = dt(new_e_year,new_e_month,new_e_day,0,0,0)
                del new_s_year,new_s_month,new_s_day
                del new_e_year,new_e_month,new_e_day
            elif lcnfd_menu_choice == "5": #execute
                
                #load+correct individual data
                ida = \
                    load_and_correct(filtered_startdate,filtered_enddate,
                        filtered_folder=filtered_folder,
                        corrections_folder=corrections_folder,
                        corr_type=corr_type)
                print "New individual data has been loaded"
                
                #bin data
                print "Now binning this data"
                bda = create_binned_set(ida,map_box,xdim,ydim)
                #dev_mean_binned = list(np.subtract(geos_VC_mean_binned,sat_VC_mean_binned))
                #NDVI should be coded in at this point
                
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
                print "Assigning country+state to indiv. data points"
                (country,state) =\
                    region_matcher(csv_lat,csv_lon,csv_country,csv_state,
                                   ida.lat,ida.lon)
                                   
                #add these to ida properly
                ida.data['country'] = d(country,'country','Country observation falls within')
                del country
                ida.data['state'  ] = d(state  ,'state'  ,'State observation falls within')
                del state
                                   
                print "Assigning country+state to binned data points"
                (country_binned,state_binned) =\
                    region_matcher(csv_lat,csv_lon,csv_country,csv_state,
                                   bda.lat,bda.lon)
                bda.data["country"] = d(country_binned,"country","Country")
                del country_binned
                bda.data["state"  ] = d(state_binned  ,"state"  ,"State"  )
                del state_binned
                print "Done. Do you wish to save these data as new pickles? Y/N"
                option = basic_menu("Done. Do you wish to save"
                                    " these data as new pickles [recommended]?",
                                    [["Y","yes"],["N","no"]],quit_option=False)
                if option == "Y":
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
                    save_pickle(bda,save_path=save_path,prefix=prefix,suffix="_binned.p")
                    bda_p = copy.deepcopy(bda)
                    current_pickle = save_path+prefix
                
                del option
                  
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
        save_pickle(bda,save_path=save_path,prefix=prefix,suffix="_binned.p")
        bda_p = copy.deepcopy(bda) #update the static store
        current_pickle = save_path+prefix
        
    elif top_level_menu_choice == "1": #Binned data
        #ubd = use binned data
        ubd_menu_choice = ""
        while ubd_menu_choice != "Z":

            ubd_menu_title = "USING BINNED DATA"
            ubd_menu_text  = [
                ["1","Plot gridded data on map"],
                ["2","Compute basic statistics"],
                ["3","Compare datasets"]
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
                                        save=do_save,save_filename=save_filename)
                    break
            
            #bda-ification done to this point
                    
            elif ubd_menu_choice == "2": #binned statisics
                stats_from_da(bda)
                
            elif ubd_menu_choice == "3": #compare_datasets 
                two_var_comparison(bda)
                     
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
                ["5","Oversample"]
                ]
            uid_menu_choice = basic_menu(uid_menu_title,
                                         uid_menu_text,
                                         quit_option=True)
            if uid_menu_choice == "1": #dots on map
                #dom = dots on map
                dom_menu_choice = ""
                while dom_menu_choice != "Z":
                    dom_menu_title = "Plotting map of binned data\n"\
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
                            plot_dots_on_map(ida.lat,ida.lon,ida.data[data_key].val,map_box,vmin=vmin,vmax=vmax,title=title,lab=unit,save=False)
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
            
            #ida-ification done to this point    
            elif uid_menu_choice == "5": #oversampler
                #print len(ida.lat)
                #raw_input("halt-->")
                os_ida = copy.deepcopy(ida)
                os_data = oversample(os_ida)
                #print len(ida.lat)
                #raw_input("halt-->")

        
               
        
       
