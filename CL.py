#!/usr/bin/env python

from CL_satpp import *
from CL_corrections import *
import os
from datetime import datetime as dt

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
south = 2.
north = 38.
west = 65.
east = 100.
compass = [north,south,east,west]
xdim = 0.3125
ydim = 0.25
startdate = dt(2014,03,01, 0, 0, 0)
enddate   = dt(2014,12,31,23,59,59)
type_geo_select = "latlon"
geo_selection = [north,south,east,west]        
current_pickle = "NONE SELECTED"
current_emfiles = ["NONE SELECTED","NONE SELECTED"]   

if initial_choice == "1":
    print "Loading default options"       
    current_pickle =   '/group_workspaces/jasmin/geoschem/'\
                       'local_users/lsurl/sat_pp/'\
                       '20140301-20141231'
    current_emfiles = ['/group_workspaces/jasmin/geoschem'\
                       '/local_users/lsurl/runs' \
                       '/geosfp_025x03125_nochem_2014/',
                       '_20140301-20141231']
    (indv_data_all,binn_data_all) = \
        load_new_pickles_all(current_pickle,verbose=True,NDVI=False)
    #(NDVI_lat,NDVI_lon,NDVI,NDVI_year,NDVI_month) = NDVI_data_all
    (lat,lon,
        sat_VC,sat_DVC,
        geos_VC,time,
        country,state,
        index_map) = indv_data_all
    (lat_binned,lon_binned,
        sat_VC_mean_binned,sat_VC_stdev_binned,
        sat_DVC_mean_binned,
        geos_VC_mean_binned,geos_VC_stdev_binned,
        country_binned,state_binned) = binn_data_all
    dev_mean_binned = list(np.subtract(geos_VC_mean_binned,sat_VC_mean_binned))
    

#main loop
top_level_menu_choice = ""
while top_level_menu_choice != "Z": #loop unless ordered to quit

    top_level_menu_title = "TOP LEVEL MENU\n" \
               "Current pickles loaded: %s\n" \
               "Current emissions files are: %sXXX%s" \
               %(current_pickle,current_emfiles[0],current_emfiles[1])
    top_level_menu_text = [
            ["1","Use binned data"],
            ["2","Use individual observations data"],
            ["G","Geographically select data to use"],
            ["C","Generate new dataset from filtered data and pacific corrections"],
            ["S","Save current dataset as new pickle"],
            ["P","Change pickle"],
            ["R","Reload pickle"],
            ["E","Change current emissions files"],   
            ["X","Inspect and change global options"]
        ]
    top_level_menu_choice = basic_menu(top_level_menu_title,
                                       top_level_menu_text)
    
    if top_level_menu_choice == "P": #change pickle
        (changed,current_pickle) = change_pickle(current_pickle)
        if changed:
            #reloading is time-consuming so only do if change made
            (indv_data_all,binn_data_all) \
                = load_new_pickles_all(current_pickle,verbose=True,NDVI=False)
            #(NDVI_lat,NDVI_lon,NDVI,NDVI_year,NDVI_month) \
            #    = NDVI_data_all
            (lat,lon,
                sat_VC,sat_DVC,
                geos_VC,time,
                country,state,
                index_map) = indv_data_all
            (lat_binned,lon_binned,
             sat_VC_mean_binned,sat_VC_stdev_binned,
             sat_DVC_mean_binned,
             geos_VC_mean_binned,geos_VC_stdev_binned,
             country_binned,state_binned) = binn_data_all
            dev_mean_binned = list(np.subtract(geos_VC_mean_binned,sat_VC_mean_binned))
        
    elif top_level_menu_choice == "E": #change emissions
        (changed,current_emfiles) = change_emfiles(current_emfiles)
        #currently no loading of the emisisons files are done here
    
    elif top_level_menu_choice == "X": #change options
        datesbefore = (startdate,enddate)
        (north,south,east,west,
         xdim,ydim,
         startdate,enddate) = global_options(
                             north,south,east,west,
                             xdim,ydim,
                             startdate,enddate)
        #if time selection has changed, reselect data.
        if datesbefore != (startdate,enddate): 
            print "Start and/or end time has changed. Selecting data..."                
            (lat,lon,
             sat_VC,sat_DVC,
             geos_VC,time,
             country,state,
             index_map) = time_select(
                startdate,enddate,time,
                lat,lon,sat_VC,sat_DVC,geos_VC,time,
                country,state,index_map)
            
    elif top_level_menu_choice == "R": #reload pickle
        #re-unpack data sets to clear any previous geographic selection
        print "Reloading..."
        (lat,lon,sat_VC,sat_DVC,geos_VC,time,country,state,index_map) = indv_data_all
        (lat_binned,lon_binned,
             sat_VC_mean_binned,sat_VC_stdev_binned,
             sat_DVC_mean_binned,
             geos_VC_mean_binned,geos_VC_stdev_binned,
             country_binned,state_binned) = binn_data_all
        dev_mean_binned = list(np.subtract(geos_VC_mean_binned,sat_VC_mean_binned))
        
        
    elif top_level_menu_choice == "G": #geographic selection option
        before_geo = (type_geo_select,geo_selection)
        (type_geo_select,geo_selection) = geo_select_menu(type_geo_select,geo_selection)
        if before_geo != (type_geo_select,geo_selection): #if geo selection has changed
            print "Updating geographic selection..."
            if type_geo_select == "country":
                (lat,lon,sat_VC,sat_DVC,geos_VC,
                 time,country,state,index_map) =\
                    geo_select_regional(country,geo_selection,
                    lat,lon,sat_VC,sat_DVC,geos_VC,
                    time,country,state,index_map)
                (lat_binned,lon_binned,
                 sat_VC_mean_binned,sat_VC_stdev_binned,
                 sat_DVC_mean_binned,
                 geos_VC_mean_binned,geos_VC_stdev_binned,
                 countries_binned,states_binned) = geo_select_regional(
                    country_binned,geo_selection,
                    lat_binned,lon_binned,
                    sat_VC_mean_binned,sat_VC_stdev_binned,
                    sat_DVC_mean_binned,
                    geos_VC_mean_binned,geos_VC_stdev_binned,
                    country_binned,state_binned)
                dev_mean_binned = list(np.subtract(geos_VC_mean_binned,sat_VC_mean_binned))
            elif type_geo_select == "state":
                (lat,lon,sat_VC,sat_DVC,geos_VC,
                 time,country,state,index_map) =\
                    geo_select_regional(state,geo_selection,
                    lat,lon,sat_VC,sat_DVC,geos_VC,
                    time,country,state,index_map)
                (lat_binned,lon_binned,
                 sat_VC_mean_binned,sat_VC_stdev_binned,
                 sat_DVC_mean_binned,
                 geos_VC_mean_binned,geos_VC_stdev_binned,
                 countries_binned,states_binned) = geo_select_regional(
                    state_binned,geo_selection,
                    lat_binned,lon_binned,
                    sat_VC_mean_binned,sat_VC_stdev_binned,
                    sat_DVC_mean_binned,
                    geos_VC_mean_binned,geos_VC_stdev_binned,
                    country_binned,state_binned)
                dev_mean_binned = list(np.subtract(geos_VC_mean_binned,sat_VC_mean_binned))
            elif type_geo_select == "rectangle":                   
                (lat,lon,sat_VC,sat_DVC,geos_VC,
                 time,country,state,index_map) =\
                    geo_select_rectangle(lat,lon,geo_selection,
                    lat,lon,sat_VC,sat_DVC,geos_VC,
                    time,country,state,index_map)
                (lat_binned,lon_binned,
                 sat_VC_mean_binned,sat_VC_stdev_binned,
                 sat_DVC_mean_binned,
                 geos_VC_mean_binned,geos_VC_stdev_binned,
                 countries_binned,states_binned) = geo_select_regional(
                    lat_binned,lon_binned,geo_selection,
                    lat_binned,lon_binned,
                    sat_VC_mean_binned,sat_VC_stdev_binned,
                    sat_DVC_mean_binned,
                    geos_VC_mean_binned,geos_VC_stdev_binned,
                    country_binned,state_binned)
                dev_mean_binned = list(np.subtract(geos_VC_mean_binned,sat_VC_mean_binned))
            elif type_geo_select == "circle":
                (lat,lon,sat_VC,sat_DVC,geos_VC,
                 time,country,state,index_map) =\
                    geo_select_circle(lat,lon,geo_selection,
                    lat,lon,sat_VC,sat_DVC,geos_VC,
                    time,country,state,index_map)
                (lat_binned,lon_binned,
                 sat_VC_mean_binned,sat_VC_stdev_binned,
                 sat_DVC_mean_binned,
                 geos_VC_mean_binned,geos_VC_stdev_binned,
                 countries_binned,states_binned) = geo_select_circle(           
                    lat_binned,lon_binned,geo_selection,
                    lat_binned,lon_binned,
                    sat_VC_mean_binned,sat_VC_stdev_binned,
                    sat_DVC_mean_binned,
                    geos_VC_mean_binned,geos_VC_stdev_binned,
                    country_binned,state_binned)
                dev_mean_binned = list(np.subtract(geos_VC_mean_binned,sat_VC_mean_binned))
    
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
                    corr_type == "advanced"
                elif corr_type == "advanced":
                    corr_type == "basic"
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
                (ULN,lat,lon,time,geos_VC,sat_VC,sat_DVC,AMF) = \
                    load_and_correct(filtered_startdate,filtered_enddate,
                        filtered_folder=filtered_folder,
                        corrections_folder=corrections_folder,
                        corr_type=corr_type)
                print "New individual data has been loaded"
                
                #bin data
                print "Now binning this data"
                (lat_binned,lon_binned,xedges,yedges,index_map,
                 sat_VC_mean_binned,sat_VC_stdev_binned,sat_VC_count_binned,
                 sat_DVC_mean_binned,
                 geos_VC_mean_binned,geos_VC_stdev_binned) = \
                 create_binned_set(lat,lon,geos_VC,sat_VC,sat_DVC,
                                   north,south,east,west,xdim,ydim)
                dev_mean_binned = list(np.subtract(geos_VC_mean_binned,sat_VC_mean_binned))
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
                                   lat,lon)
                print "Assigning country+state to binned data points"
                (country_binned,state_binned) =\
                    region_matcher(csv_lat,csv_lon,csv_country,csv_state,
                                   lat_binned,lon_binned)
                print "Done. Do you wish to save these data as new pickles? Y/N"
                option = basic_menu("Done. Do you wish to save"
                                    " these data as new pickles [recommended]?",
                                    [["Y","yes"],["N","no"]],quit_option=False)
                if option == "Y":
                    current_pickle = save_pickles(lat,lon,sat_VC,sat_DVC,geos_VC,time,
                                                  country,state,index_map,
                                                  lat_binned,lon_binned,
                                                  sat_VC_mean_binned,sat_VC_stdev_binned,
                                                  sat_DVC_mean_binned,
                                                  geos_VC_mean_binned,geos_VC_stdev_binned,
                                                  None,
                                                  country_binned,state_binned,index_map)
                
                del option
                  
    elif top_level_menu_choice == "S":
        current_pickle = save_pickles(lat,lon,sat_VC,sat_DVC,geos_VC,time,
                                      country,state,index_map,
                                      lat_binned,lon_binned,
                                      sat_VC_mean_binned,sat_VC_stdev_binned,
                                      sat_DVC_mean_binned,
                                      geos_VC_mean_binned,geos_VC_stdev_binned,
                                      None,
                                      country_binned,state_binned,index_map)
        
    elif top_level_menu_choice == "1": #Binned data
        quit_binned_data_menu = False
        while quit_binned_data_menu == False:
            os.system('clear')
            use_binned_data_menu_choice = use_binned_data_menu()
                   
            if use_binned_data_menu_choice == "1": #Binned map
                quit_binned_map_menu = False
                while quit_binned_map_menu == False:                    
                    os.system('clear')
                    binned_map_menu_choice = binned_map_menu()
                    if   binned_map_menu_choice == "1": #mean satellite observations
                        (map_preplot_menu_choice,title,vmin,vmax,units,save_filename) = map_preplot_menu("Mean satellite observation",vmin=0.,vmax=3.e16,units="molec.cm-2")
                        dataset = sat_VC_mean_binned
                    elif binned_map_menu_choice == "2": #uncertainty in satellite observations
                        (map_preplot_menu_choice,title,vmin,vmax,units,save_filename) = map_preplot_menu("Uncertainty in mean satellite observation",vmin=0.,vmax=3.e16,units="molec.cm-2")
                        dataset = sat_VC_stdev_binned
                    elif binned_map_menu_choice == "3": #mean modelled values
                        (map_preplot_menu_choice,title,vmin,vmax,units,save_filename) = map_preplot_menu("Mean modelled columns",vmin=0.,vmax=3.e16,units="molec.cm-2")
                        dataset = geos_VC_mean_binned
                    elif binned_map_menu_choice == "4": #standard deviation in modelled values
                        (map_preplot_menu_choice,title,vmin,vmax,units,save_filename) = map_preplot_menu("Standard deviation of modelled columns",vmin=0.,vmax=3.e16,units="molec.cm-2")
                        dataset = geos_VC_stdev_binned
                    elif binned_map_menu_choice == "5": #mean NDVI value
                        print "NDVI option not yet available"
                        #(map_preplot_menu_choice,title,vmin,vmax,units,save_filename) = map_preplot_menu("Mean NDVI value",vmin=0.,vmax=1.,units="(unitless)")
                        #dataset = NDVI_mean_binned
                    elif binned_map_menu_choice == "6": #sigma deviation of model from observations
                        (map_preplot_menu_choice,title,vmin,vmax,units,save_filename) = map_preplot_menu("Deviation of model from observations",vmin=-2.5,vmax=2.5,units="sigmas")
                        dataset = dev_mean_binned
                    elif binned_map_menu_choice == "Z": #go up a menu
                         quit_binned_map_menu = True
                         
                    if map_preplot_menu_choice == "": #allow for quick no-entrying
                        map_preplot_menu_choice = "P"
                                                    
                    if binned_map_menu_choice != "Z" :                                               
                        #At this point binned_map_preplot_menu_choice will be either Z (up menu), P (plot) or S (save)                    
                        if   map_preplot_menu_choice.upper() == "S": #saving the figure
                            plot_grid_from_list(lat_binned,lon_binned,dataset,xdim,ydim,north,south,east,west,title=title,vmin=vmin,vmax=vmax,lab=units,save=True,save_filename=save_filename)
                        elif map_preplot_menu_choice.upper() == "P": #plotting the figure
                            plot_grid_from_list(lat_binned,lon_binned,dataset,xdim,ydim,north,south,east,west,title=title,vmin=vmin,vmax=vmax,lab=units,save=False)
                        else: #no plot
                            pass
            elif use_binned_data_menu_choice == "2": #statisics
                [dataset_choice,stat_choice] = basic_statistics_menu("binned")
                if [dataset_choice,stat_choice] != ["Z","Z"]: #unless we're quitting
                    if dataset_choice == "1": #Mean satellite observations
                        stat = calc_statistic(sat_VC_mean_binned,  stat_choice)
                    elif dataset_choice == "2": #uncertainty in satellite observations
                        stat = calc_statistic(sat_VC_stdev_binned, stat_choice)
                    elif dataset_choice == "3": #mean modelled values
                        stat = calc_statistic(geos_VC_mean_binned, stat_choice)
                    elif dataset_choice == "4": #standard deviation in modelled values
                        stat = calc_statistic(geos_VC_stdev_binned,stat_choice)
                    elif dataset_choice == "5": #mean NDVI value
                        print "NDVI option not yet available"
                        #stat = calc_statistic(NDVI_mean_binned, stat_choice)
                    elif dataset_choice == "6": #sigma deviation of model from observations
                        stat = calc_statistic(dev_mean_binned , stat_choice)
                    print stat
                    null = raw_input("Press enter to continue-->")
                    
            if use_binned_data_menu_choice == "Z": #Exit                   
                quit_binned_data_menu = True
    
    elif top_level_menu_choice == "2": #Indiv data
        use_indiv_data_menu_choice = ""
        while use_indiv_data_menu_choice.upper() != "Z":
            os.system('clear')
            use_indiv_data_menu_choice = use_indiv_data_menu()
            if use_indiv_data_menu_choice == "1": #dots on map
                dots_on_map_menu_choice = ""
                while dots_on_map_menu_choice.upper() != "Z":                    
                    os.system('clear')
                    dots_on_map_menu_choice = dots_on_map_menu()                    
                    if   dots_on_map_menu_choice == "1": #mean satellite observations
                        (map_preplot_menu_choice,title,vmin,vmax,units,save_filename) = map_preplot_menu("Satellite observation",vmin=0.,vmax=3.e16,units="molec.cm-2")
                        dataset = sat_VC
                    elif dots_on_map_menu_choice == "2": #uncertainty in satellite observations
                        (map_preplot_menu_choice,title,vmin,vmax,units,save_filename) = map_preplot_menu("Uncertainty in satellite observation",vmin=0.,vmax=3.e16,units="molec.cm-2")
                        dataset = sat_DVC
                    elif dots_on_map_menu_choice == "3": #mean modelled values
                        (map_preplot_menu_choice,title,vmin,vmax,units,save_filename) = map_preplot_menu("Modelled columns",vmin=0.,vmax=3.e16,units="molec.cm-2")
                        dataset = geos_VC
                    
                    if map_preplot_menu_choice == "": #allow for quick no-entrying
                        map_preplot_menu_choice = "P"
                    
                    if dots_on_map_menu_choice != "Z" :                                               
                        #At this point map_preplot_menu_choice will be either Z (up menu), P (plot) or S (save)                    
                        if   map_preplot_menu_choice.upper() == "S": #saving the figure
                            plot_dots_on_map(lat,lon,dataset,north,south,east,west,vmin=vmin,vmax=vmax,title=title,lab=units,save=True,save_filename=save_filename)
                        elif map_preplot_menu_choice.upper() == "P": #plotting the figure
                            plot_dots_on_map(lat,lon,dataset,north,south,east,west,vmin=vmin,vmax=vmax,title=title,lab=units,save=False)
                        else: #no plot
                            pass     
            elif use_indiv_data_menu_choice == "2": #statisics
                [dataset_choice,stat_choice] = basic_statistics_menu("indiv")
                if [dataset_choice,stat_choice] != ["Z","Z"]: #unless we're quitting
                    if dataset_choice == "1": #Mean satellite observations
                        stat = calc_statistic(sat_VC,  stat_choice)
                    elif dataset_choice == "2": #uncertainty in satellite observations
                        stat = calc_statistic(sat_DVC, stat_choice)
                    elif dataset_choice == "3": #mean modelled values
                        stat = calc_statistic(geos_VC, stat_choice)
                    print stat
                    null = raw_input("Press enter to continue-->")
            
            elif use_indiv_data_menu_choice == "3": #timeseries statistics
                [dataset_choice,stat_choice,step_days] = timeseries_statistics_menu()
                if [dataset_choice,stat_choice,step_days] != ["Z","Z",0]: #unless we're quitting
                    if dataset_choice == "1": #mean satellite observations
                        (stat_series,time_series) = time_cycle(startdate,enddate,step_days,time,sat_VC,stat_choice)
                    elif dataset_choice == "2": #uncertainty in satellite observations
                        (stat_series,time_series) = time_cycle(startdate,enddate,step_days,time,sat_DVC,stat_choice)
                    elif dataset_choice == "3": #mean satellite observations
                        (stat_series,time_series) = time_cycle(startdate,enddate,step_days,time,geos_VC,stat_choice)
                    print stat_series
                    print time_series
                    null = raw_input("Press enter to continue-->")
          
            elif use_indiv_data_menu_choice == "4": #compare datasets
                two_var_comparison(indiv_varnames,
                                   (lat,lon,sat_VC,sat_DVC,geos_VC,time,country,state))

        
               
        
       
