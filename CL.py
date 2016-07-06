#!/usr/bin/env python

from CL_satpp import *
import os
from datetime import datetime as dt

def main():
    """The over-arching top level function"""
    
    #initialisation
    initial_choice = basic_menu(
                    "INITIAL MENU",
                    [
                     ["1","Load default files and settings"],
                     ["2","Proceed without loading"]
                    ],
                    quit_option=False)
    if initial_choice == "1":
        print "Loading default options"
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
        current_pickle =   '/group_workspaces/jasmin/geoschem/'\
                           'local_users/lsurl/sat_pp/'\
                           '20140301-20141231'
        current_emfiles = ['/group_workspaces/jasmin/geoschem'\
                           '/local_users/lsurl/runs' \
                           '/geosfp_025x03125_nochem_2014/',
                           '_20140301-20141231']
        (NDVI_data_all,indv_data_all,binn_data_all) = \
            load_new_pickles_all(current_pickle,verbose=True)
        (NDVI_lat,NDVI_lon,NDVI,NDVI_year,NDVI_month) = NDVI_data_all
        (lat,lon,
            sat_VC,sat_DVC,
            geos_VC,time,
            country,state,
            index_map) = indv_data_all
        (lat_binned,lon_binned,
            sat_mean_binned,geos_mean_binned,
            sat_uncer_binned,geos_stdev_binned,
            dev_mean_binned,NDVI_mean_binned,
            countries_binned,states_binned) = binn_data_all
    elif initial_choice == "2":
        south = 0.
        north = 0.
        west = 0.
        east = 0.
        compass = [north,south,east,west]
        xdim = 0.
        ydim = 0.
        type_geo_select = "latlon"
        geo_selection = [north,south,east,west]  
        startdate = dt(2014,01,01, 0, 0, 0)
        enddate   = dt(2014,01,01, 0, 0, 0)        
        current_pickle = "NONE SELECTED"
        current_emfiles = ["NONE SELECTED","NONE SELECTED"]
        
    
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
                (NDVI_data_all,indv_data_all,binn_data_all) \
                    = load_new_pickles_all(current_pickle,verbose=True)
                (NDVI_lat,NDVI_lon,NDVI,NDVI_year,NDVI_month) \
                    = NDVI_data_all
                (lat,lon,
                    sat_VC,sat_DVC,
                    geos_VC,time,
                    country,state,
                    index_map) = indv_data_all
                (lat_binned,lon_binned,
                    sat_mean_binned,geos_mean_binned,
                    sat_uncer_binned,geos_stdev_binned,
                    dev_mean_binned,NDVI_mean_binned,
                    countries_binned,states_binned) = binn_data_all
            
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
            (lat_binned,lon_binned,sat_mean_binned,geos_mean_binned,sat_uncer_binned,geos_stdev_binned,dev_mean_binned,NDVI_mean_binned,countries_binned,states_binned) = binn_data_all
            
            
        elif top_level_menu_choice == "G": #geographic selection option
            before_geo = (type_geo_select,geo_selection)
            (type_geo_select,geo_selection) = geo_select_menu(type_geo_select,geo_selection)
            if before_geo != (type_geo_select,geo_selection): #if geo selection has changed
                print "Updating geographic selection..."
                if type_geo_select == "country":
                    (lat,lon,sat_VC,sat_DVC,geos_VC,time,country,state,index_map) = geo_select_regional(country,geo_selection,lat,lon,sat_VC,sat_DVC,geos_VC,time,country,state,index_map)
                    (lat_binned,lon_binned,sat_mean_binned,geos_mean_binned,sat_uncer_binned,geos_stdev_binned,dev_mean_binned,NDVI_mean_binned,countries_binned,states_binned) = geo_select_regional(countries_binned,geo_selection,lat_binned,lon_binned,sat_mean_binned,geos_mean_binned,sat_uncer_binned,geos_stdev_binned,dev_mean_binned,NDVI_mean_binned,countries_binned,states_binned)
                elif type_geo_select == "state":
                    (lat,lon,sat_VC,sat_DVC,geos_VC,time,country,state,index_map) = geo_select_regional(state,geo_selection,lat,lon,sat_VC,sat_DVC,geos_VC,time,country,state,index_map)
                    (lat_binned,lon_binned,sat_mean_binned,geos_mean_binned,sat_uncer_binned,geos_stdev_binned,dev_mean_binned,NDVI_mean_binned,countries_binned,states_binned) = geo_select_regional(states_binned,geo_selection,lat_binned,lon_binned,sat_mean_binned,geos_mean_binned,sat_uncer_binned,geos_stdev_binned,dev_mean_binned,NDVI_mean_binned,countries_binned,states_binned)
                elif type_geo_select == "rectangle":                   
                    (lat,lon,sat_VC,sat_DVC,geos_VC,time,country,state,index_map) = geo_select_rectangle(lat,lon,geo_selection,lat,lon,sat_VC,sat_DVC,geos_VC,time,country,state,index_map)
                    (lat_binned,lon_binned,sat_mean_binned,geos_mean_binned,sat_uncer_binned,geos_stdev_binned,dev_mean_binned,NDVI_mean_binned,countries_binned,states_binned) = geo_select_rectangle(lat_binned,lon_binned,geo_selection,lat_binned,lon_binned,sat_mean_binned,geos_mean_binned,sat_uncer_binned,geos_stdev_binned,dev_mean_binned,NDVI_mean_binned,countries_binned,states_binned)
                elif type_geo_select == "circle":
                    (lat,lon,sat_VC,sat_DVC,geos_VC,time,country,state,index_map) = geo_select_circle(lat,lon,geo_selection,lat,lon,sat_VC,sat_DVC,geos_VC,time,country,state,index_map)
                    (lat_binned,lon_binned,sat_mean_binned,geos_mean_binned,sat_uncer_binned,geos_stdev_binned,dev_mean_binned,NDVI_mean_binned,countries_binned,states_binned) = geo_select_circle(lat_binned,lon_binned,geo_selection,lat_binned,lon_binned,sat_mean_binned,geos_mean_binned,sat_uncer_binned,geos_stdev_binned,dev_mean_binned,NDVI_mean_binned,countries_binned,states_binned)
           
            
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
                            dataset = sat_mean_binned
                        elif binned_map_menu_choice == "2": #uncertainty in satellite observations
                            (map_preplot_menu_choice,title,vmin,vmax,units,save_filename) = map_preplot_menu("Uncertainty in mean satellite observation",vmin=0.,vmax=3.e16,units="molec.cm-2")
                            dataset = sat_uncer_binned
                        elif binned_map_menu_choice == "3": #mean modelled values
                            (map_preplot_menu_choice,title,vmin,vmax,units,save_filename) = map_preplot_menu("Mean modelled columns",vmin=0.,vmax=3.e16,units="molec.cm-2")
                            dataset = geos_mean_binned
                        elif binned_map_menu_choice == "4": #standard deviation in modelled values
                            (map_preplot_menu_choice,title,vmin,vmax,units,save_filename) = map_preplot_menu("Standard deviation of modelled columns",vmin=0.,vmax=3.e16,units="molec.cm-2")
                            dataset = geos_stdev_binned
                        elif binned_map_menu_choice == "5": #mean NDVI value
                            (map_preplot_menu_choice,title,vmin,vmax,units,save_filename) = map_preplot_menu("Mean NDVI value",vmin=0.,vmax=1.,units="(unitless)")
                            dataset = NDVI_mean_binned
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
                            stat = calc_statistic(sat_mean_binned,  stat_choice)
                        elif dataset_choice == "2": #uncertainty in satellite observations
                            stat = calc_statistic(sat_uncer_binned, stat_choice)
                        elif dataset_choice == "3": #mean modelled values
                            stat = calc_statistic(geos_mean_binned, stat_choice)
                        elif dataset_choice == "4": #standard deviation in modelled values
                            stat = calc_statistic(geos_stdev_binned,stat_choice)
                        elif dataset_choice == "5": #mean NDVI value
                            stat = calc_statistic(NDVI_mean_binned, stat_choice)
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
                    dots_on_map_menu_choice == ""
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
                            #At this point binned_map_preplot_menu_choice will be either Z (up menu), P (plot) or S (save)                    
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
                                                          
                             
                             
                             
main()
        
               
        
       