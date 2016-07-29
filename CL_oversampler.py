#!/usr/bin/env python
"""Utility to "oversample" satellite observations"""

from geopy.distance import great_circle
from CL_satpp import *
import numpy as np


def oversample(var_names,vartuple,lat,lon):
    """User interface for the oversampler"""
    (var,var_name) = select_a_var(var_names,vartuple,
                                  "oversampling",numbers_only=True)
    print "Oversampling regrids to a very fine grid."
    print "The finer this grid, the greater the detail that can potentially be seen"
    print "However this also comes at a higher computational cost"
    print "What will be the dimension, in degrees, of the squares of the super-fine grid?"
    print "(0.02 degrees recommended)"
    fine_dim = float(raw_input("-->"))
    print "What radius, in km, will the %s observations be spread over?" %var_name
    print "(24 km recommeneded)"
    averaging_radius = float(raw_input("-->"))
    print "Please enter the bounds of the area you wish to inspect"
    north_view = float(raw_input("NORTH-->"))
    south_view = float(raw_input("SOUTH-->"))
    east_view  = float(raw_input("EAST -->"))
    west_view  = float(raw_input("WEST -->"))

    #at the edges, we will need to sample points outside these bounds
    #Thus we define sample bounds which are slightly larger

    north_sample = north_view
    north_OK = False
    while north_OK == False:
        print north_sample
        north_sample += 0.1 #move buffer 0.1 deg away each iteration and check if alright
        north_OK = float( great_circle((north_view           ,0.5*(east_view+west_view)),
                                  (north_sample         ,0.5*(east_view+west_view))).km)\
                                  > averaging_radius*1.2 

    south_sample = south_view
    south_OK = False
    while south_OK == False:
        south_sample -= 0.1 #move buffer 0.1 deg away each iteration and check if alright
        south_OK = float( great_circle((south_view           ,0.5*(east_view+west_view)),
                                  (south_sample         ,0.5*(east_view+west_view))).km)\
                                  > averaging_radius*1.2 

    east_sample = east_view
    east_OK = False
    while east_OK == False:
        east_sample += 0.1 #move buffer 0.1 deg away each iteration and check if alright
        east_OK = float( great_circle((0.5*(north_view+south_view),east_view           ),
                                 (0.5*(north_view+south_view),east_sample         )).km)\
                                  > averaging_radius*1.2


    west_sample = west_view
    west_OK = False
    while west_OK == False:
        west_sample -= 0.1 #move buffer 0.1 deg away each iteration and check if alright
        west_OK = float( great_circle((0.5*(north_view+south_view),west_view           ),
                                 (0.5*(north_view+south_view),west_sample         )).km)\
                                  > averaging_radius*1.2    


    #now use geo_select_rectangle to restrict our datasets to just those in the sampling box
    #this saves a lot of computational time!

    sample_box = [north_sample,south_sample,east_sample,west_sample]
    (lat,lon,var) = geo_select_rectangle(lat,lon,sample_box,lat,lon,var)

    #Define the lat and lon co-ordinates of the ultra-fine grid.

    lat_fine = []
    lon_fine = []
    
    lat_step = south_view
    while lat_step <= north_view:
        lon_step = west_view
        while lon_step <= east_view:
            lat_fine.append(lat_step)
            lon_fine.append(lon_step)
            lon_step += fine_dim
        lat_step += fine_dim

    var_fine=[]
    var_fine_count=[]
    #Now, for every cell in the ultrafine grid, determine which var are within averaging_distance
    #This will likely take a while

    fine_points=len(lat_fine)
    for i in range(0,fine_points):
        clearscreen()
        print "Processing point %i of %i" %(i,fine_points)
        this_set = [var[j] for j in range(0,len(var)) if 
                                     great_circle( (lat[j],lon[j]),
                                                   (lat_fine[i],lon_fine[i])).km
                                     <= averaging_radius ]
        var_fine.append(np.nanmean(this_set))
        var_fine_count.append(notnancount(this_set))
        del this_set

    plot_dots_on_map(lat_fine,lon_fine,var_fine,
                     north_view,south_view,east_view,west_view,vmin=0.,vmax=3.0E16,
                     title="Oversampled plot",lab=var_name,
                     save=False,save_filename="nofilenamespecified")
    plot_dots_on_map(lat_fine,lon_fine,var_fine_count,
                     north_view,south_view,east_view,west_view,vmin=0,vmax=np.nanmax(var_fine_count),
                     title="Oversampled plot",lab="data count at location")








            
              
