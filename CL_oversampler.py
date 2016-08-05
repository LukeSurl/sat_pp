#!/usr/bin/env python
"""Utility to "oversample" satellite observations"""

from geopy.distance import great_circle
from CL_satpp import *
import numpy as np


class os_data:
    """A class for oversampled data"""
    
    def __init__(self,lat,lon,data,name="unspecified"):
       self.lat = lat
       self.lon = lon
       self.data = data
       self.meta = {'name':name}
    
    def __str__(self):
        return "Oversampled %s" %self.name
       
    def set_unit(self,unit):
        self.unit = unit
       
    def set_meta(self,key,value):
        #for adding meta data to this dataset
        self.meta[key] = value
    

def oversample(var_names,vartuple,lat,lon):
    """User interface for the oversampler"""
    (var,var_name) = select_a_var(var_names,vartuple,
                                  "oversampling",numbers_only=True)
                                  
    print "Please enter the bounds of the area you wish to inspect"
    north_view = float(raw_input("NORTH-->"))
    south_view = float(raw_input("SOUTH-->"))
    east_view  = float(raw_input("EAST -->"))
    west_view  = float(raw_input("WEST -->"))                              
                                  
                                  
    print "Oversampling regrids to a very fine grid."
    print "The finer this grid, the greater the detail that can potentially be seen"
    print "However this also comes at a higher computational cost"
    print "What will be the dimension, in degrees, of the squares of the super-fine grid?"
    print "(0.02 degrees recommended)"
    fine_dim = float(raw_input("-->"))
    
    os_type = basic_menu(
                "How should individual observations be ''spread out'' by the oversampling process?",
                [["1","As circles with a defined radius (in km)"]
                 ,["2","As rectangles with defined dimensions (in degrees)"]
                 #,["3","As ''super Guassians'' [as per de Graaf et al (2016)]"]
                 ],
                 quit_option=False)
    
    if os_type == "1":    
        print "What radius, in km, will the %s observations be spread over?" %var_name
        print "(24 km recommended)"
        averaging_radius = float(raw_input("-->"))
        #at the edges, we will need to sample points outside these bounds
        #Thus we define sample bounds which are slightly larger

        #we'll assume, for safety, that 1km = 0.01 degrees; a slight overestimation
        #(even at the equator), but this is useful.

        buffer_degrees = 0.01*averaging_radius
        north_sample = north_view + buffer_degrees
        south_sample = south_view - buffer_degrees
        east_sample  = east_view  + buffer_degrees
        west_sample  = west_view  - buffer_degrees
    elif os_type == "2":
        print "What is the LATITUDE dimension of these rectangles in degrees?"
        print "1 degree latitude = 111 km; 1 km = 0.0090 degrees latitude"
        averaging_rectangle_lat_dim = float(raw_input("-->"))
        
        north_sample = north_view + 0.5*averaging_rectangle_lat_dim
        south_sample = south_view - 0.5*averaging_rectangle_lat_dim
        
        print "What is that LONGITUDE dimension of these rectangles in degrees?"
        print "Length of a degree longitude depends on latitude"
        midpoint = (north_view + south_view)/2
        km_per_degree_at_midpoint = 111.*np.cos(np.deg2rad(midpoint))
        degree_per_km_at_midpoint = 1./km_per_degree_at_midpoint
        print "At latitude of %.0f degrees, 1 degree longitude = %.0f km; 1 km = %.4f degrees" \
               %(midpoint,km_per_degree_at_midpoint,degree_per_km_at_midpoint)
        averaging_rectangle_lon_dim = float(raw_input("-->"))
        
        east_sample  = east_view  + 0.5*averaging_rectangle_lon_dim
        west_sample  = west_view  - 0.5*averaging_rectangle_lon_dim 

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
    
    if os_type == "1":
        #for circle oversampling
        for i in range(0,fine_points):
            if np.remainder(i,100) == 0: #only update screen occasionally
                clearscreen()
                print "Processing point %i of %i" %(i,fine_points)
            #for each point, first cut down a square area around the central point before doing the
            #the circle
            north_square = lat_fine[i] + buffer_degrees
            south_square = lat_fine[i] - buffer_degrees
            east_square  = lon_fine[i] + buffer_degrees
            west_square  = lon_fine[i] - buffer_degrees
            square=[north_square,south_square,east_square,west_square]
            
            (square_lat,square_lon,square_var)=geo_select_rectangle(lat,lon,square,
                                                lat,lon,var)
            this_set = [square_var[j] for j in range(0,len(square_var)) if 
                                         great_circle( (square_lat[j],square_lon[j]),
                                                       (lat_fine[i],lon_fine[i])).km
                                         <= averaging_radius ]
            var_fine.append(np.nanmean(this_set))
            var_fine_count.append(notnancount(this_set))
            del this_set
    elif os_type == "2":
        #for rectangle oversampling
        half_lat = averaging_rectangle_lat_dim * 0.5
        half_lon = averaging_rectangle_lon_dim * 0.5
        for i in range(0,fine_points):
            if np.remainder(i,100) == 0: #only update screen occasionally
                clearscreen()
                print "Processing point %i of %i" %(i,fine_points)
                
            north_square = lat_fine[i] + half_lat
            south_square = lat_fine[i] - half_lat
            east_square  = lon_fine[i] + half_lon
            west_square  = lon_fine[i] - half_lon
            square=[north_square,south_square,east_square,west_square]
                
            this_set = geo_select_rectangle(lat,lon,square,var)

            var_fine.append(np.nanmean(this_set))
            var_fine_count.append(notnancount(this_set))
            del this_set

    this_os_data = os_data(lat_fine,lon_fine,var_fine,var_name)
    this_os_data.set_meta('Fine grid dimension',fine_dim)
    
    if os_type == "1":
        this_os_data.set_meta('Oversampling type','circular')
        this_os_data.set_meta('Oversampling dimensions',averaging_radius)
    elif os_type == "2":
        this_os_data.set_meta('Oversampling type','rectangular')
        this_os_data.set_meta('Oversampling dimensions',[averaging_rectangle_lat_dim,averaging_rectangle_lon_dim])

    #this_os_data.set_meta('Spatial extent'

    plot_grid_from_list(lat_fine,lon_fine,var_fine,
                        fine_dim,fine_dim,
                        north_view,south_view,east_view,west_view,
                        title="Oversampled plot",vmin=0.0e16,vmax=3.0e16,
                        lab=var_name,
                        save=False,save_filename="nofilenamespecified")

    plot_grid_from_list(lat_fine,lon_fine,var_fine_count,
                        fine_dim,fine_dim,
                        north_view,south_view,east_view,west_view,
                        title="Oversampled plot",vmin=0,vmax=np.nanmax(var_fine_count),
                        lab="data count at location",
                        save=False,save_filename="nofilenamespecified")

