#!/usr/bin/env python
"""Utility to "oversample" satellite observations"""

from geopy.distance import great_circle
from CL_satpp import *
import numpy as np
import copy
import matplotlib.path as mplPath

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
    

def oversample(ida):
    """User interface for the oversampler"""
    
    #ida = copy.deepcopy(use_ida) #avoid memory issue.
    
    while True: #
        os_menu_title = "Select dataset for oversampling:"
        [os_menu_text,os_answers_dict] = ida.make_menu()
        os_menu_choice = basic_menu(os_menu_title,
                                    os_menu_text,
                                    quit_option=True)
        if os_menu_choice == "Z":
            break
        os_key =      os_answers_dict[os_menu_choice]    

        view_box_valid = False
        while view_box_valid == False:                          
            print "Please enter the bounds of the area you wish to inspect"
            north_view = float(raw_input("NORTH-->"))
            south_view = float(raw_input("SOUTH-->"))
            east_view  = float(raw_input("EAST -->"))
            west_view  = float(raw_input("WEST -->"))
            
            view_box = box(north_view,south_view,east_view,west_view)
            view_box_valid = view_box.valid()
                                                                      
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
                 ,["3","As rectangles with defined dimensions (in degrees) NUMPY version"]
                 ],
                 quit_option=False)
    
        if os_type == "1":    
            print "What radius, in km, will the %s observations be spread over?" %ida.data[os_key].description
            print "(24 km recommended)"
            averaging_radius = float(raw_input("-->"))
            #at the edges, we will need to sample points outside these bounds
            #Thus we define sample bounds which are slightly larger

            #we'll assume, for safety, that 1km = 0.01 degrees; a slight overestimation
            #(even at the equator), but this is useful.

            buffer_degrees = 0.01*averaging_radius
            #north_sample = north_view + buffer_degrees
            #south_sample = south_view - buffer_degrees
            #east_sample  = east_view  + buffer_degrees
            #west_sample  = west_view  - buffer_degrees
            sample_box= box(
                       view_box.n+buffer_degrees,
                       view_box.s-buffer_degrees,
                       view_box.e+buffer_degrees,
                       view_box.w-buffer_degrees)
                       
        elif os_type in ["2","3"]:
            print "What is the LATITUDE dimension of these rectangles in degrees?"
            print "1 degree latitude = 111 km; 1 km = 0.0090 degrees latitude"
            averaging_rectangle_lat_dim = float(raw_input("-->"))
                        
            print "What is that LONGITUDE dimension of these rectangles in degrees?"
            print "Length of a degree longitude depends on latitude"
            midpoint = (north_view + south_view)/2
            km_per_degree_at_midpoint = 111.*np.cos(np.deg2rad(midpoint))
            degree_per_km_at_midpoint = 1./km_per_degree_at_midpoint
            print "At latitude of %.0f degrees, 1 degree longitude = %.0f km; 1 km = %.4f degrees" \
                   %(midpoint,km_per_degree_at_midpoint,degree_per_km_at_midpoint)
            averaging_rectangle_lon_dim = float(raw_input("-->"))
            
            sample_box = box(
                            view_box.n + 0.5*averaging_rectangle_lat_dim,
                            view_box.s - 0.5*averaging_rectangle_lat_dim,
                            view_box.e + 0.5*averaging_rectangle_lon_dim,
                            view_box.w - 0.5*averaging_rectangle_lon_dim) 

        #now use geo_select_rectangle to restrict our datasets to just those in the sampling box
        #this saves a lot of computational time!

        #sample_box = [north_sample,south_sample,east_sample,west_sample]
        #print sample_box.valid()
        #print len(ida.data[os_key].val)
        #ida = geo_select_rectangle_i(sample_box,ida)
        ida = geo_select_rectangle_i(sample_box,ida)
        #print ida
        #print len(ida.lat)
        #print len(ida.data[os_key].val)
        #_ = raw_input("halt-->")

        os_var =      ida.data[os_key].val
        os_var_name = ida.data[os_key].description
        
        #reduce ida to just the variable we want
        
        #reduced_ida = d_all(ida.lat,ida.lon,time=ida.time)
        #reduced_ida.data[os_key] = ida.data[os_key]
                
        #ida = copy.deepcopy(reduced_ida)
        #del reduced_ida

        #Define the lat and lon co-ordinates of the ultra-fine grid.

        lat_fine = []
        lon_fine = []
        
        lat_step = view_box.s
        while lat_step <= view_box.n:
            lon_step = view_box.w
            while lon_step <= view_box.e:
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
                fresh_ida = copy.deepcopy(ida)
                if np.remainder(i,100) == 0: #only update screen occasionally
                    clearscreen()
                    print "Processing point %i of %i" %(i,fine_points)
                #for each point, first cut down a square area around the central point before doing the
                #the circle
                square = box(
                            lat_fine[i] + buffer_degrees,
                            lat_fine[i] - buffer_degrees,
                            lon_fine[i] + buffer_degrees,
                            lon_fine[i] - buffer_degrees)
                
                square_ida = geo_select_rectangle_i(square,fresh_ida)
                this_set = [square_ida.data[os_key].val[j] 
                            for j in range(0,len(square_ida.data[os_key].val))
                            if great_circle( (square_ida.lat[j],square_ida.lon[j]),
                                             (lat_fine[i],lon_fine[i])).km
                                             <= averaging_radius ]
                var_fine.append(np.nanmean(this_set))
                var_fine_count.append(notnancount(this_set))
                del this_set
                #del square_ida
            this_os_data = os_data(lat_fine,lon_fine,var_fine,os_var_name)
            this_os_data.set_meta('Fine grid dimension',fine_dim)
                            
        elif os_type == "2":
            #for rectangle oversampling
            half_lat = averaging_rectangle_lat_dim * 0.5
            half_lon = averaging_rectangle_lon_dim * 0.5
            for i in range(0,fine_points):
                fresh_ida = copy.deepcopy(ida)
                if np.remainder(i,100) == 0: #only update screen occasionally
                    clearscreen()
                    print "Processing point %i of %i" %(i,fine_points)
                square = box(
                             lat_fine[i] + half_lat,
                             lat_fine[i] - half_lat,
                             lon_fine[i] + half_lon,
                             lon_fine[i] - half_lon)
                    
                square_ida = geo_select_rectangle_i(square,fresh_ida)
                #print len(ida.data[os_key].val)
                #print len(square_ida.data[os_key].val)

                var_fine.append(np.nanmean(square_ida.data[os_key].val))
                var_fine_count.append(notnancount(square_ida.data[os_key].val))
                #del this_set

            this_os_data = os_data(lat_fine,lon_fine,var_fine,os_var_name)
            this_os_data.set_meta('Fine grid dimension',fine_dim)
        
        elif os_type == "3":
            #for rectangle oversampling, speed version
            half_lat = averaging_rectangle_lat_dim * 0.5
            half_lon = averaging_rectangle_lon_dim * 0.5
            ida_lat_np = np.array(ida.lat)
            ida_lon_np = np.array(ida.lon)
            ida_val_np = np.array(ida.data[os_key].val)
            
            for i in range(0,fine_points):
                #fresh_ida = copy.deepcopy(ida)
                print "Processing point %i of %i" %(i,fine_points)
                
                lats_good = abs(ida_lat_np - lat_fine[i]) <= half_lat
                lons_good = abs(ida_lon_np - lon_fine[i]) <= half_lat
                
                latlongood = lats_good*lons_good
                
                vars_to_av = ida_val_np[latlongood]
                
                var_fine.append(np.nanmean(vars_to_av))
        
        if os_type in ["1","2"]:        
            this_os_data = os_data(lat_fine,lon_fine,var_fine,os_var_name)
            this_os_data.set_meta('Fine grid dimension',fine_dim)
        else:
            this_os_data = (lat_fine,lon_fine,var_fine)
            
        if os_type == "1":
            this_os_data.set_meta('Oversampling type','circular')
            this_os_data.set_meta('Oversampling dimensions',averaging_radius)
        elif os_type == "2":
            this_os_data.set_meta('Oversampling type','rectangular')
            this_os_data.set_meta('Oversampling dimensions',[averaging_rectangle_lat_dim,averaging_rectangle_lon_dim])

        #this_os_data.set_meta('Spatial extent'

        plot_grid_from_list(lat_fine,lon_fine,var_fine,
                            fine_dim,fine_dim,
                            view_box,
                            title="Oversampled plot",vmin=0.0e16,vmax=2.0e16,
                            lab=os_var_name,
                            save=True,save_filename="/home/users/lsurl/CL/%fN-%fN--%fE-%fE--oversample.png"%(south_view,north_view,west_view,east_view))

        #plot_grid_from_list(lat_fine,lon_fine,var_fine_count,
        #                    fine_dim,fine_dim,
        #                    view_box,
        #                    title="Oversampled plot",vmin=0,vmax=np.nanmax(var_fine_count),
        #                    lab="data count at location",
        #                    save=False,save_filename="nofilenamespecified")
        return(this_os_data)


def quad_as_path(lat_corners,lon_corners):
    """Returns a path object for a quadrilateral from its corners"""
    return mplPath.Path(np.array([[lat_corners[0],lon_corners[0]],
                                  [lat_corners[1],lon_corners[1]],
                                  [lat_corners[2],lon_corners[2]],
                                  [lat_corners[3],lon_corners[3]]]))
    

def area_poly(x_cor,y_cor):
    """Returns the area of a simple polygon, given the corners"""
    
    num_corners = len(x_cor)
    #close the polygon
    x_cor = np.append(x_cor,x_cor[0])
    y_cor = np.append(y_cor,y_cor[0])
    area = 0.0
    for i in range(0,num_corners):
        area += 0.5*(x_cor[i]*y_cor[i+1]-x_cor[i+1]*y_cor[i])
    
    return(area)

def overlay(ida):
    """A more sophisticated oversampler"""

    ov_menu_title = "Select dataset for overlaying:"
    [ov_menu_text,ov_answers_dict] = ida.make_menu()
    ov_menu_choice = basic_menu(ov_menu_title,
                                ov_menu_text,
                                quit_option=False)                                                           
    ov_key =      ov_answers_dict[ov_menu_choice]
    ov_var_name = ida.data[ov_key].description      
    #Define the mapping box
    view_box_valid = False
    while view_box_valid == False:                          
        print "Please enter the bounds of the area you wish to inspect"
        north_view = float(raw_input("NORTH-->"))
        south_view = float(raw_input("SOUTH-->"))
        east_view  = float(raw_input("EAST -->"))
        west_view  = float(raw_input("WEST -->"))
        
        view_box = box(north_view,south_view,east_view,west_view)
        view_box_valid = view_box.valid()
    
    #define the ultra-fine grid
    lat_fine = []
    lon_fine = []
    fine_grid = []   
    fine_dim = 0.02                
    lat_step = view_box.s
    while lat_step <= view_box.n:
        lon_step = view_box.w
        while lon_step <= view_box.e:
            lat_fine.append(lat_step)
            lon_fine.append(lon_step)
            fine_grid.append([lat_step,lon_step])
            lon_step += fine_dim
        lat_step += fine_dim
    
    fine_grid=np.array(fine_grid)
    
    #Bring up ALL relevant data    
    all_lat_corners = np.array(ida.data["lat_corners"].val)
    all_lon_corners = np.array(ida.data["lon_corners"].val)
    all_var          = np.array(ida.data[ov_key].val)
    #all_dvar         = np.array(ida.data["sat_DVC"].val)
    
    
    #work out which tiles have at least one corner in the main box.  
    num_obs = len(all_var)
    
    print("Working out relevant observations for this area") 
    keep_region = []
    for i in range(0,num_obs):
        if view_box.within(float(all_lat_corners[i][0]),float(all_lon_corners[i][0])) or\
           view_box.within(float(all_lat_corners[i][1]),float(all_lon_corners[i][1])) or\
           view_box.within(float(all_lat_corners[i][2]),float(all_lon_corners[i][2])) or\
           view_box.within(float(all_lat_corners[i][3]),float(all_lon_corners[i][3])):
            keep_region.append(True)
        else:
            keep_region.append(False)     
    keep_region = np.array(keep_region)            
    region_lat_corners = all_lat_corners[keep_region]
    region_lon_corners = all_lon_corners[keep_region]
    region_var          = all_var[keep_region]
    #region_dvc         = all_dvc[keep]
    
    num_region=len(region_var)
    print "Number of relevant observations: %i out of %i total" %(num_region,num_obs)
    #for each observation, create a mplPath object defining its spatial extent
    print "Calculating paths..."   
    obs_paths = [quad_as_path(region_lat_corners[i],region_lon_corners[i]) for i in range(0,num_region) ]
    #work out the relative weight per fine grid cell (reciprocal of area) of each
    print "Calculating areas..." 
    weights = np.array([1./area_poly(region_lat_corners[i],region_lon_corners[i]) for i in range(0,num_region) ])
    
    #Note: For now, area is computed as if 1 degree lat = 1 degree lon
    
    #For each observation, add to a cumulative total of the weights of observations, and obs*weights
    cum_wv   = np.zeros_like(np.array(lat_fine))
    cum_w   =  np.zeros_like(np.array(lat_fine))    
    for i in range(0,num_region):
        print "Obs. %i of %i" %(i+1,num_region)
        pix_in = 1.* np.array(obs_paths[i].contains_points(fine_grid))
        #print list(pix_in)
        cum_w += weights[i]*pix_in
        cum_wv += region_var[i]*weights[i]*pix_in
        
    
    var_fine = list(np.divide(cum_wv,cum_w))    
    
    plot_grid_from_list(lat_fine,lon_fine,var_fine,
                        fine_dim,fine_dim,
                        view_box,
                        title="Overlayed plot",vmin=0.0e16,vmax=2.0e16,
                        lab=ov_var_name,
                        save=True,save_filename="/home/users/lsurl/CL/%fN-%fN--%fE-%fE--overlay.png"%(south_view,north_view,west_view,east_view))
    
    return(lat_fine,lon_fine,var_fine)  

    
    
    
    
    

