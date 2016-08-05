#!/usr/bin/env python

"""generate_pickles.py - Loads and processes filtered data from text files
and saves pickes"""

import os
import math
from datetime import date
from datetime import datetime as dt
from datetime import timedelta as td
import numpy as np
from CL_satpp import *

def is_in_range(startdate,enddate,submitdate):
    """True/False: is submitdate between startdate and enddate"""
    
    if submitdate >= startdate and submitdate <= enddate:
       return True
    else:
       return False
       
def last_day_of_month(indate):
    """Given a date, returns the last of that month""" 
    if indate.month == 12:
        return indate.replace(day=31)
    return indate.replace(month=indate.month+1, day=1) - td(days=1)
    
def find_nearest(array,value):
    """Returns the nearest entry to value in array""" 
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def load_filtered_main(data_folder,startdate,enddate):
    """Given a folder and two dates, collects data for times
    between those dates into 1D lists""" 
    
    list_of_files = [] #initialise empty list.

    #Rather crudely, we use the timestamp encoded into the filenames.
    for root, dirs, files in os.walk(data_folder):
        for filename in files:
            filedate = dt(int(filename[19:23]),
                            int(filename[24:26]),
                            int(filename[26:28]))
            if is_in_range(startdate,enddate,filedate):
                 list_of_files.append(data_folder
                                      + filename[19:23] 
                                      +"/"+filename)

    #list of files will now be full paths to all files to be processed.
    #print these to screen
    print "Will load data from %i files" % len(list_of_files)
    print "First: %s" %list_of_files[0]
    print "Last : %s" %list_of_files[len(list_of_files)-1]

    #filtered_dataset will be a list containing an obs object for each point
    filtered_dataset = []
    
    #now build arrays for all of the data in range.    
    ULN     = [] #unique line number
    lat     = [] #latitude
    lon     = [] #longitude
    time    = [] #UNIX time of observation 
    geos_VC = [] #modelled vertical column 
    sat_SC  = [] #observed slant column
    sat_DSC = [] #error in observed slant column
    AMF     = [] #AMF

    for filepath in list_of_files:
        f = open(filepath,'r')
        for line in f:
            line = line.strip() #get rid of the newline character
            columns = line.split(",") #split by columns, comma separated
            if float(columns[5]) <= 0: #Kick out datapoint if AMF is <=0
                continue
            if math.isnan(float(columns[4])): #kick out datapoint if GEOS_VC is NaN
                continue
                
            this_lat = float(columns[2])
            this_lon = float(columns[3])
            this_time= float(columns[1])
            
                
            ULN.append(int(columns[0]))    
            lat.append(float(columns[2]))
            lon.append(float(columns[3]))
            time.append(float(columns[1]))
            geos_VC.append(float(columns[4]))
            sat_SC.append(float(columns[6]))
            sat_DSC.append(float(columns[7]))
            AMF.append(float(columns[5]))
    
    return(ULN,lat,lon,time,geos_VC,sat_SC,sat_DSC,AMF)
    
def get_corrections(corrections_folder,startdate,enddate) :
    """Loads the pacific corrections within a date range from a given folder"""
    
    #init empty lists
    #these lists will be filled up line-by-line
    corr_year   = [] #the year of this correction
    corr_month  = [] #the month of this correction
    corr_lat    = [] #the latitude of this correction
    corrections = [] #For each line, 60 corrections, one for each detector
    Vcorrs      = [] #An averaged-over-all 60 correction factor
    
    startyear = startdate.year
    startmonth= startdate.month
    endyear   = enddate.year
    endmonth  = enddate.month
    
    for this_year in range(startyear,endyear+1):
        for root, dirs, files in os.walk(corrections_folder
                                         + "/" + str(this_year)):
             
            for filename in files:
                #Crudely, we get the year & month for a file from its filename
                file_year  = int(filename[0:4])
                file_month = int(filename[5:7])
              
                if is_in_range(dt(startyear,startmonth,1),
                                last_day_of_month(dt(endyear,endmonth,1)),
                                dt(file_year,file_month,2)) :

                    f = open("%s%04d/%04d_%02d_corr.csv" 
                              %(corrections_folder,
                                file_year,file_year,file_month),'r')
                                
                    #run line-by-line            
                    for line in f:
                    
                        #record year, month
                        corr_year.append(file_year)
                        corr_month.append(file_month)
                        
                        #get rid of the newline character
                        line = line.strip()
                        
                        #split by columns, comma separated
                        columns = line.split(",")
                        
                        #latitude is first column
                        corr_lat.append(float(columns[0]))
                        
                        this_lat_corrs = [] #temp list for 60 corrections
                        #fill up this list with each pixel's correction
                        for pix in range(1,61):
                            this_lat_corrs.append(float(columns[pix]))
                        #add 60-member list (nested) to corrections
                        corrections.append(this_lat_corrs)
                        del this_lat_corrs #delete temp list
                        
                        #average correction to vertical columns 
                        Vcorrs.append(float(columns[61]))
                    
    return(corr_year,corr_month,corr_lat,corrections,Vcorrs)

def AMF_and_correction(ULN,lat,time,sat_SC,sat_DSC,AMF,
                       corr_year,corr_month,corr_lat,corrections):
    """Applies the AMF and pacific correction to observed slant columns. 
    Uses detector-specific corrections"""
    
    #get a value for the number of data points
    data_points = len(lat)
    
    #how many different corrections lines are there?
    corrections_lines = len(corrections)
    
    #initialise lists for vertical columns
    sat_VC = []
    sat_DVC = []
    
    #The detector (0-59) can be identified from the ULN 
    pixel = np.mod(ULN,60)
    
    #Get a list of all the different lats at which
    #corrections are calculated
    lat_levels = np.array(list(set(corr_lat)))
    
    #loop over every data point
    for i in range(0,data_points):
    
        #identify which latitude is it is correct to
        #take corrections from
        lat_to_find = find_nearest(lat_levels,lat[i])
        
        #work out the year and month of this data point
        this_date = date.fromtimestamp(time[i])
        this_year = this_date.year
        this_month = this_date.month

        #identify the correction that applies to this
        #detector (pixel[i]), this year, this month and
        #this latitude
        this_corr = [corrections[j][pixel[i]] 
                     for j in range(0,corrections_lines) 
                     if corr_year[j] == this_year 
                     and corr_month[j] == this_month 
                     and corr_lat[j] == lat_to_find]

        #now, we SHOULD have exactly 1 entry in this_corr.
        #If no correction was identified we can call this
        #a NaN. If more than 1 correction is identified, that's
        #an error worth terminating the script for.
        #Save the desired correction in this_corr_f
        if len(this_corr) == 1:
            this_corr_f = float(this_corr[0])
        elif len(this_corr) > 1:
            raise ValueError("More than one correction number" 
                             "found in AMF_and_correction."
                             "This shouldn't happen")
        else:
            this_corr_f = None
        del this_corr #delete variable before overwrite on next iteration
        
        #Assign NaNs if no valid correction was found (or it
        #was a NaN). Otherwise calculate the VC and DVC.
        if this_corr_f == None or np.isnan(this_corr_f):
            sat_VC.append(float("nan"))
            sat_DVC.append(float("nan"))
        else: 
            sat_VC.append((sat_SC[i] - this_corr_f)/AMF[i])
            sat_DVC.append(abs(sat_DSC[i]/AMF[i]))
        
        del this_corr_f #delete variable before overwrite on next iteration
    
    return(sat_VC,sat_DVC)

def BASIC_AMF_and_correction(lat,time,sat_SC,sat_DSC,
                             AMF,corr_year,corr_month,corr_lat,
                             Vcorrections):
    """Applies the AMF and a pacific correction to observed slant columns
    Unlike AMF_and_correction, these corrections are not pixel specific,
    rather the VCs are calcualted and an average correction to that (for all
    pixels) is applied"""
    
    #get a value for the number of data points
    data_points = len(lat)
    
    #how many different corrections lines are there?
    corrections_lines = len(Vcorrections)
    
    #initialise lists for vertical columns 
    sat_VC = []
    sat_DVC = []
    
    #Get a list of all the different lats at which
    #corrections are calculated
    lat_levels = np.array(list(set(corr_lat)))
    
    #loop over every data point    
    for i in range(0,data_points):
    
        #identify which latitude is it is correct to
        #take corrections from
        lat_to_find = float(find_nearest(lat_levels,lat[i]))
        
        #work out the year and month of this data point        
        this_date = date.fromtimestamp(time[i])
        this_year = this_date.year
        this_month = this_date.month

        #Assign True/False whether each criteria
        #(year, month, lat) is matched 
        year_matched     = [k == this_year for k in corr_year ]
        month_matched    = [k == this_month for k in corr_month ]
        lat_matched      = [k == lat_to_find for k in corr_lat ]
        
        #if all three match at a given posisition,
        #then that's the correction needed
        this_corr = [Vcorrections[j] for j in range(0,corrections_lines)
                                     if year_matched[j]
                                     and month_matched[j]
                                     and lat_matched[j]]
        
        #delete variables before overwrite on next iteration
        del year_matched, month_matched, lat_matched
                                      
        #now, we SHOULD have exactly 1 entry in this_corr.
        #If no correction was identified we can call this
        #a NaN. If more than 1 correction is identified, that's
        #an error worth terminating the script for.
        #Save the desired correction in this_corr_f
        if len(this_corr) == 1: #exactly one match, good
            this_corr_f = float(this_corr[0])
        elif len(this_corr) > 1: #more than one match, v bad
            raise ValueError("More than one correction number" 
                             "found in AMF_and_correction."
                             "This shouldn't happen")
        else: #0 matches, bad but not fatal.
            print "Alert: no correction value (not even a NaN)"\
                  "found for year=%i, month = %i, lat = %f"\
                   %(this_year,this_month,lat[i])
            this_corr_f = None
        if this_corr_f == None or np.isnan(this_corr_f):
            sat_VC.append(float("nan"))
            sat_DVC.append(float("nan"))
        else :
            this_sat_VC = sat_SC[i]/AMF[i] - this_corr_f    
            sat_VC.append(this_sat_VC)
            sat_DVC.append(abs(sat_DSC[i]/AMF[i]))

    return(sat_VC,sat_DVC)
    
def load_and_correct(startdate,enddate,
                     filtered_folder="NOTAVALIDPATH",
                     corrections_folder="NOTAVALIDPATH",
                     corr_type="ask"):
    """Loads filtered data, pacific correction, and AMF and applies these"""
    
    
    #check if provided path for filtered data is an actual folder
    #The default value of filtered_folder will return False
    #This is far from foolproof, but it's a
    #useful quick catch of invalid input
    path_valid = os.path.isdir(filtered_folder)
    
    #if we haven't got a valid path, ask for one.
    while path_valid == False:
        print "Need a valid path"
        print "Enter the folder in which filtered data has been saved."
        print "Do not include year subdirectories, these will be navigated automatically."
        filtered_folder = raw_input("-->")
        #Check that is a valid path.
        path_valid = os.path.isdir(filtered_folder)
        if path_valid == False:
            print "%s is not a valid path" %filtered_folder
    
    del path_valid
    
    #load the data        
    print "loading filtered data..."
    (ULN,lat,lon,time,geos_VC,sat_SC,sat_DSC,AMF) =\
        load_filtered_main(filtered_folder,startdate,enddate)

    #check if provided path for filtered data is valid
    #Same process & caveats as before...
    
    path_valid = os.path.isdir(corrections_folder)
    #if we haven't got a valid path, ask for one.
    while path_valid == False:
        print "Need a valid path"
        print "Enter the folder in which correctoins have been saved."
        print "Do not include year subdirectories, these will be navigated automatically."
        filtered_folder = raw_input("-->")
        #Check that is a valid path.
        path_valid = os.path.isdir(corrections_folder)
        if path_valid == False:
            print "%s is not a valid path" %filtered_folder
    
    del path_valid
    
    #load the corrections    
    print "loading corrections"
    (corr_year,corr_mon,corr_lat,corrections,Vcorrections) = get_corrections(corrections_folder,startdate,enddate)
    
    #if we need to prompt the user for the type of correction used
    while corr_type not in ["basic","advanced"]:
        print "Pacfic corrections"
        print "[1] Calculate VCs, apply same pacific correction to each observation at same latitude"
        print "[2] Apply detector-specific corrections to SCs, then convert to VCs"
        corr_option = raw_input("-->")
        if corr_option == "1":
            corr_type = "basic"
        elif corr_option == "2":
            corr_type = "advanced"
        del corr_option
    
    print "Applying AMFs and corrections"
    if corr_type == "advanced":
        (sat_VC,sat_DVC) = \
            AMF_and_correction(ULN,lat,time,
                               sat_SC,sat_DSC,AMF,
                               corr_year,corr_mon,corr_lat,
                               corrections)
    elif corr_type == "basic":
        (sat_VC,sat_DVC) = \
            BASIC_AMF_and_correction(lat,time,
                                     sat_SC,sat_DSC,AMF,
                                     corr_year,corr_mon,corr_lat,
                                     Vcorrections)
    
    #Build ida
    ida = build_ida(ULN,lat,lon,time,geos_VC,sat_VC,sat_DVC,AMF)
    
    #add pacific correction type to metadata where applicable
    ida.data['sat_VC'].add_meta( 'Pacific correction type',corr_type)
    ida.data['sat_DVC'].add_meta('Pacific correction type',corr_type)
                                   
    return(ida)
    
    

