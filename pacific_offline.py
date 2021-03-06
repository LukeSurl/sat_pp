#!/usr/bin/env python
"""pacific_offline.py - Daily correction factors based on reference region written to file"""

from datetime import date as date
from datetime import datetime as dt
from datetime import timedelta as td
import numpy as np
import cPickle
from os import listdir
from os.path import isfile, join
import h5py
from bpch import bpch

class correction:
    """A date, and its median slant and model numbers"""
    def __init__(self,this_date):
        self.date = this_date
        self.med_obs = None #placeholder
        self.med_geos = None #placeholder


def pacific_ol(area,dates,pacific_sat_dir,geos_dir,save_location,geosfilestr='ts_satellite_omi'):
    
    #==Options==

    #(ref_n,ref_s,ref_e,ref_w) = (40.,0.,180.,160.) #reference area.
    (ref_n,ref_s,ref_e,ref_w) = area #reference area.

    load_pickle = False #update an existing pickle (True) or create a new one (False)
    do_obs = True #work out median obs in reference sector
    do_geos = True #work out median model column in reference sector

    #dates
    first_date = date(2014,1,1)
    last_date = date(2014,12,31)
    (first_date,last_date) = dates

    #get a list of dates
    list_of_dates = []
    dateclock = first_date
    while dateclock <= last_date:
        list_of_dates.append(dateclock)
        dateclock += td(days=1)

    #Load/save loaction
    #pickle_location = "/group_workspaces/jasmin/geoschem/local_users/lsurl/CL_PP/bark_2014b.p"

    #data locations
    #pacific_data_dir = "/group_workspaces/cems2/nceo_generic/nceo_ed/OMHCHO_pacific_bark/"
    #geos_dir = "/group_workspaces/jasmin/geoschem/local_users/lsurl/runs/geosfp_2x25_tropchem_2012-2014/"

    #==Processing==

    #assign dates to files
    if do_obs:
        files_in_obs_folder = [f for f in listdir(pacific_sat_dir) if isfile(join(pacific_sat_dir, f))] #all files in folder
        files_in_obs_folder = [f for f in files_in_obs_folder if f.startswith('OMI-Aura_L2-OMHCHO_')] #only OMI files
        obs_file_dates = []
        for f in files_in_obs_folder:
            this_year = int(f[19:23])
            this_month = int(f[24:26])
            this_day= int(f[26:28])
            obs_file_dates.append(date(this_year,this_month,this_day))
            
    if do_geos:
        geosfilestr  ='ts_satellite_omi'
        print "Looking for ND51 files that begin with '%s'." %geosfilestr
        files_in_geos_folder = [f for f in listdir(geos_dir) if isfile(join(geos_dir, f))] #all files in folder
        print len(files_in_geos_folder)
        #print files_in_geos_folder
        files_in_geos_folder = [f for f in files_in_geos_folder if f.startswith(geosfilestr)] #only ND51 files
        print len(files_in_geos_folder)
        geos_file_dates = []
        for f in files_in_geos_folder:
            this_year = int(f[len(geosfilestr)+1:len(geosfilestr)+5])
            this_month = int(f[len(geosfilestr)+5:len(geosfilestr)+7])
            this_day= int(f[len(geosfilestr)+7:len(geosfilestr)+9])
            geos_file_dates.append(date(this_year,this_month,this_day))


    all_corrections = []
    #Cycle through dates
    for dateclock in list_of_dates:

        today_correction = correction(dateclock)
        
        #observation median
        if do_obs:
            today_files = [files_in_obs_folder[i] for i in range(0,len(files_in_obs_folder)) if obs_file_dates[i] == dateclock]
            today_SC = np.array([])
            for f in today_files:
                try:
                    f_HDF = h5py.File(join(pacific_sat_dir, f))
                except IOError: #dud files
                    #print "Cannot open %s" %f
                    continue
                
                #shorthands
                df = f_HDF['HDFEOS']['SWATHS']['OMI Total Column Amount HCHO']['Data Fields']
                gf = f_HDF['HDFEOS']['SWATHS']['OMI Total Column Amount HCHO']['Geolocation Fields']
                
                #get data
                this_lat = np.array(gf['Latitude']).flatten()
                this_lon = np.array(gf['Longitude']).flatten()
                this_SC =  np.array(df['ColumnAmount']).flatten()
                this_XTQF = np.array(gf['XtrackQualityFlags']).flatten()
                this_MDQF = np.array(df['MainDataQualityFlag']).flatten()

                #geo and quality filter
                keep = np.logical_and.reduce((this_lat > ref_s,
                       this_lat < ref_n,
                       this_lon > ref_w,
                       this_lon < ref_e,
                       this_XTQF == 0,
                       this_MDQF == 0,
                       this_SC != np.nan))
                this_SC = this_SC[keep]
                #print len(this_SC)
                today_SC = np.append(today_SC,this_SC)
                
            today_correction.med_obs = np.median(today_SC)
            print "Median observed slant    column for %s: %g" %(dateclock.strftime("%Y-%m-%d"),today_correction.med_obs)
        
        #geos median
        if do_geos:
            today_files = [files_in_geos_folder[i] for i in range(0,len(files_in_geos_folder)) if geos_file_dates[i] == dateclock]
            if len(today_files) != 1:
                today_correction.med_geos = np.nan
                print "Median modelled vertical column for %s: nan" %dateclock.strftime("%Y-%m-%d")
            else:        
                f = bpch(join(geos_dir, today_files[0]))
                
                this_hcho = np.array(f.variables["IJ-AVG-$_CH2O"]) * 1e-9 #convert to v/v
                box_height = np.array(f.variables['BXHGHT-$_BXHEIGHT']) * 100 #box height in cm
                air_den =    np.array(f.variables['TIME-SER_AIRDEN'])  # air density molec.cm-3
                air_amount = np.multiply(box_height,air_den) #air per grid box molec.cm-2
                this_data_col = np.sum(np.multiply(this_hcho,air_amount),axis=1) #column molec.cm-2
                
                today_correction.med_geos = np.median(this_data_col)
                print "Median modelled vertical column for %s: %g" %(dateclock.strftime("%Y-%m-%d"),today_correction.med_geos)
                                  
        all_corrections.append(today_correction)     

    print "saving..."
    cPickle.dump(all_corrections,open(save_location,"wb"))
