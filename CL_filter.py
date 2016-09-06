#!/usr/bin/env python
"""CL_filter.py -- Filtering of dataset"""

from datetime import datetime as dt
from datetime import timedelta as td
import h5py
import os
from os import listdir
from os.path import isfile, join
import pickle
import numpy as np

from CL_satpp import box,d,d_all,add_slash

#Tools to filter satellite observational data.

#This script will look up output from the AMF calculator
#and also inspect the HDF5 files directly.

def OMI_from_HDF(HDF_file):
    """Returns various sets of data from an OMI HCHO file as lists"""
     
    #shorthands
    df = HDF_file['HDFEOS']['SWATHS']['OMI Total Column Amount HCHO']['Data Fields']
    gf = HDF_file['HDFEOS']['SWATHS']['OMI Total Column Amount HCHO']['Geolocation Fields']
    
    (n_scans,n_tracks) = gf['Latitude'].shape #shape of the data
    
    #get simple location data 
    lat = list(np.array(gf['Latitude']).flatten())
    lon = list(np.array(gf['Longitude']).flatten())
    
    #time is a little fiddly
    time = []
    for scan in range(0,n_scans):
        [year,month,day,hour,minute,second] = list(gf['TimeUTC'][scan])
        for track in range(0,n_tracks):
            time.append(dt(year,month,day,hour,minute,second))
    
    #scan and row/track identifiers
    scan_n_list = []
    track_n_list = []
    for scan in range(0,n_scans):
        for track in range(0,n_tracks):
            scan_n_list.append(scan)
            track_n_list.append(track)
    
    scan_n  = d(scan_n_list,"Scan",description="Scan number")
    track_n = d(track_n_list,"Row",description="Row number")
            
    #angles of observation    
    SAA = d(
            list(np.array(gf['SolarAzimuthAngle']).flatten()),
            "SAA",description="Solar Azimuth Angle")
    SAA.unit = u"\u00b0" #degree sign        
    SZA = d(
            list(np.array(gf['SolarZenithAngle']).flatten()),
            "SZA",description="Solar Zenith Angle")
    SZA.unit = u"\u00b0" #degree sign        
    VAA = d(
            list(np.array(gf['ViewingAzimuthAngle']).flatten()),
            "VAA",description="Viewing Azimuth Angle")
    VAA.unit = u"\u00b0" #degree sign
    VZA = d(
            list(np.array(gf['ViewingZenithAngle']).flatten()),
            "VZA",description="Viewing Zenith Angle")
    VZA.unit = u"\u00b0" #degree sign
    
    
    #HCHO measurements
        
    #NASA files, unhelpfully, do not report the slant columns.
    #Instead, they report the Vertical Column worked out from their AMF
    #However, they also report their AMF, so we can easily undo the AMF
    #and get the slant columns
    #Also treat uncertainty the same wat
    CA  = np.array(df['ColumnAmount']).flatten()
    DCA= np.array(df['ColumnUncertainty']).flatten()
    NASA_AMF = np.array(df['AirMassFactor']).flatten()
    
    SC = d(
           list(np.multiply(CA,NASA_AMF)),
           "Slant Column",description="Slant Column Amount")
    SC.unit = "molec/cm2"
    DSC= d(
           list(np.multiply(DCA,NASA_AMF)),
           "Slant Column Uncertainty",description="Slant Column Uncertainty")
    DSC.unit = "molec/cm2"
    
    #it may also be useful to get the final VC here, even though we'll
    #calculate our own
    
    NASA_VC = d(
                list(np.array(df['ReferenceSectorCorrectedVerticalColumn'])),
                "Vertical Column (NASA)",
                description="Vertical Column (as calulated in OMHCHO)")
    NASA_VC.unit = "molec/cm2"
    
    #Cloud info
    #Cloud infomation is taken from OMCLDO2 product
    cloud_fraction = d(
                       list(np.array(df['AMFCloudFraction']).flatten()),
                       "Cloud fraction",
                       description="Cloud fraction from OMCLDO2")
    cloud_pressure = d(
                       list(np.array(df['AMFCloudPressure']).flatten()),
                       "Cloud top pressure",
                       description="Cloud-top pressure from OMCLDO2")
    cloud_pressure.unit = "hPa"
                       
    
    #Quality flags
    main_data_quality_flag = d(
                                list(np.array(df['MainDataQualityFlag']).flatten()),
                                "Main data quality flag")
    Xtrack_flag = d(
                    list(np.array(gf['XtrackQualityFlags']).flatten()),
                    "Xtrack quality flag")
    
    #now construct a d_all object
    
    OMI_all = d_all(lat,lon,time)
    OMI_all.add_data([scan_n,track_n,
                    SAA,SZA,VAA,VZA,
                    SC,DSC,NASA_VC,
                    cloud_fraction,cloud_pressure,
                    main_data_quality_flag,Xtrack_flag])    
    return(OMI_all)
        
def merge_OMI(all_OMIs):
    """Merges multiple OMI_all files together"""
    for OMI in all_OMIs:
        if OMI == all_OMIs[0]:
            #if first set
            merged = OMI
        else:
            merged.lat.extend(OMI.lat)
            merged.lon.extend(OMI.lon)
            merged.time.extend(OMI.time)
            for dataset in merged.data:
                merged.data[dataset].val.extend(OMI.data[dataset].val)
    return(merged)
    
def write_for_AMF(OMI_all,write_directory,write_name):
    """writes a text file for the AMF calculator to read"""
    
    if not(write_directory.endswith("/")):
        write_directory = write_directory + "/"
    
    len_all = len(OMI_all.lat)
    outfile = open(write_directory + write_name,'w')
    delimiter=","
    for line in range(0,len_all):
        outfile.write("%i,%i,%i,%f,%f,%g,%g,%g,%g\n" \
                        %(line,
                          OMI_all.data["Scan"].val[line],
                          OMI_all.data["Row"].val[line],
                          OMI_all.lat[line],
                          OMI_all.lon[line],
                          OMI_all.data["SZA"].val[line],
                          OMI_all.data["VZA"].val[line],
                          OMI_all.data["Cloud fraction"].val[line],
                          OMI_all.data["Cloud top pressure"].val[line])
                     )
    outfile.close()

def check_HDFs_in_directory(directory):
    """For a given directory, returns list of good HDF files and their dates"""
    
    #print directory
    directory = add_slash(directory)
    #get a list of all files.
    files = [directory+f for f in listdir(directory) if isfile(join(directory, f))] 
    
    #keep only those ending in .he5
    files = [f for f in files if f.endswith('.he5') ]
    
    #if the file is less than 10K in size, it contains no useful data    
    files = [f for f in files if os.stat(f).st_size > 10240]
    
    #now get the dates
    dates_of_files = []
    for f in files:
        try:
            f_HDF = h5py.File(f)
        except IOError:
            #sometimes some files can't be read (corrupted?)
            print "%s cannot be read" %f
            files.remove(f)
            continue
        granule_date = f_HDF['HDFEOS']['SWATHS']['OMI Total Column Amount HCHO']['Geolocation Fields']['TimeUTC'][0]
        this_date = dt(granule_date[0],granule_date[1],granule_date[2],0,0,0)
        dates_of_files.append(this_date)
        
    return(files,dates_of_files)
    
def assemble_sat_data(sat_directory,start_date,end_date,OMI_save,AMF_txt_save):
    """Initial reading of satellite data"""
    (HDF_files,dates_of_files) = check_HDFs_in_directory(sat_directory)
    
    i_to_keep = [i for i in range(0,len(dates_of_files)) 
                 if dates_of_files[i] >= start_date and dates_of_files[i] < end_date ]

    HDF_files = [HDF_files[i] for i in i_to_keep]
    dates_of_files = [dates_of_files[i] for i in i_to_keep]
    
    
    #do process for each day
    have_data = [] #an empty dictionary that will be populated 
    date_clock = start_date
    while date_clock < end_date:
        print "Assembling data for "+ date_clock.strftime("%Y-%m-%d")
        this_day_HDF_files = [HDF_files[i] for i in range(0,len(HDF_files))
                          if dates_of_files[i] == date_clock ]
        if this_day_HDF_files == []: #if no data for this day
            have_data.append([date_clock,False])
        else:
            have_data.append([date_clock, True])
        #print this_day_HDF_files
        these_OMIs = []
        for H in this_day_HDF_files:
            these_OMIs.append(OMI_from_HDF(h5py.File(H)))
        #print these_OMIs
        this_day_OMI = merge_OMI(these_OMIs)
        del these_OMIs
        #save day's OMIs as pickle
        pickle_name = add_slash(OMI_save)+date_clock.strftime("%Y-%m-%d")+"_OMI.p"
        pickle.dump(this_day_OMI,open(pickle_name,"wb"))
        
        #write an AMF input .csv file
        AMF_write_name = date_clock.strftime("%Y-%m-%d")+"_for_AMF.csv"
        write_for_AMF(this_day_OMI,add_slash(AMF_txt_save),AMF_write_name)
       
        date_clock += td(days=1)
    
    return(have_data)

def YYYYMMDD_sub(date,string):
    """substitutes YYYY, MM, DD with appropriate values from a date"""
    YYYY = date.strftime("%Y")
    MM =   date.strftime("%m")
    DD =   date.strftime("%d")
    string = string.replace("YYYY",YYYY)
    string = string.replace("MM",MM)
    string = string.replace("DD",DD)
    return(string)
    
def write_bsub(AMF_txt_save_dir,date_dos):
    """writes a shell script that can be executed to do AMF processing in parallel"""
    #fixed variables
    AMF_rundir = '/home/users/lsurl/amf581CL'
    GEOS_ND51_files = '/group_workspaces/jasmin/geoschem/local_users/lsurl/runs/geosfp_025x03125_tropchem_in_2014/ts_omi.YYYYMMDD.bpch'
    AMF_inputs = add_slash(AMF_txt_save_dir) + "YYYY-MM-DD_for_AMF.csv"
    bsub_stuff = 'bsub -q lotus -n 1 -o /home/users/lsurl/CL/templog.log -J amf-YYYYMMDD -W 2:00'
    output_location = '/home/users/lsurl/amf581CL/output/YYYYMMDD_amfs.hcho'
    outfile = open("batch_jobs.sh",'w')
    outfile.write("#!/bin/bash\n")
    outfile.write("cd "+AMF_rundir+"\n")
    
    for date_do in date_dos:
        if date_do[1] == False:
            continue
        this_date = date_do[0]
        outfile.write(YYYYMMDD_sub(this_date,bsub_stuff)) #bsub commands
        outfile.write(' ')
        outfile.write(add_slash(AMF_rundir)+"amf.run") #path to amf.run
        outfile.write(' ')
        outfile.write(YYYYMMDD_sub(this_date,AMF_inputs)) #observation data
        outfile.write(' ')
        outfile.write(YYYYMMDD_sub(this_date,GEOS_ND51_files)) #GOES ND51
        outfile.write(' ')
        outfile.write(YYYYMMDD_sub(this_date,output_location)) #where to write output
        outfile.write(' ')
        outfile.write(YYYYMMDD_sub(this_date,'YYYY MM DD')) #date stamp
        outfile.write('\n') #new line
    
    outfile.close()
         
