#!/usr/bin/env python
"""CL_preprocess.py -- Python utility for pre-processing"""

from datetime import datetime as dt
from datetime import timedelta as td
import h5py
import os
from os import listdir
from os.path import isfile, join
import pickle
import cPickle
import numpy as np
import copy

from CL_satpp import box,d,d_all,add_slash,clearscreen,basic_menu,geo_select_rectangle_i
from CL_run_management import run_select,new_run
#from CL_omi_matchup import match_BRUG

#Tools to filter satellite observational data.

#This script will look up output from the AMF calculator
#and also inspect the HDF5 files directly.

def preprocessing(base_dir=None,run_name=None):
    """control script for the pre-processing steps"""
    
    if base_dir==None:
        #set to current dir if none supplied
        base_dir=os.path.dirname(os.path.realpath(__file__))
    
    while True: #break to exit
        clearscreen()
        print "Pre-processing"
        print "Current top-level directory: %s" %base_dir
        print "Currently selected run     : %s" %run_name
        print "[D] Change top-level directory"       
        if run_name == None:
            print "[R] Select a run in this directory"
        else:
            print "[R] Select a different run in this directory"
        print "[N] Set up new run" 
        if base_dir != None and run_name != None:
            print "[1] Download new 'main' satellite data"
            print "[2] Download new 'pacific' satellite data"
            print "[3] Extract data from main satellite files and generate input for AMF calculator"
            print "[4] Run AMF calculator"
            print "[5] Apply pacific correction to slant columns"
            print "[6] Read AMF output and calculate VCs"
            print "[7] Filter for various critera"
            print "[8] Get equivalent BRUG data"
            print "[A] Do steps 3-7"
        print "[M] Merge some pickles"
        print "[Z] Go to previous menu"
        option = raw_input("-->").upper()
        
        if option == "Z": #quit
            break
        
        if option =="D": #change top-level directory
            base_dir = raw_input("Enter new top-level directory -->")
            if not os.path.isdir(base_dir):
                print "Not a valid path"
                base_dir = None
                run_name = None
                _ = raw_input("Press enter to continue")
            else:
                run_name = None
                continue
            
        if option == "R":
            if base_dir == None:
                print "Need to select a valid top-level directory first"
                _ = raw_input("Press enter to continue-->")
                continue
            run_name = run_select(base_dir)
            continue
                
        if option == "N":
            if base_dir == None:
                print "Need to select a valid top-level directory first"
                _ = raw_input("Press enter to continue-->")
                continue
            (valid,new_run_name)=new_run(base_dir=base_dir)
            if valid:
                run_name=new_run_name
            else:
                _ = raw_input("Run setup failed. Press enter to continue-->")
            continue
            
        sat_main    = add_slash(os.path.join(base_dir,run_name,"sat_main/"))
        sat_pacific = add_slash(os.path.join(base_dir,run_name,"sat_pacific/"))
        geos_main   = add_slash(os.path.join(base_dir,run_name,"geos_main/"))
        geos_pacific= add_slash(os.path.join(base_dir,run_name,"geos_pacific/"))
        amf_input   = add_slash(os.path.join(base_dir,run_name,"amf_input/"))
        amf_output  = add_slash(os.path.join(base_dir,run_name,"amf_output/"))
        p1          = add_slash(os.path.join(base_dir,run_name,"p1/"))
        p2          = add_slash(os.path.join(base_dir,run_name,"p2/"))
        p3          = add_slash(os.path.join(base_dir,run_name,"p3/"))
        p4          = add_slash(os.path.join(base_dir,run_name,"p4/"))
        pp          = add_slash(os.path.join(base_dir,run_name,"pp/"))
        date_list   = pickle.load(open(os.path.join(base_dir,run_name,"date_list.p"),"rb"))
                    
        if option == "1":
            download_sat(base_dir,run_name,main_pacific="main")
            continue
        
        if option == "2":
            download_sat(base_dir,run_name,main_pacific="pacific")
        
        
        if option == "3":
            date_dos = assemble_sat_data(sat_main,date_list,p1,amf_input)
            
            print "Enter AMF calculator directory"
            amf_run_dir = os.path.join(base_dir,"amf581/")
            print "Press enter to use %s" %amf_run_dir          
            amf_run_dir_in = raw_input("-->")
            if amf_run_dir_in !=  "":
                amf_run_dir = amf_run_dir_in
                      
            write_bsub(date_dos,amf_run_dir,amf_input,amf_output,geos_main)
            continue
        
        if option == "4":
            print 'chmod +x %sbatch_jobs.sh' %amf_input
            os.system('chmod +x %sbatch_jobs.sh' %amf_input)
            print '%sbatch_jobs.sh' %amf_input
            os.system('%sbatch_jobs.sh' %amf_input)
            print "Jobs have been submitted"
            print "Calculation of the AMFs will probably take a few hours"
            print "Output will be written to %s" %amf_output
            print "This script is not able to detect when this process is complete"
            print "Check the ''bjobs'' UNIX command to see when processing is complete"
            _ = raw_input("Press enter to continue-->")
            continue
        
        if option == "5": #apply pacific correction
            print "Loading pickled data"
            in_p1s = load_daily_pickles(p1)
            print "Calculating corrections"
            p2s = pacific_correction(in_p1s,base_dir,run_name)
            print "Saving pickle"
            save_daily_pickles(p2s,p2,"slant_corrected_unfiltered")
            
        if option == "6":
            read_in_AMF(amf_output,p2,p3,"vertical_corrected_unfiltered",date_list)
            
        if option == "7":
            filtering(p3,p4)
        
        if option == "8":
            this_box = box(40.,0.,105.,60.)
            if not(this_box.valid):
                raise IOError("Not a valid box!")
            BRUG = get_BRUG(date_list,this_box)
            cPickle.dump(BRUG,open(join(p4,"BRUG.p"),"wb"))
        
        if option == "M":
            merge_pickles()

def download_sat(base_dir,run_name,main_pacific):
    clearscreen()
    print "Data cannot be automatically downloaded using this script"
    print "Downloading instructions for %s satellite data:" %main_pacific
    print "1: Go to http://mirador.gsfc.nasa.gov/"
    if main_pacific=="main":
        print "2: Use the keyword 'omhcho', and select the desired spatial and temporal bounds for the area/time of interest"
    elif main_pacific=="pacific":
        print "2: Use the keyword 'omhcho', choose the spatial bounds (-90,-180),(90,-120), and the desired temporal bounds"
    print "3: Click 'Search GES-DISC'"
    print "4: Select the top check-box and click 'Add Selected Files To Cart'"
    print "5: Click 'Continue to Cart'"
    print "6: Click 'Checkout'"
    if main_pacific == "main":
        download_location = os.path.join(base_dir, run_name,"sat_main/")
    elif main_pacific=="pacific":
        download_location = os.path.join(base_dir, run_name,"sat_pacific/")
    print "7: Follow the intructions on the page to download the files (metadata optional) to %s" %download_location
    _ = raw_input("Press enter to continue-->")


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
    
    scan_n  = d(scan_n_list,"scan",description="Scan number")
    track_n = d(track_n_list,"row",description="Row number")
            
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
           list(CA),
           "SC",description="Slant Column Amount")
    SC.unit = "molec/cm2"
    DSC= d(
           list(DCA),
           "DSC",description="Slant Column Uncertainty")
    DSC.unit = "molec/cm2"
    
    #it may also be useful to get the final VC here, even though we'll
    #calculate our own
    
    NASA_VC = d(
                list(np.array(df['ReferenceSectorCorrectedVerticalColumn']).flatten()),
                "Vertical Column (NASA)",
                description="Vertical Column (as calulated in OMHCHO)")
    NASA_VC.unit = "molec/cm2"
    
    #Cloud info
    #Cloud infomation is taken from OMCLDO2 product
    cloud_fraction = d(
                       list(np.array(df['AMFCloudFraction']).flatten()),
                       "cloud fraction",
                       description="Cloud fraction from OMCLDO2")
    cloud_pressure = d(
                       list(np.array(df['AMFCloudPressure']).flatten()),
                       "cloud top pressure",
                       description="Cloud-top pressure from OMCLDO2")
    cloud_pressure.unit = "hPa"
                       
    
    #Quality flags
    main_data_quality_flag = d(
                                list(np.array(df['MainDataQualityFlag']).flatten()),
                                "main data quality flag")
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

def merge_pickles():
    i = 0
    pickle_in_l=[]
    while True:
        i+=1
        pickle_in = raw_input("Full path to pickle number %i. Blank entry will finalise list.-->"%i)
        if pickle_in == "":
            break
        pickle_in_l.append(pickle_in)
    
    to_merge = []    
    for i in range(0,len(pickle_in_l)):
        to_merge.append(cPickle.load(open(pickle_in_l[i],"rb")))
           
    merged =merge_OMI(to_merge,do_daylist=False)
    pickle_out = raw_input("Enter full path to save merged pickle to-->")
    cPickle.dump(merged,open(pickle_out,"wb"))
        
        
def merge_OMI(all_OMIs,do_daylist=True):
    """Merges multiple OMI_all files together"""
    if do_daylist:
        daylist = []
    for OMI in all_OMIs:
        if do_daylist:
            daylist.append(OMI.meta["day"])
        if OMI == all_OMIs[0]:
            #if first set
            merged = OMI
        else:
            merged.lat.extend(OMI.lat)
            merged.lon.extend(OMI.lon)
            try:
                merged.time.extend(OMI.time)
            except AttributeError:
                #OMI does not have to have a time field
                pass
                
            for dataset in merged.data:
                merged.data[dataset].val.extend(OMI.data[dataset].val)
    if do_daylist:
        del merged.meta["day"]
        all_days = list(set(daylist))
        all_days.sort()
        merged.meta["days"] = all_days
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
                          OMI_all.data["scan"].val[line],
                          OMI_all.data["row"].val[line],
                          OMI_all.lat[line],
                          OMI_all.lon[line],
                          OMI_all.data["SZA"].val[line],
                          OMI_all.data["VZA"].val[line],
                          OMI_all.data["cloud fraction"].val[line],
                          OMI_all.data["cloud top pressure"].val[line])
                     )
    outfile.close()

def check_HDFs_in_directory(directory,HDF_type="OMHCHO"):
    """For a given directory, returns list of good HDF files and their dates"""
    
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
        if HDF_type=="OMHCHO":
            granule_date = f_HDF['HDFEOS']['SWATHS']['OMI Total Column Amount HCHO']['Geolocation Fields']['TimeUTC'][0]
            this_date = dt(granule_date[0],granule_date[1],granule_date[2],0,0,0)
        elif HDF_type=="BRUG":
            gf = f_HDF['HDFEOS']['SWATHS']['Geolocation Fields']
            this_date = dt(int(gf['Year'][0]),int(gf['Month'][0]),int(gf['Day'][0]),0,0,0)
            del gf
        
        dates_of_files.append(this_date)
        
    return(files,dates_of_files)
    
def assemble_sat_data(sat_directory,date_list,OMI_save,AMF_txt_save):
    """Initial reading of satellite data"""
    (HDF_files,dates_of_files) = check_HDFs_in_directory(sat_directory)
    
    start_date = min(date_list)
    end_date = max(date_list)
    
    i_to_keep = [i for i in range(0,len(dates_of_files)) 
                 if dates_of_files[i] in date_list ]

    HDF_files = [HDF_files[i] for i in i_to_keep]
    dates_of_files = [dates_of_files[i] for i in i_to_keep]
    
    print "%i HDF files found in %s" %(len(HDF_files),sat_directory)
    
    if start_date==None:
        start_date = min(dates_of_files)
    if end_date==None:
        end_date = max(dates_of_files)
    
    #do process for each day
    have_data = [] #an empty dictionary that will be populated 
    date_clock = start_date
    while date_clock <= end_date:
        print "Assembling data for "+ date_clock.strftime("%Y-%m-%d")
        this_day_HDF_files = [HDF_files[i] for i in range(0,len(HDF_files))
                          if dates_of_files[i] == date_clock ]
        print "Found %i HDF files for this day" %len(this_day_HDF_files)
        if this_day_HDF_files == []: #if no data for this day
            have_data.append([date_clock,False])
            date_clock += td(days=1)
            continue
        else:
            have_data.append([date_clock, True])
        #print this_day_HDF_files
        these_OMIs = []
        for H in this_day_HDF_files:
            these_OMIs.append(OMI_from_HDF(h5py.File(H)))
        #print these_OMIs
        this_day_OMI = merge_OMI(these_OMIs,do_daylist=False)
        del these_OMIs
        
        #assign day
        this_day_OMI.meta["day"] = date_clock
        
        #assign ULN
        this_day_OMI.add_data(d(range(0,len(this_day_OMI.lat)),"ULN",description="Unique Line Number"))
        
        #save day's OMIs as pickle        
        pickle_name = add_slash(OMI_save)+date_clock.strftime("%Y-%m-%d")+"_OMI.p"
        print "Saving data as pickle: %s" %pickle_name
        cPickle.dump(this_day_OMI,open(pickle_name,"wb"))
        
        #write an AMF input .csv file        
        AMF_write_name = date_clock.strftime("%Y-%m-%d")+"_for_AMF.csv"
        print "Writing AMF input: %s" %AMF_write_name 
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
    
def write_bsub(date_dos,AMF_rundir,AMF_input,AMF_output,geos_main):
    """writes a shell script that can be executed to do AMF processing in parallel"""
    #fixed variables
    
    GEOS_ND51_files = '%sts_omi.YYYYMMDD.bpch' %geos_main
    AMF_inputs = add_slash(AMF_input) + "YYYY-MM-DD_for_AMF.csv"
    bsub_stuff = 'bsub -q short-serial -o %s/log-YYYYMMDD.log -J amf-YYYYMMDD -W 2:00' %AMF_output
    output_location = '%sYYYYMMDD_amfs.hcho' %AMF_output
    outfile = open(add_slash(AMF_input)+"batch_jobs.sh",'w')
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
         
def read_in_AMF(AMF_dir,pickle_dir_2,
                pickle_dir_3,save_pickle_name,
                date_list):
    """Reads in AMFs and calculates VCs"""
    #INPUTS
    #AMF_dir = folder for output of AMF calculator
    #pickle_dir_2 = folder with pickles containing SC data
    #pickle_dir_3 = folder to which to save pickle with VCs
    #date_list = list of dates to consider
    
    #add slashes to paths if not already supplied.
    AMF_dir=add_slash(AMF_dir)
    pickle_dir_2=add_slash(pickle_dir_2)
    pickle_dir_3=add_slash(pickle_dir_3)
    
    #add .p to filename if not already supplied.
    if not(save_pickle_name.endswith(".p")):
        save_pickle_name = save_pickle_name + ".p"
    
    print_missing_dates=True
    
    #loop through days
    daily_new_OMIs = []
    for date_clock in date_list :
        print "Processing %s..." %date_clock.strftime("%Y-%m-%d")
        #open this day's pickle
        pickle2_path = pickle_dir_2 + date_clock.strftime("%Y-%m-%d")+"_slant_corrected_unfiltered.p"
        try:
            OMI_in = pickle.load(open(pickle2_path,"rb"))
        except IOError:
            #file doesn't exist
            if print_missing_dates:
                print "%s does not exist" %pickle2_path           
            continue #go to next date
        
        #open this day's AMF output
        AMF_path = AMF_dir + date_clock.strftime("%Y%m%d")+"_amfs.hcho"
        try:
            AMFs_in = open(AMF_path,"rb")
        except IOError:
            #file doesn't exist
            if print_missing_dates:
                print "%s does not exist" %AMF_path            
            continue #go to next date
        
        #loop through lines in AMFs_in
        a_ULN=[]
        a_AMF=[]
        a_GEOS=[]       
        for line in AMFs_in:
            line = line.strip() #get rid of the newline character
            columns = line.split(",") #split up by columns, space separated.
            a_ULN.append(int(columns[1]))
            a_AMF.append(float(columns[2]))
            a_GEOS.append(float(columns[3]))           
        AMFs_in.close() #close file
        
        len_OMI=len(OMI_in.data['ULN'].val)
        len_a  =len(a_ULN)
        
        AMF_list = []
        sat_VC_list = []
        sat_DVC_list = []
        geos_VC_list = []
        #print OMI_in.data['ULN'].val[100:200]
        #print a_ULN[100:200]
        for j in range(0,len_OMI):
            if OMI_in.data['ULN'].val[j] == a_ULN[j]:
                #simple match
                if a_GEOS[j] <= 0.:
                    geos_VC_list.append(np.nan)
                else:                
                    geos_VC_list.append(a_GEOS[j])
                if np.isnan(a_AMF[j]) or a_AMF[j] <= 0.:
                    AMF_list.append(np.nan)
                    sat_VC_list.append(np.nan)
                    sat_DVC_list.append(np.nan)
                else:
                    AMF_list.append(a_AMF[j])
                    if OMI_in.data['SC'].val[j] == np.nan:
                        sat_VC_list.append(np.nan)
                        sat_DVC_list.append(np.nan)
                    else:
                        sat_VC_list.append(OMI_in.data['SC'].val[j]/a_AMF[j])
                        sat_DVC_list.append(OMI_in.data['DSC'].val[j]/a_AMF[j])
            else:
                #ULNs not perfectly aligned.
                print "ULNs not aligned between %s and %s" %(
                    pickle2_path,AMF_path)
                print "You need to do more coding Luke"
                raise IOError
        
        OMI_in.add_data(d(AMF_list,"AMF","Air mass factor"))
        OMI_in.add_data(d(sat_VC_list,"sat_VC","Observed vertical HCHO column"))
        OMI_in.data["sat_VC"].unit = 'molec/cm2'
        OMI_in.add_data(d(sat_DVC_list,"sat_DVC","Error in observed vertical HCHO column"))
        OMI_in.data["sat_DVC"].unit = 'molec/cm2'
        OMI_in.add_data(d(geos_VC_list,"geos_VC","Modelled vertical HCHO column"))
        OMI_in.data["geos_VC"].unit = 'molec/cm2'
        
        daily_new_OMIs.append(copy.deepcopy(OMI_in))        
        del OMI_in
                
        
            
    #now merge these OMIs
    print "Merging into one file..."
    OMI_all  = merge_OMI(daily_new_OMIs)
    
    #now save this
    print "Saving file..."
    pickle_save_path = pickle_dir_3 + save_pickle_name
    cPickle.dump(OMI_all,open(pickle_save_path,"wb"))
    print "AMFs retrived and VCs calculated. Saved to %s" %pickle_save_path
    
    return(OMI_all)


def pacific_correction(in_p1s,base_dir,run_name):
    """Corrects the data using pacific measurements"""
    
    out_p2s = []
    
    sat_pacific_dir = os.path.join(base_dir,run_name,"sat_pacific")  
    
    reference_sector_full = box(90,-90,-120,-180)
    if reference_sector_full.valid == False:
        raise ValueError("reference_sector_full incorrectly defined")
    reference_sector_eq = box(15,-15,-120,-180)
    if reference_sector_eq.valid   == False:
        raise ValueError("reference_sector_eq incorrectly defined")
        
    for p1 in in_p1s:
        date_clock = p1.meta["day"]
        print "Calculating pacific correction for %s" %date_clock.strftime("%Y-%m-%d")
        pacific_data = get_pacific(date_clock,sat_pacific_dir)
        if pacific_data == False:
            #if there were no files
            all_output.append(None) #empty object
            del pacific_data
            date_clock += td(days=1)
            continue
        
        
        
        #pacific_data is a d_all object holding pacific measurements for this day
        
        #do basic main data quality filtering
        for i in range(0,len(pacific_data.lat)):
            if pacific_data.data["main data quality flag"].val[i] != 0:
                pacific_data.data["SC"].val[i] = np.nan
            if pacific_data.data["Xtrack quality flag"].val[i] != 0:
                pacific_data.data["SC"].val[i] = np.nan
                
        good_i = np.isfinite(pacific_data.data["SC"].val)        
        #strictly restrict to designated area 
        
        print "pixels in pacific data downloaded: %i" %len(pacific_data.lat)
        
        #geo_select_rectangle_i(reference_sector_full,pacific_data)
        pacific_data = geo_select_rectangle_i(reference_sector_full,pacific_data)

        
        print "pixels in pacific data, full reference sector: %i" %len(pacific_data.lat)

        #No filtering of this data

        #make a pacific_data_equatorial object    
        pacific_data_eq = copy.deepcopy(pacific_data)
        pacific_data_eq = geo_select_rectangle_i(reference_sector_eq,pacific_data_eq)
        #geo_select_rectangle_i(reference_sector_eq,pacific_data_eq)
        
        print "pixels in pacific data, equatorial reference sector: %i" %len(pacific_data_eq.lat)
        
        good_i = np.isfinite(pacific_data.data["SC"].val)    

        #get row-dependant medians of SC
        print "Calcululating row-dependant medians in equatorial reference sector"
        row_dependant_median = {}
        for row in range(0,60):
            these_SCs = np.array([pacific_data_eq.data["SC"].val[i] 
                         for i in range(0,len(pacific_data_eq.data["SC"].val))
                         if pacific_data_eq.data["row"].val[i] == row
                         and pacific_data_eq.data["SC"].val[i] != np.nan ])
                         
                         
            row_dependant_median[row] = np.nanmedian(these_SCs)
        
        print "Row medians in reference column:"
        print row_dependant_median
        
        #Correct the offset in pacific_data
        
        for i in range(0,len(pacific_data.data["SC"].val)):
            pacific_data.data["SC"].val[i] -= row_dependant_median[pacific_data.data["row"].val[i]]

            
            
        good_SC = np.array(pacific_data.data["SC"].val)[good_i]
        good_lat = np.array(pacific_data.lat)[good_i]
        #determine polynomial for latitudinal dependance
        print "Calculating polynominal latitudinal dependance"
        poly = np.polyfit(good_lat,good_SC,4)
        #deg being 4 is arbitrary choice
        #could try and weight fitting points by uncertainties...
        
        print "Polynominal: %gy^4 + %gy^3 + %gy^2 + %gy + %g"\
               %(poly[0],poly[1],poly[2],poly[3],poly[4]) 
        print "Applying corrections"

        for i in range(0,len(p1.data["SC"].val)):
            p1.data["SC"].val[i] \
            -= row_dependant_median[p1.data["row"].val[i]] \
               + np.polyval(poly,p1.lat[i])
        
        out_p2s.append(p1)
    
    return(out_p2s)
        
def get_pacific(this_date,sat_pacific_dir):
    """Loads pacific data from he5 files"""
    
    start_date = this_date
    end_date = this_date + td(days=1)
    
    (HDF_files,dates_of_files) = check_HDFs_in_directory(sat_pacific_dir)
        

    i_to_keep = [i for i in range(0,len(dates_of_files)) 
                 if dates_of_files[i] >= start_date 
                 and dates_of_files[i] < end_date ]

    HDF_files = [HDF_files[i] for i in i_to_keep]
    dates_of_files = [dates_of_files[i] for i in i_to_keep]
       
    #do process for each day    
    
    print "Assembling pacific data for "+ this_date.strftime("%Y-%m-%d")
    this_day_HDF_files = [HDF_files[i] for i in range(0,len(HDF_files))
                      if dates_of_files[i] == this_date ]
    if this_day_HDF_files == []: #if no data for this day
        return False
    
    
    these_OMIs = []
    for H in this_day_HDF_files:
        these_OMIs.append(OMI_from_HDF(h5py.File(H)))
    
    this_day_OMI = merge_OMI(these_OMIs,do_daylist=False)
    del these_OMIs
    
    return(this_day_OMI)

def load_daily_pickles(directory):
    """Loads all of the pickles saved in a folder and returns them as a list"""
    #Each pickle must contain a single object
    
    directory = add_slash(directory)
    #get a list of all files.
    files = [directory+f for f in listdir(directory) if isfile(join(directory, f))] 
    
    #keep only those ending in .p
    files = [f for f in files if f.endswith('.p') ]
    
    pickle_list = [] 
    for f in files:
        pickle_list.append(cPickle.load(open(f,"rb")))
        
    return pickle_list
    
def save_daily_pickles(p_in,directory,suffix):
    """Saves a list of objects each corresponding to a different date"""
    for p in p_in:
        this_date = p.meta["day"]
        pickle_name = add_slash(directory)+this_date.strftime("%Y-%m-%d")+"_"+suffix+".p"
        cPickle.dump(p,open(pickle_name,"wb")) 
        
def filtering(p3_dir,p4_dir):
    """Filters out data which fails various quality standards"""
    
    #load_the_pickle
    p3_files = [add_slash(p3_dir)+f
              for f in listdir(add_slash(p3_dir))
              if isfile(join(p3_dir, f))]
    
    if len(p3_files) > 1:
        print "More than one pickle in %s" %p3_dir
        print "Choose one file:"
        for i in range(0,len(p3_files)):
            print "[%i] %s" %(i,p3_files[i])
        choice = int(raw_input("-->"))
        p3_file = p3_files[choice]
        
    else:
        p3_file = p3_files[0]       
    
    f = p3_file.replace(add_slash(p3_dir),"")
              
    print "Opening pickled data in %s" %p3_file
    
    ida = cPickle.load(open(p3_file,"rb"))
   
    #let the user choose default or custom filtering option
    filter_option = basic_menu("Default or custom filtering?",
                               [["1","Default filtering"],
                                ["2","Custom filtering"]],
                                quit_option=True)
    if filter_option == "Z":
        return
    elif filter_option == "1":
        print "Default filtering chosen"
        #Main data quality flag
        print "Removing points where Main Data Quality Flag =! 0"
        print "Removing points where Xtrack Quality Flag != 0"
        print "Removing points where Solar Zenith Angle > 70 deg"
        print "Removing points where Cloud fraction > 0.4"
        print "Clearing points where AMF or pacific correction failed"
        print "Removing points where vertical column errors larger than 3 times column"
        filter_i = []
        num_pre_filter = len(ida.lat)
        for i in range(0,num_pre_filter):
            if ida.data["main data quality flag"].val[i] != 0 :
                filter_i.append(False)
                continue
            if ida.data["Xtrack quality flag"   ].val[i] != 0 :
                filter_i.append(False)
                continue
            if ida.data["SZA"                   ].val[i] > 70.:
                filter_i.append(False)
                continue
            if ida.data["cloud fraction"        ].val[i] > 0.4 :
                filter_i.append(False)
                continue
            if np.isnan(ida.data["sat_VC"       ].val[i]):
                #AMF or pacific correction failed here
                filter_i.append(False)
                continue
            if ida.data["sat_DVC"].val[i] > ida.data["sat_VC"].val[i]*3:
                filter_i.append(False)
                continue
            #if we get this far, this data point has passed
            filter_i.append(True)
        #filter this data
        ida.filter_all(filter_i)

        #add meta information about filtering
        ida.meta["Filtering criteria"] = \
                    "main data quaility flag == 0\n"+\
                    "Xtrack quaility flag == 0\n"+\
                    "Solar zenity angle <= 70.\n"+\
                    "Cloud fraction <= 0.4\n"+\
                    "DVC < VC*3"
        ida.meta["Filtering stats"] = \
            {"Total datapoints":num_pre_filter,
             "Retained datapoints":len(ida.lat),
             "Removed datapoints":num_pre_filter-len(ida.lat)}
        
        print ida.meta["Filtering stats"]     
        
        save_path = join(p4_dir,"vertical_corrected_filtered.p")
        print "Saving as %s" %save_path
        cPickle.dump(ida,open(join(save_path),"wb"))
        _ = raw_input("Press enter to continue")
        return
    elif filter_option == "2":
        print "I haven't programmed this yet"
        return   

def get_BRUG(date_list,this_box):
    """For a particular pickle, brings in other OMI data from a particular directory"""
    
    #match_dir = raw_input("In which directory can BRUG data be found?-->")
    match_dir = '/group_workspaces/cems2/nceo_generic/nceo_ed/BE_H2CO/OMI/YYYY/MM/DD/'
    #ida is a d_all object containing all the HCHO data. No data has been removed for filtering.
    
    #get list of days covered by ida
    
    
    #find all the HDF files
    HDF_files = []
    dates_of_files = []
    for today in date_list:
        today_match_dir = YYYYMMDD_sub(today,match_dir)
        print "Checking for HDF files in %s" %today_match_dir
        try:
            (today_HDF_files,today_dates_of_files) = check_HDFs_in_directory(today_match_dir,HDF_type="BRUG")
            HDF_files.extend(today_HDF_files)
            dates_of_files.extend(today_dates_of_files)
        except OSError:
            print "%s cannot be accessed (may not exist)"%today_match_dir
            pass
    
    start_date = min(date_list)
    end_date = max(date_list)
    
    i_to_keep = [i for i in range(0,len(dates_of_files)) 
                 if dates_of_files[i] in date_list ]

    HDF_files = [HDF_files[i] for i in i_to_keep]
    dates_of_files = [dates_of_files[i] for i in i_to_keep]
    
    #extract information from BRUG
    print "extracting information from BRUG"
    BRUG_collection = []
    for H in HDF_files:
        print "Extracting from %s" %H
        this_BRUG = OMI_BRUG_from_HDF(h5py.File(H))
        #print (this_BRUG.lat)
        this_BRUG_nest = geo_select_rectangle_i(this_box,this_BRUG)
        BRUG_collection.append(copy.deepcopy(this_BRUG_nest))
        del this_BRUG, this_BRUG_nest
    print "Merging..."
    BRUG = merge_OMI(BRUG_collection)
    
    #Directly associating this data seems like a bust.
    #Let's just return it and see if we can play with the data somehow.
    
    return(BRUG)
            
def TAI93_to_UTC(TAI93):
    """converts a TAI93 timestamp to a datetime object"""
    #Doesn't account for leap seconds stuff, so +/- a few seconds may be necessary
    tai_epoch = dt(1993,1,1,0,0,0)
    return(tai_epoch+td(seconds=TAI93))


def OMI_BRUG_from_HDF(HDF_file):
    """Returns various sets of data from an OMI BRUG file as lists"""
     
    #shorthands
    df = HDF_file['HDFEOS']['SWATHS']['Data Fields']
    gf = HDF_file['HDFEOS']['SWATHS']['Geolocation Fields']
    
    (n_scans,n_tracks) = df['H2CO Columns']["VCD"].shape #shape of the data
    
    #get simple location data 
    #there are 5 lats for each data point in this file, just get the last one for centre
    lat_5 = list(np.array(gf['Latitudes']).flatten())
    lat=lat_5[4::5]

    lon_5 = list(np.array(gf['Longitudes']).flatten())
    lon=lon_5[4::5]

    #Now get the corners:
    lat_corners = []
    lon_corners = []
    for i in range(0,len(lat)):
        lat_corners.append([lat_5[i*5+0],lat_5[i*5+1],lat_5[i*5+2],lat_5[i*5+3]])
        lon_corners.append([lon_5[i*5+0],lon_5[i*5+1],lon_5[i*5+2],lon_5[i*5+3]])
    
    
    LAT_cor = d(lat_corners,"lat_corners",description="Pixel corner latitudes")
    LAT_cor.unit = u"\u00b0" #degree sign
    LON_cor = d(lon_corners,"lon_corners",description="Pixel corner longitudes")
    LON_cor.unit = u"\u00b0" #degree sign
    
    #scan and row/track identifiers
    scan_n_list = []
    track_n_list = []

    for scan in range(0,n_scans):
        for track in range(0,n_tracks):
            scan_n_list.append(scan)
            track_n_list.append(track)
    
    scan_n  = d(scan_n_list,"scan",description="Scan number")
    track_n = d(track_n_list,"row",description="Row number")
    
    #print "time bit start..."
    time_1 = list(np.array(gf['Time']).flatten())
    time = []
    for t in range(0,len(time_1)):
        time.append(TAI93_to_UTC(int(time_1[t])))
    
    #print "time bit end..."  
             
    #AMF    
    AMF = d(
            list(np.array(df['H2CO AMF']['AMF']).flatten()),
            "BRUG_AMF",description="AMF (BRUG)")
    
    #VC    
    VC = d(
            list(np.array(df['H2CO Columns']['VCD']).flatten()),
            "BRUG_VC",description="Vertical column amount (BRUG)")
    VC.unit = "molec/cm2"
    
    #DVC
    DVC = d(
            list(np.array(df['H2CO Errors']['VCD0E_S']).flatten()),
            "BRUG_DVC",description="Vertical column uncertainty (BRUG)")
    DVC.unit = "molec/cm2"
    
    #Normalised Slant
    
    NSC = d(
            list(np.array(df['H2CO Columns']['SCD3']).flatten()),
            "BRUG_DVC",description="Normalised sland column (BRUG)")
    NSC.unit = "molec/cm2"
    
    #QF
    QF = d(
            list(np.array(df['H2CO Errors']['QF']).flatten()),
            "BRUG_QF",description="Quality Flag (BRUG)")
    QF.unit = "molec/cm2"
    
    #day of file
    file_date = dt(int(gf['Year'][0]),int(gf['Month'][0]),int(gf['Day'][0]),0,0,0)
    
      
    OMI_BRUG_all = d_all(lat,lon,time=time)
    OMI_BRUG_all.add_data([scan_n,track_n,AMF,VC,DVC,QF])
    OMI_BRUG_all.add_data([LAT_cor,LON_cor])
    OMI_BRUG_all.meta['day']=file_date
    
    return(OMI_BRUG_all)                        
      
preprocessing(base_dir="/group_workspaces/jasmin/geoschem/local_users/lsurl/CL_PP")         
