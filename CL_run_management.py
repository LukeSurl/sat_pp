#!/usr/bin/env python
"""CL_run_management.py - scripts to record and report where data is saved"""

import os
import numpy as np
from CL_satpp import add_slash,basic_menu,petc
from datetime import datetime as dt
from datetime import timedelta as td
import pickle
import cPickle

def change_dates(run_dir):
    #define dates
    print "What is the first date to be considered in this run?"
    start_year = input("YEAR-->")
    start_month = input("MONTH-->")
    start_day = input("DAY-->")
    start_date = dt(start_year,start_month,start_day,0,0,0)
    print "What is the last date to be considered in this run?"
    end_year = input("YEAR-->")
    end_month = input("MONTH-->")
    end_day = input("DAY-->")
    end_date = dt(end_year,end_month,end_day,0,0,0)
    
    date_list = []
    date_clock = start_date
    while date_clock <= end_date:
        date_list.append(date_clock)
        date_clock+=td(days=1)
    cPickle.dump(date_list,open(run_dir+"/date_list.p","wb"))

def new_run(base_dir=None,run_name=None,sat_main=None,sat_pacific=None,geos_main=None,geos_pacific=None):
    """A script to create the space required for a completely new run"""
    
    print "Creating directories for a new run"
    if run_name == None:
        run_name=raw_input("Enter name of new run -->")
            
    print "Creating directories for run %s" %run_name
    if base_dir==None:
        current_dir = os.path.dirname(os.path.realpath(__file__))
        print "Run directory will be subdirectory of %s" %current_dir
        print "If different location required, enter it here, otherwise press enter to use the above"
        option = raw_input("-->")
        if option == "":
            base_dir = add_slash(current_dir) + run_name
        else:
            base_dir = add_slash(option) + run_name
        del option
    elif base_dir=="current":
        current_dir = os.path.dirname(os.path.realpath(__file__))
        base_dir = add_slash(current_dir) + run_name
    else:
        base_dir = base_dir
    
    
    run_dir = add_slash(os.path.join(base_dir,run_name))
    
    
    #check that directory does not already exist
    if os.path.exists(run_dir):
        print "Error, %s already exists" %run_dir
        return(False,run_name)
    else:
        os.mkdir(run_dir)
    
    #define dates
    change_dates(run_dir)
    
    #species = "HCHO"
    
    #choose species
    species=""
    while species not in ["HCHO","NO2"]:
        species = raw_input("Enter species:\n[HCHO]\n[NO2]\n-> ").upper()
    species_file = open(run_dir+"species.txt","wb")
    species_file.write(species)
    species_file.close()   
            
    #create directories and subdirectories
    
    #Main satellite data
    if sat_main==None:
        sat_main = raw_input("Enter directory containing main satellite data \n(symbolic link will be created)-->")    
    os.symlink(sat_main,run_dir+"sat_main")
    print "Symbolic link to %s created at %s" %(sat_main,base_dir+"sat_main/")
    
    #Pacific satellite data
    if sat_pacific==None:
        sat_pacific = raw_input("Enter directory containing satellite data covering pacific calibration region \n(symbolic link will be created)-->")
    os.symlink(sat_pacific,run_dir+"sat_pacific")
    print "Symbolic link to %s created at %s" %(sat_pacific,base_dir+"sat_pacific/")
    
    #GEOS outout main
    if geos_main==None:
        geos_main = raw_input("Enter directory containing the main GEOS Chem output \n(symbolic link will be created)-->")
    os.symlink(geos_main,run_dir+"geos_main")
    print "Symbolic link to %s created at %s" %(geos_main,base_dir+"geos_main/")
    
    #GEOS pacific
    if geos_pacific==None:
        geos_pacific = raw_input("Enter directory containing the GEOS Chem output"+\
                                       "covering the Pacific.\n"+\
                                       "This will probably be global run output. \n(symbolic link will be created)-->")
    os.symlink(geos_pacific,run_dir+"geos_pacific")
    print "Symbolic link to %s created at %s" %(geos_main,base_dir+"geos_main/")
    
    #AMF input
    os.mkdir(run_dir+"amf_input/")
    print "Directory %s created" %(run_dir+"amf_input/")
    
    #AMF output
    os.mkdir(run_dir+"amf_output/")
    print "Directory %s created" %(run_dir+"amf_output/")
    
    #Slant column pickles
    os.mkdir(run_dir+"p1/")
    print "Directory %s created" %(run_dir+"p1/")
  
    #Vertical column pickles, unfiltered, uncorrected
    os.mkdir(run_dir+"p2/")
    print "Directory %s created" %(run_dir+"p2/")
    
    os.mkdir(run_dir+"p3/")
    print "Directory %s created" %(run_dir+"p3/")
    
    #Vertical column pickles, filtered, corrected
    os.mkdir(run_dir+"p4/")
    print "Directory %s created" %(run_dir+"p4/")
    
    #Post-processing folder
    os.mkdir(run_dir+"pp/")
    print "Directory %s created" %(run_dir+"pp/")
    
    return(True,run_name,species)
                                              
def run_select(base_dir):
    """Asks the user to pick an existing run"""
    
    #get all subdirectories of base_dir
    subdirs = [name for name in os.listdir(base_dir)
            if os.path.isdir(os.path.join(base_dir, name))]
    
    #restrict to those with the correct subdirectories    
    subdirs = [subdir for subdir in subdirs
               if  os.path.isdir(os.path.join(base_dir,subdir,"p1"))
               and os.path.isdir(os.path.join(base_dir,subdir,"p2"))
               and os.path.isdir(os.path.join(base_dir,subdir,"p3"))
               and os.path.isdir(os.path.join(base_dir,subdir,"p4"))
               and os.path.isdir(os.path.join(base_dir,subdir,"pp"))
               ]
    menu_text = []
    for i in range(0,len(subdirs)):
        menu_text.append([str(i),subdirs[i]])
    
    if menu_text == []:
        print "There are no valid runs in %s" %base_dir
        _ = raw_input("Press enter to continue-->")
        return None
        
    choice = basic_menu("Choose an existing run:",menu_text)
    if choice == "Z":
        return None
    else:
        return subdirs[int(choice)]           

def read_species(run_dir):
    """Reads which species species.txt is pointing to. Assumes HCHO if no info"""
    
    #check species.txt exists
    if not os.path.exists(os.path.join(run_dir,"species.txt")):
        #if it doesn exist, assume HCHO and write the file
        switch_species(run_dir,"HCHO")
        return("HCHO")
    else:
        species_file = open(os.path.join(run_dir,"species.txt"),"rb")
        species = species_file.read()
        if species not in ["HCHO","NO2"]:
            print "Error with read_species. species.txt file contents: %s.\n I do not recognise this as HCHO or NO2.\n I will assume HCHO."
            petc()
            species_file.close()
            switch_species(run_dir,"HCHO") #write a good HCHO file
            return ("HCHO")
        else:
            return(species)

def switch_species(run_dir,new_species):
    """Switches the defined species for this run"""
    
    species_file_path = os.path.join(run_dir,"species.txt")
    #delete prior species.txt if exists
    try:
        os.remove(species_file_path)
    except OSError:
        pass
    
    #write a new one
    species_file = open(species_file_path,"wb")
    species_file.write(new_species)
    species_file.close()   
    
def clone_run(current_run,main_script):
    """A script to create a run by partially copying another run"""
    
    #WRITE THIS
    
    
    
