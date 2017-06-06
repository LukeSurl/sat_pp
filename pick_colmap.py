#!/usr/bin/env python

"""pick_colmap.py - A command-line interface for the brewer2mpl package"""


#usage:
#Place this script in the working directory.
# 
#import pick_colmap
#colmap = pick_colmap.pick_colmap(returntype="matplotlib")
#
#The interface will start and "colmap" will then be a matplotlib colormap object.
 

try:
    import brewer2mpl
except ImportError:
    print "pick_colmap.py needs to be in a (virtual) environment with brewer2mpl installed to work!"
    sys.exit()

import os

def free_colorbar(colmap,vmin=0.,vmax=1.):
    import matplotlib.pyplot as plt
    import matplotlib.colorbar as colorbar
    import matplotlib.cm as cm
    import matplotlib.colors as colors

    cm1 = colors.LinearSegmentedColormap.from_list("MyCmapName",colmap.mpl_colors)
    cnorm = colors.Normalize(vmin=vmin,vmax=vmax)
    cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
    cpick.set_array([])
    plt.colorbar(cpick)
    plt.show()


def pick_colmap(preview=False,returntype="normal"):
    os.system('cls' if os.name == 'nt' else 'clear') #clear the screen
    print "==COLOR MAP PICKER - refer to http://colorbrewer2.org for visuals==="
    map_type = ""
    while map_type not in ["S","D","Q"]:
        map_type = raw_input("Choose map type: (S)equential, (D)iverging, (Q)ualitative --> ").upper()
    
    if map_type == "S":
        map_type_full = "Sequential"
        options = ["Blues","BuGn","BuPu","GnBu","Greens",
                    "Greys","OrRd","Oranges","PuBu","PuBuGn",
                    "PuRd","Purples","RdPu","Reds","YlGn",
                    "YlGnBu","YlOrBr","YlOrRd"]
    elif map_type == "D":
        map_type_full = "Diverging"
        options = ["BrBG","PRGn","PiYG","PuOr","RdBu","RdGy","RdYlBu","RdYlGn","Spectral"]
    elif map_type == "Q":
        map_type_full = "Qualitative"
        options = ["Accent","Dark2","Paired","Pastel1","Pastel2","Set1","Set2","Set3"]
    
    map_choice = ""
    print "Pick a color map from the following:"
    for o in options:
        print "> %s" %o
    while map_choice not in options:
        map_choice = raw_input("--> ")
    
    print "Colour bar %s selected" %(map_type+" - "+map_choice)
    reverse = ""
    while reverse not in ["Y","N"]:
        reverse = raw_input("Reverse the color map? Y/N -- >").upper()
    if reverse=="Y":
        rev=True
    else:
        rev=False
    for n in range(12,7,-1): #get maximum number of colors
        try:
            bmap = brewer2mpl.get_map(map_choice, map_type_full, n, reverse=rev)
            break
        except ValueError:
            pass
    
    if preview:
        free_colorbar(bmap)
        
    #return options. See https://pypi.python.org/pypi/brewer2mpl/1.4 for other things you can do with this object    
    if returntype=="normal":
        return(bmap)
    elif returntype=="matplotlib":
        return(bmap.mpl_colormap)
    
        
        
