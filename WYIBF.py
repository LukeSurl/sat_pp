#!/usr/bin/env python
"""A function to compute the Williamson-York Iterative Bivariate Fit"""

"""This is a reverse engineering of the spreadsheet with is the supplementary
material to Cantrell (2008) [Atmos. Chem. Phys., 8, 5477-5487, 2008]."""

__author__      = "Luke Surl"
__contact__     = "L.Surl@ed.ac.uk"

import numpy as np
from scipy.optimize import minimize

def WYIBF(x,y,dx,dy,err_type="var",iterate=True):
    #x = x-axis variable, 1D variable
    #y = y-axis variable, 1D variable
    #dx = error in each entry x, 1D variable
    #dy = error in each entry y, 1D varaible
    #--x, y, dx, dy must all be the same length---
    #err_type = can be either "stdev", "var", "weight", tell the code what dx & dy are. 
    
    #Input validation
    #number of data points
    num_points = len(x)
    if num_points != len(y) or num_points != len(dx) or num_points != len(dy):
        raise ValueError('in call to WYIBF, x, y, dx, dy are not all the same length')
        
    #err_type
    if err_type not in ["stdev","var","w"]:
        raise ValueError('in call to WYIBF, err_type must be either "stdev", "var", or "w"') 
    
    #convert stdevs to vars
    if err_type=="stdev":
        dx = [this_dx**2 for this_dx in dx]
        dy = [this_dy**2 for this_dy in dy]
        err_type="var"
    
    #compute standard linear fit
    par = np.polyfit(x, y, 1, full=True)
    #m_standard=par[0][0]
    m_standard = 1.0
    #b_standard=par[0][1]
    b_standard = 0.0
    
    
    print "First guess m = ",m_standard
    
    #use code letters as in spreadsheet,                   
    if err_type in ["var"]:
        err_code = "S"
    else :
        err_code = "W"
        
    #strip out points where variance in both x and y are zero
    if err_code == "S":
        x_fil = [ x[i] for i in range(0,num_points) if not(dx[i] == 0. and dy[i] == 0.)]
        y_fil = [ y[i] for i in range(0,num_points) if not(dx[i] == 0. and dy[i] == 0.)]
        dx_fil = [ dx[i] for i in range(0,num_points) if not(dx[i] == 0. and dy[i] == 0.)]
        dy_fil = [ dy[i] for i in range(0,num_points) if not(dx[i] == 0. and dy[i] == 0.)]
    
    if iterate:
        #solve to get best value for m
        #res = minimize(WYIBF_minimizer,m_standard,args=(x_fil,y_fil,dx_fil,dy_fil,err_code),method='Nelder-Mead',options={'xtol': 0.01, 'disp': True, 'maxiter': 10})
        res = minimize(WYIBF_minimizer,m_standard,args=(x_fil,y_fil,dx_fil,dy_fil,err_code),method='Nelder-Mead',options={'xtol': 0.01, 'disp': False, 'maxiter': 15})
        #run the iterator one last time with this best value to get the final results for all 5 output variables.
        (m,b,merr,berr,gof) = WYIBF_iterator(res.x,x_fil,y_fil,dx_fil,dy_fil,err_code)
    else:
        (m,b,merr,berr,gof) = WYIBF_iterator(m_standard,x_fil,y_fil,dx_fil,dy_fil,err_code)
    
    
    return (m,b,merr,berr,gof)
    
    
def WYIBF_minimizer(guess,x,y,dx,dy,err_code) :
    #print "iterating..."
    (m,b,merr,berr,gof) = WYIBF_iterator(guess,x,y,dx,dy,err_code)
    #print "CURRRENT: y = (",m," +/- ",merr,")*x + (",b," +/- ",berr,")    Goodness of fit: ",gof 
    diff = guess-m #this is the function the iterator will try to minimise
    return abs(diff)


def WYIBF_iterator(guess,x,y,dx,dy,err_code):
    
    #create some empty arrays of the same length as the input
    w         = np.zeros_like(x)
    wx        = np.zeros_like(x)
    wy        = np.zeros_like(x)
    B         = np.zeros_like(x)
    wBU       = np.zeros_like(x)
    wBV       = np.zeros_like(x)
    Bw        = np.zeros_like(x)
    wBminBbar2= np.zeros_like(x)
    S         = np.zeros_like(x)
    
    #count the number of points
    num_points = len(x)
    
    #Calculate w, wx, wy (cols GHI in spreadsheet)  
    for i in range(0,num_points):
        if err_code == "W":
            w[i]       = 1./((guess*guess/dx[i])+(1./dy[i]))
        elif err_code == "S":
            w[i]       = 1./((guess*guess*dx[i])+(dy[i]  ))
        else :
            raise ValueError('in call to WYIBF_iterator, err_code must be "W" or "S"')
        wx[i]      = w[i]*x[i]
        wy[i]      = w[i]*y[i]
        
    xbar = sum(wx)/sum(w) #H27
    ybar = sum(wy)/sum(w) #I27
    
    #Calculate B, wBU, wBV Bw (cols JKLM)
    for i in range(0,num_points):
        if err_code == "W":
            B[i]       = w[i]*((x[i]-xbar)/dy[i]+guess*(y[i]-ybar)/dx[i])
        elif err_code == "S":
            B[i]       = w[i]*((x[i]-xbar)*dy[i]+guess*(y[i]-ybar)*dx[i])
        else :
            raise ValueError('in call to WYIBF_iterator, err_code must be "W" or "S"')    
        wBU[i]     = w[i]*B[i]*(y[i]-ybar)
        wBV[i]     = w[i]*B[i]*(x[i]-xbar)
        Bw[i]      = B[i]*w[i]
        
    Bbar = sum(Bw)/sum(w) #M27
    
    #Calculate w(B-Bbar)^2 (col N)
    for i in range(0,num_points):    
        wBminBbar2[i] = w[i]*(B[i]-Bbar)*(B[i]-Bbar)
        
    m_out = sum(wBU)/sum(wBV) #G1
    b_out = ybar - m_out*xbar #G2
    
    #Calculate S (col P)
    for i in range(0,num_points):
        S[i] = w[i]*(y[i]-m_out*x[i]-b_out)*(y[i]-m_out*x[i]-b_out)
    
    gof_out  = sum(S)/float(num_points-2) #G5
    
    merr_out = np.sqrt(1./sum(wBminBbar2))*np.sqrt(gof_out) #G3
    berr_out = np.sqrt(1./sum(w)+(xbar+Bbar)*(xbar+Bbar)*1./sum(wBminBbar2))*np.sqrt(gof_out) #G4
    
    #print m_out

    return(m_out,b_out,merr_out,berr_out,gof_out)
    
    
        
            
            
            
            
        
        
            
        
    
    
