# -*- coding: utf-8 -*-
"""
@author: zadavi
"""
import wget
import os
import numpy as np
import pandas as pd
import os
from scipy.interpolate import griddata
import random
import calendar
#This script download VMF files from VMF server.Here, only operational VMF1 data is downloading. 
# VMF1_OP should be changed to VMF1_FC to download forecating type.
def mjd_to_ymd(mjd):
    jd = mjd + 2400000.5
    jd_int = int(jd)
    fraction = jd - jd_int

    # Calculate year, month, and day
    a = int((jd_int - 1867216.25) / 36524.25)
    b = jd_int + 1 + a - int(a / 4)
    c = b + 1524
    d = int((c - 122.1) / 365.25)
    e = int(365.25 * d)
    f = int((c - e) / 30.6001)

    day = int(c - e - int(30.6001 * f) + fraction)
    month = f - 1 - 12 * int(f / 14)
    year = d - 4715 - int((7 + month) / 10)
    return year,month,day



def download_VMF(year=None,month=None,day=None,directory_vmf=None): 

   LISTt=['00','06','12','18']
   for lstnum in range(0,len(LISTt)):
         try:
            filename='VMFG_'+str(year)+str(month).zfill(2)+str(day).zfill(2)+'.H'+LISTt[lstnum] 
            url="https://vmf.geo.tuwien.ac.at/trop_products/GRID/2.5x2/VMF1/VMF1_OP/"+str(year)+'/'+filename
            if os.path.exists(filename):
               os.remove(filename) # if exist, remove it directly
            wget.download(url,directory_vmf)
         except:
             pass
         
# =============================================================================                   

def read_vmf1_grid(mjd = None,ell = None): 
# =============================================================================
#     Function read_vmf1_grid searches for grid file close in time and computes
#     mapping coefficients ah, aw for selected station

#     Input parameters:
#     -------------------     
#     mjd   ... Modified Julian Date
#     ell   ... ellipsoidal station coordinates lat [degrees], lon [degrees] ,hell [m] 'Should be in list form e.g 'ell=[42.72,56.13,570]'' 
#     File Converted to python script by Zohreh Adavi zohreh.adavi@tuwien.ac.at, 2023-02-20
# =============================================================================
    
# Closest 6 hour period
    mjd=np.array(mjd)
    mjd_6h = np.round(mjd * 4) / 4

# List of relevant 6 hour intervals
    mjd_6h_list = np.unique(mjd_6h)

# Read vmf1 grid of period before and after current date
# (needed for temporal interpolation)
    
    if len(mjd_6h_list) == 1:
        mjd_6h_list = np.array([mjd_6h_list - 0.25,mjd_6h_list,mjd_6h_list + 0.25])
    
# Initialise variable
    ah_ell =[]
    aw_ell =[]
    
    zhd_ell =[]
    zwd_ell = []
    for i in range(0,len(mjd_6h_list)):
        # Convert mjd to date
        yr,mn,dy = mjd_to_ymd(np.rint(np.amin(mjd_6h_list[i])))
        # Hour of the day
        pwd=os.getcwd()
        directory_vmf_down=pwd+"\\Mapping_Fcn\\vmf1\\"
        hh = (mjd_6h_list[i]- np.fix(mjd_6h_list[i])) * 24
        if int(hh)==-6:

         #download_VMF(int(yr),int(mn),int(dy-1),directory_vmf_down)
         hour=18
         day_new=int(dy-1)
        else:
         #download_VMF(int(yr),int(mn),int(dy),directory_vmf_down)
         hour=int(hh)
         day_new=int(dy)
                  
        if day_new==0:
           mn=mn-1
           _,day_new=calendar.monthrange(yr, mn)
        
        # Open grid file with mapping coefficients
        directory_vmf=pwd+"\\Mapping_Fcn\\vmf1\\VMFG_"+str(int(yr))+str(int(mn)).zfill(2)+str(day_new).zfill(2)+'.H'+str(hour).zfill(2)
        
        if not os.path.exists(directory_vmf):
           download_VMF(int(yr),int(mn),int(day_new),directory_vmf_down)
        
        file=open(directory_vmf)
        
        ###############################
        pos = 0
        line0='! Comment:'
        cur_line = file.readline()
        Columns=['lat', 'lon', 'ah', 'aw', 'zhd', 'zwd']
        while not cur_line.startswith(line0):
               pos = pos+1
               cur_line = file.readline()
        
        file.close()  
        file=open(directory_vmf)
        dat = pd.read_csv(file ,skiprows=[i for i in range(0,pos)],sep='\s+',skipfooter =0, skip_blank_lines=True,header = 0, index_col= False ,names=Columns)
        file.close()    

        
        # Restore grid data
        lat = np.array(dat.iloc[:,0])

        lon = np.array(dat.iloc[:,1])
        
        ah = np.array(dat.iloc[:,2])
        
        aw = np.array(dat.iloc[:,3])
        
        zhd = np.array(dat.iloc[:,4])
        
        zwd = np.array(dat.iloc[:,5])
        # Interpolate data to station coordinates
        ah_ell.append(griddata((lat,lon),ah,(ell[0],ell[1])))
        
        aw_ell.append(griddata((lat,lon),aw,(ell[0],ell[1])))
        
        zhd_ell.append(griddata((lat,lon),zhd,(ell[0],ell[1])))
        
        zwd_ell.append(griddata((lat,lon),zwd,(ell[0],ell[1])))
        
        
    mjdd=pd.DataFrame(mjd_6h_list)
    
    ah_ell_new=pd.DataFrame(ah_ell)
    aw_ell_new=pd.DataFrame(aw_ell)
    zhd_ell_new=pd.DataFrame(zhd_ell)
    zwd_ell_new=pd.DataFrame(zwd_ell)
    # Interpolate and store final coefficients to return
    ah_vec=np.interp(np.transpose(mjd),mjdd.iloc[:,0], ah_ell_new.iloc[:,0])
    aw_vec=np.interp(np.transpose(mjd),mjdd.iloc[:,0], aw_ell_new.iloc[:,0])
    
    zhd_vec=np.interp(np.transpose(mjd),mjdd.iloc[:,0], zhd_ell_new.iloc[:,0])
    zwd_vec=np.interp(np.transpose(mjd),mjdd.iloc[:,0], zwd_ell_new.iloc[:,0])
    
    return ah_vec,aw_vec,zhd_vec,zwd_vec
   

# =============================================================================
def vmf1_ht(dmjd = None,dlat = None,dlon=None,ht = None,zd = None,ah= None,aw= None): 
# =============================================================================
#     !!! This is the version with height correction !!!
#     !!! It has to be used with the grid !!!
#
#     This subroutine determines the VMF1 (Vienna Mapping Functions 1)
#     Reference: Boehm, J., B. Werl, H. Schuh (2006),
#     Troposphere mapping functions for GPS and very long baseline interferometry
#     from European Centre for Medium-Range Weather Forecasts operational analysis data,
#     J. Geoph. Res., Vol. 111, B02406, doi:10.1029/2005JB003629.
#
#     Please mind that the coefficients in this paper are wrong. The corrected version of
#     the paper can be found at:
#     http://ggosatm.hg.tuwien.ac.at/DOCS/PAPERS/2006Boehm_etal_VMF1.pdf
#
#     input data
#     ----------
#     ah:   hydrostatic coefficient a (http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/)
#     aw:   wet coefficient a         (http://ggosatm.hg.tuwien.ac.at/DELAY/GRID/)
#     dmjd: modified julian date
#     dlat: ellipsoidal latitude in radians
#     dlon: ellipsoidal longitude in radians
#     ht:   ellipsoidal height in meter
#     zd:   zenith distance in radians
#
#     output data
#     -----------
#     vmf1h: hydrostatic mapping function
#     vmf1w: wet mapping function
#
#     Johannes Boehm, 2005 October 2
#      Rev 2011 July 21: latitude -> ellipsoidal latitude

#     implicit double precision (a-h,o-z)
    
#     pi = 3.14159265359d0
    
#     reference day is 28 January
#     this is taken from Niell (1996) to be consistent
#     File Converted to python script by Zohreh Adavi zohreh.adavi@tuwien.ac.at, 2023-02-20
# =============================================================================
    # ell=[np.degrees(dlat),np.degrees(dlon),ht]
    # ah,aw,zhd,zwd=read_vmf1_grid(dmjd,ell)
    doy = dmjd - 44239.0 + 1 - 28
    bh = 0.0029
    c0h = 0.062
    if (dlat < 0):
        phh = np.pi
        c11h = 0.007
        c10h = 0.002
    else:
        phh = 0.0
        c11h = 0.005
        c10h = 0.001
    
    ch = c0h + ((np.cos(doy / 365.25 * 2.0 * np.pi + phh) + 1.0) * c11h / 2.0 + c10h) * (1.0 - np.cos(dlat))
    sine = np.sin(np.pi / 2.0 - zd)
    beta = bh / (sine + ch)
    gamma = ah / (sine + beta)
    topcon = (1.0 + ah / (1.0 + bh / (1.0 + ch)))
    vmf1h = topcon / (sine + gamma)
    # C  height correction for hydrotatic part [Niell, 1996]
    a_ht = 2.53e-05
    b_ht = 0.00549
    c_ht = 0.00114
    hs_km = ht / 1000.0
    beta = b_ht / (sine + c_ht)
    gamma = a_ht / (sine + beta)
    topcon = (1.0 + a_ht / (1.0 + b_ht / (1.0 + c_ht)))
    ht_corr_coef = 1.0 / sine - topcon / (sine + gamma)
    ht_corr = ht_corr_coef*hs_km
    vmf1h = vmf1h + ht_corr
    bw = 0.00146
    cw = 0.04391
    beta = bw / (sine + cw)
    gamma = aw / (sine + beta)
    topcon = (1.0 + aw / (1.0 + bw / (1.0 + cw)))
    vmf1w = topcon / (sine + gamma)
    return vmf1h,vmf1w

# =============================================================================
def vmf1(dmjd = None,dlat = None,dlon=None,ht = None,zd = None,ah= None,aw= None):
# =============================================================================
#   This subroutine determines the VMF1 (Vienna Mapping Functions 1) for specific sites.
#   Reference: Boehm, J., B. Werl, H. Schuh (2006),
#   Troposphere mapping functions for GPS and very long baseline interferometry
#   from European Centre for Medium-Range Weather Forecasts operational analysis data,
#   J. Geoph. Res., Vol. 111, B02406, doi:10.1029/2005JB003629.
#
#   Please mind that the coefficients in this paper are wrong. The corrected version of
#   the paper can be found at:
#   http://vmf.geo.tuwien.ac.at/documentation/Boehm#20et#20al.,#202006a#20#20#20(VMF1).pdf
#
#
#    input data
#    ----------
#    ah:   hydrostatic coefficient a (http://vmf.geo.tuwien.ac.at/trop_products/)
#    aw:   wet coefficient a         (http://vmf.geo.tuwien.ac.at/trop_products/)
#    dmjd: modified julian date
#    dlat: ellipsoidal latitude in radians
#    dlon: ellipsoidal longitude in radians
#    ht:   ellipsoidal height in meter
#    zd:   zenith distance in radians
#
#    output data
#    -----------
#    vmf1h: hydrostatic mapping function
#    vmf1w: wet mapping function
#
#    Johannes Boehm, 2005 October 2
#    rev 2011 July 21: latitude -> ellipsoidal latitude
#
 
#    implicit double precision (a-h,o-z)    
#    pi = 3.14159265359e0  
#    reference day is 28 January
#    this is taken from Niell (1996) to be consistent
#    File Converted to python script by Zohreh Adavi zohreh.adavi@tuwien.ac.at, 2023-02-20
# =============================================================================
    # ell=[np.degrees(dlat),np.degrees(dlon),ht]
    # ah,aw,zhd,zwd=read_vmf1_grid(dmjd,ell)   
    
    doy = dmjd - 44239.0 + 1 - 28
    bh = 0.0029
    c0h = 0.062
    if (dlat < 0):
        phh = np.pi
        c11h = 0.007
        c10h = 0.002
    else:
        phh = 0.0
        c11h = 0.005
        c10h = 0.001
    
    ch = c0h + ((np.cos(doy / 365.25 * 2.0 * np.pi + phh) + 1.0) * c11h / 2.0 + c10h) * (1.0 - np.cos(dlat))
    sine = np.sin(np.pi / 2.0 - zd)
    beta = bh / (sine + ch)
    gamma = ah / (sine + beta)
    topcon = (1.0 + ah / (1.0 + bh / (1.0 + ch)))
    vmf1h = topcon / (sine + gamma)
    bw = 0.00146
    cw = 0.04391
    beta = bw / (sine + cw)
    gamma = aw / (sine + beta)
    topcon = (1.0 + aw / (1.0 + bw / (1.0 + cw)))
    vmf1w = topcon / (sine + gamma)
    return vmf1h,vmf1w