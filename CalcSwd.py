# -*- coding: utf-8 -*-
"""
@author: zadavi
"""    


import numpy as np
import pandas as pd
from GPT3 import gpt3_5_fast
from VMF import read_vmf1_grid,vmf1_ht
from GMF_GradinetMF import gmf,gradient_MF
# File created by Zohreh Adavi zohreh.adavi@tuwien.ac.at, February 2023
# =============================================================================

def ZHD_Saastamoinen(mjd = None,h_ell = None,h_ortho=None,lat = None,lon = None):
    
#     Input parameters:
#     -------------------    
#     mjd:   modified Julian date (scalar, only one epoch per call is possible)
#     lat:   ellipsoidal latitude in radians [-pi/2:+pi/2] (vector)
#     lon:   longitude in radians [-pi:pi] or [0:2pi] (vector)
#     h_ell: ellipsoidal height in m (vector)
#     it:    case 1: no time variation but static quantities
#     case 0: with time variation (annual and semiannual terms)

    P,T,dT,Tm,e,ah,aw,la,undu,Gn_h,Ge_h,Gn_w,Ge_w= gpt3_5_fast(mjd,lat,lon,h_ell,it =0)
    ZHD=pd.DataFrame(np.zeros((len(P),len(P))))
    
    ZHD = (0.0022768*P) / (1 - (0.00266*np.cos(2*lat)) -(0.00000028*h_ortho))
    
    return ZHD
###############################################################################
# def ZHD_VMF(mjd = None,ell = None):
    
# #     Input parameters:
# #     -------------------     
# #     mjd   ... Modified Julian Date
# #     ell   ... ellipsoidal station coordinates lat [degrees], lon [degrees] ,hell [m] 'Should be in list form e.g 'ell=[42.72,56.13,570]'' 

#     ah,aw,ZHD,ZWD=read_vmf1_grid(mjd,ell)
    
#     return ZHD
    
###############################################################################
#ZHD calculated from GPT model and Saastamoinen Function zhd_GPT
#ZHD derived from VMF1 grid zhd_VMF
def CalcZWD(mjd = None,ZTD = None,h_ell= None,h_ortho=None,lat = None,lon = None,Type_zhd=None,zhd=None):
    
#     Input parameters:
#     -------------------    
#     mjd:   modified Julian date (scalar, only one epoch per call is possible)
#     lat:   ellipsoidal latitude in radians [-pi/2:+pi/2] (vector)
#     lon:   longitude in radians [-pi:pi] or [0:2pi] (vector)
#     h_ell: ellipsoidal height in m (vector)
#     ZTD: Zenith tropospheric delay
#     Type_zhd: 
#     zhd_GPT: ZHD calculated from GPT model and Saastamoinen 
#     zhd_VMF: ZHD derived from VMF1 grid

    if (Type_zhd=='zhd_GPT'):
        ZHD = ZHD_Saastamoinen(mjd,h_ell,h_ortho,lat,lon)
    else:
        ZHD=zhd
        
    
    ZWD = ZTD - ZHD
    return ZWD
###############################################################################
def CalcSWDgrad(mjd = None,ZTD = None,dlat = None,dlon = None,dhgt = None,h_ortho=None,gn = None,ge = None,AZ = None,ELV = None,Type_MF=None,Type_zhd=None,grad=None): 

    #     Input parameters:
#     -------------------    
#     mjd:   modified Julian date (scalar, only one epoch per call is possible)
#     lat:   ellipsoidal latitude in radians [-pi/2:+pi/2] (vector)
#     lon:   longitude in radians [-pi:pi] or [0:2pi] (vector)
#     dhgt: ellipsoidal height in m (vector)
#     ZWD: Zenith wet delay
#     gn & ge: Horizontal gradients
#     az & ELV: Azimuth and elevation angles
#     Type_MF:
#     GMF or VMF mapping function
#     Type_zhd: 
#     zhd_GPT: ZHD calculated from GPT model and Saastamoinen 
#     zhd_VMF: ZHD derived from VMF1 grid
    zhd=0
    if (Type_zhd=='zhd_VMF' or Type_MF=='VMF'):
        ell=[np.degrees(dlat),np.degrees(dlon),dhgt]
        ah,aw,zhd,zwd=read_vmf1_grid(mjd,ell)
        

    ZWD=CalcZWD(mjd ,ZTD,dhgt,h_ortho,dlat,dlon,Type_zhd,zhd)


    zd = np.radians(90 - ELV)    
    
    if (Type_MF=='GMF'):
       
       mf_h,mf_w = gmf(mjd,dlat,dlon,dhgt,zd)
    else:
       mf_h,mf_w = vmf1_ht(mjd,dlat,dlon,dhgt,zd,ah,aw)
    
    
    MF_G=gradient_MF(ELV)
    Grad_part=np.multiply(gn,np.cos(np.radians(AZ)))+np.multiply(ge,np.sin(np.radians(AZ)))
    if grad=='True':
       SWD=ZWD*np.array(mf_w)+np.multiply(np.array(MF_G),np.array(Grad_part))
    else:
       SWD=ZWD*np.array(mf_w)
    return SWD,ZWD,mf_w,MF_G
 ###############################################################################
 def CalcSTD(mjd = None,dlat = None,dlon = None,dhgt = None,h_ortho=None,AZ = None,ELV = None,Type_MF=None,Type_zhd=None): 

    #     Input parameters:
#     -------------------    
#     mjd:   modified Julian date (scalar, only one epoch per call is possible)
#     lat:   ellipsoidal latitude in radians [-pi/2:+pi/2] (vector)
#     lon:   longitude in radians [-pi:pi] or [0:2pi] (vector)
#     dhgt: ellipsoidal height in m (vector)
#     ZWD: Zenith wet delay
#     gn & ge: Horizontal gradients
#     az & ELV: Azimuth and elevation angles
#     Type_MF:
#     GMF or VMF mapping function
#     Type_zhd: 
#     zhd_GPT: ZHD calculated from GPT model and Saastamoinen 
#     zhd_VMF: ZHD derived from VMF1 grid
    zhd=0
    if (Type_zhd=='zhd_VMF' or Type_MF=='VMF'):
        ell=[np.degrees(dlat),np.degrees(dlon),dhgt]
        ah,aw,zhd,zwd=read_vmf1_grid(mjd,ell)
        



    zd = np.radians(90 - ELV)    
    
    if (Type_MF=='GMF'):
       
       mf_h,mf_w = gmf(mjd,dlat,dlon,dhgt,zd)
    else:
       mf_h,mf_w = vmf1_ht(mjd,dlat,dlon,dhgt,zd,ah,aw)
    
    

    STD=zwd*np.array(mf_w)+zhd*np.array(mf_h)
    ZTD=zhd+zwd
    return STD,ZTD