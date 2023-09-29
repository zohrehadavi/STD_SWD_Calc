# -*- coding: utf-8 -*-
"""
@author: zadavi
"""

from Read_Write_TRP import Readwrite_TRO_ZTD
from satellite_read_int_elv_az import Orbit_ELV_AZ
from datetime import date
from GPT3 import gpt3_5_fast
from VMF import read_vmf1_grid,vmf1_ht
from GMF_GradinetMF import gmf,gradient_MF
from CalcSwd import CalcSWDgrad

from satellite_read_int_elv_az import *
from Read_Write_TRP import *
from GPT3 import *
from VMF import *
from GMF_GradinetMF  import *
from CalcSwd import *

import pandas as pd
import numpy as np
from astropy.time import Time
import pygeodesy
from pygeodesy.ellipsoidalKarney import LatLon
ginterpolator = pygeodesy.GeoidKarney("./geoids/egm2008-5.pgm")


def write_Sinex_SWD(year=None,month=None,day=None,folder_orb=None,email=None,Name_Datasource_orb=None,Type_orb=None,Hour=0,folder_source_tro=None,file_name_tro=None,folder_dest_snx=None,fle_out_snx_ztd=None,fle_out_snx_swd=None,Type_MF=None,Type_zhd=None,elv=None,Type_data_source=None,New_rate=None,Sat_cons=None,int_need_TRO=None,grad=None):
#     Input parameters:
#     -------------------    
#    folder_orb: The saving path for orbit file
#    email: Registered email for CDDIS
#    Name_Datasource_orb: IGS, GFZ, CODE, GRG, ESA 'For ultra rapid before week 2038 please enter 'igv' if the 'Type_data_source' is not 'MGEX'
#    Type_orb: Final or Ultra
#    Hour: In the case of Ultra, the hour can be selected; otherwise, the orbit file for hour 0,6,12,18 will download and then merge them (it applies for ultra rapid file 'igv')
#    file_source_tro: The path of TRO files
#    file_name_tro: The name of TRO file
#    file_dest_snx: Destination of created ztd or swd SINEX files
#    fle_out_snx_ztd: File name for SINEX ztd
#    fle_out_snx_swd: File name for SINEX swd
#     Type_MF:
#     GMF or VMF1 mapping function
#     Type_zhd: 
#     zhd_GPT: ZHD calculated from GPT model and Saastamoinen 
#     zhd_VMF: ZHD derived from VMF1 grid
#    elv: Cut of angle (defult:5 [deg])
#    New_rate: new rate for interpolating data (ORB and ZTD) and creating finer dataset (defult: 150 [sec] #2.5 minutes)
#    Type_data_source: If you want to download from MGEX please set option as 'MGEX', otherwise it will download from the main directory of CDDIS
#    Sat_cons        : Which GNSS constellation should be used (IT must be consistent with the ORB data source (for instance ['G','R','E'])) 
#    int_need_TRO: 1 if interpolation is needed and 0 otherwise (default=0)
#    File created by Zohreh Adavi zohreh.adavi@tuwien.ac.at, February 2023
# =============================================================================

  columns_snx_swd=['SITE','lon', 'lat', 'h', 'Xsta', 'Ysta', 'Zsta', 'Xsat', 'Ysat', 'Zsat', 'PRN', 'AZ', 'ELV', 'hour', 'minute', 'second','mjd', 'SWD','ZTD']
  doy=date(year, month, day).timetuple().tm_yday
  
  TRO=Readwrite_TRO_ZTD(folder_source_tro,file_name_tro,folder_dest_snx,fle_out_snx_ztd,doy,int_need_TRO,New_rate)
  Sta_name=TRO['SITE'].unique()
  
  snx_swd=pd.DataFrame(np.zeros((0,len(columns_snx_swd))),columns=columns_snx_swd)
  df_final=[]
  
  for it in range(0,len(Sta_name)):
    jloc0=TRO.loc[(TRO['SITE']==Sta_name[it])]
    jloc0.reset_index(drop=True, inplace=True)
    
    lon=jloc0.iloc[0,1]
    lat=jloc0.iloc[0,2]
    h=jloc0.iloc[0,3]
    Xo=jloc0.iloc[0,4]
    Yo=jloc0.iloc[0,5]
    Zo=jloc0.iloc[0,6]
    

    
    
    single_position=LatLon(lat, lon)
    Ngeo= ginterpolator(single_position)
    h_ortho=h-Ngeo
    
    ORB=Orbit_ELV_AZ(year,month,day,folder_orb,email,Name_Datasource_orb,Hour,float(Xo),float(Yo),float(Zo),elv,Type_orb,Type_data_source,New_rate,Sat_cons)

    for jt in range(0,len(jloc0)):
      ZTD=jloc0.iloc[jt,10]/1000
      ge=jloc0.iloc[jt,11]/1000
      gn=jloc0.iloc[jt,12]/1000
        
      jloc1=ORB.loc[(ORB['hour'] ==jloc0.iloc[jt,7]) & (ORB['minute'] ==jloc0.iloc[jt,8]) & (ORB['second'] == jloc0.iloc[jt,9]) & ((ORB['ELV'])>elv)]
      jloc1.reset_index(drop=True, inplace=True)

      if len(jloc1)>0:
        t = Time('{}-{}-{}T{}:{}:{}'.format(year,month, day, int(jloc0.iloc[jt,7]), int(jloc0.iloc[jt,8]), int(jloc0.iloc[jt,9])), format='isot', scale='utc')
        mjd=t.mjd
         
        ELV=jloc1.loc[:,'ELV']
        AZ =jloc1.loc[:,'AZ']
        SWD,ZWD,mf_w,MF_G=CalcSWDgrad(mjd,float(ZTD),np.radians(float(lat)),np.radians(float(lon)),float(h),h_ortho,float(gn),float(ge),AZ,ELV,Type_MF,Type_zhd,grad)
        df_new=pd.DataFrame(np.zeros((len(ELV),len(columns_snx_swd))),columns=columns_snx_swd)  
         
        df_new.loc[:,'SITE']=jloc0.iloc[0,0]
        df_new.loc[:,'lon']=lon
        df_new.loc[:,'lat']=lat
        df_new.loc[:,'h']=h
         
        df_new.loc[:,'Xsta']=float(Xo)
        df_new.loc[:,'Ysta']=float(Yo)
        df_new.loc[:,'Zsta']=float(Zo)
         
        df_new.loc[:,'Xsat']=jloc1.loc[:,'X']
        df_new.loc[:,'Ysat']=jloc1.loc[:,'Y']
        df_new.loc[:,'Zsat']=jloc1.loc[:,'Z']
            
        df_new.loc[:,'PRN']=jloc1.loc[:,'SAT']
        df_new.loc[:,'AZ']=jloc1.loc[:,'AZ']
        df_new.loc[:,'ELV']=jloc1.loc[:,'ELV']
         
        df_new.loc[:,'hour']=jloc0.iloc[jt,7]
        df_new.loc[:,'minute']=jloc0.iloc[jt,8]
        df_new.loc[:,'second']=jloc0.iloc[jt,9]
        df_new.loc[:,'mjd']=mjd
        
        df_new.loc[:,'SWD']=np.transpose(SWD)
        df_new.loc[:,'ZTD']=float(ZTD)
        df_final.append(df_new) 
        del df_new
    del jloc1
  del jloc0,ORB
  snx_swd=pd.concat(df_final)     
#########################################################
  output_path=folder_dest_snx+'\\'+fle_out_snx_swd
# finaldf.to_csv(output_path, columns = False, index = False)
  with open(output_path, 'a') as f:
        df_string = snx_swd.to_string(header=False, index=False)
        f.write(df_string)
        
  return snx_swd
        