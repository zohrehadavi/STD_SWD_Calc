# -*- coding: utf-8 -*-
"""
@author: zadavi
"""
import pandas as pd
import numpy as np
import time
import pyproj
import linecache
from astropy.time import Time

transformer = pyproj.Transformer.from_crs(
    {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
    {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
    )

def Readwrite_TRO_ZTD(file_source=None,file_name=None,file_dest=None,fle_out=None,doy=None,int_need=None,New_rate=None):
# -------------Input Parameters
# file_source: The directory where TRO files are located
# file_name  : The full name of TRO file
# file_dest  : Where the created file should be saved
# fle_out    : The name of output file
# doy        : Day of year
# int_need   : 1 if interpolation is needed and 0 otherwise (default=0)
# New_rate   : New interpolation interval [second] (default=150 s #2.5 minutes)

# File created by Zohreh Adavi zohreh.adavi@tuwien.ac.at, February 2023
# =============================================================================    

####SAMPLING TROP
 line0='SAMPLING TROP'
 directory_trp=file_source+'/'+file_name
 with open(directory_trp) as f:
          pos = 0
          cur_line = f.readline()
          index_pos=cur_line.find(line0)
          while (index_pos<0):
                pos = pos+1
                cur_line = f.readline()
                index_pos=cur_line.find(line0)
 f.close()
 

 particular_line = linecache.getline(directory_trp, pos+1)
 sample_tro=int(particular_line[50:53])
  
####Reading CRD part
 Columnsmain0=['SITE','lon','lat','h','X','Y','Z']
 Columns_crd=['SITE','PT','SOLN','T','STA_X','STA_Y','STA_Z','SYSTEM','REMRK']
 line0='+TROP/STA_COORDINATES'
 appended_Data0 = []  
 directory_trp=file_source+'/'+file_name
 with open(directory_trp) as f:
          pos = 0
          cur_line = f.readline()
          index_pos=cur_line.find(line0)
          
          while index_pos<0:
                pos = pos+1
                cur_line = f.readline()
                index_pos=cur_line.find(line0)

 f.close()
 
 ##
 line0='-TROP/STA_COORDINATES'
 with open(directory_trp) as f:
          pose = 0
          cur_line = f.readline()
          index_pos=cur_line.find(line0)
          
          while index_pos<0:
                pose = pose+1
                cur_line = f.readline()
                index_pos=cur_line.find(line0)

 f.close()
 
 
 ##
 
 fle=open(directory_trp)
 file_crd = fle.readlines()[pos+2:pose] 
 fle.close()  

 df = pd.DataFrame({'Sta_inf': file_crd})
 df[['SITE','PT','SOLN','T','STA_X','STA_Y','STA_Z','SYSTEM','REMRK']] = df['Sta_inf'].str.split(expand=True)
 df = df.drop('Sta_inf', axis=1)
    
 Mainfle=pd.DataFrame(np.zeros((len(df),len(Columnsmain0))),columns=Columnsmain0) 
 lon1, lat1, alt1 = transformer.transform(df['STA_X'],df['STA_Y'],df['STA_Z'],radians=False) 
 Mainfle.iloc[:,0]=df['SITE']
 Mainfle.iloc[:,1]=lon1
 Mainfle.iloc[:,2]=lat1
 Mainfle.iloc[:,3]=alt1
 Mainfle.iloc[:,4]=df['STA_X']
 Mainfle.iloc[:,5]=df['STA_Y']
 Mainfle.iloc[:,6]=df['STA_Z']
 
 appended_Data0.append(Mainfle)
 del df
 del fle
 del Mainfle
 finaldf0 = pd.concat(appended_Data0)
####Reading TRP part
 Columnsmain=['SITE','Year','DoY','hour','minute','second','ZTD','GN','GE']
 Columns_trp=['SITE','EPOCH','TROTOT','STDDEVt','TGNTOT','STDDEVn','TGETOT','STDDEVe']
 line='+TROP/SOLUTION' 
 directory_trp=file_source+'/'+file_name
 with open(directory_trp) as f:
          pos = 0
          cur_line = f.readline()
         
          while not cur_line.startswith(line):
                pos = pos+1
                cur_line = f.readline()
          f.seek(pos)
          f.close() 
 fle=open(directory_trp)                      
 df = pd.read_csv(fle ,skiprows=[i for i in range(0,pos+1)],sep='\s+',skipfooter = 2, skip_blank_lines=True,header = 0, index_col= False ,names=Columns_trp)
 fle.close()     
 
 
 
 
 
 
 EPOCHall=df['EPOCH'].str[7:12].astype(int)
 DoY=df['EPOCH'].str[3:6].astype(int)
 Year=df['EPOCH'].str[0:2]
 
 finaldf1=pd.DataFrame(np.zeros((0,len(Columnsmain))),columns=Columnsmain)  
 jt=0
 for it in range(0,len(EPOCHall)):
     
     if DoY[it]==doy:
        convert = time.strftime("%H:%M:%S", time.gmtime(EPOCHall[it]))
     

        finaldf1.loc[jt,'SITE']=df.loc[it,'SITE']
        finaldf1.loc[jt,'Year']=int(str(20)+Year[it])
        finaldf1.loc[jt,'DoY']=DoY[it]
        finaldf1.loc[jt,'hour']=convert[0:2]
        finaldf1.loc[jt,'minute']=convert[3:5]
        finaldf1.loc[jt,'second']=convert[6:8]
        finaldf1.loc[jt,'ZTD']=df.loc[it,'TROTOT']
        finaldf1.loc[jt,'GN']=df.loc[it,'TGETOT']
        finaldf1.loc[jt,'GE']=df.loc[it,'TGNTOT']
        jt=jt+1
 
 
 del df
 del fle
########Concating CRD to TRP File
 year=int(finaldf1['Year'].unique())
 Columnsmain_new=['SITE','lon','lat','h','X','Y','Z','hour','minute','second','ZTD','GE','GN']
 
 appended_Data = [] 
 for it in range(0,len(finaldf0)):
    jloc0=finaldf1.loc[(finaldf1['SITE'] ==finaldf0.iloc[it,0])]
    jloc0.reset_index(drop=True, inplace=True)
    Mainfle=pd.DataFrame(np.zeros((len(jloc0),len(Columnsmain_new))),columns=Columnsmain_new)
    Mainfle['SITE']= jloc0['SITE']
    Mainfle['lon']= finaldf0.iloc[it,1]
    Mainfle['lat']= finaldf0.iloc[it,2]  
    Mainfle['h']= finaldf0.iloc[it,3]
    Mainfle['X']= finaldf0.iloc[it,4]
    Mainfle['Y']= finaldf0.iloc[it,5]
    Mainfle['Z']= finaldf0.iloc[it,6]
    
    Mainfle['hour']= jloc0['hour']
    Mainfle['minute']= jloc0['minute']
    Mainfle['second']= jloc0['second']
    Mainfle['ZTD']= jloc0['ZTD']
    Mainfle['GE']= jloc0['GE']
    Mainfle['GN']= jloc0['GN']

    Mainfle.reset_index(drop=True, inplace=True)
    appended_Data.append(Mainfle)
    del jloc0,Mainfle

 finaldf_main = pd.concat(appended_Data)
 
########Interpolation if needed
  #defult: 2.5 minutes
  
 Int_fle=pd.DataFrame(np.zeros((0,len(Columnsmain_new))),columns=Columnsmain_new) 
 nt=0 
 Sta_name=finaldf_main['SITE'].unique()

#######Interpolation of ZTD, GN, and GE
 for it in range(0,len(Sta_name)):
         jloc0=finaldf_main.loc[(finaldf_main['SITE']==Sta_name[it])]
         jloc0.reset_index(drop=True, inplace=True)
         MJDd=[]
         interval_data=[]
         for jt in range(0,len(jloc0)):
             t = Time('{}:{}:{}:{}'.format(year, doy, int(jloc0.loc[jt,'hour']), int(jloc0.loc[jt,'minute']), int(jloc0.loc[jt,'second'])), format='yday', scale='utc')
             MJDd.append(t.mjd)
             interval_data.append((int(jloc0.loc[jt,'hour'])*3600)+(int(jloc0.loc[jt,'minute'])*60)+(int(jloc0.loc[jt,'second'])))
             
         
         intervald=np.array(interval_data)
         dmjd= pd.DataFrame({'values': MJDd})
         del MJDd
         del interval_data
         
         #mjd for 2.5 minutes
         
         convert_par=(New_rate/(3600*24))
         mjd_new=pd.DataFrame(np.arange(float(dmjd.min()),float(dmjd.max()),convert_par,dtype=np.float64))
         interval = range(intervald.min(),intervald.max()+New_rate,New_rate)
         
         ZTD_interp=np.interp(mjd_new.iloc[:,0],dmjd.iloc[:,0], jloc0['ZTD'])
         GE_interp =np.interp(mjd_new.iloc[:,0],dmjd.iloc[:,0], jloc0['GE'])
         GN_interp =np.interp(mjd_new.iloc[:,0],dmjd.iloc[:,0], jloc0['GN'])
         
         
         
         step_fix=int(sample_tro/New_rate)
         lt=0
         def_xs=[]
         for jt in range(0,len(ZTD_interp),step_fix):
             #if len(jloc0)>lt:
              def_xs.append(ZTD_interp[jt]-jloc0.loc[lt,'ZTD'])
              ZTD_interp[jt]=jloc0.loc[lt,'ZTD']
              GE_interp[jt]=jloc0.loc[lt,'GE']
              GN_interp[jt]=jloc0.loc[lt,'GN']
              lt=lt+1
             
         
         for kt in range(0,len(ZTD_interp)):
                convert = time.strftime("%H:%M:%S", time.gmtime(interval[kt]))
                Int_fle.loc[nt,'SITE']=Sta_name[it]
                Int_fle.loc[nt,'lon']=jloc0.iloc[0,1]
                Int_fle.loc[nt,'lat']=jloc0.iloc[0,2]
                Int_fle.loc[nt,'h']=jloc0.iloc[0,3]
                Int_fle.loc[nt,'X']=jloc0.iloc[0,4]
                Int_fle.loc[nt,'Y']=jloc0.iloc[0,5]
                Int_fle.loc[nt,'Z']=jloc0.iloc[0,6]
                
                Int_fle.loc[nt,'hour']=convert[0:2]
                Int_fle.loc[nt,'minute']=convert[3:5]
                Int_fle.loc[nt,'second']=convert[6:8]
                
                Int_fle.loc[nt,'ZTD']=ZTD_interp[kt]
                Int_fle.loc[nt,'GE']=GE_interp[kt]  
                Int_fle.loc[nt,'GN']=GN_interp[kt]  
             
                nt=nt+1
         del ZTD_interp,GE_interp,GN_interp

#########################################################
 if int_need>0:
    finaldf=Int_fle
 else:
    finaldf=finaldf_main
    
 output_path=file_dest+fle_out
 # finaldf.to_csv(output_path, columns = False, index = False)
 with open(output_path, 'a') as f:
    df_string = finaldf.to_string(header=False, index=False)
    f.write(df_string)
    
 return finaldf