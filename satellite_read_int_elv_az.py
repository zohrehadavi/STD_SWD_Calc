# -*- coding: utf-8 -*-
"""
@author: zadavi
"""

#Name_Datasource: e.g: week<2238  igr & week >= 2238 GRG0OPSFIN_
from ftplib import FTP_TLS
import datetime
from datetime import date
import gnsscal
import gzip
import os
import unlzw3
from pathlib import Path
import csv
import time
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# File created by Zohreh Adavi zohreh.adavi@tuwien.ac.at, February 2023
# =============================================================================

def merge_files(files_to_merge, output_file):
    merged_data = pd.DataFrame()

# Loop through the files and merge them
    for file in files_to_merge:
    # Read the data from the file
       if file !=files_to_merge[len(files_to_merge)-1]:
           data = pd.read_csv(file,header = None, index_col= False, skip_blank_lines=True,skipfooter = 1)
       else:
           data = pd.read_csv(file,header = None, index_col= False, skip_blank_lines=True,skipfooter = 0)

    # Remove the header from all files except the first one
       if file != files_to_merge[0]:
           data0 = data.loc[22:,0]
       else:
           data0 = data
           
    
    # Merge the data into the merged_data dataframe
           merged_data = pd.concat([merged_data, data0])
           del data,data0
    merged_data.to_csv(output_file, header = None,index=False)
    
###############################################################################
               
def Download_Ultra_Orbit(year=None,month=None,day=None,folder=None,email=None,Name_Datasource=None,Hour=None):

    
    if (Name_Datasource=='igv'):
       try:
          fle_destinations=[]
          if Hour==None:
             Hour=[0,6,12,18,0]
             for it in range(0,len(Hour)):
                if it!=len(Hour)-1:
                   DoY=date(year, month, day).timetuple().tm_yday
                   dayt = datetime.date(year, month, day)
                   D=gnsscal.date2gpswd(dayt)
                   destination = folder +'/'+ "IGR" +str(D[0])+str(D[1])+'.SP3'
                else:
                   DoY=date(year, month, day+1).timetuple().tm_yday
                   dayt = datetime.date(year, month, day+1)
                   D=gnsscal.date2gpswd(dayt)
                   
                directory = "/glonass/products/"+str(D[0])
                file_name =Name_Datasource.lower()+str(D[0])+str(D[1])+'_'+str(Hour[it]).zfill(2)+'.sp3.Z'
                ftps = FTP_TLS(host = 'gdc.cddis.eosdis.nasa.gov')
                ftps.login(user='anonymous', passwd=email)
                ftps.prot_p()
                ftps.nlst()
                ftps.cwd(directory)
                ftps.retrbinary("RETR " + file_name, open(file_name, 'wb').write)
                ftps.quit()
             #####
                uncompressed_data = unlzw3.unlzw(Path(file_name))
                data = uncompressed_data.decode('utf-8').splitlines()
                with open(file_name[0:15], "w", newline="") as csv_file:
                 writer = csv.writer(csv_file)
                 for line in data:
                     # making cells by using , delimiters
                     writer.writerow(line.split(","))
             ######
                destination0 = folder +'/'+ "IGR" +str(D[0])+str(D[1])+'_'+str(Hour[it]).zfill(2)+'.SP3'
                source = folder +'/'+ file_name[0:15]
                os.rename(source, destination0)
                os.remove(file_name)
                fle_destinations.append(destination0)
             
            
             merge_files(fle_destinations, destination)
             #remove downloaded files
             for filename in fle_destinations:
       
                 file_path = filename
                 if os.path.isfile(file_path):
                         os.remove(file_path)
          else:
              
             DoY=date(year, month, day+1).timetuple().tm_yday
             dayt = datetime.date(year, month, day+1)
             D=gnsscal.date2gpswd(dayt)
             
             destination = folder +'/'+ "IGR" +str(D[0])+str(D[1])+'.SP3'
             directory = "/glonass/products/"+str(D[0])
             file_name =Name_Datasource.lower()+str(D[0])+str(D[1])+'_'+str(Hour).zfill(2)+'.sp3.Z'
             ftps = FTP_TLS(host = 'gdc.cddis.eosdis.nasa.gov')
             ftps.login(user='anonymous', passwd=email)
             ftps.prot_p()
             ftps.nlst()
             ftps.cwd(directory)
             ftps.retrbinary("RETR " + file_name, open(file_name, 'wb').write)
             ftps.quit()
          #####
             uncompressed_data = unlzw3.unlzw(Path(file_name))
             data = uncompressed_data.decode('utf-8').splitlines()
             with open(file_name[0:15], "w", newline="") as csv_file:
              writer = csv.writer(csv_file)
              for line in data:
                  # making cells by using , delimiters
                  writer.writerow(line.split(","))
          ######
             source = folder +'/'+ file_name[0:15]
             os.rename(source, destination)
             os.remove(file_name)
       except:
         pass

    else:
        try:
           directory = "/gnss/products/"+str(D[0])
         
           if Name_Datasource=='IGS':
               file_name=Name_Datasource+'0OPSULT_'+str(year)+str(DoY)+'0000_02D_15M_ORB.SP3.gz'
           else:
               file_name=Name_Datasource+'0OPSULT_'+str(year)+str(DoY)+'0000_02D_05M_ORB.SP3.gz'
           
           ftps = FTP_TLS(host = 'gdc.cddis.eosdis.nasa.gov')
           ftps.login(user='anonymous', passwd=email)
           ftps.prot_p()
           ftps.nlst()
           ftps.cwd(directory)
           ftps.retrbinary("RETR " + file_name, open(file_name, 'wb').write)
           ftps.quit()
           #####
           input = gzip.GzipFile(folder+'/'+file_name, 'rb')
           s= input.read()
           input.close()

           output = open(folder+'/'+file_name[0:38], 'wb')
           output.write(s)
           output.close()
           ######
           destination = folder + "IGR" +str(D[0])+str(D[1])+'.SP3'
           source = folder +'/'+ file_name[0:38]
           os.rename(source, destination)
           
        except:
          pass
    return destination
################################################################################
import os
import shutil

def move_and_rename_file(source_file, destination_folder, new_filename):
    # Move the file to the destination folder
    shutil.move(source_file, destination_folder)

    # Rename the file
    old_filename = os.path.basename(source_file)
    old_path = os.path.join(destination_folder, old_filename)
    new_path = os.path.join(destination_folder, new_filename)
    os.rename(old_path, new_path)
################################################################################
#Type_orb='Final' or 'Ultra'
def Int_Orbit(year=None,month=None,day=None,folder=None,email=None,Name_Datasource=None,Hour=None,Type_orb=None,Type_data_source=None,new_int_obs=None,Sat_cons=None):
 
 if Type_data_source=='MGEX':
      destination=Download_MGEX_Orbit(year,month,day,folder,email,Name_Datasource,Type_orb)
 else:
    if (Type_orb=='Final'):
       destination=Download_Final_Orbit(year,month,day,folder,email,Name_Datasource)
    else:
       destination=Download_Ultra_Orbit(year,month,day,folder,email,Name_Datasource,Hour)
                  
 Columnsmain0=['SAT','Year','Month','Day','hour','minute','second','X','Y','Z']
 line0='/* PCV:IGS'
 with open(destination) as f:
        cur_line = f.readline()
        cur_line = f.readline()
        int_obs=int(cur_line[25:29])
        cur_line = f.readline()
        NUMSAT=int(cur_line[4:7])
        f.close()
        
 with open(destination) as f:
        pos = 0
        cur_line = f.readline()
        
        while not cur_line.startswith(line0):
              pos = pos+1
              cur_line = f.readline()
        f.close()
        
 if ((Name_Datasource=='GFZ') or (Name_Datasource=='JAX')):
     pos=pos+3
 elif ((Name_Datasource=='WUM') and (Type_orb!='Final')):
    pos=pos+0
 elif ((Name_Datasource=='WUM') and (Type_orb=='Final')):
     pos=pos+2
 elif (Name_Datasource=='GRG'):
     pos=pos+2
 else:
     pos=pos+0
  
 if ((Name_Datasource=='igv') and (Name_Datasource=='CODE_orbits')):
    const_footer=0
 else:
    const_footer=1
    
 df = pd.read_csv(destination ,skiprows=[i for i in range(0,pos+1)],sep='\s+',skipfooter = const_footer, skip_blank_lines=True,header = None, index_col= False,names=['a','b','c','d','e','f','g','h','i','j','k'] )

 Mainfle0=pd.DataFrame(np.zeros((0,len(Columnsmain0))),columns=Columnsmain0)
 kt=0

 for it in range(0,len(df),NUMSAT+1):
    time_acq=df.iloc[it,0:7]
    for jt in range(it+1,it+NUMSAT+1):
        if jt<len(df):
           Mainfle0.loc[kt,'SAT']=df.iloc[jt,0][1:5]
           Mainfle0.loc[kt,'Year']=time_acq[1]
           Mainfle0.loc[kt,'Month']=time_acq[2]
           Mainfle0.loc[kt,'Day']=time_acq[3]
           Mainfle0.loc[kt,'hour']=time_acq[4]
           Mainfle0.loc[kt,'minute']=time_acq[5]
           Mainfle0.loc[kt,'second']=time_acq[6]
           Mainfle0.loc[kt,'X']=df.iloc[jt,1]
           Mainfle0.loc[kt,'Y']=df.iloc[jt,2]
           Mainfle0.loc[kt,'Z']=df.iloc[jt,3]
           kt=kt+1

 del df
 os.remove(destination) 

################ 
 Mainfle=Mainfle0.where(Mainfle0['Day'] == day)
 Mainfle = Mainfle.dropna()
 Mainfle.reset_index(drop = True, inplace = True)
 
 Mainfle.drop_duplicates(subset=['SAT','Day','hour','minute','second'], keep='last', inplace=True)
 Mainfle.reset_index(drop = True, inplace = True)
##########Interpolation for every 2.5 minutes
 end_time=int(int(Mainfle.iloc[-1,4])*3600+int(Mainfle.iloc[-1,5])*60+int(float(Mainfle.iloc[-1,6])))
 fitTime = range(0,end_time+int_obs,int_obs)



 sorted_Mainfle = Mainfle.sort_values(by=['SAT','hour','minute','second'], ascending=[True,True,True,True])
 sorted_Mainfle.reset_index(drop = True, inplace = True)

  #2.5 minutes
 interval = range(0,end_time+new_int_obs,new_int_obs)


 Int_fle=pd.DataFrame(np.zeros((0,len(Columnsmain0))),columns=Columnsmain0) 
 nt=0 
 mt=0 
 for kt in range(0,NUMSAT):
    
    df_new=sorted_Mainfle.loc[mt:len(fitTime)+mt-1,:]
    df_new.reset_index(drop = True, inplace = True)
    fitx=np.polyfit(fitTime,sorted_Mainfle.loc[mt:mt+len(fitTime)-1,'X'].astype('float'),deg=16)
    fity=np.polyfit(fitTime,sorted_Mainfle.loc[mt:mt+len(fitTime)-1,'Y'].astype('float'),deg=16)
    fitz=np.polyfit(fitTime,sorted_Mainfle.loc[mt:mt+len(fitTime)-1,'Z'].astype('float'),deg=16)
    
    x_interp = np.polyval(fitx, interval)
    y_interp = np.polyval(fity, interval)
    z_interp = np.polyval(fitz, interval)
    
    step_fix=int(int_obs/new_int_obs)
    lt=0
    def_xs=[]
    for jt in range(0,len(x_interp),step_fix):
        x_interp[jt]=df_new.loc[lt,'X']
        y_interp[jt]=df_new.loc[lt,'Y']
        z_interp[jt]=df_new.loc[lt,'Z']
        lt=lt+1
        
    for it in range(0,len(x_interp)):
           convert = time.strftime("%H:%M:%S", time.gmtime(interval[it]))
           Int_fle.loc[nt,'SAT']=df_new.iloc[0,0]
           Int_fle.loc[nt,'Year']=df_new.iloc[0,1]
           Int_fle.loc[nt,'Month']=df_new.iloc[0,2]
           Int_fle.loc[nt,'Day']=df_new.iloc[0,3]
           Int_fle.loc[nt,'hour']=convert[0:2]
           Int_fle.loc[nt,'minute']=convert[3:5]
           Int_fle.loc[nt,'second']=convert[6:8]
           Int_fle.loc[nt,'X']=x_interp[it]
           Int_fle.loc[nt,'Y']=y_interp[it]
           Int_fle.loc[nt,'Z']=z_interp[it]  
        
           nt=nt+1
    
    mt=mt+len(fitTime)
 ############### 
    Int_fle_sat=[]
    for jt in range(0,len(Sat_cons)):
     df_sat=Int_fle.where(Int_fle['SAT'].astype(str).str[0]==Sat_cons[jt])
     Int_fle_sat.append(df_sat)
     del df_sat
  
 Int_fle_final0=pd.concat(Int_fle_sat)
 Int_fle_final = Int_fle_final0.dropna()
 Int_fle_final.reset_index(drop = True, inplace = True) 
 del Mainfle,Mainfle0,sorted_Mainfle   
 return Int_fle_final

##############################################################################
import pyproj
transformer = pyproj.Transformer.from_crs(
    {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
    {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
    )

def Orbit_ELV_AZ(year=None,month=None,day=None,folder=None,email=None,Name_Datasource=None,Hour=None,Xo=None,Yo=None,Zo=None,elv=None,Type_orb=None,Type_data_source=None,new_int_obs=None,Sat_cons=None):

    Columnsmain=['SAT','Year','Month','Day','hour','minute','second','X','Y','Z','ELV','AZ']
    df0=Int_Orbit(year,month,day,folder,email,Name_Datasource,Hour,Type_orb,Type_data_source,new_int_obs,Sat_cons)
    lambdaO, phiO, h = transformer.transform(Xo,Yo,Zo,radians=True) 
    Xs=df0.loc[:,'X']*1000
    Ys=df0.loc[:,'Y']*1000
    Zs=df0.loc[:,'Z']*1000
    ## *******************Transformation To ECEF*****************************
    r = np.sqrt((Xs - Xo) ** 2 + (Ys - Yo) ** 2 + (Zs - Zo) ** 2)
    dx = (Xs - Xo) / r
    dy = (Ys - Yo) / r
    dz = (Zs - Zo) / r
    E = (np.multiply((- np.sin(lambdaO)),(dx))) + (np.multiply((np.cos(lambdaO)),(dy)))
    N = (np.multiply((np.multiply(- np.sin(phiO),np.cos(lambdaO))),(dx))) + (np.multiply(np.multiply(- np.sin(phiO),np.sin(lambdaO)),(dy))) + (np.multiply(np.cos(phiO),(dz)))
    U = (np.multiply((np.multiply(np.cos(phiO),np.cos(lambdaO))),(dx))) + (np.multiply((np.multiply(np.cos(phiO),np.sin(lambdaO))),(dy))) + (np.multiply((np.sin(phiO)),(dz)))
    ## ********************Elevation and Azimuth Computation*****************
    elev = np.degrees(np.arctan(U / (np.sqrt((E) ** 2 + (N) ** 2))))
    az = np.degrees(np.arctan2(E,N))
    
    df_new=pd.DataFrame(np.zeros((0,len(Columnsmain))),columns=Columnsmain) 

    
    for it in range(0,len(elev)):
        if (elev[it]>=elv) and (az[it]>=0):
            df_new.loc[it,'SAT']=df0.iloc[it,0]
            df_new.loc[it,'Year']=df0.iloc[it,1]
            df_new.loc[it,'Month']=df0.iloc[it,2]
            df_new.loc[it,'Day']=df0.iloc[it,3]
            df_new.loc[it,'hour']=df0.iloc[it,4]
            df_new.loc[it,'minute']=df0.iloc[it,5]
            df_new.loc[it,'second']=df0.iloc[it,6]
            df_new.loc[it,'X']=df0.iloc[it,7]
            df_new.loc[it,'Y']=df0.iloc[it,8]
            df_new.loc[it,'Z']=df0.iloc[it,9]
            df_new.loc[it,'ELV']=elev[it] 
            df_new.loc[it,'AZ']=az[it] 
            
        elif (elev[it]>=elv) and (az[it]<0):

           df_new.loc[it,'SAT']=df0.iloc[it,0]
           df_new.loc[it,'Year']=df0.iloc[it,1]
           df_new.loc[it,'Month']=df0.iloc[it,2]
           df_new.loc[it,'Day']=df0.iloc[it,3]
           df_new.loc[it,'hour']=df0.iloc[it,4]
           df_new.loc[it,'minute']=df0.iloc[it,5]
           df_new.loc[it,'second']=df0.iloc[it,6]
           df_new.loc[it,'X']=df0.iloc[it,7]
           df_new.loc[it,'Y']=df0.iloc[it,8]
           df_new.loc[it,'Z']=df0.iloc[it,9]
           df_new.loc[it,'ELV']=elev[it] 
           df_new.loc[it,'AZ']=az[it]+360
           
        else:
            df_new.loc[it,'SAT']=df0.iloc[it,0]
            df_new.loc[it,'Year']=df0.iloc[it,1]
            df_new.loc[it,'Month']=df0.iloc[it,2]
            df_new.loc[it,'Day']=df0.iloc[it,3]
            df_new.loc[it,'hour']=df0.iloc[it,4]
            df_new.loc[it,'minute']=df0.iloc[it,5]
            df_new.loc[it,'second']=df0.iloc[it,6]
            df_new.loc[it,'X']=df0.iloc[it,7]
            df_new.loc[it,'Y']=df0.iloc[it,8]
            df_new.loc[it,'Z']=df0.iloc[it,9]
            df_new.loc[it,'ELV']=0 
            df_new.loc[it,'AZ']=0
    
    
    return df_new