#/usr/bin/env python

#
# Loading Python modules:
#
import ftplib
import numpy as np
import numpy 
from netCDF4 import Dataset
import csv
import glob 
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import os
from mpl_toolkits.basemap import Basemap
from scipy.ndimage.filters import gaussian_filter
from wrf import ALL_TIMES
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
cartopy_ylim, latlon_coords)
from wrf import getvar, interplevel, interpz3d
import matplotlib.colors as mcolors
import pygrib
import cfgrib
from datetime import datetime
from datetime import datetime, timedelta
from scipy.interpolate import griddata
import pickle


#
# Defining some functions:
#

#
# Defining function to obtain initial and end date for each cyclone:
#
def get_cyclone_dates(folder_wrf_path, cyc_num):

    # Listing sorted WRF files:
    file_name = folder_wrf_path+'medCYC_{:04d}_d02*_00_00'.format(cyc_num)
    files = np.sort(glob.glob(file_name))

    # Initialize variables to store the earliest and latest dates
    earliest_date = None
    latest_date = None

    # Loop through the filenames
    for file_name in files:
        # Extract the date part from the filename
        parts = file_name.split('_')
        date_str = parts[4] + ' ' + parts[5]  + parts[6] + parts[7] # Combine date and time parts
        # Convert the string to a datetime object
        date = datetime.strptime(date_str, '%Y-%m-%d %H%M%S')

        # Update earliest and latest dates
        if earliest_date is None or date < earliest_date:
            earliest_date = date
        if latest_date is None or date > latest_date:
            latest_date = date

    # Print the results
    print("Earliest date:", earliest_date)
    print("Latest date:", latest_date)

    return files, earliest_date, latest_date


#
# Defining function to create hourly dates from the first initial date to the final date:
#
def generate_dates(earliest_date, latest_date):
    # Define the start and end dates
    start_date = earliest_date
    end_date = latest_date

    # Initialize the current date as the start date
    current_date = start_date

    # List to hold all the dates
    date_list = []

    # Loop through from start date to end date
    while current_date <= end_date:
        date_list.append(current_date)
        # Increment the current date by one hour
        current_date += timedelta(hours=1)
    
    return(date_list)


#
# Defining function to download Copernicus data:
#
def downloading_copernicus_data(start_date,end_date,cyc_num):
    start_year    = f"{start_date.year:04d}"
    start_month   = f"{start_date.month:02d}"
    start_day     = f"{start_date.day:02d}"
    start_hour    = f"{start_date.hour:02d}"
    start_minute  = "00"
    start_seconds = "00"
    #
    end_year    = f"{end_date.year:04d}"
    end_month   = f"{end_date.month:02d}"
    end_day     = f"{end_date.day:02d}"
    end_hour    = f"{end_date.hour:02d}"
    end_minute  = "00"
    end_seconds = "00"
    #
    start_date_format = f"{start_year}-{start_month}-{start_day}T{start_hour}:{start_minute}:{start_seconds}"
    end_date_format = f"{end_year}-{end_month}-{end_day}T{end_hour}:{end_minute}:{end_seconds}"
    print("")
    print("Downloading data: Initial date: " +start_date_format)
    print("Downloading data: End date: " +end_date_format)
    filename = f"cmems_obs-wind_CYC_{cyc_num:04d}.nc"
    #
    # WIND_GLO_PHY_L4_MY_012 data is only available from 11 Jan 2007 to 21 Dec 2023
    #
    import datetime as datetime2
    if start_date < datetime2.datetime(2007, 1, 11) or end_date > datetime2.datetime(2023, 12, 22):    
        print("Copernincus Data is only available from 11 Jan 2007 onwards.")
        print("Copernicus Data for cyclone #"+str(cyc_num)+" will be not obtained.")
    else:
        #
        # Downloading data:
        #
        copernicusmarine.subset(
            dataset_id = "cmems_obs-wind_glo_phy_my_l4_0.125deg_PT1H",
            variables = ["eastward_wind","northward_wind"],
            start_datetime = start_date_format,
            end_datetime = end_date_format,
            minimum_longitude = -10,
            maximum_longitude = 40,
            minimum_latitude = 25,
            maximum_latitude = 46,
            output_filename = filename,
            force_download=True
        )   


# Function to download files from FTP server
def download_files(ftp, remote_path, local_path):
    try:
        os.makedirs(local_path, exist_ok=True)
        filenames = ftp.nlst(remote_path)
        for filename in filenames:
            local_filename = os.path.join(local_path, os.path.basename(filename))
            with open(local_filename, 'wb') as local_file:
                ftp.retrbinary(f'RETR {filename}', local_file.write)
        return True
    except Exception as e:
        with open(LOG_FILE, 'a') as log:
            log.write(f"Error downloading {remote_path}: {e}\n")
        return False


#
# Defining function to download EUMETSAT FTP data:
#
def downloading_eumetsat_ftp_data(start_date,end_date,cyc_num):
    start_year    = f"{start_date.year:04d}"
    start_month   = f"{start_date.month:02d}"
    start_day     = f"{start_date.day:02d}"
    start_hour    = f"{start_date.hour:02d}"
    start_minute  = "00"
    start_seconds = "00"
    #
    end_year    = f"{end_date.year:04d}"
    end_month   = f"{end_date.month:02d}"
    end_day     = f"{end_date.day:02d}"
    end_hour    = f"{end_date.hour:02d}"
    end_minute  = "00"
    end_seconds = "00"
    #
    start_YMD_date    = start_date.strftime("%Y-%m-%d")
    end_YMD_date      = f"{end_year}{end_month}{end_day}"
    start_date_format = f"{start_year}-{start_month}-{start_day}T{start_hour}:{start_minute}:{start_seconds}"
    end_date_format   = f"{end_year}-{end_month}-{end_day}T{end_hour}:{end_minute}:{end_seconds}"
    print("Downloading data: Initial date: " +start_date_format)
    print("Downloading data: End date: " +end_date_format)
    filename = f"cmems_obs-wind_CYC_{cyc_num:04d}.nc"
    #
    #
    import datetime as datetime2
    if start_date < datetime2.datetime(2010, 7, 1) or end_date > datetime2.datetime(2021, 12, 1):
        print("EUMETSAT Data is only available from 1 July 2010 onwards.")
        print("EUMETSAT Data for cyclone #"+str(cyc_num)+" will be not obtained.")
    else:
        #
        # Downloading data:
        #
        print("Downloading data ...")
        current_date = start_date
        while current_date <= end_date:
            day_of_the_year = current_date.timetuple().tm_yday
            year            = current_date.year
            month           = current_date.month
            day             = current_date.day
            print("day of the year and year-month-day: " + str(day_of_the_year) +", " +str(year)+"-"+str(month)+"-"+str(day))
            remote_path = f"{REMOTE_DIR}/{year}/{day_of_the_year:03d}"
            local_path = os.path.join(LOCAL_DIR, f"{year}/CYC_{cyc_num}/{day_of_the_year:03d}")
            with ftplib.FTP(FTP_SERVER) as ftp:
                ftp.login(USERNAME, PASSWORD)
                success = download_files(ftp, remote_path, local_path)

            with open(LOG_FILE, 'a') as log:
                if success:
                    log.write(f"Download for {current_date} (day {day_of_the_year}) completed successfully at {datetime.now()}\n")
                else:
                    log.write(f"Download for {current_date} (day {day_of_the_year}) failed at {datetime.now()}\n")

            #
            current_date += timedelta(days=1)





#
# Start Downloading Copernincus data for each cyclone:
#
#cyc_tot = 201    # Entire cyclones database
cyc_tot = 5

FTP_SERVER = "eftp.ifremer.fr"
USERNAME = "ig1fdaf"
PASSWORD = "CERSAT-Data"
REMOTE_DIR = "/provider/knmi/satellite/l2b/metop-a/ascat/12.5km/data"
LOCAL_DIR = "/hpcperm/bidc/Project_Spgrflao/KNMI_EUMETSAT_Obs/EUMETSAT_Obs"
LOG_FILE = os.path.join(LOCAL_DIR, "download.log")


for cyc_num in range(1,cyc_tot):
    print ('') 
    print('CYCLONE NUMBER: ', cyc_num)
    #LOG_FILE = os.path.join(LOCAL_DIR, "download_cyc_{:04d}.log".format(cyc_num)
 
    # Folder WRF path:
    folder_wrf_path = '/scratch/bidc/DCC/medTRACK_{:04d}/'.format(cyc_num)

    # Select range of times for each Cyclone simulation:
    [files, earliest_date,latest_date] = get_cyclone_dates(folder_wrf_path, cyc_num)

    # Generate Dates every hour:
    date_list = generate_dates(earliest_date,latest_date)
    #print(date_list) 
    # Donwloading Copernicus Data:
    downloading_eumetsat_ftp_data(earliest_date,latest_date,cyc_num)





#
# End Python Script
#
