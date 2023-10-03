# STD_SWD_Calc
Python package for creating GNSS STD/SWD as well as modelled (GPT/VMF) ZTD
Install all required packages indicated in the requirements.
All required inputs shown in 'input_data.py', can be changed to the desired ones
You can use all of the subroutines individually, but to calculate the SWD SINEX file, a sample 'Main_program.py' can be run to see how it will work.


---------------#Requirements-------------

pandas      # pip install pandas

numpy       # pip install numpy

pyproj      # pip install pyproj

astropy     # pip install astropy

gnsscal     # pip install gnsscal

unlzw3      # pip install unlzw3

wget        # pip install wget

# Notes:

CDDIS provides the orbit file for download. Make an account there (https://urs.earthdata.nasa.gov/) if you don't have one and enter your registered email address in the 'email' field in the inputs part.


Sat_cons=['G','R'] -> Satellite constellation  Sat_cons=['G','R','E'],...



If you use this package on GitHub, please cite this paper as a reference:

Adavi, Z., Weber, R. & Rohm, W. Pre-analysis of GNSS tomography solution using the concept of spread of model resolution matrix. J Geod 96, 27 (2022). https://doi.org/10.1007/s00190-022-01620-1
