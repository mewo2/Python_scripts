from __future__ import division

import os
import time as time # for reading in a timer
import numpy as np # maths functions (arrays etc.)
import math
from matplotlib import pyplot as plt # for ploting
from scipy import signal # for convolution function
from scipy import ndimage # for resampling image
from scipy.fftpack import fft2
from matplotlib import cm # colour mapping
import matplotlib as matplotlib
import matplotlib.pyplot as plt
import copy as cp
from numpy import *
import datetime
import random
from scipy import ndimage
from scipy import misc
from scipy.stats import rankdata
from time import gmtime, strftime
import pylab as pl
from glob import glob
from osgeo import gdal, gdalconst # for reading in raster
from osgeo.gdalconst import * # for reading in raster
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import sys
sys.path.insert(0, 'C:/Users/ggwillc/Desktop/Python_scripts/functions')
import FFT_functions 
reload(FFT_functions)
import FFT_filter_functions 
reload(FFT_filter_functions)
import raster_functions 
reload(raster_functions)
import Filter_functions 
reload(Filter_functions)

startTime = time.time()

# Register driver
#gdal.AllRegister() #<-- useful ont *ly if reading in 
driver = gdal.GetDriverByName('ENVI') ## http://www.gdal.org/formats_list.html
driver.Register()

#for file_name in glob("helheim_222a_sample_MAX_ROI_2_butterworth_1_low_pass_fft_output_butterworth_1_filter_50_percent.bin"):
#for file_name in glob("helheim_222a_sample_ELEVATION_ROI_2_butterworth_1_low_pass_fft_main_butterworth_1_filter_50_percent.bin"):
#for file_name in glob("*ELEV*output*butter*1*50_percent.bin"):
#for file_name in glob("*MAX*output*butter*1*50_percent.bin"):

roi_number = 4
#file = "helheim_222a_sample_ELEVATION_ROI_%i_no_filtering_fft_output_no_filtering.bin" %roi_number ### MAGNITUDE SURFACE FROM ALs SCRIPT
#file = "helheim_222a_sample_ELEVATION_ROI_%i_no_filtering_fft_main_no_filtering.bin" %roi_number  ### ELEVATION SURFACE
file = "222a.helheim_post_0.5m.bin"
for file_name in glob(file):

	snip_file_name = file_name.split('.')[0]

	full_snip_1 = file_name.split('_sample')[0]
	print full_snip_1
		
	temp_snip_2 = file_name.split('_pass')[0]
	
	if "MAX" in temp_snip_2:
		full_snip_2 = temp_snip_2.split('MAX_')[1]
		print full_snip_2
	elif "ELEVATION" in temp_snip_2:
		full_snip_2 = temp_snip_2.split('ELEVATION_')[1]
		print full_snip_2
	else:
		"Name snipping error!!"
		os._exit(1)
	
	inds, cols, rows, bands, originX, originY, pixelWidth, pixelHeight, image_array, image_array_name = raster_functions.ENVI_raster_binary_to_2d_array(file_name)	
	
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	print 'CHECK OUTPUT DIRECTORY EXISTS'
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

	opath = r'C:\Users\ggwillc\Desktop\FFT\FFT_OUTPUTS_AL_ROUTINE'

	if os.path.isdir(opath):
		print "output_path exists"	
	else:
		print "output_path DOESN'T exist...\n"
		os.makedirs(opath)
		print "...but it does now"		
		
	### CALCULATE FFT USING SCIPY
	
	fft_python_complex = np.fft.fft2(image_array)
	fft_python_mag = abs(fft_python_complex)
	print "fft_python_complex data type:"
	print fft_python_complex.dtype
	print "fft_python_mag data type:"
	print fft_python_mag.dtype
	print "max: %f" %fft_python_mag.max()
	
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~'
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~'
	print 'Processing starts here....'
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~'
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~'

	## Get rolled FFT image object
	FFT_surface = fft_python_mag
	freq = 300
	post = 0.5
	input_x, input_y, magnitude = FFT_functions.magnitude_2D_RETURN(FFT_surface, freq)
		
	#FFT_functions.plot_FFT_2D_interactive_z_AXIS_FREQ(FFT_surface, opath, "current_input", 300)
	brown_noise = FFT_functions.brown_noise_surface(FFT_surface, 300)
	#FFT_functions.plot_brown_noise(brown_noise, 300)
			
	#~~~~~~~~~~
	
	## Image trimming variables:
	## Cuts off 10 pixels at the end of the y and x axis to account for the Gibbs effect
	## Also ignores a cut off zone as defined around the centre of the full rolled FFT image - currently ASSUMES A SQUARE REGION
	gibbs_allowance= 20 # pixels
	cutoff_zone = 40 # pixels
	mod_length = freq - gibbs_allowance
	
	## Subsample quadrants (require top and bottom left - 1 and 3)
	## Test which quadrant has greater maximum value
	
	#######################
	test_for_quadrant_max = 1
	if(test_for_quadrant_max == 1):
		Q1_title = "Q1 (top left)"
		Q1 = FFT_filter_functions.subsample_quadrant(np.log(brown_noise),0,300,0,300,opath,Q1_title )
		GIBBS_TRIMMED_Q1 = Q1[0:mod_length,0:mod_length].copy()
		ny, nx = Q1.shape
		filtered_quart_GIBBS_TRIMMED_Q1 = FFT_filter_functions.FFT_gaussian_filter_values_SIMPLE(post, nx, ny, GIBBS_TRIMMED_Q1, opath, 25)
		x_axis_limit_Q1 = len(filtered_quart_GIBBS_TRIMMED_Q1) - cutoff_zone
		y_axis_limit_Q1 = len(filtered_quart_GIBBS_TRIMMED_Q1) - cutoff_zone 
		gibbs_trimmed_and_cutoff_area_max_Q1 = filtered_quart_GIBBS_TRIMMED_Q1[0:y_axis_limit_Q1,0:x_axis_limit_Q1].max()
	
		Q3_title = "Q3 (bottom left)"
		Q3 = FFT_filter_functions.subsample_quadrant(np.log(brown_noise),300,600,0,300,opath,Q3_title)
		GIBBS_TRIMMED_Q3 = Q3[gibbs_allowance:freq,gibbs_allowance:freq].copy()
		ny, nx = Q3.shape
		filtered_quart_GIBBS_TRIMMED_Q3 = FFT_filter_functions.FFT_gaussian_filter_values_SIMPLE(post, nx, ny, GIBBS_TRIMMED_Q3, opath, 25)
		x_axis_limit_Q3 = len(filtered_quart_GIBBS_TRIMMED_Q3) - cutoff_zone
		gibbs_trimmed_and_cutoff_area_max_Q3 = filtered_quart_GIBBS_TRIMMED_Q3[cutoff_zone:len(filtered_quart_GIBBS_TRIMMED_Q3),0:x_axis_limit_Q3].max()
	
		# corner = 3 # 1 = top left, 3 = bottom left
		if(gibbs_trimmed_and_cutoff_area_max_Q1 > gibbs_trimmed_and_cutoff_area_max_Q3):
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			print "Top left corner will be used"
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			corner = 1
		elif(gibbs_trimmed_and_cutoff_area_max_Q1 < gibbs_trimmed_and_cutoff_area_max_Q3):
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			print "Bottom left corner will be used"
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			corner = 3
		else:
			print "Both quarters have same maximum value"
			print "Program now closing"
			os._exit(1)
		
		# Define mesh coordinate arrays
		if(corner == 1):
			title = "Q1 (top left)"
		#	Q1 = FFT_filter_functions.subsample_quadrant(np.log(brown_noise),0,300,0,300,opath,title )
			ny, nx = Q1.shape
		elif(corner == 3):
			title = "Q3 (bottom left)"
		#	Q3 = FFT_filter_functions.subsample_quadrant(np.log(brown_noise),300,600,0,300,opath,title)
			ny, nx = Q3.shape
			
		x = np.linspace(0,nx-1,nx)
		y = np.linspace(0,ny-1,ny)
			
		if(corner == 1):
			
			for i in range(len(filtered_quart_GIBBS_TRIMMED_Q1)):
				for j in range(len(filtered_quart_GIBBS_TRIMMED_Q1[i])):
			
					if(i > (0 + cutoff_zone)) and (j < (300 - cutoff_zone)):
						if(filtered_quart_GIBBS_TRIMMED_Q1[i,j] == filtered_quart_GIBBS_TRIMMED_Q1.max()):
							print filtered_quart_GIBBS_TRIMMED_Q1[i,j]
							print "y = %.2f" %y[i]
							filtered_quart_GIBBS_TRIMMED_max_value_y = y[i]
							print "x = %.2f" %x[j]
							filtered_quart_GIBBS_TRIMMED_max_value_x = x[j]
							max_position_orientation = FFT_filter_functions.FFT_max_value_POINT_bearings_NORTH_FRQ(filtered_quart_GIBBS_TRIMMED_max_value_y, filtered_quart_GIBBS_TRIMMED_max_value_x, freq, 1)
							dist_frq = FFT_filter_functions.peak_distance(filtered_quart_GIBBS_TRIMMED_max_value_x,filtered_quart_GIBBS_TRIMMED_max_value_y,input_x,input_y,freq,corner)
							print "max position distance (frq): %f" %dist_frq
							#print "max position distance (px): %f" %dist_px
							#print "max position distance (m): %f" %(dist_px*post)
							print "max position orientation: %f degreesN" %(max_position_orientation)			
							print "max position orientation (rotated -90degN): %f degreesN" %(max_position_orientation - 90.)	
					else:
						#print "Avoiding central cutoff zone...."
						#print "[%f.0,%f.0]" %(j,i)
						pass
			
			plt.imshow(np.log(filtered_quart_GIBBS_TRIMMED_Q1)), plt.colorbar()
			plt.gca().format_coord = FFT_functions.formatter_FRQ(filtered_quart_GIBBS_TRIMMED_Q1, freq)
			#plt.show()		
			
						
		elif(corner == 3):
			
			for i in range(len(filtered_quart_GIBBS_TRIMMED_Q3)):
				for j in range(len(filtered_quart_GIBBS_TRIMMED_Q3[i])):
				
					if(i < (cutoff_zone)) and (j > (300 - cutoff_zone)):
						#print "[%i,%i]" %(i,j)
						pass
					else:
						if(filtered_quart_GIBBS_TRIMMED_Q3[i,j] == filtered_quart_GIBBS_TRIMMED_Q3.max()):
							print filtered_quart_GIBBS_TRIMMED_Q3[i,j]
							print "y = %.2f" %y[i]
							filtered_quart_GIBBS_TRIMMED_Q3_max_value_y = y[i]
							print "x = %.2f" %x[j]
							filtered_quart_GIBBS_TRIMMED_Q3_max_value_x = x[j]
							max_position_orientation = FFT_filter_functions.FFT_max_value_POINT_bearings_NORTH_FRQ(filtered_quart_GIBBS_TRIMMED_Q3_max_value_y, filtered_quart_GIBBS_TRIMMED_Q3_max_value_x, freq, corner)
							dist_frq = FFT_filter_functions.peak_distance(filtered_quart_GIBBS_TRIMMED_Q3_max_value_x,filtered_quart_GIBBS_TRIMMED_Q3_max_value_y,input_x,input_y,freq,corner)
							
							print "max position distance (frq): %f" %dist_frq
							#print "max position distance (px): %f" %dist_px
							#print "max position distance (m): %f" %(dist_px*post)
							print "max position orientation: %f degreesN" %(max_position_orientation)			
							print "max position orientation (rotated -90degN): %f degreesN" %(max_position_orientation - 90.)	
						
			plt.imshow(np.log(filtered_quart_GIBBS_TRIMMED_Q3)), plt.colorbar()
			plt.gca().format_coord = FFT_functions.formatter_FRQ(filtered_quart_GIBBS_TRIMMED_Q3, freq)
			#plt.show()		
			
endTime = time.time()
print "Script was running for: %.2f seconds to run " %(endTime - startTime)