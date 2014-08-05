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
sys.path.insert(0, '/home/staff/ggwillc/Desktop/Python_scripts/functions')
#sys.path.insert(0, 'C:/Users/ggwillc/Desktop/Python_scripts-master/functions')
import FFT_functions as FFT_functions
import FFT_filter_functions as FFT_filter_functions
import raster_functions as raster_functions	

startTime = time.time()

# Register driver
#gdal.AllRegister() #<-- useful ont *ly if reading in 
driver = gdal.GetDriverByName('ENVI') ## http://www.gdal.org/formats_list.html
driver.Register()

#for file_name in glob("helheim_222a_sample_MAX_ROI_2_butterworth_1_low_pass_fft_output_butterworth_1_filter_50_percent.bin"):
#for file_name in glob("helheim_222a_sample_ELEVATION_ROI_2_butterworth_1_low_pass_fft_main_butterworth_1_filter_50_percent.bin"):
#for file_name in glob("*ELEV*output*butter*1*50_percent.bin"):
#for file_name in glob("*MAX*output*butter*1*50_percent.bin"):

#for file_name in glob("helheim_222a_sample_ELEVATION_ROI_2_tophat_high_pass_fft_output_tophat_filter_1_tenth_image.bin"):
for file_name in glob("C:\Users\ggwillc\Desktop\FFT\FFT_surfaces_binary\helheim_222a_sample_ELEVATION_ROI_2_tophat_high_pass_fft_output_tophat_filter_1_tenth_image.bin"):

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

	#opath = r'/home/staff/ggwillc/Desktop/FFT_image_tests/'
	#opath = r'/geog/data/sirius/epsilon/ggwillc/FFT_2D/Helheim/AL_FFT_outputs/small_roi/max_filter/'
	opath = r'C:\Users\ggwillc\Desktop\FFT\FFT_OUTPUTS'

	if os.path.isdir(opath):
		print "output_path exists"	
	else:
		print "output_path DOESN'T exist...\n"
		os.makedirs(opath)
		print "...but it does now"		

	print '~~~~~~~~~~~~~~~~~~~~~~~~~~'
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~'
	print 'Processing starts here....'
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~'
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~'

	## Get rolled FFT image object
	FFT_surface = image_array
	plot_title = "Rolled FFT"
	freq = 300
	post = 0.5
	
	input_x, input_y, magnitude = FFT_functions.magnitude_2D_RETURN(FFT_surface, freq)
	
	#~~~~~~
	plt.imshow(magnitude), plt.colorbar()
	plt.title("Magnitude surface - basic imshow output") 
	output_filename_temp = opath + 'MAGNITUDE_basic_imshow_output.png'  
	plt.savefig(output_filename_temp)
	plt.clf()
	
	os._exit(1)
	#~~~~~~
	
	FFT_functions.plot_FFT_2D_interactive_z_AXIS_FREQ(FFT_surface, opath, "current_input", 300)
	brown_noise = FFT_functions.brown_noise_surface(FFT_surface, 300)
	FFT_functions.plot_brown_noise(brown_noise, 300)
	
	#max_filtered_FFT = FFT_filter_functions.FFT_max_filter_values(freq, post, input_x, input_y, magnitude, 50)
	#FFT_filter_functions.FFT_max_filter_values(freq, post, input_x, input_y, brown_noise, 50)
		
	#max_filtered_FFT_brown_noise_flattened = FFT_filter_functions.FFT_max_filter_values_SIMPLE(freq, post, input_x, input_y, brown_noise, 50)
		
	max_filtered_FFT = FFT_filter_functions.FFT_max_filter_values_SIMPLE(freq, post, input_x, input_y, brown_noise, opath, 25)
	gibbs_removed_filtered_surface = FFT_filter_functions.Stifle_Gibbs_effect_on_FFT_max_filter(max_filtered_FFT, freq, opath)
	gibbs_removed_filtered_surface_NO_ZEROS = FFT_filter_functions.remove_zeros(gibbs_removed_filtered_surface, opath)
		
	ONLY_MAX_FFT = FFT_filter_functions.zero_all_but_max(gibbs_removed_filtered_surface, opath)
	'''
	xy_list_1, xy_list_2, xy_list_3, xy_list_4 = FFT_filter_functions.FFT_max_filtered_Max_positions_and_values(ONLY_MAX_FFT, freq, post, input_x, input_y, opath)
	FFT_filter_functions.FFT_max_streak_vector_calculation(xy_list_1, xy_list_2, xy_list_3, xy_list_4)
	'''
	
	aaron_rank = FFT_filter_functions.get_order_array_aaron(gibbs_removed_filtered_surface,opath)
	#aaron_rank = FFT_filter_functions.get_order_array_aaron(ONLY_MAX_FFT, opath)
	print aaron_rank.max()
	print aaron_rank.min()
	
	#FFT_filter_functions.return_value_of_given_rank_position_SINGLE(gibbs_removed_filtered_surface_NO_ZEROS, aaron_rank, 1)	
	#FFT_filter_functions.return_value_of_given_rank_position_MULTIPLE(gibbs_removed_filtered_surface_NO_ZEROS, aaron_rank, 10)	
	
	# Save ranked surfaces (individually set)
	
	rank_pos = 1022
	ranked_surface_value_int, ranked_surface_value_log = FFT_filter_functions.return_value_of_given_rank_position_SINGLE_IMAGE_RETURN_PART_1(gibbs_removed_filtered_surface_NO_ZEROS, aaron_rank, rank_pos)
	rank_output_SINGLE = FFT_filter_functions.zero_all_but_rank_value_IMAGE_RETURN_PART_2(gibbs_removed_filtered_surface_NO_ZEROS, opath, ranked_surface_value_log, rank_pos, 1)
	
	# Save ranked surfaces (calculated for multiple positions)
	'''
	rank_pos_max = 7
	rank_pos_min = 1
	for i in range(rank_pos_max):
		ranked_surface_value_int, ranked_surface_value_log = FFT_filter_functions.return_value_of_given_rank_position_SINGLE_IMAGE_RETURN_PART_1(gibbs_removed_filtered_surface_NO_ZEROS, aaron_rank, rank_pos_min)
		FFT_filter_functions.zero_all_but_rank_value_IMAGE_RETURN_PART_2(gibbs_removed_filtered_surface_NO_ZEROS, opath, ranked_surface_value_log, rank_pos_min, 1)
		rank_pos_min += 1
	'''
	
	#print 'Clear variables' 
	band = None
	image_array = None
	
