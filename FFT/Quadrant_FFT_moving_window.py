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
sys.path.insert(0, 'C:/Users/ggwillc/Documents/GitHub/Python_scripts/functions')
import FFT_functions 
reload(FFT_functions)
import FFT_filter_functions 
reload(FFT_filter_functions)
import raster_functions 
reload(raster_functions)
import Filter_functions 
reload(Filter_functions)

startTime = time.time()

driver = gdal.GetDriverByName('ENVI') ## http://www.gdal.org/formats_list.html
driver.Register()

file = "222a.helheim_post_0.5m_FFT_SUBSAMPLE.bin"
for file_name in glob(file):

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
	
	######## MOVING WINDOW STARTS HERE
	nywind,nxwind = image_array.shape
	kernel_size = 1001 ## FFT is faster for a window to the power of 2 OR a prime number
	stepsize = 200
	freq = int(nxwind/2)
	post = 0.5
	
	###########################
	## SET MOVING WINDOW UP  ##
	###########################

	if(kernel_size % 2 != 0 and kernel_size >= 3):
		window = np.ones([kernel_size,kernel_size])
	elif(kernel_size % 2 == 0 or kernel_size < 3):
		print "kernel is even - it needs to be odd and at least of a value of 3"
		os._exit(1)

	number_of_output_variables = 2 # spacing and orientation
	imgout = [np.zeros((1 + nxwind // stepsize, 1 + nywind // stepsize)) for i in xrange(number_of_output_variables)]
	window_border_indent = int((kernel_size - 1)/2 )
	window_border_indent_end_X = int(nxwind - window_border_indent)
	window_border_indent_end_Y = int(nywind - window_border_indent) 
	
	#print "window_border_indent: %f"  %window_border_indent
	#print "window_border_indent_end_X: %f"  %window_border_indent_end_X
	#print "window_border_indent_end_Y: %f"  %window_border_indent_end_Y
	
	window_iteration = 0
	for ii in range(window_border_indent, window_border_indent_end_X, stepsize):
		for jj in range(window_border_indent, window_border_indent_end_Y, stepsize):
			
			window_iteration+=1
			print "Loop iteration: %i" %window_iteration
						
			imin=max(0,ii-window_border_indent) # gets the maximum of either 0 or i-kernel_size/2...
			imax=min(nxwind-1,ii+window_border_indent)+1
			jmin=max(0,jj-window_border_indent)
			jmax=min(nywind-1,jj+window_border_indent)+1
			print imin, imax, jmin, jmax
		
			calc_moving_window_size_x = imax - imin 
			calc_moving_window_size_y = jmax- jmin

			# create subsample
			subsample_elevation_surface = image_array[imin:imax,jmin:jmax]
			subsample_nx,subsample_ny = subsample_elevation_surface.shape	
			
			fft_complex = np.fft.fft2(subsample_elevation_surface)
			FFT_surface = fft_mag = abs(fft_complex) ## magnitude surface
			
			##FFT_functions.plot_FFT_2D_interactive_z_AXIS_FREQ(FFT_surface, opath, "current_input", 300)
			nyqvist_frq = ((subsample_nx-1)//2) ## FFT image can't present more frequencies than are equal to half the image size
			print nyqvist_frq
			brown_noise = FFT_functions.brown_noise_surface(FFT_surface, nyqvist_frq)
			#FFT_functions.plot_brown_noise(brown_noise, 300)
	
			## Image trimming variables:
			## Cuts off 10 pixels at the end of the y and x axis to account for the Gibbs effect
			## Also ignores a cut off zone as defined around the centre of the full rolled FFT image - currently ASSUMES A SQUARE REGION
			gibbs_allowance= 50 # pixels
			cutoff_zone = 40 # pixels
			#mod_length = freq - gibbs_allowance
			mod_length = (kernel_size/2) - gibbs_allowance
					
			## Subsample quadrants (require top and bottom left - 1 and 3)
			## Test which quadrant has greater maximum value
			test_for_quadrant_max = 1
			if(test_for_quadrant_max == 1):
				Q1_title = "Q1 (top left)"
				#Q1 = FFT_filter_functions.subsample_quadrant(np.log(brown_noise),0,300,0,300,opath,Q1_title )
				Q1 = FFT_filter_functions.subsample_quadrant(np.log(brown_noise),0,kernel_size//2,0,kernel_size//2,opath,Q1_title )
				GIBBS_TRIMMED_Q1 = Q1[0:mod_length,0:mod_length].copy()
				ny, nx = Q1.shape
				filtered_quart_GIBBS_TRIMMED_Q1 = FFT_filter_functions.FFT_gaussian_filter_values_SIMPLE(post, nx, ny, GIBBS_TRIMMED_Q1, opath, 25)
				x_axis_limit_Q1 = len(filtered_quart_GIBBS_TRIMMED_Q1) - cutoff_zone
				y_axis_limit_Q1 = len(filtered_quart_GIBBS_TRIMMED_Q1) - cutoff_zone 
				gibbs_trimmed_and_cutoff_area_max_Q1 = filtered_quart_GIBBS_TRIMMED_Q1[0:y_axis_limit_Q1,0:x_axis_limit_Q1].max()
												
				Q3_title = "Q3 (bottom left)"
				#Q3 = FFT_filter_functions.subsample_quadrant(np.log(brown_noise),500,1000,0,500,opath,Q3_title)
				Q3 = FFT_filter_functions.subsample_quadrant(np.log(brown_noise),kernel_size//2,kernel_size,0,kernel_size//2,opath,Q3_title)
				GIBBS_TRIMMED_Q3 = Q3[gibbs_allowance:nyqvist_frq,gibbs_allowance:nyqvist_frq].copy()
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
					title = "Q1 (top left)"
					ny, nx = Q1.shape
				elif(gibbs_trimmed_and_cutoff_area_max_Q1 < gibbs_trimmed_and_cutoff_area_max_Q3):
					print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
					print "Bottom left corner will be used"
					print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
					corner = 3
					title = "Q3 (bottom left)"
					ny, nx = Q3.shape
				else:
					print "Both quarters have same maximum value"
					print "Program now closing"
					#os._exit(1)
					
					#print " BUT FOR NOW JUST USE CORNER VALUE AS 1"
					#corner = 1
					#title = "Q1 (top left)"
					#ny, nx = Q1.shape
					
					print " BUT FOR NOW JUST USE CORNER VALUE AS 3"
					corner = 3
					title = "Q3 (bottom left)"
					ny, nx = Q3.shape
					
				x = np.linspace(0,nx-1,nx)
				y = np.linspace(0,ny-1,ny)
					
				if(corner == 1):
					
					print "The origin of this quadrant is [%i,%i]" %( nx , ny )
					gibbs_y, gibbs_x = filtered_quart_GIBBS_TRIMMED_Q1.shape
					
					for i in range(len(filtered_quart_GIBBS_TRIMMED_Q1)):
						for j in range(len(filtered_quart_GIBBS_TRIMMED_Q1[i])):
					
							#if(i > (0 + cutoff_zone)) and (j < (300 - cutoff_zone)):
							# if(i > (0 + cutoff_zone)) and (j < (kernel_size/2 - cutoff_zone)):
								# if(filtered_quart_GIBBS_TRIMMED_Q1[i,j] == filtered_quart_GIBBS_TRIMMED_Q1.max()):
									# do some stuff
								# else:
									# #print "Avoiding central cutoff zone...."
									# #print "[%f.0,%f.0]" %(j,i)
									# pass
							if(i > (gibbs_y - cutoff_zone)) and (j > gibbs_x - cutoff_zone):
								pass
							else:
								if(filtered_quart_GIBBS_TRIMMED_Q1[i,j] == filtered_quart_GIBBS_TRIMMED_Q1.max()):
									print filtered_quart_GIBBS_TRIMMED_Q1[i,j]
									print "y = %.2f" %y[i]
									filtered_quart_GIBBS_TRIMMED_Q1_max_value_y = y[i]
									print "x = %.2f" %x[j]
									filtered_quart_GIBBS_TRIMMED_Q1_max_value_x = x[j]
									max_position_orientation = FFT_filter_functions.FFT_max_value_POINT_bearings_NORTH_FRQ(filtered_quart_GIBBS_TRIMMED_Q1_max_value_y, filtered_quart_GIBBS_TRIMMED_Q1_max_value_x, nyqvist_frq, corner, nx, ny)
									dist_frq = FFT_filter_functions.peak_distance(filtered_quart_GIBBS_TRIMMED_Q1_max_value_x,filtered_quart_GIBBS_TRIMMED_Q1_max_value_y,calc_moving_window_size_x,calc_moving_window_size_y,nyqvist_frq,corner)
									print "max position distance (frq): %f" %dist_frq
									#print "max position distance (px): %f" %dist_px
									#print "max position distance (m): %f" %(dist_px*post)
									print "max position orientation: %f degreesN" %(max_position_orientation)			
									print "max position orientation (rotated -90degN): %f degreesN" %(max_position_orientation - 90.)	
												
					plt.imshow(np.log(filtered_quart_GIBBS_TRIMMED_Q1)), plt.colorbar()
					plt.title("Corner 1: [%i,%i]" %(ii, jj))
					plt.gca().format_coord = FFT_functions.formatter_FRQ(filtered_quart_GIBBS_TRIMMED_Q1, nyqvist_frq)
					plt.show()		
					
					# populate spacing array
					# populate orientation array
												
				elif(corner == 3):
					
					print "The origin of this quadrant is [%i,%i]" %( nx, 0)
					gibbs_y, gibbs_x = filtered_quart_GIBBS_TRIMMED_Q1.shape
					
					for i in range(len(filtered_quart_GIBBS_TRIMMED_Q3)):
						for j in range(len(filtered_quart_GIBBS_TRIMMED_Q3[i])):
						
							#if(i < (cutoff_zone)) and (j > (300 - cutoff_zone)):
							#if(i < (cutoff_zone)) and (j > (kernel_size/2 - cutoff_zone)):
							if((i < cutoff_zone) and (j > (gibbs_x - cutoff_zone))):
								pass
							else:
								if(filtered_quart_GIBBS_TRIMMED_Q3[i,j] == filtered_quart_GIBBS_TRIMMED_Q3.max()):
									print filtered_quart_GIBBS_TRIMMED_Q3[i,j]
									print "y = %.2f" %y[i]
									filtered_quart_GIBBS_TRIMMED_Q3_max_value_y = y[i]
									print "x = %.2f" %x[j]
									filtered_quart_GIBBS_TRIMMED_Q3_max_value_x = x[j]
									max_position_orientation = FFT_filter_functions.FFT_max_value_POINT_bearings_NORTH_FRQ(filtered_quart_GIBBS_TRIMMED_Q3_max_value_y, filtered_quart_GIBBS_TRIMMED_Q3_max_value_x, nyqvist_frq, corner, nx, 0)
									dist_frq = FFT_filter_functions.peak_distance(filtered_quart_GIBBS_TRIMMED_Q3_max_value_x,filtered_quart_GIBBS_TRIMMED_Q3_max_value_y,calc_moving_window_size_x,calc_moving_window_size_y,nyqvist_frq,corner)
									
									print "max position distance (frq): %f" %dist_frq
									#print "max position distance (px): %f" %dist_px
									#print "max position distance (m): %f" %(dist_px*post)
									print "max position orientation: %f degreesN" %(max_position_orientation)			
									print "max position orientation (rotated -90degN): %f degreesN" %(max_position_orientation - 90.)	
								
					plt.imshow(np.log(filtered_quart_GIBBS_TRIMMED_Q3)), plt.colorbar()
					plt.title("Corner 3: [%i,%i]" %(ii, jj))
					plt.gca().format_coord = FFT_functions.formatter_FRQ(filtered_quart_GIBBS_TRIMMED_Q3, nyqvist_frq)
					plt.show()		
					
					# populate spacing array
					# populate orientation array
						
				print "End of loop %i" %window_iteration
				#sys.exit()

endTime = time.time()
run_time = endTime - startTime 
print "Run time: %.4f seconds" %run_time