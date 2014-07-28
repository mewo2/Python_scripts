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
import FFT_functions as FFT_functions


def sorting_FFT_values_UNTESTED_METHODS(dem_maximum_filter, magnitude):
	
		def get_order_array(a):
			b = np.empty(a.shape, dtype=int)
			for k, row in enumerate(a):
				b[k] = rankdata(-row, method='dense') - 1
			return b
				
		#ordered_array = get_order_array(dem_maximum_filter)
				
		def get_order_array_2(a):
			#a_idx = np.argsort(a, axis=-1)[:, ::-1]
			a_idx = np.argsort(a, axis=1)[:, ::1]
			a_sorted = a[np.arange(a.shape[0])[:, None], a_idx]
			#a_diff = np.zeros_like(a_sorted, dtype=np.bool)
			a_diff = np.zeros(a_sorted.shape, dtype=np.bool)
			a_diff[:, 1:] = a_sorted[:, :-1] != a_sorted[:, 1:]
			a_sorted_ranks = np.cumsum(a_diff, axis=1)
			a_ranks = a_sorted_ranks[np.arange(a_sorted_ranks.shape[0])[:, None], np.argsort(a_idx, axis=1)]
			return a_ranks

		ordered_array = get_order_array_2(dem_maximum_filter)	
				
		def get_order_x_Akavall(x):
			unique_x = np.unique(x) ## gets all of the unique values in the list (i.e. not duplicates)
			step_1 = np.argsort(unique_x)[::-1] # sorts the unique list (in reverse)
			temp_dict = dict(zip(unique_x, step_1)) # zips together the unique values (1 col) and their sort positions (col 2) creating a dictionary 
			print temp_dict
			return np.vectorize(temp_dict.get)(x) 
			
			
		def get_order_array_Akavall(x):
			new_array = np.empty(x.shape, dtype=np.int)
			#print "working"
			for i in xrange(x.shape[0]):
				new_array[i] = get_order_x_Akavall(x[i])
				#print "running"
			return new_array	
		
		#ordered_array = get_order_array_Akavall(dem_maximum_filter)	
		
		print "ordered_array:"
		print ordered_array
		print "ordered_array shape:"
		print ordered_array.shape
		print "ordered_array min:"
		print ordered_array.min()
		print "ordered_array max:"
		print ordered_array.max()
		print "dem_maximum_filter at 'ordered_array' max (99):"		
		print dem_maximum_filter[ordered_array == ordered_array.max()]
		print "ordered_array at 'maximum filtered' max (589344.0):"		
		print ordered_array[dem_maximum_filter == dem_maximum_filter.max()]
		print "dem_maximum_filter max:"
		print dem_maximum_filter.max()
		print "dem_maximum_filter dtype:"
		print dem_maximum_filter.dtype
		print dem_maximum_filter.shape
		
		print "Central magnitude:"
		print np.log(magnitude[300,300])
		print "Central max_filter:"
		print np.log(dem_maximum_filter[300,300])
		print "Central ordered array:"
		print ordered_array[300,300]
		
		plt.subplot(3,1,1)
		plt.imshow(np.log(magnitude)),plt.colorbar()
		
		plt.subplot(3,1,2)
		plt.imshow(np.log(dem_maximum_filter)),plt.colorbar()

		plt.subplot(3,1,3)		
		plt.imshow(ordered_array),plt.colorbar()
		plt.show()
		
		os._exit(1)

# Essentially calc. angle of non-right triangle ABC
# C is the angle between the vertical axis(for a geolocated image this is the N line...) at the origin and the grid coordinates of the peak
# c(side length) runs from the peak position (x,y : pos_x_coord,pos_y_coord) to the top of the N arrow (x,y : frq,frq*2)
# a (side length) runs from the origin (x,y : frq,frq) to the peak position (x,y : pos_x_coord,pos_y_coord)
# b (side length) runs from the origin (x,y : frq,frq) to the top of the "north line" (x,y : frq,frq*2)
def FFT_max_value_POINT_bearings_NORTH_FRQ(pos_x_coord,pos_y_coord,frq):
	frq_float = frq + 0.0
	origin_x = frq_float
	origin_y = frq_float
	north_x = frq_float
	north_y = frq_float*2 	
	
	'''
	print "pos_x_coord: %f" %pos_x_coord
	print "pos_y_coord: %f" %pos_y_coord
	print "frq_float: %f" %frq_float
	print "origin_x: %f" %origin_x
	print "origin_y: %f" %origin_y
	print "north_x: %f" %north_x
	print "north_y: %f" %north_y
	'''
	if(pos_y_coord == origin_y and pos_x_coord == origin_x):
		C_deg = 0.0
		return C_deg
	elif((pos_y_coord < origin_y and pos_x_coord >= origin_x) or (pos_x_coord >= origin_x and pos_y_coord >= origin_y)):
		a = math.sqrt(((origin_x - pos_x_coord)**2)+((origin_y - pos_y_coord)**2))
		b = frq_float
		c = math.sqrt(((north_x - pos_x_coord)**2)+((north_y - pos_y_coord)**2))
		#cosC = (((a**2) + (b**2)) - (c**2)) // (2*((a)*(b)))
		a2 = a**2
		b2 = b**2
		c2 = c**2
		ab_2 = 2*(a*b)
		
		cosC_1 = (a2 + b2) - c2
		cosC_2 = cosC_1/ab_2
		
		C_rad = math.acos(cosC_2)
		C_deg = math.degrees(C_rad)
		#C_deg = math.degrees(math.acos(cosC))
		
		'''
		print "a: %f" %a
		print "b: %f" %b
		print "c: %f" %c
		print "a2: %.50f" %a2
		print "b2: %.50f" %b2
		print "c2: %.50f" %c2
		print "ab_2: %f" %ab_2
		print "cosC_1: %f" %cosC_1
		print "cosC_2: %f" %cosC_2
		print "C_rad: %f" %C_rad
		print "C_deg: %f" %C_deg
		'''
		
		if(C_deg > 180.):
			print "Fail inside bearing calculation - angle size too great...."
			os._exit(1)
		else:
			return C_deg
	elif((pos_y_coord < origin_y and pos_x_coord < origin_x) or (pos_y_coord >= origin_y and pos_x_coord < origin_x)):
		a = math.sqrt(((origin_x - pos_x_coord)**2)+((origin_y - pos_y_coord)**2))
		b = frq_float
		c = math.sqrt(((north_x - pos_x_coord)**2)+((north_y - pos_y_coord)**2))
		#cosC = (((a**2) + (b**2)) - (c**2)) // (2*((a)*(b)))
		a2 = a**2
		b2 = b**2
		c2 = c**2
		ab_2 = 2*(a*b)
		
		cosC_1 = (a2 + b2) - c2
		cosC_2 = cosC_1/ab_2
		
		C_rad = math.acos(cosC_2)
		C_deg = 360. - math.degrees(C_rad)
		#C_deg = math.degrees(math.acos(cosC))
		
		'''
		print "a: %f" %a
		print "b: %f" %b
		print "c: %f" %c
		print "a2: %.50f" %a2
		print "b2: %.50f" %b2
		print "c2: %.50f" %c2
		print "ab_2: %f" %ab_2
		print "cosC_1: %f" %cosC_1
		print "cosC_2: %f" %cosC_2
		print "C_rad: %f" %C_rad
		print "C_deg: %f" %C_deg
		'''
		
		if(C_deg < 180.):
			print "Fail inside bearing calculation - angle size too great...."
			os._exit(1)
		else:
			return C_deg
	
## FFT_max_filter_values is an older function that maximum filtered the surface, got the xy coords of the maximum pixels and then calculated their bearings (degN) relative to the centre of the image (assuming "up" is north) - outputs were then written to a .txt file - this has been reworked now so is OLD
def FFT_max_filter_values_OLD(freq, post, input_x, input_y, magnitude, kernel=50):
		
		print "magnitude.shape"
		print magnitude.shape
		print "magnitude.dtype"
		print magnitude.dtype
		print "magnitude.max"
		print magnitude.max()
		print "magnitude.min"
		print magnitude.min()
		
		print "input_x: %f" %input_x
		print "input_y: %f" %input_y

		print '~~~~~~~~~~~~~~'
		print 'Maximum filter'
		print '~~~~~~~~~~~~~~'

		print "Calculating...."
		method = "_maximum"
		dem_maximum_filter = ndimage.filters.maximum_filter(magnitude,size=(kernel,kernel),mode='reflect')
		filter_max = dem_maximum_filter.max()
		#dem_maximum_filter[300,300] = 0.0
		'''
		print "magnitude.shape"	
		print dem_maximum_filter.shape
		print "magnitude.dtype"
		print dem_maximum_filter.dtype
		print "Max filter: max"
		
		print filter_max
		print "Max filter: min"
		print dem_maximum_filter.min()
		'''
	
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		print ' Maximum value coordinates '
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		
		b = dem_maximum_filter == (filter_max)
		coords =  np.column_stack(np.where(b))
		#print coords.shape
		max_coord_x = coords[:,0]
		max_coord_y = coords[:,1] 
		
				
		print "############"
		print "     X      " 
		print "############"
		print max_coord_x
		print max_coord_x.shape
		print max_coord_x.dtype

		print "############"
		print "     Y      " 
		print "############"
		print max_coord_y
		print max_coord_y.shape
		print max_coord_y.dtype

		print "##############"
		print " Co-ordinates " 
		print "##############"
		
		#sorting_FFT_values(dem_maximum_filter, magnitude)
		
		accum_dist_frq = 0.0
		accum_dist_px = 0.0

		file_name = "FFT_test.txt"
		f = open( file_name, 'w')
		f.write("mag,dist.(frq),dist.(px),dist.(m),bearing (degreesN)\n")
		
		# peak of interest (position in maximum array)
		for POI in range(len(max_coord_x)):
		
			x_co = max_coord_x[POI]
			y_co = max_coord_y[POI]
			print "x_co: %f" %x_co
			print "y_co: %f" %y_co
	
			#mag = magnitude[x_co][y_co] # WRONG : this would get the unfiltered value from the position of the condition required for the max value coords of the filtered surface
			mag = dem_maximum_filter[x_co][y_co]
	
			#print 'Get x and y positions as distances from central origin of [0,0]'
			#x_co_frq = max_coord_x - freq/2
			#y_co_frq = max_coord_y - freq/2
			
			if (x_co>freq):
				x_co_frq = x_co - freq
			elif(x_co == 300.0):
				x_co_freq = 0.1
			else:
				x_co_frq = (freq - x_co)
			
			if (y_co>freq):
				y_co_frq = y_co - freq
			elif(y_co == 300.0):
				y_co_freq = 0.1
			else:
				y_co_frq = (freq - y_co)
				
			#x_co_frq = x_co - freq/2
			#y_co_frq = y_co - freq/2
			print "x_co_frq: %f" %x_co_frq
			print "y_co_frq: %f" %y_co_frq
			
			print "Freq: %f" %freq	
			
			'''
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			print "Calc. diagnonal distance of x,y at position 1 from origin"
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			'''
				
			### In frequencies
			x_co_frq_pos = x_co_frq #[POI]
			y_co_frq_pos = y_co_frq #[POI]
			print "x_co_frq_pos: %f" %x_co_frq_pos
			print "y_co_frq_pos: %f" %y_co_frq_pos
			xy_dist_frq = math.sqrt(x_co_frq_pos**2 + y_co_frq_pos**2)	# relative to the centre of the FFT which is [0,0]
			print "xy_dist_frq: %f" %xy_dist_frq
			
			### In pixels
			x_co_px_pos = input_x/x_co_frq_pos
			y_co_px_pos = input_y/y_co_frq_pos
			print "x_co_px_pos: %f" %x_co_px_pos
			print "y_co_px_pos: %f" %y_co_px_pos
			xy_dist_px = math.sqrt(x_co_px_pos**2 + y_co_px_pos**2)	# relative to the centre of the FFT which is [0,0]
			print "xy_dist_px: %f" %xy_dist_px
			
			'''
			print "Position: %i" %POI
			print "Magnitude at position:"
			print mag
	
			print "Distance of position:"
			print "%f (frq.)" %(xy_dist_frq)
			print "%f (px)" %(xy_dist_px)
			print "%f (m)" %(xy_dist_px*post)
			'''
			
			#accum_dist_frq += xy_dist_frq
			#accum_dist_px += xy_dist_px
			
			print "#################"
			print "Calculate bearing"
			print "#################"
						
			bearing = FFT_max_value_bearings_NORTH_FRQ(x_co,y_co,freq)
			
			instance_dist_frq = xy_dist_frq
			instance_dist_px = xy_dist_px
			
			#mean_dist_frq = len(max_coord_x)/accum_dist_frq
			#mean_dist_px = len(max_coord_x)/accum_dist_px
			f.write("%f,%f,%f,%f,%f" %(mag,instance_dist_frq,instance_dist_px,instance_dist_px*post,bearing))
			f.write("\n")
	
		## Calc. max magnitude mean point distance
		'''	
		mean_dist_frq = accum_dist_frq/len(max_coord_x)
		mean_dist_px = accum_dist_px/len(max_coord_x)
		
		print "Mean spacing of maximum FFT 'spike':"
		print "%f (frq.)" %(mean_dist_frq) 
		print "%f (px)" %(mean_dist_px)
		print "%f (m)" %(mean_dist_px*post)
		'''
		plt.imshow(np.log(dem_maximum_filter)), plt.colorbar()
		plt.show()
		plt.clf()
		
		#file_name = "FFT_test.txt"
		#f = open( file_name, 'w')
		#f.write("mag,dist.(frq),dist.(px),dist.(m)\n")
		#f.write("%f,%f,%f,%f" %(mag,mean_dist_frq,mean_dist_px,mean_dist_px*post))
		#f.write("\n")
		f.close()
		
		print "Output file written: %s" %file_name
		
		return dem_maximum_filter

# FFT_max_filter_values_SIMPLE creates a maximum filtered FFT surface and returns the array
def FFT_max_filter_values_SIMPLE(freq, post, input_x, input_y, magnitude, opath, kernel=50):
		
		print "magnitude.shape"
		print magnitude.shape
		print "magnitude.dtype"
		print magnitude.dtype
		print "magnitude.max"
		print magnitude.max()
		print "magnitude.min"
		print magnitude.min()
		
		
		print "input_x: %f" %input_x
		print "input_y: %f" %input_y

		#os._exit(1)
		
		print '~~~~~~~~~~~~~~'
		print 'Maximum filter'
		print '~~~~~~~~~~~~~~'

		print "Calculating...."
		method = "_maximum"
		dem_maximum_filter = ndimage.filters.maximum_filter(magnitude,size=(kernel,kernel),mode='reflect')
		filter_max = dem_maximum_filter.max()
		
		plt.imshow(np.log(dem_maximum_filter)), plt.colorbar()
		plt.title("Maximum surface")
		output_filename = opath + 'maximum_surface_DEVELOPMENT_TEST.png' 
		plt.savefig(output_filename)
		#plt.show()
		
		plt.clf()
		
		return dem_maximum_filter

# Stifles the apparent gibbs effect that is clear on the maximum filter surfaces - limits on the x and y axis to be stifled are hardwired in the function
def Stifle_Gibbs_effect_on_FFT_max_filter(max_filter_surface, freq, opath):
		#  if coordinates are between 250 and 350 in x, make value NaN
		#  if coordinates are between 250 and 350 in y, make value NaN
		
		gibbs_removed_filtered_surface = max_filter_surface
		
		a = 0
		b = 0
		
		lower_limit = 250
		upper_limit = 350
			
		print "lower_limit = %i" %lower_limit
		print "upper_limit = %i" %upper_limit
			
		print "#######################"
		print "######## ROWS #########"
		print "#######################" 
		for i in range(len(gibbs_removed_filtered_surface)-1):
			a+=1
			#print "%i" %a
					
			if(a >= lower_limit and a <=  upper_limit):
				gibbs_removed_filtered_surface[:,a] = 1.0
				#print "value changed!"
			else:
				wtf =0
				#print "value not changed!"
						
		print "#######################"
		print "######## COLS #########"
		print "#######################" 
		for i in range(len(gibbs_removed_filtered_surface[0])-1):
			b+=1
			#print "%i" %b
			
			if(b >= lower_limit and b <=  upper_limit):
				gibbs_removed_filtered_surface[b,:] = 1.0
				#print "value changed!"
			#else:
				wtf =0
				#print "value not changed!"
				
		
		plt.imshow(np.log(gibbs_removed_filtered_surface)), plt.colorbar()
		plt.title("Gibbs effect suppressed")
		#plt.show()
		
		output_filename = opath + 'gibbs_effect_stifled_DEVELOPMENT_TEST.png' 
		plt.savefig(output_filename)
		plt.clf()
		
		return gibbs_removed_filtered_surface

# Outputs the rank position of each pixel values in an array (1 can be the biggest or smallest)
def get_order_array_aaron(a, opath):
	values, inverse = np.unique(a, return_inverse=True)
	sort_values = np.arange(values.size)
	
	zero_smallest = sort_values[inverse].reshape(a.shape)
	one_biggest = sort_values.max() + 1 - zero_smallest
	
	plt.imshow(one_biggest), plt.colorbar()
	plt.title("Ranked FFT surface")
	output_filename = opath + '\gibbs_effect_removed_RANKED_DEVELOPMENT_TEST_AARON.png' 
	plt.savefig(output_filename)
	plt.clf()
	
	print output_filename
	
	return one_biggest

# Another routine (part 1) that ranks pixel positions - image development NOT WORKING
def get_order_part_1_of_2(x):
    unique_x = np.unique(x)
    step_1 = np.argsort(unique_x)[::-1]
    temp_dict = dict(zip(unique_x, step_1))
    return np.vectorize(temp_dict.get)(x)
	
# Another routine (part 2) that ranks pixel positions - image development NOT WORKING	
def get_order_array_part_2_of_2(surface_to_be_ranked, opath):
    ranked_array = np.empty(surface_to_be_ranked.shape, dtype=np.int)
    for i in xrange(surface_to_be_ranked.shape[0]):
        ranked_array[i] = get_order_part_1_of_2(surface_to_be_ranked[i])
    
	print "Max rank number:"
	print ranked_array.max()
	print "Min rank number:"
	print ranked_array.min()
	print ranked_array
	
	a = np.log(surface_to_be_ranked[ranked_array == ranked_array.min()])
	print a.max()
	print np.log(surface_to_be_ranked.max())
	
	plt.imshow(ranked_array), plt.colorbar()
	plt.title("Ranked gibbs removed FFT surface")
	output_filename = opath + 'gibbs_effect_removed_RANKED_DEVELOPMENT_TEST.png' 
	plt.savefig(output_filename)
	plt.clf()
	
	return ranked_array

	
def FFT_max_filtered_Max_positions_and_values_and_bearings_OLD(maximum_filtered_fft, freq, post, input_x, input_y, opath):
		
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		print ' Maximum value coordinates '
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		
		filter_max = maximum_filtered_fft.max() ### THIS WORKS
		
				## @ 16/07/14 20:02
				## now need to get coordinates of all distances of filter_max in maximum_filtered_fft - use numpy mesh to create x and y and then pull off the x and y positions in those arrays where the value of maximum_filtered_fft ==  filter_max
				## Those values can then be appended to the arrays of max_coord_x and  max_coord_y
				
		b = maximum_filtered_fft == (filter_max) ### THIS IS THE PROBLEM - LOOK AT USING NUMPY MESH??
		coords =  np.column_stack(np.where(b))  ### THIS IS THE PROBLEM - LOOK AT USING NUMPY MESH??
		
		#print coords.shape
		max_coord_x = coords[:,0]
		max_coord_y = coords[:,1] 
		debug = 1
		
		if(debug == 1):
			print "b.shape:"
			print b.shape
			print "coords.shape:"
			print coords.shape
			print "filter_max (log):"
			print np.log(filter_max)
			os._exit(1)
		else:
			wtf = 0
				
		print "############"
		print "     X      " 
		print "############"
		print max_coord_x
		print max_coord_x.shape
		print max_coord_x.dtype

		print "############"
		print "     Y      " 
		print "############"
		print max_coord_y
		print max_coord_y.shape
		print max_coord_y.dtype

		print "##############"
		print " Co-ordinates " 
		print "##############"
		
		accum_dist_frq = 0.0
		accum_dist_px = 0.0

		file_name = "%s/FFT_test.txt" %opath
		f = open( file_name, 'w')
		f.write("mag,dist.(frq),dist.(px),dist.(m),bearing (degreesN), x(grid_pos), y(grid_pos)\n")
		
		# peak of interest (position in maximum array)
		for POI in range(len(max_coord_x)):
		
			x_co = max_coord_x[POI]
			y_co = max_coord_y[POI]
			'''
			print "x_co: %f" %x_co
			print "y_co: %f" %y_co
			'''
			#mag = magnitude[x_co][y_co] # WRONG : this would get the unfiltered value from the position of the condition required for the max value coords of the filtered surface
			mag = maximum_filtered_fft[x_co][y_co]
	
			#print 'Get x and y positions as distances from central origin of [0,0]'
			#x_co_frq = max_coord_x - freq/2
			#y_co_frq = max_coord_y - freq/2
			
			if (x_co>freq):
				x_co_frq = x_co - freq
			elif(x_co == 300.0):
				x_co_freq = 0.1
			else:
				x_co_frq = (freq - x_co)
			
			if (y_co>freq):
				y_co_frq = y_co - freq
			elif(y_co == 300.0):
				y_co_freq = 0.1
			else:
				y_co_frq = (freq - y_co)
				
			#x_co_frq = x_co - freq/2
			#y_co_frq = y_co - freq/2
			'''
			print "x_co_frq: %f" %x_co_frq
			print "y_co_frq: %f" %y_co_frq
			
			print "Freq: %f" %freq	
			'''
			'''
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			print "Calc. diagnonal distance of x,y at position 1 from origin"
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
			'''
				
			### In frequencies
			x_co_frq_pos = x_co_frq #[POI]
			y_co_frq_pos = y_co_frq #[POI]
			
			#print "x_co_frq_pos: %f" %x_co_frq_pos
			#print "y_co_frq_pos: %f" %y_co_frq_pos
			
			xy_dist_frq = math.sqrt(x_co_frq_pos**2 + y_co_frq_pos**2)	# relative to the centre of the FFT which is [0,0]
			
			#print "xy_dist_frq: %f" %xy_dist_frq
			
			### In pixels
			x_co_px_pos = input_x/x_co_frq_pos
			y_co_px_pos = input_y/y_co_frq_pos
			#print "x_co_px_pos: %f" %x_co_px_pos
			#print "y_co_px_pos: %f" %y_co_px_pos
			xy_dist_px = math.sqrt(x_co_px_pos**2 + y_co_px_pos**2)	# relative to the centre of the FFT which is [0,0]
			#print "xy_dist_px: %f" %xy_dist_px
			
			'''
			print "Position: %i" %POI
			print "Magnitude at position:"
			print mag
	
			print "Distance of position:"
			print "%f (frq.)" %(xy_dist_frq)
			print "%f (px)" %(xy_dist_px)
			print "%f (m)" %(xy_dist_px*post)
			'''
			
			#accum_dist_frq += xy_dist_frq
			#accum_dist_px += xy_dist_px
			
			#print "#################"
			#print "Calculate bearing"
			#print "#################"
						
			bearing = FFT_max_value_POINT_bearings_NORTH_FRQ(x_co,y_co,freq)
			
			instance_dist_frq = xy_dist_frq
			instance_dist_px = xy_dist_px
			
			#mean_dist_frq = len(max_coord_x)/accum_dist_frq
			#mean_dist_px = len(max_coord_x)/accum_dist_px
			f.write("%f,%f,%f,%f,%f,%f,%f" %(mag,instance_dist_frq,instance_dist_px,instance_dist_px*post,bearing,x_co,y_co))
			f.write("\n")
	
		## Calc. max magnitude mean point distance
		'''	
		mean_dist_frq = accum_dist_frq/len(max_coord_x)
		mean_dist_px = accum_dist_px/len(max_coord_x)
		
		print "Mean spacing of maximum FFT 'spike':"
		print "%f (frq.)" %(mean_dist_frq) 
		print "%f (px)" %(mean_dist_px)
		print "%f (m)" %(mean_dist_px*post)
		'''
		plt.imshow(np.log(maximum_filtered_fft)), plt.colorbar()
		plt.show()
		plt.clf()
		#file_name = "FFT_test.txt"
		#f = open( file_name, 'w')
		#f.write("mag,dist.(frq),dist.(px),dist.(m)\n")
		#f.write("%f,%f,%f,%f" %(mag,mean_dist_frq,mean_dist_px,mean_dist_px*post))
		#f.write("\n")
		f.close()
		
		print "Output file written: %s" %file_name
	
# Where a string of maximum values occur in a quarter if a rolled FFT matrix, this function calculates the most extreme vectors
# ASSUMES x and y are corrected to represent pixel positions with origin (0,0) in the bottom left hand corner
# NB/ NEEDS FIXING SO THAT IT CALCS VECTORS ISOLATED TO FFT MAXIMUM PIXEL STREAKS..... i.e. doesn't pick up other random maximum values

# Following functions get lists of the x,y positions of the max values located in the different quarters of the 2D FFT array
def list_max_x_and_y_coords_quad_1(x, y, xy_list_1, freq):
	
	#if(x < freq and y > freq):
	if(x < freq and y < freq):
		xy_list_1.append([x,y])
	
	#print "list development factory 1...."
	return xy_list_1

def list_max_x_and_y_coords_quad_2(x, y, xy_list_2, freq):
	
	#if(x > freq and y > freq):
	if(x > freq and y < freq):
		xy_list_2.append([x,y])
	
	#print "list development factory 2...."
	return xy_list_2

def list_max_x_and_y_coords_quad_3(x, y, xy_list_3, freq):
	
	#if(x < freq and y < freq):
	if(x < freq and y > freq):
		xy_list_3.append([x,y])
	
	#print "list development factory 3...."
	return xy_list_3

def list_max_x_and_y_coords_quad_4(x, y, xy_list_4, freq):
	
	#if(x > freq and y < freq):
	if(x > freq and y > freq):
		xy_list_4.append([x,y])
	
	#print "list development factory 4...."
	return xy_list_4	

# Calculates maximum values and positions (takes in the maximum filtered FFT surface)
def FFT_max_filtered_Max_positions_and_values(maximum_filtered_fft, freq, post, input_x, input_y, opath):
	
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	print ' Maximum value coordinates '
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		
	filter_max = np.log(maximum_filtered_fft.max())
	
	print "filter_max value (log):"
	print filter_max
	filter_max = np.int_(filter_max)
	print filter_max
		
	nx,ny = maximum_filtered_fft.shape
		
	x = np.linspace(0,nx-1,nx)
	y = np.linspace(0,ny-1,ny)

	file_name = "%s/FFT_max_xyz_test.txt" %opath
	f = open( file_name, 'w')
	f.write("mag,x_pos,y_pos,bearing_degN\n")
		
	max_filtered_log = np.log(maximum_filtered_fft)
	
	list_1 = []
	list_2 = []
	list_3 = []
	list_4 = []
		
	for i in range(len(maximum_filtered_fft)):
		for j in range(len(maximum_filtered_fft[i])):
				
			#if(maximum_filtered_fft[i,j] == filter_max):
			if(max_filtered_log[i,j] >= filter_max):
			#if(max_filtered_log[i,j] >= 18.0):
				#print "[%i,%i]" %(x[i],y[j]) 
				#a.append([x[i],y[j]]) 
					
				mag = np.log(maximum_filtered_fft[i,j])
				x_at_max = x[i]
				y_at_max = y[j]
					
				## This has to happen as the indexing of the FFT array (pixel positions) has a 
				## bottom left origin of (0,600) and NOT (0,0) as normal
				if(y_at_max >= 300):
					y_at_max_for_bearing_calc = y_at_max - freq
				elif(y_at_max < 300):
					y_at_max_for_bearing_calc = y_at_max + freq
					
				bearing = FFT_max_value_POINT_bearings_NORTH_FRQ(x_at_max,y_at_max_for_bearing_calc,freq)
				
				'''
				## This uses "corrected" bearings i.e. origin is 0,0 at bottom left (not 0,600)
				## Creates a lift for each quarter (1 = top left, 2 = top right, 3 = bottom left, 4 = bottom right)
				xy_list_1 = np.array(list_max_x_and_y_coords_quad_1(x_at_max, y_at_max_for_bearing_calc, list_1, freq))
				xy_list_2 = np.array(list_max_x_and_y_coords_quad_2(x_at_max, y_at_max_for_bearing_calc, list_2, freq))
				xy_list_3 = np.array(list_max_x_and_y_coords_quad_3(x_at_max, y_at_max_for_bearing_calc, list_3, freq))
				xy_list_4 = np.array(list_max_x_and_y_coords_quad_4(x_at_max, y_at_max_for_bearing_calc, list_4, freq))
				'''
				
				
				## This uses "uncorrected" bearings i.e. origin is 0,600 at bottom left (not 0,0)
				## Creates a lift for each quarter (1 = top left, 2 = top right, 3 = bottom left, 4 = bottom right)
				xy_list_1 = np.array(list_max_x_and_y_coords_quad_1(x_at_max, y_at_max, list_1, freq))
				xy_list_2 = np.array(list_max_x_and_y_coords_quad_2(x_at_max, y_at_max, list_2, freq))
				xy_list_3 = np.array(list_max_x_and_y_coords_quad_3(x_at_max, y_at_max, list_3, freq))
				xy_list_4 = np.array(list_max_x_and_y_coords_quad_4(x_at_max, y_at_max, list_4, freq))
				
				
				f.write("%f,%f,%f,%f" %(mag,x_at_max,y_at_max,bearing))
				f.write("\n")
							
	f.close()
	
	#print xy_list_1
	print xy_list_2.shape
	#return xy_list_1
	#return xy_list
	return xy_list_1, xy_list_2, xy_list_3, xy_list_4
	
def FFT_max_streak_vector_calculation(xy_list_1, xy_list_2, xy_list_3, xy_list_4):

	# This needs to get the coordinate pair that is furthest from the origin - individual 
	# max/min x and y values will unlikely be representative of a position of interest
	# TRIAL 1: Use minimum x and y at the equivalent position - y will likely have two values for this position if it is a pixel corner....
	def streak_vectors(xy_list):
		x_min = xy_list[:,0].min() 
		#y_min = xy_list[:,1].min() 
		y_min = xy_list[x_min,1]
		x_max = xy_list[:,0].max()
		#y_max = xy_list[:,1].max()
		y_max = xy_list[x_max,1]
		
		print "Min vector for list:"
		print "[%f, %f]" %(x_min,y_min)
		print "Max vector for list:"	
		print "[%f, %f]" %(x_max,y_max)
		return x_min, y_min, x_max, y_max
	
	def vector_direction(x_min, y_min, x_max, y_max, maxmin=1):
		print "Vector direction being calculated...."
		print "x_min"
		print x_min
		print "y_min"
		print y_min
		print "x_max"
		print x_max
		print "y_max"
		print y_max
		tan_of_angle = ((y_min - y_max)//(x_min - x_max))
		vector_rad = math.atan(tan_of_angle)
		
		vector_deg = math.degrees(vector_rad)
		print "Vector direction (degN):"
		print vector_deg
		return vector_rad, vector_deg
				
	if(xy_list_1 != []):
		print "List 1 has values"
		x_min_1, y_min_1, x_max_1, y_max_1 = streak_vectors(xy_list_1)	
		xy_vector_direction_quad_1_RAD, xy_vector_direction_quad_1_DEG = vector_direction(x_min_1, y_min_1, x_max_1, y_max_1)
		print "xy_vector_direction_quad_1 (radians): %f" %(xy_vector_direction_quad_1_RAD)
		print "xy_vector_direction_quad_1 (degrees): %f" %(xy_vector_direction_quad_1_DEG)
	else:
		print "List 1 has no values"
	
	if(xy_list_2 != []):
		print "List 2 has values"
		x_min_2, y_min_2, x_max_2, y_max_2 = streak_vectors(xy_list_2)
		xy_vector_direction_quad_2_RAD, xy_vector_direction_quad_2_DEG = vector_direction(x_min_2, y_min_2, x_max_2, y_max_2)
		print "xy_vector_direction_quad_2 (radians): %f" %(xy_vector_direction_quad_2_RAD)
		print "xy_vector_direction_quad_2 (degrees): %f" %(xy_vector_direction_quad_2_DEG)
	else:
		print "List 2 has no values"
	
	if(xy_list_3 != []):
		print "List 3 has values"
		x_min_3, y_min_3, x_max_3, y_max_3 = streak_vectors(xy_list_3)
		xy_vector_direction_quad_3_RAD, xy_vector_direction_quad_3_DEG = vector_direction(x_min_3, y_min_3, x_max_3, y_max_3)
		print "xy_vector_direction_quad_3 (radians): %f" %(xy_vector_direction_quad_3_RAD)
		print "xy_vector_direction_quad_3 (degrees): %f" %(xy_vector_direction_quad_3_DEG)
	else:
		print "List 3 has no values"
	
	if(xy_list_4 != []):
		print "List 4 has values"
		x_min_4, y_min_4, x_max_4, y_max_4 = streak_vectors(xy_list_4)
		xy_vector_direction_quad_4_RAD, xy_vector_direction_quad_4_DEG = vector_direction(x_min_4, y_min_4, x_max_4, y_max_4)
		print "xy_vector_direction_quad_4 (radians): %f" %(xy_vector_direction_quad_4_RAD)
		print "xy_vector_direction_quad_4 (degrees): %f" %(xy_vector_direction_quad_4_DEG)
	else:
		print "List 4 has no values"
	
# Removes any zeros present in an incoming array - and converts them to very small values (prevents any log routines failing)		
def remove_zeros(FFT_surface, opath, replacement_value=0.0000001):

	for i in range(len(FFT_surface)):
		for j in range(len(FFT_surface[i])):
			if(FFT_surface[i,j] == 0):
				FFT_surface[i,j] = replacement_value
			elif(FFT_surface[i,j] == 0.0):
				FFT_surface[i,j] = replacement_value
			else:
				FFT_surface[i,j] = FFT_surface[i,j]
	
	plt.imshow(np.log(FFT_surface)), plt.colorbar()
	plt.title("Gibbs effect suppressed - no zeros")
	#plt.show()
		
	output_filename = opath + 'gibbs_effect_stifled_no_zeros_DEVELOPMENT_TEST.png' 
	plt.savefig(output_filename)
	plt.clf()
	return FFT_surface
	
def zero_all_but_max(FFT_surface, opath, replacement_value=0.0000001):
	filter_max = np.log(FFT_surface.max())
	print "filter_max value (log):"
	print filter_max
	#filter_max = np.int_(filter_max)
	print filter_max
	
	log_FFT = np.log(FFT_surface)
	
	for i in range(len(FFT_surface)):
		for j in range(len(FFT_surface[i])):
			if(log_FFT[i,j] < (filter_max - 0.5)):
			#if(log_FFT[i,j] < filter_max):
				log_FFT[i,j] = replacement_value
			else:
				log_FFT[i,j] = log_FFT[i,j]
	
	#plt.imshow(np.log(FFT_surface)), plt.colorbar()
	plt.imshow(log_FFT), plt.colorbar()
	plt.title("Gibbs effect suppressed - only max")
		
	output_filename = opath + 'gibbs_effect_stifled_ONLY_MAX_DEVELOPMENT_TEST.png' 
	plt.savefig(output_filename)
	
	#plt.show()
	plt.clf()
	return log_FFT


def return_value_of_given_rank_position_SINGLE(input_FFT_surface, aaron_rank, rank_position=7):
	
	value = input_FFT_surface[aaron_rank == rank_position]
	value_int = np.int_(value)
	value_log = np.log(value)
	print "Rank position %i value (int):" %rank_position
	print value_int[0]
	print "Rank position %i value (log):" %rank_position
	print value_log[0]
	
		
def return_value_of_given_rank_position_MULTIPLE(input_FFT_surface, aaron_rank, number_of_rank_positions=7):
	
	rank = 1
	for i in range(number_of_rank_positions):
		value = input_FFT_surface[aaron_rank == rank]
		value_int = np.int_(value)
		value_log = np.log(value)
		print "Rank position %i value (int):" %rank
		print value_int[0]
		print "Rank position %i value (log):" %rank
		print value_log[0]
		rank += 1
		#print "Rank position value = %f" %value
	
## Both:
## return_value_of_given_rank_position_SINGLE_IMAGE_RETURN_PART_1
## zero_all_but_rank_value_IMAGE_RETURN_PART_1
## are used to take in a complete ranked surface (from the aaron_rank method), calculate the value for a given rank position and then output a matrix showing the position ## of these values with all other values reduced to a fixed value of 0.0000001 - it must be stated if a log surface is required
def return_value_of_given_rank_position_SINGLE_IMAGE_RETURN_PART_1(input_FFT_surface, aaron_rank, rank_position=7):
	
	value = input_FFT_surface[aaron_rank == rank_position]
	value_int = np.int_(value)
	value_log = np.log(value)
	print "Rank position %i value (int):" %rank_position
	print value_int[0]
	print "Rank position %i value (log):" %rank_position
	print value_log[0]
	
	return value_int[0], value_log[0]
		
def zero_all_but_rank_value_IMAGE_RETURN_PART_2(FFT_surface, opath, value_to_return, rank_position, log = 1):
	
	replacement_value=0.0000001
		
	if(log == 0):
		print "Non-log value used to filter surface"
	elif(log == 1):
		print "Log value used to filter surface"
		FFT_surface = np.log(FFT_surface)
	
	for i in range(len(FFT_surface)):
		for j in range(len(FFT_surface[i])):
			'''
			## This displays all values > value_to_return and is useful for deciding on which range of FFT maximum spikes to use....
			if(FFT_surface[i,j] < (value_to_return)):
				FFT_surface[i,j] = replacement_value
			else:
				FFT_surface[i,j] = FFT_surface[i,j]
			'''
			## This only displays values == value_to_return and is useful for specific rank point positions
			if(FFT_surface[i,j] != (value_to_return)):
				FFT_surface[i,j] = replacement_value
			else:
				FFT_surface[i,j] = FFT_surface[i,j]
	
	#plt.imshow(np.log(FFT_surface)), plt.colorbar()
	plt.imshow(FFT_surface), plt.colorbar()
	plt.title("FFT_surface - only rank position %s displayed" %(rank_position)) 
		
	output_filename = opath + 'rank_position_%i_fft_surface.png'  %(rank_position)
	plt.savefig(output_filename)
	
	#plt.show()
	plt.clf()
	return FFT_surface