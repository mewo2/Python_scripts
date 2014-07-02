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


def sorting_FFT_values(dem_maximum_filter, magnitude):
	
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		print '  Sort max filter surface  '
		print '    Calc. top 10 values 	  '
		print '~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		
		'''
		dem_maximum_filter_int = dem_maximum_filter.astype(int)
		indices = dem_maximum_filter_int.ravel().argsort()
		print dem_maximum_filter_int.ravel()[indices[-10:]]
		'''
		
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
				
		def get_order_array_aaron(a):
			values, inverse = np.unique(a, return_inverse=True)
			sort_values = np.sort(values.size)
			return sort_values[inverse].reshape(a.shape)
				
		#ordered_array = get_order_array_aaron(dem_maximum_filter)	
		
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
		############

# Essentially calc. angle of non-right triangle ABC
# C is the angle between the vertical axis(for a geolocated image this is the N line...) at the origin and the grid coordinates of the peak
# c(side length) runs from the peak position (x,y : pos_x_coord,pos_y_coord) to the top of the N arrow (x,y : frq,frq*2)
# a (side length) runs from the origin (x,y : frq,frq) to the peak position (x,y : pos_x_coord,pos_y_coord)
# b (side length) runs from the origin (x,y : frq,frq) to the top of the "north line" (x,y : frq,frq*2)
		
def FFT_max_value_bearings_NORTH_FRQ(pos_x_coord,pos_y_coord,frq):
	print "pos_x_coord: %f" %pos_x_coord
	print "pos_y_coord: %f" %pos_y_coord
	frq_float = frq + 0.0
	origin_x = frq_float
	origin_y = frq_float
	north_x = frq_float
	north_y = frq_float*2 	
	
	print "frq_float: %f" %frq_float
	print "origin_x: %f" %origin_x
	print "origin_y: %f" %origin_y
	print "north_x: %f" %north_x
	print "north_y: %f" %north_y
	
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
	
			
def FFT_max_filter_values(freq, post, input_x, input_y, magnitude, kernel=50):
		
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
		
		#file_name = "FFT_test.txt"
		#f = open( file_name, 'w')
		#f.write("mag,dist.(frq),dist.(px),dist.(m)\n")
		#f.write("%f,%f,%f,%f" %(mag,mean_dist_frq,mean_dist_px,mean_dist_px*post))
		#f.write("\n")
		f.close()
		
		print "Output file written: %s" %file_name
