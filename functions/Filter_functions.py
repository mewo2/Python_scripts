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

# As described here: http://homepages.inf.ed.ac.uk/rbf/HIPR2/csmooth.htm
# Using an odd sized kernel, a given point in an array greater or smaller 
# than the surrounding points will be converted to the neighbourhood maximum or minimum

# A copy of the entire array is made so that i,j of the array can be temporarily made to NaN  - can't just copy the moving window area (which would be better) as i,j of main array will not be i,j of the moving window region array
# A simple "if value at position in moving window == value in main array" test as this will fail if values are duplicated

def conservative_smooth_inefficient(array2D, kernel_size = 3):
		
	stepsize = 1	
	if(kernel_size % 2 != 0 and kernel_size >= 3):
		window = np.ones([kernel_size,kernel_size])
	elif(kernel_size % 2 == 0 or kernel_size < 3):
		print "kernel is even - it needs to be odd and at least of a value of 3"
		os._exit(1)
		
	nxwind, nywind = array2D.shape
		
	for i in range(0, nxwind, stepsize):
		for j in range(0, nywind, stepsize):
			
		# CALCULATE MAX AND MIN RANGES OF ROWS AND COLS THAT CAN BE ACCESSED BY THE WINDOW
			imin=max(0,i-((kernel_size-1)/2)) # gets the maximum of either 0 or i-kernel_size/2...
			imax=min(nxwind-1,i+((kernel_size-1)/2))+1
			jmin=max(0,j-((kernel_size-1)/2))
			jmax=min(nywind-1,j+((kernel_size-1)/2))+1
		
			array2D_temp = array2D.copy() # THIS IS THE MOST INEFFICIENT PART OF THE CODE
			array2D_temp[i,j] = np.nan
			
			data_wind=array2D_temp[imin:imax,jmin:jmax]
			print data_wind
			
			centre_value = array2D[i,j]
			#print "centre_value: %f" %centre_value
		
			max_value = np.nanmax(data_wind) 
			min_value = np.nanmin(data_wind) 
			
			if(centre_value > max_value):
				centre_value = max_value
			elif(centre_value < min_value):
				centre_value = min_value
			else:
				centre_value = centre_value
			
			#print "max_value: %f" %max_value	
			#print "min_value: %f" %min_value
			#print "New centre_value: %f" %centre_value
			
			## Append new centre value to output array
			array2D[i,j] = centre_value			
			
	return array2D
	
# Removed the moving window print out as it wasn;t required (compare with "conservative_smooth_inefficient()")
def conservative_smooth_efficient(array2D, kernel_size = 3):
	stepsize = 1    
	if(kernel_size % 2 != 0 and kernel_size >= 3):
		window = np.ones([kernel_size,kernel_size])
	elif(kernel_size % 2 == 0 or kernel_size < 3):
		print "kernel is even - it needs to be odd and at least of a value of 3"
		os._exit(1)
   
	nxwind, nywind = array2D.shape
	print array2D
   
	for i in range(0, nxwind, stepsize):
		for j in range(0, nywind, stepsize):
			imin=max(0,i-((kernel_size-1)/2)) 
			imax=min(nxwind-1,i+((kernel_size-1)/2))+1
			jmin=max(0,j-((kernel_size-1)/2))
			jmax=min(nywind-1,j+((kernel_size-1)/2))+1
			centre_value = array2D[i,j]
			array2D[i,j] = np.nan
			max_value = np.nanmax(array2D[imin:imax,jmin:jmax]) 
			min_value = np.nanmin(array2D[imin:imax,jmin:jmax]) 
			if(centre_value > max_value):
				centre_value = max_value
			elif(centre_value < min_value):
				centre_value = min_value
			else:
				centre_value = centre_value
			## Append new centre value to output array
			array2D[i,j] = centre_value      
	print "~~~~"
	print array2D
	
	return array2D