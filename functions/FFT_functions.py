'''
FFT functions
'''

from __future__ import division

import sys
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
from time import gmtime, strftime
import pylab as pl
from glob import glob

from osgeo import gdal, gdalconst # for reading in raster
from osgeo.gdalconst import * # for reading in raster

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


def FFT2_processing(image_array):
	print "Calculating FFT..."
	FFT2_output = fft2(image_array)
	return FFT2_output
	
	
# Creates a function which alters the behaviour of the matplotlib function "format_coord"
# see http://matplotlib.org/examples/api/image_zcoord.html
# ONLY WORKS FOR POSITIVE AXES - otherwise alter the format coord if statement 
def formatter(mat):
	numrows,numcols = mat.shape
	def format_coord(x,y):
		col = int(x+0.5)
		row = int(y+0.5)
		if col>=0 and col<numcols and row>=0 and row<numrows:
			z = mat[row,col]
			return 'x=%1.4f, y=%1.4f, z=%1.4f'%(x,y,z)
		else:
			return  'x=%1.4f, y=%1.4f' %(x,y)
	return format_coord


# Creates a function which alters the behaviour of the matplotlib function "format_coord"
# see http://matplotlib.org/examples/api/image_zcoord.html
# Deals with negative frequencies along the axis (otherwise use "formatter()")
def formatter_FRQ(mat, frq):
	numrows,numcols = mat.shape
	neg_frq = frq - (frq*2)
	def format_coord(x,y):
		col = int(x+0.5)
		row = int(y+0.5)
		if col>=neg_frq and col<numcols and row>=neg_frq and row<numrows:
			z = mat[row,col]
			return 'x=%1.4f, y=%1.4f, z=%1.4f'%(x,y,z)
		else:
			return  'x=%1.4f, y=%1.4f' %(x,y)
	return format_coord

	
# Plots the binary and uses the "formatter" function to enable interactive viewing of x,y and z values 
# Option to alter frequenciy resolution by changing the "frq" value	
def plot_FFT_2D_interactive_z_AXIS_FREQ(FFT2_output, opath, snip_file_name, frq=50):
	plt.clf()
	magnitude = np.absolute(FFT2_output) # gives magnitude component of FFT_output
	magnitude[0,0] = 1 # magnitude [0,0] is constant - ignore and make 1
	x, y = magnitude.shape
	magnitude = np.roll(magnitude, x//2, 0) # shifts whole image to middle of axis (x//2)
	magnitude = np.roll(magnitude, y//2, 1)
	maxfreq = frq
	magnitude = magnitude[x//2 - maxfreq: x//2 + maxfreq, y//2 - maxfreq: y//2 + maxfreq]  # cuts down to within 50 frequencies
	
	#plt.imshow(np.log(magnitude)),plt.colorbar()
	plt.imshow(np.log(magnitude), extent=(-maxfreq, maxfreq, -maxfreq, maxfreq)),plt.colorbar() # gives x and y as distances from the origin
	
	time_stamp = strftime("%H.%M.%S")
	output_filename = opath + '%s_2D.png' %(snip_file_name)
	plt.gca().format_coord = formatter_FRQ(magnitude, frq)
	#plt.savefig(output_filename)
	print "Frequency about to be displayed: %i" %frq
	plt.xlabel("Freq. distance from origin")
	plt.ylabel("Freq. distance from origin")
	plt.show()		


def plot_FFT_2D_interactive_z_AXIS_WAVELENGTHS(FFT2_output, opath, snip_file_name, frq=50):
	plt.clf()
	magnitude = np.absolute(FFT2_output) # gives magnitude component of FFT_output
	magnitude[0,0] = 1 # magnitude [0,0] is constant - ignore and make 1
	x, y = magnitude.shape
	magnitude = np.roll(magnitude, x//2, 0) # shifts whole image to middle of axis (x//2)
	magnitude = np.roll(magnitude, y//2, 1)
	maxfreq = frq
	magnitude = magnitude[x//2 - maxfreq: x//2 + maxfreq, y//2 - maxfreq: y//2 + maxfreq]  # cuts down to within 50 frequencies
	
	#plt.imshow(np.log(magnitude)),plt.colorbar()
	plt.imshow(np.log(magnitude), extent=(-maxfreq/x, maxfreq/x, -maxfreq/y, maxfreq/y)),plt.colorbar() # gives x and y as distances from the origin
	
	time_stamp = strftime("%H.%M.%S")
	output_filename = opath + '%s_2D.png' %(snip_file_name)
	plt.gca().format_coord = formatter_FRQ(magnitude, frq)
	#plt.savefig(output_filename)
	print "Frequency about to be displayed: %i" %frq
	plt.xlabel("Wavelength distance from origin (m)")
	plt.ylabel("Wavelength distance from origin (m)")
	plt.show()	

############

# Frequency units not corrected from centre (i.e. not zero at centre)
def plot_FFT_2D(FFT2_output, opath):
	plt.clf()
	magnitude = np.absolute(FFT2_output) # gives magnitude component of FFT_output
	magnitude[0,0] = 1 # magnitude [0,0] is constant - ignore and make 1
	x, y = magnitude.shape
	magnitude = np.roll(magnitude, x//2, 0) # shifts whole image to middle of axis (x//2)
	magnitude = np.roll(magnitude, y//2, 1)
	maxfreq = 50
	magnitude = magnitude[x//2 - maxfreq: x//2 + maxfreq, y//2 - maxfreq: y//2 + maxfreq]  # cuts down to within 50 frequencies
	plt.imshow(np.log(magnitude)),plt.colorbar()
	time_stamp = strftime("%H.%M.%S")
	output_filename = opath + '%s_2D.png' %(snip_file_name)
	plt.gca().format_coord = format_coord
	plt.savefig(output_filename)
	plt.show()	

# Frequency scale corrected from centre	(0 at centre)
def plot_FFT_2D_axis_frequency(FFT2_output, opath, snip_file_name, frq=50):
	plt.clf()
	magnitude = np.absolute(FFT2_output) # gives magnitude component of FFT_output
	magnitude[0,0] = 1 # magnitude [0,0] is constant - ignore and make 1
	x, y = magnitude.shape
	magnitude = np.roll(magnitude, x//2, 0) # shifts whole image to middle of axis (x//2)
	magnitude = np.roll(magnitude, y//2, 1)
	maxfreq = frq
	magnitude = magnitude[x//2 - maxfreq: x//2 + maxfreq, y//2 - maxfreq: y//2 + maxfreq]  # cuts down to within 50 frequencies
	
	print "x: %i" %x
	print "y: %i" %y
	
	fig = plt.figure()
			
	plt.imshow(np.log(magnitude), extent=(-maxfreq, maxfreq, -maxfreq, maxfreq)),plt.colorbar() # gives x and y as distances from the origin
	
	time_stamp = strftime("%H.%M.%S")
	output_filename = opath + '%s_2D_AXIS_FREQ.png' %(snip_file_name)
	
	plot_title = "%s" %(snip_file_name)
	plt.title(plot_title)
	plt.xlabel("Freq. distance from origin")
	plt.ylabel("Freq. distance from origin")
			
	plt.savefig(output_filename)
	print "Magnitude maximum: %f" %(magnitude.max())
	#return magnitude
	
# Wavelength units used instead of frequency units
def plot_FFT_2D_axis_wavelength(FFT2_output, opath, snip_file_name):
	plt.clf()
	magnitude = np.absolute(FFT2_output) # gives magnitude component of FFT_output
	magnitude[0,0] = 1 # magnitude [0,0] is constant - ignore and make 1
	x, y = magnitude.shape
	magnitude = np.roll(magnitude, x//2, 0) # shifts whole image to middle of axis (x//2)
	magnitude = np.roll(magnitude, y//2, 1)
	maxfreq = 50
	magnitude = magnitude[x//2 - maxfreq: x//2 + maxfreq, y//2 - maxfreq: y//2 + maxfreq]  # cuts down to within 50 frequencies
	#fig = plt.figure()
	
	fig = plt.figure()
		
	#plt.imshow(np.log(magnitude), extent=(x/-maxfreq, x/maxfreq, y/-maxfreq, y/maxfreq)),plt.colorbar() # gives x and y as distances from the origin
	plt.imshow(np.log(magnitude), extent=(-maxfreq, maxfreq, -maxfreq, maxfreq)),plt.colorbar() # gives x and y as distances from the origin
	
	#ticks=np.arange(x/-maxfreq, x/maxfreq, 10)
	#plt.xticks(ticks)
	
	plt.xticks(range(x/-maxfreq, x/maxfreq, 10))
	#plt.yticks(range(y/-maxfreq, y/maxfreq, 5))
										
	time_stamp = strftime("%H.%M.%S")
	output_filename = opath + '%s_2D_AXIS_WAVELENGTHS.png' %(snip_file_name)
	plot_title = "%s" %(snip_file_name)
	plt.title(plot_title)
	plt.xlabel("Wavelength distance from origin (m)")
	plt.ylabel("Wavelength distance from origin (m)")
	plt.savefig(output_filename)

	
def plot_FFT_2D_filter_size_half(FFT2_output, opath, snip_file_name):
	plt.clf()
	magnitude = np.absolute(FFT2_output) # gives magnitude component of FFT_output
	magnitude[0,0] = 1 # magnitude [0,0] is constant - ignore and make 1
	x, y = magnitude.shape
	magnitude = np.roll(magnitude, x//2, 0) # shifts whole image to middle of axis (x//2)
	magnitude = np.roll(magnitude, y//2, 1)
	maxfreq = 50
	magnitude = magnitude[x//2 - maxfreq: x//2 + maxfreq, y//2 - maxfreq: y//2 + maxfreq]  # cuts down to within 50 frequencies
	fig = plt.figure()
	plt.imshow(np.log(magnitude)),plt.colorbar()
	time_stamp = strftime("%H.%M.%S")
	output_filename = opath + '%s_2D_filter_50pc.png' %(snip_file_name)
	plot_title = "%s" %(snip_file_name)
	plt.title(plot_title)
	plt.savefig(output_filename)
	#return FFT_plot
	
	
def plot_FFT_2D_filter_size_quarter(FFT2_output, opath, snip_file_name):
	plt.clf()
	magnitude = np.absolute(FFT2_output) # gives magnitude component of FFT_output
	magnitude[0,0] = 1 # magnitude [0,0] is constant - ignore and make 1
	x, y = magnitude.shape
	magnitude = np.roll(magnitude, x//2, 0) # shifts whole image to middle of axis (x//2)
	magnitude = np.roll(magnitude, y//2, 1)
	maxfreq = 50
	magnitude = magnitude[x//2 - maxfreq: x//2 + maxfreq, y//2 - maxfreq: y//2 + maxfreq]  # cuts down to within 50 frequencies
	fig = plt.figure()
	plt.imshow(np.log(magnitude)),plt.colorbar()
	time_stamp = strftime("%H.%M.%S")
	output_filename = opath + '%s_2D_filter_25pc.png' %(snip_file_name)
	plot_title = "%s" %(snip_file_name)
	plt.title(plot_title)
	plt.savefig(output_filename)
	#return FFT_plot
	
## REQUIRES NEW VERSION OF MATPLOTLIB		
def plot_FFT_3D(FFT2_output, opath):
	plt.clf()
	magnitude = np.absolute(FFT2_output) # gives magnitude component of FFT_output
	magnitude[0,0] = 1 # magnitude [0,0] is constant - ignore and make 1
	x, y = magnitude.shape
	magnitude = np.roll(magnitude, x//2, 0) # shifts whole image to middle of axis (x//2)
	magnitude = np.roll(magnitude, y//2, 1)
	maxfreq = 50
	magnitude = magnitude[x//2 - maxfreq: x//2 + maxfreq, y//2 - maxfreq: y//2 + maxfreq]  # cuts down to within 50 frequencies

	fig = plt.figure()
	ax = fig.gca(projection='3d')
	x, y = np.meshgrid(x, y)
	z = magnitude
	surface = ax.plot_surface(x,y,z,rstride=1,cstride=1,cmap=cm.coolwarm, linewidth=0, antialised=False)
	fig.colorbar(surface, shrink=0.5, aspect=5)

	output_filename = opath + '%s_3D.png' %(snip_file_name)
	plt.savefig(output_filename)
	
########

def brown_noise_surface(fft_surface, frq=50):
	
	magnitude = np.absolute(fft_surface) # gives magnitude component of FFT_output
	a = magnitude.dtype
	length = magnitude.size
	magnitude[0,0] = 1 # magnitude [0,0] is constant - ignore and make 1
	x, y = magnitude.shape
	magnitude = np.roll(magnitude, x//2, 0) # shifts whole image to middle of axis (x//2)
	magnitude = np.roll(magnitude, y//2, 1)
	magnitude = magnitude[x//2 - frq: x//2 + frq, y//2 - frq: y//2 + frq]  # cuts down to within 50 frequencies
	
	brown_noise = magnitude
	
	x_frq = frq*2
	y_frq = frq*2
		
	## i_coord and j_coord are the coordinates at a given point in the array relative to the centre of the image		
	for i in range(x_frq): 
		for j in range(y_frq):
	
	#for i in range(50): 
		#for j in range(50):
			
			print "Creating brown noise ratio surface"
			
			if (i>frq):
				i_coord = i - frq
			else:
				i_coord = (frq - i)
			
			if (j>frq):
				j_coord = j - frq
			else:
				j_coord = (frq - j)
			
			#temp[i,j] = magnitude[i,j]/i_coord**2 + j_coord**2
			pos = brown_noise[i,j]
			brown_noise[i,j] = pos * (i_coord**2 + j_coord**2)  # this is equivalent to multiplying by frequency squared
																# output is ratio of fft to brown noise fft
			#print pos
			#print brown_noise[i,j]
			#print "i: %f" %(i)
			#print "j: %f" %(j)
			#print "i (coord): %f" %(i_coord)
			#print "j (coord): %f" %(j_coord)
			#print "i^2 (coord): %f" %(i_coord**2)
			#print "j^2 (coord): %f" %(j_coord**2)
			
			#brown_noise[i,j]/(i_coord**2 + j_coord**2)
	
	print "Brown noise maximum: %f" %(brown_noise.max())
	return brown_noise

def plot_brown_noise(brown_noise, frq):
	
	x, y = brown_noise.shape
	fig = plt.figure()
	plt.imshow(np.log(brown_noise), extent=(-frq/x, frq/x, -frq/y, frq/y)),plt.colorbar()
	plt.title('Brown noise (rolled as FFT)')
	plt.gca().format_coord = formatter_FRQ(brown_noise, frq)
	plt.show()	

#def noise_fft_subtraction():

	
	
