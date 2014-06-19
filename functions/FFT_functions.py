'''
FFT functions
'''

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
	
	
def plot_FFT_2D_axis_frequency(FFT2_output, opath, snip_file_name):
	plt.clf()
	magnitude = np.absolute(FFT2_output) # gives magnitude component of FFT_output
	magnitude[0,0] = 1 # magnitude [0,0] is constant - ignore and make 1
	x, y = magnitude.shape
	magnitude = np.roll(magnitude, x//2, 0) # shifts whole image to middle of axis (x//2)
	magnitude = np.roll(magnitude, y//2, 1)
	maxfreq = 50
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
	##return FFT_plot
	
	
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
