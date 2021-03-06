from __future__ import division

import numpy as np
from matplotlib import pyplot as plt 
import datetime
import time as time
from time import gmtime, strftime
import math


plt.clf()

# read in file
data_path = r'/geog/data/sirius/epsilon/ggwillc/FFT_1D/transects'
transect_number = 3
site = 'KNS'
winsize = 121
transect = "kns_201_median_difference_%ikernel_TRANSECT_%i_ASCII.txt" %(winsize, transect_number)
output = r'/geog/data/sirius/epsilon/ggwillc/FFT_1D/%s_FFT_1D_plots/kernel_%ipx' %(site, winsize)
pixel_size = 0.5

transect_location = "%s/%s" %(data_path, transect)
print transect_location

f_np = np.loadtxt(transect_location, skiprows=3, usecols = (0,1))

x = f_np [:,0]
y = f_np [:,1]

fft_1d = np.absolute(np.fft.fft(y))
fft_1d_snip = fft_1d[1:100] # value 2 - 101 (value 1 is noise)
wavelengths = len(x)/np.arange(1,100)

fig = plt.figure()
ax = fig.add_subplot(111) # add a subplot in position 111 (rows/cols/position) to figure called "ax"

# plot in freq. space
def frequency_plot(fft_1d_snip):
	x_label = 'Index position'
	y_label = 'Magnitude'
	plot_type = "frequency_space"
	ax.plot(fft_1d_snip) 
	return x_label, y_label, plot_type, ax
	
# plot wavelengths (nb/  in pixels)
def wavelength_plot (wavelengths,fft_1d_snip):
	x_label = 'Wavelength'
	y_label = 'Magnitude'
	plot_type = "freq_vs_wavelengths"
	ax.plot(wavelengths,fft_1d_snip) 
	return x_label, y_label, plot_type, ax
	
def magnitude_metres(fft_1d_snip, array_length):
	mag_m = fft_1d_snip/math.sqrt(array_length)
	return mag_m
	
def wavelength_metres(pixel_size, wavelengths):
	wavelength_m = wavelengths*pixel_size
	return wavelength_m

def FFT_plot(site, ax, transect_number, plot_type, output_path, x_label, y_label, magnitude_units='px', wavelength_units='px'):
	#plt.clf()
	time_stamp = strftime("%H.%M.%S")
	name = "%s_transect_%i_%s_%s.pdf" %(site, transect_number, plot_type, time_stamp) 
	#ax.set_title("Helheim_transect_%i_FFT" %(transect_number))
	output_name = "%s/%s" %(output_path, name)
	plt.xlabel("%s (%s)" %(x_label, wavelength_units))
	plt.ylabel("%s (%s)" %(y_label, magnitude_units))
	plt.savefig(output_name) 
	#plt.show() 
	

###
#~~~~~~~~~~~~~~~~~~~~~~~
###

## Units in pixels (wavelength plot)
'''
x_label, y_label, plot_type, ax = wavelength_plot(wavelengths,fft_1d_snip)
FFT_plot(ax, transect_number, plot_type, output, x_label, y_label)

## Units in metres (wavelength plot)
'''
magnitude_m = magnitude_metres(fft_1d_snip, len(x))
wavelength_m = wavelength_metres(pixel_size, wavelengths)
x_label, y_label, plot_type, ax = wavelength_plot(wavelength_m,magnitude_m)
FFT_plot(site, ax, transect_number, plot_type, output, x_label, y_label, 'm', 'm')


###
#~~~~~~~~~~~~~~~~~~~~~~~
###

## Units in pixels (frequency plot)
#x_label, y_label, plot_type, ax = frequency_plot(fft_1d_snip)
#output = r'/geog/data/sirius/epsilon/ggwillc/FFT_1D/helheim_FFT_1D_plots'
#FFT_plot(ax, transect_number, plot_type, output, x_label, y_label)

## Units in metres (frequency plot)
#x_label, y_label, plot_type, ax = frequency_plot(fft_1d_snip)
#output = r'/geog/data/sirius/epsilon/ggwillc/FFT_1D/helheim_FFT_1D_plots'
#FFT_plot(ax, transect_number, plot_type, output, x_label, y_label)

print "FFT 1D COMPLETE"

