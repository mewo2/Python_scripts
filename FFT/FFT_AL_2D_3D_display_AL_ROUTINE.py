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
from time import gmtime, strftime
import pylab as pl
from glob import glob
from osgeo import gdal, gdalconst # for reading in raster
from osgeo.gdalconst import * # for reading in raster
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

#http://stackoverflow.com/questions/4383571/importing-files-from-different-folder-in-python
import sys
#sys.path.insert(0, '/home/staff/ggwillc/Desktop/Python_scripts/functions')
sys.path.insert(0, 'C:/Users/ggwillc/Desktop/Python_scripts/functions')
import FFT_functions as FFT_functions
import raster_functions as raster_functions	
	
#import home.staff.ggwillc.Desktop.Python_scripts.functions.FFT_functions

# start timing
startTime = time.time()

## Set file location
roi_number = 4
file = "helheim_222a_sample_ELEVATION_ROI_%i_no_filtering_fft_output_no_filtering.bin" %roi_number
#for file_name in glob("*ELEV*output*butter*1*50_percent.bin"):
#for file_name in glob("*MAX*output*butter*1*50_percent.bin"):
#for file_name in glob("helheim_222a_sample_ELEVATION_ROI_2_tophat_high_pass_fft_output_tophat_filter_1_tenth_image.bin"):
#for file_name in glob("*50_percent.bin"):
#for file_name in glob("*25_percent.bin"):

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
	
	inds, cols, rows, bands, originX, originY, pixelWidth, pixelHeight, image_array, image_array_name = raster_functions.ENVI_raster_binary_to_2d_array(file_name)
			
	print '~~~~~~~~~~~~~~' 
	print 'Subsample 2D array'
	print '~~~~~~~~~~~~~~'
	#image_array_subsample_DATA = image_array[7119:7219, 7219:7319] ##  [i1:i2, j1:j2]
	#image_array_subsample_DATA = image_array[3230:4230, 6117:7117] ## HELHEIM
	#image_array_subsample_DATA = image_array[5117:8665, 1209:6028] ## HELHEIM
	#print type(image_array_subsample_DATA)
	#print shape(image_array_subsample_DATA)

	print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	print 'CHECK OUTPUT DIRECTORY EXISTS'
	print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

	#opath = r'/geog/data/sirius/epsilon/ggwillc/vario_outputs/subsample_tests/'
	#opath = r'/geog/data/sirius/epsilon/ggwillc/vario_outputs/subsample_tests/winsize_%i_stepsize_%i_random_hits_%i/' %(winsize, stepsize, nsample_input)
	#opath = r'/geog/data/sirius/epsilon/ggwillc/FFT2/tests2/'
	
	#opath = r'/geog/data/sirius/epsilon/ggwillc/FFT_2D/Helheim/AL_FFT_outputs/high_pass/heavy_crevasse_ROI/FFT_surface_plots/'
	#opath = r'/geog/data/sirius/epsilon/ggwillc/FFT_2D/Helheim/AL_FFT_outputs/high_pass/lesser_crevasse_1_ROI/FFT_surface_plots/'
	#opath = r'/geog/data/sirius/epsilon/ggwillc/FFT_2D/Helheim/AL_FFT_outputs/high_pass/lesser_crevasse_2_ROI/FFT_surface_plots/'
	
	#opath = r'/home/staff/ggwillc/Desktop/FFT_image_tests/'
	#opath = r'/home/staff/ggwillc/Desktop/'
	#opath = r'C:/Users/ggwillc/Desktop/FFT/FFT_OUTPUTS/'
	
	#opath = r'C:/Users/ggwillc/Desktop/FFT/FFT_OUTPUTS_AL_ROUTINE/'
	opath = r'C:/Users/ggwillc/Desktop/'
	
	if os.path.isdir(opath):
		print "output_path exists"	
	else:
		print "output_path DOESN'T exist...\n"
		os.makedirs(opath)
		print "...but it does now"		

	FFT_surface = image_array

	print "brown noise!!!!!!!"
	
	#plot_title = "%s %s" %(full_snip_1, full_snip_2)
	plot_title = "ROI %i | cols: %i | rows: %i\n [FREQ IMAGE ROLLED]" %(roi_number, cols, rows)
	print plot_title
	
	
	#brown_noise = FFT_functions.brown_noise_surface(FFT_surface, 300)
	##FFT_functions.plot_brown_noise(brown_noise, 350)
	#FFT_functions.save_brown_noise(brown_noise, 300, opath, plot_title, snip_file_name)
	FFT_functions.plot_FFT_2D_axis_frequency(FFT_surface, opath, snip_file_name, plot_title, 300)

#	FFT_functions.plot_FFT_2D_axis_wavelength(FFT_surface, opath, snip_file_name)

	## OLD
	##FFT_functions.plot_FFT_2D(FFT_surface, opath, snip_file_name)
	##FFT_functions.plot_FFT_2D_filter_size_half(FFT_surface, opath, snip_file_name)
	##FFT_functions.plot_FFT_2D_filter_size_quarter(FFT_surface, opath, snip_file_name)
	#FFT_functions.plot_FFT_3D(FFT_surface, opath)

	print '~~~~~~~~~~~~~~'
	print 'Clear variables'
	print '~~~~~~~~~~~~~~'
 
	band = None
	image_array = None
