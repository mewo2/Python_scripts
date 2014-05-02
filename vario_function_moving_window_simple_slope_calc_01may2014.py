from __future__ import division

import os
import time as time # for reading in a timer
import numpy as np # maths functions (arrays etc.)
import math
from matplotlib import pyplot as plt # for ploting
from scipy import signal # for convolution function
from scipy import ndimage # for resampling image
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
from osgeo import gdal, gdalconst # for reading in raster
from osgeo.gdalconst import * # for reading in raster

# start timing
startTime = time.time()

variance = []
lag = []

# start timing
startTime = time.time()

# Register driver
#gdal.AllRegister() #<-- useful only if reading in 
driver = gdal.GetDriverByName('ENVI') ## http://www.gdal.org/formats_list.html
driver.Register()

# Clear any previous plots
plt.clf()

# Set file location
#file_name = r"/geog/data/altair/epsilon/ggwillc/AL_ARSF_GRNLND_2013/LiDAR/201a/post_0.5/bin/dem_median_filter_kernel_121_crevasse_surface"
#file_name = r"/geog/data/sirius/epsilon/ggwillc/Maximum_surface_filtering/Helheim/222/HELHEIM_222a_dem_maximum_filter_kernel_239_20_percent_reduction_crevasse_surface_20_percent_reduction" 
file_name = r"/geog/data/arcturus/epsilon/ggwillc/Maximum_surface_filtering/Helheim/222/HELHEIM_222a_dem_maximum_filter_kernel_239_20_percent_reduction_crevasse_surface_20_percent_reduction_SUBSAMPLE"

# open file
inds = gdal.Open(file_name, GA_ReadOnly)

if inds is None:
	print "Really sorry Sir but I couldn't open this blasted file: " + file_name
	print '\nPerhaps you need an ENVI .hdr file? If so, just open the binary up in ENVI and one will be created for you!'
	os._exit(1)
else:
	print "%s opened successfully" %file_name
	
print '~~~~~~~~~~~~~~'
print 'Get image size'
print '~~~~~~~~~~~~~~'
cols = inds.RasterXSize
rows = inds.RasterYSize
bands = inds.RasterCount

print "columns: %i" %cols
print "rows: %i" %rows
print "bands: %i" %bands

print '~~~~~~~~~~~~~~'
print 'Get georeference information'
print '~~~~~~~~~~~~~~'
geotransform = inds.GetGeoTransform()
originX = geotransform[0]
originY = geotransform[3]
pixelWidth = geotransform[1]
pixelHeight = geotransform[5]

print "origin x: %i" %originX
print "origin y: %i" %originY
print "width: %2.2f" %pixelWidth
print "height: %2.2f" %pixelHeight

print '~~~~~~~~~~~~~~' 
print 'Convert image to 2D array'
print '~~~~~~~~~~~~~~'

band = inds.GetRasterBand(1)
image_array = band.ReadAsArray(0, 0, cols, rows)
image_array_name = file_name
print type(image_array)
print shape(image_array)

print '~~~~~~~~~~~~~~' 
print 'Subsample 2D array'
print '~~~~~~~~~~~~~~'
#image_array_subsample_DATA = image_array[7119:7219, 7219:7319] ##  [i1:i2, j1:j2]
#image_array_subsample_DATA = image_array[3230:4230, 6117:7117] ## HELHEIM
image_array_subsample_DATA = image_array[5117:8665, 1209:6028] ## HELHEIM
print type(image_array_subsample_DATA)
print shape(image_array_subsample_DATA)

print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
print 'CHECK OUTPUT DIRECTORY EXISTS'
print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

winsize = 479
stepsize = 1000
nsample_input = 1000000

opath = r'/geog/data/sirius/epsilon/ggwillc/vario_outputs/subsample_tests/'
opath = r'/geog/data/sirius/epsilon/ggwillc/vario_outputs/subsample_tests/winsize_%i_stepsize_%i_random_hits_%i/' %(winsize, stepsize, nsample_input)

if os.path.isdir(opath):
	print "output_path exists"	
else:
	print "output_path DOESN'T exist...\n"
	os.makedirs(opath)
	print "...but it does now"
	

def variogram(image_array_subsample, bins):
	var_accum = np.zeros(len(bins))
	var_accum_counter = np.zeros(len(bins))
	var_mean = np.zeros(len(bins))

# Loop array and accumulate values according to bin - integrate counter for each accumulation instance --> creates one list of bins and total variance
	isize, jsize = image_array_subsample.shape
	for i in xrange(isize):
		for j in xrange(jsize):
			value = image_array_subsample[i,j]
			img_diff_abs = (image_array_subsample - value)**2 # calc difference between every position in array and the value (done in one calculation)
			for i2 in xrange(i, isize):
				for j2 in xrange(jsize):
					lag_val = ((i2 - i) ** 2 + (j2 - j) ** 2) ** .5
					bin_No = int((lag_val/bin_size)-1)
					var_accum[bin_No] += img_diff_abs[i2,j2] # get the value of img_diff_abs at position [i2,j2]
					var_accum_counter[bin_No] += 1
					
	# At end of array, divide accumulations by individual counters to give mean variance and plot
	var_mean = var_accum / np.maximum(var_accum_counter, 1)
	return var_mean


def random_variogram(image_array_subsample, bins, nsamples=10000):
	var_accum = np.zeros(len(bins))
	var_accum_counter = np.zeros(len(bins))
	var_mean = np.zeros(len(bins))

	isize, jsize = image_array_subsample.shape
	#print "isize: %i" %(isize)
	#print "jsize: %i" %(jsize)
	for k in xrange(nsamples):
		i = np.random.randint(isize)
		j = np.random.randint(jsize)
		i2 = np.random.randint(isize)
		j2 = np.random.randint(jsize)
		
		lag_val = ((i2 - i) ** 2 + (j2 - j) ** 2) ** .5
		bin_No = int((lag_val/bin_size)-1)
		var_accum[bin_No] += (image_array_subsample[i,j] - image_array_subsample[i2,j2])**2
		var_accum_counter[bin_No] += 1
			
	var_mean = var_accum / np.maximum(var_accum_counter, 1)
	return var_mean


def vario_stats(bins, var_mean):
	diffs = np.diff(var_mean) # calc difference between values in the array
	signdiffs = np.diff(np.sign(diffs)) # get the sign of the differences
	maxes = np.where(signdiffs < 0)[0] + 1 # if diff between values for the position after that in which we are in is negative, then we have a maximum
	mins = np.where(signdiffs > 0)[0] + 1 # if diff between values for the position after that in which we are in is positive, then we have a minimum
	try:
		max1 = maxes[0]
		min1 = mins[mins > max1][0]
		max_var = var_mean[max1]
		max_lag = bins[max1]
		min_var = var_mean[min1]
		min_lag = bins[min1]
	except IndexError:
		max_var = -9999.0
		max_lag = -9999.0
		min_var = -9999.0
		min_lag = -9999.0
		return 0, 0, 0, 0
	# pond, mindest, p1, p2
	assert max_var >= min_var, "Max %f, min %f" % (max_var, min_var)
	assert min_lag >= max_lag, "Max %f, min %f" % (max_lag, min_lag)
	return max_var, min_lag - max_lag, (max_var - min_var) / (min_lag - max_lag), (max_var - min_var) / max_var


def plot_variogram(window_iteration, bins, variance_mean, ii, jj):
	plt.clf()
	plt.plot(bins,variance_mean)
	#plt.title('Vario-plot: array of length %i' %(len(image_array_subsample)))
	plt.title('Vario-plot %i | pos (%i, %i)' %(window_iteration,ii,jj))
	#plt.title('pos (%f, %f)' %(ii,jj))
	plt.xlabel('Lag (binned)')
	plt.ylabel('Variance')
	time_stamp = strftime("%H.%M.%S")
	#vario_plot =  r'/home/staff/ggwillc/Desktop/variofunction_test_outputs/'+ 'test_vario_plot_%s.png' %(time_stamp)
	vario_plot_windows =  r'/geog/data/sirius/epsilon/ggwillc/vario_outputs/'+ 'vario_plot_RANDOM_%s_window_iteration_%s.png' %(time_stamp,window_iteration)
	plt.savefig(vario_plot_windows)
	#plt.show()

## TEST ARRAY (IGNORE ALL OF THE ABOVE APART FROM THE IMPORTS AND USE THIS IF YOU FANCY)
'''
image_array_subsample_DATA = pl.rand(10,10) #np.ndarray(shape = (30,30)) * random.randrange(10,80,1) # random values between 0 and 80
indNan=np.isnan(image_array_subsample_DATA)
image_array_subsample_DATA[indNan]=10.
'''
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET MOVING WINDOW UP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

#winsize = 479#119#239

if winsize > 1 and winsize%2 != 0:
	print "WINDOW SIZE VALID"
else:
	print "WINDOW SIZE INVALID - MUST BE AN ODD NUMBER AND > 1"

nxwind, nywind = image_array_subsample_DATA.shape

#~~~~~~~~~~~~
# DEFINE BINS
#~~~~~~~~~~~~

start_bin_value = 0.0

# Calculate maximum lag possible (OF MOVING WINDOW - assumes window has equal dimensions in x and y)
max_lag = math.sqrt(2)*winsize

end_bin_value = max_lag
print 'max_lag: %.4f' %(max_lag)
bin_no = 40.0
bin_int = int(bin_no)
print 'bin_no: %d' %(bin_no)
bin_size = (max_lag/bin_no)
print 'bin_size: %.4f' %(bin_size)

# Returns evenly spaced values within a given interval
# numpy.digitize() will "bin" the values; i.e. x_binned[i] will give the bin number of value in i-th index
bins = np.arange(start_bin_value, end_bin_value, bin_size)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IMPLEMENT LOOP FOR MOVING WINDOW TO WORK THROUGH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nvars = 4 # number of output variables
# datout=np.ones((nvars,nxwind,nywind))

#stepsize = 1000

imgout = [np.zeros((1 + nxwind // stepsize, 1 + nywind // stepsize)) for i in xrange(nvars)]
window_iteration=0
for ii in range(0, nxwind, stepsize):
	print "ii: %i" %(ii)
	for jj in range(0, nywind, stepsize):
		
		#print "ii: %i" %(ii)
		#print "jj: %i" %(jj)
		
		# CALCULATE MAX AND MIN RANGES OF ROWS AND COLS THAT CAN BE ACCESSED BY THE WINDOW
		imin=max(0,ii-winsize/2) # gets the maximum of either 0 or ii-winsize/2...
		imax=min(nxwind-1,ii+winsize/2)+1
		jmin=max(0,jj-winsize/2)
		jmax=min(nywind-1,jj+winsize/2)+1
		
		'''
		print ii,imin,imax,jj,jmin,jmax # shows row and column ranges being accessed
		print "imin: %i" %(imin)
		print "imax: %i" %(imax)
		print "jmin: %i" %(jmin)
		print "jmax: %i" %(jmax)
		'''
		
		### Cell ii,jj is at top left of array
		#datwind=image_array_subsample_DATA[ii:winsize*(ii+1),jj:winsize*(jj+1)]
		### Cell ii,jj is at centre of array
		datwind=image_array_subsample_DATA[imin:imax,jmin:jmax] # IF CELLS ARE OUTSIDE OF THE ARRAY, THEY TAKE THE VALUES OF THE OTHER END OF THE ARRAY
		##print datwind.max(),datwind.min() 
		
		'''
		print "Max value in data window: %f" %(datwind.max())
		print "Min value in data window: %f" %(datwind.min())
		'''
		
		window_iteration += 1
		
		#var_mean = variogram(datwind, bins)
		#nsample_input = 1000000
		var_mean = random_variogram(datwind, bins, nsample_input)
		vtuple = vario_stats(bins, var_mean)
		
		# print "v1 (pond): %f" %(v1)
		# print "v2 (mindest): %f" %(v2)
		# print "v3 (p1): %f" %(v3)
		# print "v4 (p2): %f" %(v4)
		# pond_LIST.append(v1)
		# mindest_LIST.append(v2)
		# p1_LIST.append(v3)
		# p2_LIST.append(v4)
		
		#pond_array[ii][jj] = v1
		
		'''
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		print "MOVING WINDOW ARRAY OUTPUT CALCULATION"
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		'''
			
		#print "pond array test value (same position as in big data set): %i" %(array_test)
		
		for vv in range(nvars):
			imgout[vv][ii // stepsize, jj // stepsize] = vtuple[vv]
		
		plot_variogram(window_iteration, bins, var_mean, ii, jj)
		
# SAVE pond, mindest, p1 and p2 arrays (values given to centre cell of moving window)
for vv in range(nvars):
	fig = plt.figure()
	plt.imshow(imgout[vv]),plt.colorbar()
	
	if vv == 0:
		fig.suptitle('Pond (centre cell)')
	elif vv == 1:
		fig.suptitle('Mindest (centre cell)')
	elif vv == 2:
		fig.suptitle('p1 (centre cell)')
	elif vv == 3:
		fig.suptitle('p2 (centre cell)')

	time_stamp = strftime("%H.%M.%S")	
	image_output_centre_cell =  opath + str(vv+1) + '_centre_cell_output_%s_RANDOM_winsize_%i_stepsize_%i_nsamples_%i.png' %(time_stamp, winsize, stepsize, nsample_input)
	#plt.show()
	fig.savefig(image_output_centre_cell)

# SAVE pond, mindest, p1 and p2 arrays (values averaged for cells across moving window)

'''
for vv in range(nvars):			

	v = "v%i" %(vv+1)
	variable_value = v
	
	fig = plt.figure()
	
	if vv is 0:
		for i in range(len(pond_array)):
			for j in range(len(pond_array)):
				if pond_array[i][j] == -9999.:
					pond_array[i][j] = 0
				else:
					pond_array[i][j] = pond_array[i][j] 
					
		plt.imshow(pond_array),plt.colorbar()
		fig.suptitle('Pond (moving window mean)')
	elif vv is 1:
		plt.imshow(mindest_array),plt.colorbar()
		fig.suptitle('Mindest (moving window mean)')
	elif vv is 2:
		plt.imshow(p1_array),plt.colorbar()
		fig.suptitle('p1 (moving window mean)')
	elif vv is 3:
		plt.imshow(p2_array),plt.colorbar()
		fig.suptitle('p2 (moving window mean)')
	
	time_stamp = strftime("%H.%M.%S")
	image_output_moving_window = r'/geog/data/sirius/epsilon/ggwillc/vario_outputs/' + variable_value + '_moving_window_mean_output_%s_ASSERT_samples_%i.png' %(time_stamp, nsample_input)
	#plt.show()
	fig.savefig(image_output_moving_window)
	
'''
##~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Accumulate all value windows and calculate mean surface
# Add conditional to prevent averaging where one or more values has a value of -9999
##~~~~~~~~~~~~~~~~~~~~~~~~~~~

print "Cell hits: %i" %(window_iteration)
