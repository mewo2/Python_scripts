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
file_name = r"/geog/data/altair/epsilon/ggwillc/AL_ARSF_GRNLND_2013/LiDAR/201a/post_0.5/bin/dem_median_filter_kernel_121_crevasse_surface"

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
image_array_subsample_DATA = image_array[7119:7219, 7219:7319]

def test_function(image_array_subsample, window_iteration):

	var_accum = np.zeros(len(bins))
	var_accum_counter = np.zeros(len(bins))
	var_mean = np.zeros(len(bins))

# Loop array and accumulate values according to bin - integrate counter for each accumulation instance --> creates one list of bins and total variance
	isize, jsize = image_array_subsample.shape
	for i in xrange(isize):
		for j in xrange(jsize):
			value = image_array_subsample[i,j]
			img_diff_abs = np.abs(image_array_subsample - value)
			for i2 in xrange(i, isize):
				for j2 in xrange(jsize):
					lag_val = ((i2 - i) ** 2 + (j2 - j) ** 2) ** .5
					bin_No = int((lag_val/bin_size)-1)
					#print "bin_No: %i" %(bin_No)
					####lag.append(lag_val)
										
					var_accum[bin_No] += img_diff_abs[i2,j2]
					var_accum_counter[bin_No] += 1
					#print "counter for bin 1: %i" %(var_accum_counter[5])
					####variance.append(var_value_ABSOLUTE)
					
					# At end of array, divide accumulations by individual counters to give mean variance and plot
	
	var_mean = var_accum / np.maximum(var_accum_counter, 1)

	#### plot vario-function

	# make both arrays are of same type
	bin_np = np.array(bins)
	'''
	### PLOT VARIO-FUNCTION
	plt.clf()
	plt.plot(bin_np,variance_mean_np)
	#plt.title('Vario-plot: array of length %i' %(len(image_array_subsample)))
	plt.title('Vario-plot %i | Array length: %i ' %(window_iteration, len(image_array_subsample)))
	plt.xlabel('Lag (binned)')
	plt.ylabel('Variance')
	time_stamp = strftime("%H.%M.%S")
	#vario_plot =  r'/home/staff/ggwillc/Desktop/variofunction_test_outputs/'+ 'test_vario_plot_%s.png' %(time_stamp)
	vario_plot_windows =  r'/geog/data/sirius/epsilon/ggwillc/vario_outputs/'+ 'vario_plot_NEW_%s_window_iteration_%s.png' %(time_stamp,window_iteration)
	plt.savefig(vario_plot_windows)
	#plt.show()
	'''
	# CALCULATE SLOPES BETWEEN POINTS
	
	diffs = var_mean[1:] - var_mean[:-1]
	signdiffs = np.diff(np.sign(diffs))
	maxes = np.where(signdiffs > 0)[0] + 1
	mins = np.where(signdiffs < 0)[0] + 1
	try:
		max1 = maxes[0]
		min1 = mins[mins > max1][0]
		max_var = var_mean[max1]
		max_lag = bin_np[max1]
		min_var = var_mean[min1]
		min_lag = bin_np[min1]
	except IndexError:
		max_var = -9999.0
		max_lag = -9999.0
		min_var = -9999.0
		min_lag = -9999.0
		return 0, 0, 0, 0
	print "Real values!"
	return max_var, min_lag - max_lag, (max_var - min_var) / (min_lag - max_lag), (max_var - min_var) / max_var
		
## TEST ARRAY (IGNORE ALL OF THE ABOVE APART FROM THE IMPORTS AND USE THIS IF YOU FANCY)
'''
image_array_subsample_DATA = pl.rand(10,10) #np.ndarray(shape = (30,30)) * random.randrange(10,80,1) # random values between 0 and 80
indNan=np.isnan(image_array_subsample_DATA)
image_array_subsample_DATA[indNan]=10.
'''
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SET MOVING WINDOW UP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

winsize = 60

if winsize > 1 and winsize%2 != 0:
	print "WINDOW SIZE VALID"
else:
	print "WINDOW SIZE INVALID - MUST BE AN ODD NUMBER AND > 1"
	
nywind=len(image_array_subsample_DATA)#)sizeofdatay/winsize
nxwind=len(image_array_subsample_DATA)#sizeofdatax/winsize

#~~~~~~~~~~~~
# DEFINE BINS
#~~~~~~~~~~~~

start_bin_value = 0.0
## Calculate maximum lag possible (assumes window has equal dimensions in x and y)
##max_lag = math.sqrt(((((len(image_array_subsample_DATA)) - 1) - 0)*(((len(image_array_subsample_DATA)) - 1) - 0)) + ((((len(image_array_subsample_DATA)) - 1) - 0)*(((len(image_array_subsample_DATA)) - 1) - 0)))

# Calculate maximum lag possible (OF MOVING WINDOW - assumes window has equal dimensions in x and y)
#max_lag = math.sqrt((((winsize - 1) - 0)*((winsize - 1) - 0)) + (((winsize - 1) - 0)*((winsize - 1) - 0)))
max_lag = math.sqrt(((winsize-0)*(winsize-0))+((winsize-0)*(winsize-0)))

end_bin_value = max_lag
print 'max_lag: %.4f' %(max_lag)
bin_no = 20.0
bin_int = int(bin_no)
print 'bin_no: %d' %(bin_no)
bin_size = (max_lag/bin_no)
print 'bin_size: %.4f' %(bin_size)

# Returns evenly spaced values within a given interval
# numpy.digitize() will "bin" the values; i.e. x_binned[i] will give the bin number of value in i-th index
bins = np.arange(start_bin_value, end_bin_value, bin_size)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAKE ARRAYS OF -9999 FOR POND, MINDEST, P1 AND P2 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pond_array = np.ones([len(image_array_subsample_DATA),len(image_array_subsample_DATA)])
pond_array = pond_array - 10000
mindest_array = np.ones([len(image_array_subsample_DATA),len(image_array_subsample_DATA)])
mindest_array = mindest_array - 10000
p1_array = np.ones([len(image_array_subsample_DATA),len(image_array_subsample_DATA)])
p1_array = p1_array - 10000
p2_array = np.ones([len(image_array_subsample_DATA),len(image_array_subsample_DATA)])
p2_array = p2_array - 10000

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAKE EMPTY LISTS FOR POND, MINDEST, P1 AND P2 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pond_LIST = []
mindest_LIST = []
p1_LIST = []
p2_LIST = []

pond_check_LIST = []
mindest_check_LIST = []
p1_check_LIST = []
p2_check_LIST = []

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# IMPLEMENT LOOP FOR MOVING WINDOW TO WORK THROUGH
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nvars = 4 # number of output variables
# datout=np.ones((nvars,nxwind,nywind))
imgout = [np.zeros((nxwind, nywind)) for i in xrange(nvars)]
window_iteration=0

for ii in range(nxwind):
	for jj in range(nywind):
		print "###########################################"
		print "###########################################"
		print "~~~~~~~~~~~~~~~~NEW WINDOW~~~~~~~~~~~~~~~~~"
		print "###########################################"
		print "###########################################"
		
		print "ii: %i" %(ii)
		print "jj: %i" %(jj)
		
		# CALCULATE MAX AND MIN RANGES OF ROWS AND COLS THAT CAN BE ACCESSED BY THE WINDOW
		imin=max(0,ii-winsize/2) # gets the maximum of either 0 or ii-winsize/2...
		imax=min(nxwind-1,ii+winsize/2)+1
		jmin=max(0,jj-winsize/2)
		jmax=min(nywind-1,jj+winsize/2)+1
		
		print ii,imin,imax,jj,jmin,jmax # shows row and column ranges being accessed
		print "imin: %i" %(imin)
		print "imax: %i" %(imax)
		print "jmin: %i" %(jmin)
		print "jmax: %i" %(jmax)
		
		### Cell ii,jj is at top left of array
		#datwind=image_array_subsample_DATA[ii:winsize*(ii+1),jj:winsize*(jj+1)]
		### Cell ii,jj is at centre of array
		datwind=image_array_subsample_DATA[imin:imax,jmin:jmax] # IF CELLS ARE OUTSIDE OF THE ARRAY, THEY TAKE THE VALUES OF THE OTHER END OF THE ARRAY
		##print datwind.max(),datwind.min() 
		print "Max value in data window: %f" %(datwind.max())
		print "Min value in data window: %f" %(datwind.min())
		
		window_iteration += 1
		vtuple = test_function(datwind, window_iteration)
		
		# print "v1 (pond): %f" %(v1)
		# print "v2 (mindest): %f" %(v2)
		# print "v3 (p1): %f" %(v3)
		# print "v4 (p2): %f" %(v4)
		# pond_LIST.append(v1)
		# mindest_LIST.append(v2)
		# p1_LIST.append(v3)
		# p2_LIST.append(v4)
		
		#pond_array[ii][jj] = v1
		
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		print "MOVING WINDOW ARRAY OUTPUT CALCULATION"
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
		
		i_diff = imax-imin
		print "idiff: %f" %(i_diff)
		j_diff = jmax-jmin
		print "jdiff: %f" %(j_diff)
		
		array_test = pond_array[ii][jj]
		print "pond array test value (same position as in big data set): %i" %(array_test)
		
		for vv in range(nvars):
			imgout[vv][ii, jj] = vtuple[vv]
		print "###########################################"
		print "###########################################"
		print "###########################################"

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
	image_output_centre_cell = r'/geog/data/sirius/epsilon/ggwillc/vario_outputs/' + str(vv+1) + '_centre_cell_output_%s.png' %(time_stamp)
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
	image_output_moving_window = r'/geog/data/sirius/epsilon/ggwillc/vario_outputs/' + variable_value + '_moving_window_mean_output_%s.png' %(time_stamp)
	#plt.show()
	fig.savefig(image_output_moving_window)
	
'''
##~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Accumulate all value windows and calculate mean surface
# Add conditional to prevent averaging where one or more values has a value of -9999
##~~~~~~~~~~~~~~~~~~~~~~~~~~~

print "Cell hits: %i" %(window_iteration)
