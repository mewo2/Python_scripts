import os
import time as time # for reading in a timer
import numpy as np # maths functions (arrays etc.)
from osgeo import gdal, gdalconst # for reading in raster
from osgeo.gdalconst import * # for reading in raster
from matplotlib import pyplot as plt # for ploting
from scipy import signal # for convolution function
from scipy import ndimage # for resampling image
from matplotlib import cm # colour mapping
import matplotlib as matplotlib
import matplotlib.pyplot as plt
import copy as cp
from numpy import *

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
#file_name = r"/geog/data/sirius/epsilon/ggwillc/Helheim/helheim_lidar_sorting/222a_lidar/bin/222a.helheim_post_0.5m.bin"
#file_name = r"/geog/data/sirius/epsilon/ggwillc/Helheim/helheim_lidar_sorting/223-_lidar/bin/223-.helheim_post_0.5m.bin"

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

# Set pixel offset.....
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
#image_array_subsample = image_array[4790:4910, 7000:7120]
#image_array_subsample = image_array[4790:4850, 7000:7060]
image_array_subsample = image_array[4790:4850, 7000:7060]


##for i in range (len(image_array)):
##for j in range (len(image_array[i])):
#temp_window_size = 10 # 50 gives a 25m window size
for i in range (len(image_array_subsample)):
	for j in range (len(image_array_subsample[i])):
		pos_now_row = i
		pos_now_col = j
		value = image_array_subsample[i][j]
		print "%i,%i" %(pos_now_row,pos_now_col) # position
		print "value = %f" %(value)

		rows = len(image_array_subsample) # gets no. rows (count starts at 1)
		cols = len(image_array_subsample[0]) # gets no. cols (count starts at 1)
		print "Number of rows (size starting at 1): %d" %(rows)
		print "Number of columns (size starting at 1): %d" %(cols)

		rows2 = len(image_array_subsample)-1 # gets no. rows in terms of position (count starts at 0)
		cols2 = len(image_array_subsample[0])-1 # gets no. cols in terms of position (count starts at 0)
		print "Number of rows (in terms of index starting at 0): %d" %(rows2)
		print "Number of cols (in terms of index starting at 0): %d" %(cols2)

		for i2 in range (len(image_array_subsample)):
			for j2 in range (len(image_array_subsample[i2])):
					
				pos_now_row_2 = i2
				pos_now_col_2 = j2
					
				var_value = value - image_array_subsample[i2][j2]
				var_value_ABSOLUTE = math.fabs(var_value)
				variance.append(var_value_ABSOLUTE)
				
				lag_val_pre = ((pos_now_row_2 - pos_now_row)*(pos_now_row_2 - pos_now_row))+((pos_now_col_2 - pos_now_col)*(pos_now_col_2 - pos_now_col))
				lag_val = math.sqrt(lag_val_pre)
				lag.append(lag_val)

#### plot initial graph
#plt.clf()
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.plot(lag,variance)
#plt.xlabel('Values')
#plt.ylabel('Frequency')
#test_plot_1 =  r'/home/staff/ggwillc/Desktop/variofunction_test_outputs/'+ 'test_plot_1.png'
#plt.savefig(test_plot_1)
#plt.show()

#### Loop through lists, set bins and create binned lists

for i in range (len(lag)):
	
	#print out list values
	variance_value = variance[i]
	lag_value = lag[i]
	print 'lag: %.2f variance: %.2f' %(lag_value, variance_value)

#define bins
start_bin_value = 0.0
max_lag = max(lag)
end_bin_value = max_lag
print 'max_lag: %.4f' %(max_lag)
bin_no = 20.0
print 'bin_no: %d' %(bin_no)
bin_size = (max_lag/bin_no)
print 'bin_size: %.4f' %(bin_size)

# Returns evenly spaced values within a given interval
# numpy.digitize() will "bin" the values; i.e. x_binned[i] will give the bin number of value in i-th index
# converts variance list into a numpy array
bins = np.arange(start_bin_value, end_bin_value, bin_size)
x_binned = np.digitize(lag,bins)#,right=False) <-- deprecated from version on UNIX at Swansea?
y_numpyarray = np.array(variance)

# creates a numpy array of mean y values per bin
# There is quite a bit of jugglery here.
# x_binned==i gives a boolean of size x_binned
# y_numpyArray[x_binned==i] returns the elements of y_numpyArray where the boolean is true
# The len() check is to make sure that mean() is not called for an empty array (which results in NAN)

# **************************************************************************************************#
# ONLY WORKS FOR NEWER VERSION OF PYTHON - OTHERWISE USE HARDWIRED VERSION BELOW TO POPULATE Y_MEANS 
# **************************************************************************************************#

#y_means = np.array([
#	y_numpyarray[x_binned == i].mean()
#	if len(y_numpyarray[x_binned == i]) > 0
#	else 0
#	for i in range(1, len(bins)+1)])
# remove the +1 if you want to create the array "binnedData" - NOTE: this will invalidate "variance_mean_value" in the vario-function plotting due to an out_of bounds error - this is related to the way "binnedData" calculates the bin_start and bin_end values (the i+1 operator) - "i" in terms of the bin parameters is positioned towards the smaller value e.g. min<bin<max - in this instance the min == i and max == i+1

y_means = np.array([
	(y_numpyarray[x_binned == 1].mean()),
	(y_numpyarray[x_binned == 2].mean()),
	(y_numpyarray[x_binned == 3].mean()),
	(y_numpyarray[x_binned == 4].mean()),
	(y_numpyarray[x_binned == 5].mean()),
	(y_numpyarray[x_binned == 6].mean()),
	(y_numpyarray[x_binned == 7].mean()),
	(y_numpyarray[x_binned == 8].mean()),
	(y_numpyarray[x_binned == 9].mean()),
	(y_numpyarray[x_binned == 10].mean()),
	(y_numpyarray[x_binned == 11].mean()),
	(y_numpyarray[x_binned == 12].mean()),
	(y_numpyarray[x_binned == 13].mean()),
	(y_numpyarray[x_binned == 14].mean()),
	(y_numpyarray[x_binned == 15].mean()),
	(y_numpyarray[x_binned == 16].mean()),
	(y_numpyarray[x_binned == 17].mean()),
	(y_numpyarray[x_binned == 18].mean()),
	(y_numpyarray[x_binned == 19].mean()),
	(y_numpyarray[x_binned == 20].mean())])

variance_mean_value = [(y_means[i]) for i in range(len(y_means))]
print 'number of mean values: %i' %(len(variance_mean_value))

## Loop through lag and bin arrays to cross check
for i in range(len(lag)):
	lag_value = lag[i]
	bin_value = x_binned[i]
	var_value = variance[i]
	print 'lag: %.2f, bin: %.2f, variance: %.2f' %(lag_value, bin_value, var_value)

#### plot vario-function

# make both arrays are of same type
variance_mean_np = np.array(variance_mean_value)
#variance_mean_np = np.array(y_means)
bin_np = np.array(bins)
#bin_mid_point_np = np.array(bin_mid_point)
plt.clf()
plt.plot(bin_np,variance_mean_np)
plt.xlabel('Lag (binned)')
plt.ylabel('Variance')
test_vario_plot =  r'/home/staff/ggwillc/Desktop/variofunction_test_outputs/'+ 'test_vario_plot.png'
plt.savefig(test_vario_plot)
#plt.show()

########### WORKS TO HERE 15:17 09/04/14 #############

# CALCULATE SLOPES BETWEEN POINTS
# MAX 1 IS WHERE  THE FIRST + -> - SLOPE CHANGE OCCURS
# MAX 2 IS WHERE THE FIRST - -> +  SLOPE CHANGE OCCURS

binvar_unsorted = zip(bin_np,variance_mean_np) # combine lag and var lists in to 2 column list
binvar_sorted = sorted(binvar_unsorted, key=lambda value: value[0]) # sort combined column list by lag (1st col[0])

####
## Calc slope between values
####
#for i in range(len(lagvar_sorted)):

slope_gradient_list = []
for i in range(10):
	row_value = binvar_sorted[i][0]
	col_value = binvar_sorted[i][1]
	print "row value: %f" %(row_value)
	print "column value: %f" %(col_value)

	# don't calculate this for i == [0]
	if i > 0 & i < 15:
		print "Calculating slope between current and previous point"
		O = binvar_sorted[i][1] - binvar_sorted[i-1][1]
		A = binvar_sorted[i][0] - binvar_sorted[i-1][0]
		print "O value: %f" %(O)
		print "A value: %f" %(A)

		if binvar_sorted[i][0] == 0:
			print "adjacent == 0 therefore no slope angle as same point...."
		else:
			tan_x = O/A
			x_rad  = math.atan(tan_x)
			print "tan_x value: %f" %(tan_x)
			print "x_rad value: %f" %(x_rad)
			slope_gradient_list.append(x_rad)
	else:
		print "First point in binned lag/variance list therefore no slope to calculate"
		slope_gradient_list.append(0.0)

####
# Check for max AND min values
# NB/ slope of [i] is that between [i-1] and [i] - so to check for a slope up to or down to [i], you check slope at position [i] in array - not [i-1]!!!
####

max_counter = 0
min_counter = 0

max_number = []
max_lag = []
max_variance = []

min_number = []
min_lag = []
min_variance = []


#########################

for i in range(len(slope_gradient_list)):
#for i in range(10):
	print i
	print slope_gradient_list[i]

	try: # TESTS IF STILL WITHIN COLUMN (J) LIMITS
		next_cell = slope_gradient_list[i+1]
		next_cell_test = 1
	except IndexError: # CATCHES THE OUTOFBOUNDS ERROR
		next_cell_test = 0

	if next_cell_test == 1:
		if slope_gradient_list[i] > 0.0 and slope_gradient_list[i+1] <= 0.0:
			print "Maximum"
			var_max = binvar_sorted[i][1]
			lag_max = binvar_sorted[i][0]
			max_counter += 1
			print "Max %i: Lag = %f, Var = %f" %(max_counter,lag_max, var_max)
			max_number.append(max_counter)
			max_lag.append(lag_max)
			max_variance.append(var_max)
		if slope_gradient_list[i] < 0.0 and slope_gradient_list[i+1] > 0.0:
			print "Minimum"
			var_min = binvar_sorted[i][1]
			lag_min = binvar_sorted[i][0]
			min_counter += 1
			print "Min %i: Lag = %f, Var = %f" %(min_counter,lag_min,var_min)
			min_number.append(min_counter)
			min_lag.append(lag_min)
			min_variance.append(var_min)
		elif slope_gradient_list[i+1] > 0.0:
			print "Not Maximum: gradient after point is positive"
			print "slope_gradient_list[i+1] = %f" %(slope_gradient_list[i+1])
		elif slope_gradient_list[i] < 0.0:
			print "i: %i" %(i)
			print "slope_gradient_list[i]: %f" %(slope_gradient_list[i])
			print "slope_gradient_list[i-1]: %f" %(slope_gradient_list[i-1])
			print "slope_gradient_list[i+1]: %f" %(slope_gradient_list[i+1])
			print "slope_gradient_list[i-1]: %f" %(slope_gradient_list[i-1])
			print "Not Maximum: gradient leading up to point is negative"
		elif slope_gradient_list[i] == 0.0:
			print "No gradient leading up to point - plateau in variance"
		else:
			print "Not Maximum"
		
	else:
		print "out of bounds - min/max values here are at the greatest lag of the vario-function"
		if slope_gradient_list[i] > 0.0:
			print "Maximum"
			var_max = binvar_sorted[i][1]
			lag_max = binvar_sorted[i][0]
			max_counter += 1
			print "Max %i: Lag = %f, Var = %f" %(max_counter,lag_max, var_max)
			max_number.append(max_counter)
			max_lag.append(lag_max)
			max_variance.append(var_max)
		elif slope_gradient_list[i] < 0.0:
			print "Minimum"
			var_min = binvar_sorted[i][1]
			lag_min = binvar_sorted[i][0]
			min_counter += 1
			print "Min %i: Lag = %f, Var = %f" %(min_counter,lag_min,var_min)
			min_number.append(min_counter)
			min_lag.append(lag_min)
			min_variance.append(var_min)
		else:
			print "Plateau leading to the last point of the vario-function - no max/min value identified"
			
print "max_counter: %i" %(max_counter)
print "min_counter: %i" %(min_counter)
max_list = zip(max_number,max_lag,max_variance)
min_list = zip(min_number,min_lag,min_variance)

##############

####
## Establish variables
####

print "VARIO-FUNCTION CLASSIFICATION: CALCULATING FIRST ORDER VARIABLES"

##pond
print "pond = %f" %(max(variance_mean_np))

try:
	min_lag[0]
	empty_min_list = 0
except Exception:
	empty_min_list = 1
	
try:
	max_lag[0]
	empty_max_list = 0
except Exception:
	empty_max_list = 1
	
if empty_min_list == 0 and empty_max_list == 0:
	##mindest
	mindest = min_lag[0]  - max_lag[0]
	if mindest > 0:
		print "mindest = %f" %(mindest)
	elif mindest < 0:
		print "min_1 precedes max_1 - will use min_2_lag instead"
		mindest = min_lag[1] - max_lag[0] 
	elif mindest == 0:
		print "min_1 does not lag after max_1 - will use min_2_lag instead"
		mindest = min_lag[1] - max_lag[0] 

	##p1
	p1 = (max_variance[0] - min_variance[0]) / (min_lag[0] - max_lag[0])
	print "p1 = %f" %(p1)

	##p2
	p2 = (max_variance[0] - min_variance[0]) / max_variance[0]
	print "p2 = %f" %(p2)

elif empty_min_list == 0 and empty_max_list == 1:
	print "WARNING: mindest, p1 and p2 variables not calculated - not enough maximum values"
	print "CHECK DATA - WIDEN WINDOW DIMENSIONS?"
elif empty_min_list == 1 and empty_max_list == 0:
	print "WARNING: mindest, p1 and p2 variables not calculated - not enough minimum values"
	print "CHECK DATA"
	print "CHECK DATA - WIDEN WINDOW DIMENSIONS?"
elif empty_min_list == 1 and empty_max_list == 1:
	print "WARNING: mindest, p1 and p2 variables not calculated - no maximum or minimum values available from vario-function"
	print "CHECK DATA - WIDEN WINDOW DIMENSIONS?"