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

# Set file location
#file_name = r"/geog/data/altair/epsilon/ggwillc/AL_ARSF_GRNLND_2013/LiDAR/201a/post_0.5/bin/201a.wholeGlacier.0.5.bin" # the r escapes numbers (and special characters like capital letters) in the pathname
#file_name = r"/geog/data/altair/epsilon/ggwillc/AL_ARSF_GRNLND_2013/LiDAR/203a/post_0.5/bin/203a.wholeGlacier.0.5.bin" # the r escapes numbers (and special characters like capital letters) in the pathname
file_name = r"/geog/data/sirius/epsilon/ggwillc/Helheim/helheim_lidar_sorting/222a_lidar/bin/222a.helheim_post_0.5m.bin.snip.100px"

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
print image_array.shape

##~~~This works to display ~~~#
#print "%s resampled by a factor of 0.25" % (image_array_name)
#dem_025_refactor = ndimage.zoom(image_array, 0.25, order=0) ## essentially refactors the array by a factor the order value denotes the interpolation algorithm used
#print dem_025_refactor
#display = plt.imshow(dem_025_refactor)
#plt.title(image_array_name + ' (resampled)')
#plt.colorbar(display)
#image_output_original =  r'/home/staff/ggwillc/Desktop/filtering_output_images/'+ 'HELHEIM_223_original_surface.png'
#plt.savefig(image_output_original)
#plt.show()

temp_window_size = 50 # gives a 25m window size

#for i in range (len(image_array)):
#	for j in range (len(image_array[i])):		
for i in range (temp_window_size):
	for j in range (temp_window_size):		
		
		pos_now_row = i
		pos_now_col = j
		value = image_array[i][j]
		print "%i,%i" %(pos_now_row,pos_now_col) # position
		print "value = %f" %(value)

		rows = len(image_array) # gets no. rows (count starts at 1)
		cols = len(image_array[0]) # gets no. cols (count starts at 1)
		print "Number of rows (size starting at 1): %d" %(rows)
		print "Number of columns (size starting at 1): %d" %(cols)

		rows2 = len(image_array)-1 # gets no. rows in terms of position (count starts at 0)
		cols2 = len(image_array[0])-1 # gets no. cols in terms of position (count starts at 0)
		print "Number of rows (in terms of index starting at 0): %d" %(rows2)
		print "Number of cols (in terms of index starting at 0): %d" %(cols2)

		for i2 in range (len(image_array)):
			for j2 in range (len(image_array[i2])):
					
				pos_now_row_2 = i2
				pos_now_col_2 = j2
					
				var_value = value - image_array[i2][j2]
				var_value_ABSOLUTE = math.fabs(var_value)
				variance.append(var_value_ABSOLUTE)
				
				lag_val_pre = ((pos_now_row_2 - pos_now_row)*(pos_now_row_2 - pos_now_row))+((pos_now_col_2 - pos_now_col)*(pos_now_col_2 - pos_now_col))
				lag_val = math.sqrt(lag_val_pre)
				lag.append(lag_val)
		
#### plot initial graph	

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(lag,variance)
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.show()
  
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
bin_no = 10.0
print 'bin_no: %d' %(bin_no)
bin_size = max_lag/bin_no
print 'bin_size: %.4f' %(bin_size)

# Returns evenly spaced values within a given interval
# numpy.digitize() will "bin" the values; i.e. x_binned[i] will give the bin number of value in i-th index
# converts variance list into a numpy array
bins = np.arange(start_bin_value, end_bin_value, bin_size) 
x_binned = np.digitize(lag,bins,right=False) 
y_numpyarray = np.array(variance) 

# creates a numpy array of mean y values per bin
# There is quite a bit of jugglery here.
# x_binned==i gives a boolean of size x_binned
# y_numpyArray[x_binned==i] returns the elements of y_numpyArray where the boolean is true
# The len() check is to make sure that mean() is not called for an empty array (which results in NAN)

y_means = np.array([
	y_numpyarray[x_binned == i].mean() 
	if len(y_numpyarray[x_binned == i]) > 0 
	else 0
	for i in range(1, len(bins)+1)]) 

## creates a three column array: |bin start|bin end|mean y| - created to the length of y_means	
#binnedData = [(bins[i], bins[i + 1], y_means[i]) for i in range(len(y_means))] 
#print 'bin_start | bin_end | mean variance'
#print binnedData

variance_mean_value = [(y_means[i]) for i in range(len(y_means))] 
print 'number of mean values: %i' %(len(variance_mean_value))
#bin_mid_point = [(bins[i] + bin_size) for i in range(len(y_means))] 

## Loop through lag and bin arrays to cross check
#for i in range(len(lag)):
#	lag_value = lag[i]
#	bin_value = x_binned[i]
#	var_value = variance[i]
#	print 'lag: %.2f	bin: %.2f	variance: %.2f' %(lag_value, bin_value, var_value)

#### plot vario-function

# make both arrays are of same type
variance_mean_np = np.array(variance_mean_value)
bin_np = np.array(bins)
#bin_mid_point_np = np.array(bin_mid_point)

plt.plot(bin_np,variance_mean_np)
plt.xlabel('Lag (binned)')
plt.ylabel('Variance')
plt.show()

# Define max_1 and min_1 variables
# max_1 = max velocity for smallest lag
# min_1 = smallest variance where lag cannot be < lag of max_1
max_1 = max(variance_mean_np) 

# NEED TO CALCULATE SLOPES BETWEEN POINTS
# MAX 1 IS WHERE  THE FIRST + -> - SLOPE CHANGE OCCURS
# MAX 2 IS WHERE THE FIRST - -> +  SLOPE CHANGE OCCURS
# SHOULD BE ABLE TO ITERATE FROM THERE
# CREATE NEW SLOPE ANGLE LIST ANDUSE TO GET VARIANCE AND LAG VALUES

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

for i in range(len(slope_gradient_list)):
#for i in range(10):
	print i
	print slope_gradient_list[i]

	try: # TESTS IF STILL WITHIN COLUMN (J) LIMITS
		next_cell = slope_gradient_list[i+1]
		next_cell_test = 1
	except IndexError: # CATCHES THE OUTOFBOUNDS ERROR
		next_cell_test = 0
		
	if 	next_cell_test == 1:
		if slope_gradient_list[i] > 0.0 and slope_gradient_list[i+1] <= 0.0:
			print "Maximum"
			var_max = binvar_sorted[i][1]
			lag_max = binvar_sorted[i][0]
			max_counter += 1
			print "Max %i: Lag = %f		Var = %f" %(max_counter,lag_max, var_max)
			max_number.append(max_counter)
			max_lag.append(lag_max)
			max_variance.append(var_max)
		if slope_gradient_list[i] < 0.0 and slope_gradient_list[i+1] > 0.0:
			print "Minimum"
			var_min = binvar_sorted[i][1]
			lag_min = binvar_sorted[i][0]
			min_counter += 1
			print "Min %i: Lag = %f		Var = %f" %(min_counter,lag_min,var_min)
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
			print "slope_gradient_list[i-1] = %f" %(slope_gradient_list[i-1])
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
			print "Max %i: Lag = %f		Var = %f" %(max_counter,lag_max, var_max)
			max_number.append(max_counter)
			max_lag.append(lag_max)
			max_variance.append(var_max)
		elif slope_gradient_list[i] < 0.0:
			print "Minimum"
			var_min = binvar_sorted[i][1]
			lag_min = binvar_sorted[i][0]
			min_counter += 1
			print "Min %i: Lag = %f		Var = %f" %(min_counter,lag_min,var_min)
			min_number.append(min_counter)
			min_lag.append(lag_min)
			min_variance.append(var_min)
		else:
			print "Plateau leading to the last point of the vario-function - no max/min value identified"
			
print "max_counter: %i"	%(max_counter)
print "min_counter: %i"	%(min_counter)
max_list = zip(max_number,max_lag,max_variance)
min_list = zip(min_number,min_lag,min_variance)

####
## Establish variables
####

print "VARIO-FUNCTION CLASSIFICATION: CALCULATING FIRST ORDER VARIABLES"

##pond
print "pond = %f" %(max(variance_mean_np))

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

print "VARIO-FUNCTION CLASSIFICATION: VARIABLE ARRAYS CALCULATING"

#pond_array = cp.copy(image_array_copy) # creates new array equal to size of initial image_array
pond_array = np.zeros((temp_window_size,temp_window_size))
pond_array.fill(pond)
#mindest_array = cp.copy(image_array_copy) # creates new array equal to size of initial image_array
mindest_array = np.zeros((temp_window_size,temp_window_size))
mindest_array.fill(mindest)
#p1_array = cp.copy(image_array_copy)  # creates new array equal to size of initial image_array
p1_array = np.zeros((temp_window_size,temp_window_size))
p1_array.fill(p1)
#p2_array = cp.copy(image_array_copy)  # creates new array equal to size of initial image_array
p2_array = np.zeros((temp_window_size,temp_window_size))
p2_array.fill(p2)

print "Pond array:"
print pond_array

print "Mindest array:"
print mindest_array

print "p1 array:"
print p1_array

print "p2 array:"
print p2_array
print "VARIO-FUNCTION CLASSIFICATION: VARIABLE ARRAYS COMPLETE"
