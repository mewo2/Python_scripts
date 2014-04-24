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
#image_array_subsample_DATA = image_array[4790:4910, 7000:7120]
#image_array_subsample_DATA = image_array[4790:4850, 7000:7060]
#image_array_subsample_DATA = image_array[4790:4910, 7000:7120]
image_array_subsample_DATA = image_array[7119:7489, 7489:7689]

def test_function(image_array_subsample, window_iteration):

	var_accum = np.zeros(len(bins))
	var_accum_counter = np.zeros(len(bins))
	var_mean = np.zeros(len(bins))

# Loop array and accumulate values according to bin - integrate counter for each accumulation instance --> creates one list of bins and total variance

	for i in range(len(image_array_subsample)):
		for j in range(len(image_array_subsample[i])):
			pos_now_row = i
			pos_now_col = j
			value = image_array_subsample[i][j]
			#print "%i,%i" %(pos_now_row,pos_now_col) # position
			#print "value = %f" %(value)

			rows = len(image_array_subsample) # gets no. rows (count starts at 1)
			cols = len(image_array_subsample[0]) # gets no. cols (count starts at 1)
			#print "Number of rows (size starting at 1): %d" %(rows)
			#print "Number of columns (size starting at 1): %d" %(cols)

			rows2 = len(image_array_subsample)-1 # gets no. rows in terms of position (count starts at 0)
			cols2 = len(image_array_subsample[0])-1 # gets no. cols in terms of position (count starts at 0)
			#print "Number of rows (in terms of index starting at 0): %d" %(rows2)
			#print "Number of cols (in terms of index starting at 0): %d" %(cols2)

			for i2 in range (len(image_array_subsample)):
				for j2 in range (len(image_array_subsample[i2])):
				
					try: # TESTS IF STILL WITHIN COLUMN (J) LIMITS
						instance_i2 = image_array_subsample[i2]
						instance_j2 = image_array_subsample[j2]
						instance_test = 1
					except IndexError: # CATCHES THE OUTOFBOUNDS ERROR
						instance_test = 0
					
					# IF VALUE IS INSIDE OF WINDOW THEN DO SOMETHING
					# IF VALUE IS OUTSIDE OF WINDOW THEN DO SOMETHING
					
					if instance_test == 1:
						pos_now_row_2 = i2
						pos_now_col_2 = j2
						
						lag_val_pre = ((pos_now_row_2 - pos_now_row)*(pos_now_row_2 - pos_now_row))+((pos_now_col_2 - pos_now_col)*(pos_now_col_2 - pos_now_col))
						lag_val = math.sqrt(lag_val_pre)
						bin_No = int((lag_val/bin_size)-1)
						#print "bin_No: %i" %(bin_No)
						####lag.append(lag_val)
											
						var_value = value - image_array_subsample[i2][j2]
						var_value_ABSOLUTE = math.fabs(var_value)
						var_accum[bin_No] += var_value_ABSOLUTE
						var_accum_counter[bin_No] += 1
						#print "counter for bin 1: %i" %(var_accum_counter[5])
						####variance.append(var_value_ABSOLUTE)
						
						# At end of array, divide accumulations by individual counters to give mean variance and plot
					
					elif instance_test == 1: 
						print "OUT OF BOUNDS REACHED - SORT OUT A BOUNDARY CONDITION"
						
	for i in range(len(bins)):
		var_mean[i] = var_accum[i]/var_accum_counter[i]

	#### plot vario-function

	# make both arrays are of same type
	variance_mean_np = np.array(var_mean)
	bin_np = np.array(bins)

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

	# CALCULATE SLOPES BETWEEN POINTS
	
	binvar_unsorted = zip(bin_np,variance_mean_np) # combine lag and var lists in to 2 column list
	binvar_sorted = sorted(binvar_unsorted, key=lambda value: value[0]) # sort combined column list by lag (1st col[0])
	print "length of binvar_sorted = %f" %(len(binvar_sorted))

	####
	## Calc slope between values
	####
	#for i in range(len(lagvar_sorted)):

	slope_gradient_list = []
	for i in range(bin_int):
		row_value = binvar_sorted[i][0]
		col_value = binvar_sorted[i][1]
		#print "row value: %f" %(row_value)
		#print "column value: %f" %(col_value)

		# don't calculate this for i == [0]
		if i > 0 and i < bin_int:
			#print "Calculating slope between current and previous point"
			O = binvar_sorted[i][1] - binvar_sorted[i-1][1]
			A = binvar_sorted[i][0] - binvar_sorted[i-1][0]
			#print "O value: %f" %(O)
			#print "A value: %f" %(A)

			if binvar_sorted[i][0] == 0:
				print "adjacent == 0 therefore no slope angle as same point...."
			else:
				tan_x = O/A
				x_rad  = math.atan(tan_x)
				#print "tan_x value: %f" %(tan_x)
				#print "x_rad value: %f" %(x_rad)
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
	p1=-9999.
	p2=-9999.
	mindest=-9999.
	pond=-9999.

	##pond
	pond = max(variance_mean_np)
	print "pond = %f" %(max(variance_mean_np))
	
	###################
	# CHECK LISTS EXIST
	if not min_lag:
		empty_min_list = 0
	else:
		empty_min_list = 1
			
	if not max_lag:
		empty_max_list = 0
	else:
		empty_max_list = 1
	###################
	
	###################
	# CHECK POSITION IN min_lag LIST EXISTS
	try:
		min_lag[1]
		min_lag_presence = 1
	except:
		min_lag_presence = 0
	###################
	
	if empty_min_list == 1 and empty_max_list == 1:
		##mindest
		mindest = min_lag[0]  - max_lag[0]
		if mindest > 0:
			print "mindest = %f" %(mindest)
		elif mindest < 0 and min_lag_presence == 1:
			print "min_1 precedes max_1 - will use min_2_lag instead"
			mindest = min_lag[1] - max_lag[0] 
		elif mindest == 0 and min_lag_presence == 1:
			print "min_1 does not lag after max_1 - will use min_2_lag instead"
			mindest = min_lag[1] - max_lag[0] 
		elif mindest == 0 and min_lag_presence == 0:
			print "min_1 does not lag after max_1 - would use min_2_lag instead but value DOES NOT EXIST"
			mindest = mindest 
	
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
	
	return(pond,mindest,p1,p2)

'''	
## TEST ARRAY (IGNORE ALL OF THE ABOVE APART FROM THE IMPORTS AND USE THIS IF YOU FANCY)

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
	plt.imshow(datout[vv]),plt.colorbar()
	
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