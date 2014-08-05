import matplotlib as matplotlib
import numpy as np
import matplotlib.pyplot as plt
from numpy import *

variance = []
lag = []

image_array = array([(10.0,12.0,38.0),(5.0,19.0,13.0),(2.0,41.0,23.0),(4.0,12.0,1.0)]) # equal dimensions
#image_array = array([(10,12,38),(5,19,13,26),(2,41,23),(4,12,1)]) # unequal dimensions

for i in range (len(image_array)):
	for j in range (len(image_array[i])):		
		
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

# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# ax1.plot(lag,variance)
# plt.xlabel('Values')
# plt.ylabel('Frequency')
# plt.show()
  
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
	y_numpyarray[x_binned == i].mean() # calcs mean of values in y array when i == to the bin number 
	if len(y_numpyarray[x_binned == i]) > 0 
    else 0
    for i in range(1, len(bins)+1)]) # remove the +1 if you want to create the array "binnedData" - NOTE: this will invalidate "variance_mean_value" in the vario-function plotting due to an out_of bounds error - this is related to the way "binnedData" calculates the bin_start and bin_end values (the i+1 operator) - "i" in terms of the bin parameters is positioned towards the smaller value e.g. min<bin<max - in this instance the min == i and max == i+1

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
