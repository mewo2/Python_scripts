#####
## numpy - arrays etc.
#####

import numpy as np

# Print a 1D array of values 0 - 19
a1 = arange(20)
a1
# Print a 2D array of values 0 - 19 comprising 4 rows and 5 cols
a2 = arange(20).reshape(4,5)
a2
# Print a 1D array of 20 1's
b1 = np.ones(20)
b1
# Print a 2D array of 20 1's comprising 4 rows and 5 cols
b2 = np.ones(20).reshape(4,5)
b2

# Print a 2D array of 1's (20 rows and 20 cols)
c1 = np.ones((20,20))
# Print a 2D array of of values ranging 0 - 19
c2 = np.ones((20,20)) * np.arange(20) #(np.ones((20,20)).arange(20) won't work as it calls for a method np.ones().arange() which doesn't exist
#array([[  0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 11.,  12.,  13.,  14.,  15.,  16.,  17.,  18.,  19.],			#row 0
 #      [  0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10., 11.,  12.,  13.,  14.,  15.,  16.,  17.,  18.,  19.],
	#   .....
     #  [  0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,         11.,  12.,  13.,  14.,  15.,  16.,  17.,  18.,  19.]])	#row 19
	 
####
# Numpy arrays can be used instead of lists also....
# see: http://stackoverflow.com/questions/568962/how-do-i-create-an-empty-array-matrix-in-numpy
# NB_1/ You need to know how long the list will be before it can be used however
# NB_2/ You probably also need a counter so that you can assign the array position to which data is to be stored
####

# Create a np.ndarray full of random values
import random
import numpy as np
from numpy import *

image_array_subsample = np.ndarray(shape = (20,20)) * random.random() # random values - could include NaN instances
#or
image_array_subsample = np.ndarray(shape = (20,20)) * random.randrange(0,80) # random values between 0 and 80
#see https://docs.python.org/2/library/random.html

~~~OR~~~
# Create 10x10 array of random numbers between 0.0 and 1.0
z = np.random.rand(10,10)

