####
#Python  Moving windows
####

# Example here uses scipy.ndimage.filters.generic_filter (see http://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.filters.generic_filter.html)

#Following boundary conditions are currently supported (these are set using mode='condition'):
#  “nearest” 	Use the value at the boundary 	[1 2 3]->[1 1 2 3 3]
#  “wrap” 	Periodically replicate the array 	[1 2 3]->[3 1 2 3 1]
#  “reflect” 	Reflect the array at the boundary 	[1 2 3]->[1 1 2 3 3]
#  “constant” 	Use a constant value, default is 0.0 	[1 2 3]->[0 1 2 3 0]
	
	
# Working example

import numpy as np
from scipy import ndimage

window = np.ones((3,3))
# array([[ 1.,  1.,  1.],
#       [ 1.,  1.,  1.],
#       [ 1.,  1.,  1.]])

image = np.array([(1,2,3,4),(5,6,7,8),(9,10,11,12),(13,14,15,16),(17,18,19,20)])
#array([[ 1,  2,  3,  4],
#       [ 5,  6,  7,  8],
#       [ 9, 10, 11, 12],
#       [13, 14, 15, 16],
#       [17, 18, 19, 20]])

def test_function(x):
	return(x+x).sum()
	
output_image = ndimage.generic_filter(image, test_function, footprint=window, mode ='reflect')
#array([[ 48,  60,  78,  90],
       #[ 96, 108, 126, 138],
       #[168, 180, 198, 210],
       #[240, 252, 270, 282],
       #[288, 300, 318, 330]])
	   
# The array position in 'image' takes the centre of the 3x3 moving window
# The window moves along every position in 'image'
# The function takes in all values within the window as a 1D array, does whatever to them and then outputs them as float

# Taking position image[0][0] , considering the boundary condition (reflect), the function takes in a 1D array of [1,1,2,1,2,5,5,6] 
# The function then adds this array to itself, giving an output of 48 which is then used to populate output_image[0][0]
# This 1D array is equivalent to the 2D form of:
#[1, 1, 2]
#[1, 1, 2]
#[5, 5, 6]
# Notice the effect of the boundary condition:
#[1, 1, 2]
#[1,
#[5,
# Without the boundary you are left with the values of the array and space at the edges:
#[       ]
#[   1, 2]
#[   5, 6]