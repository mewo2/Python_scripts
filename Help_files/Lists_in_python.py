###
# All about lists in python
###

from numpy import *

a = [1,2,3,4,5,6,7,8,9,10] # defines a 1D array
a.append(987) # appends 987 to array a
print a.count(987) # prints out the number of instances of 987
print a.index(987) # prints out the position of 987 in the array

# define empty lists

variance = []
lag = []
var = 0.5

# image image_array is a 2D array
# loop through image_array and for every cell, compare it to every other cell calculating variance and lag, inserting all values into lists "variance" and "lag"

image_array = array([(1,2,3),(4,5,6)])

 for i in range (len(image_array2)):
     for j in range (len(image_array2[i])):
         variance.append(image_array2[i][j])
         lag_value = image_array2[i][j] - var
         lag.append(lag_value)

print variance
print lag
		 
# loop through the lists

for i in range (len(lag)):
    print lag[i]
	
	# loop through multiple lists of same length
for i in range (len(lag)):
    print lag[i]
    print variance[i]	
	
	# loop through any lists as long as same length (lag set as the length but variance read)
for i in range (len(lag)):
    print variance[i]