import numpy as np
import random

# x and y in this instance are the other way than you might think...
# It works rows/cols (not cols/rows...)
# "0" on the left hand axis also starts at the top of the array!

z = np.random.rand(3,4)
ny, nx = z.shape

# Creates an array starting at 0, ending at value of 2, with the number of instances equal to nx (or ny) (i.e. if nx = 3, start = 0 and end = 2, you will have values of [0, 1, 2]

#x = np.linspace(0,2,nx)
#y = np.linspace(0,2,ny)

x = np.linspace(0,nx-1,nx)
y = np.linspace(0,ny-1,ny)

# Reverse the direction of the y array
y_reversed = y[::-1]

# Get value at position in z (call the y axis first)
z[y[1],x[0]]
z[y[1]x[3]]


## Then loop z 2d array to get position of max
## i relates to the y axis
for i in range(len(z)):
    for j in range(len(z[i])):
        if(z[i,j] == z.max()):
            print z[i,j]
            print y_reversed[i]
			print x[j]