import os
import numpy as np # maths functions (arrays etc.)
from matplotlib import pyplot as plt # for ploting

# Create 2 vectors
x = (1,2,3,4,5,6,7,8,9,10)
y = (1*1,2*2,3*3,4*4,5*5,6*6,7*7,8*8,9*9,10*10)

# Create a figure (easily used for multiple figures)
fig = plt.figure()
ax = fig.add_subplot(111) # add a subplot in position 111 (rows/cols/position) to figure called "ax"
ax.plot(x,y) # create xy plot using predefined vectors
ax.set_title("Example plot")
name = "test" # set file name to save
output = r'/home/staff/ggwillc/Desktop/filtering_output_images/' + name + '.png'
# set output directory
plt.savefig(output) # save image
plt.show() # display image on screen
