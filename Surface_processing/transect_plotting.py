import numpy as np
from matplotlib import pyplot as plt # for ploting

# read in file
f_np = np.loadtxt('trn2_gaussian_sigma_120px.txt', skiprows=3, usecols = (0,1))

x = f_np [:,0]
y = f_np [:,1]


fig = plt.figure()
ax = fig.add_subplot(111) # add a subplot in position 111 (rows/cols/position) to figure called "ax"
ax.plot(x,y) 
ax.set_title("Transect")
#name = "test" # set file name to save
#output = r'/home/staff/ggwillc/Desktop/'
#plt.savefig(output) 
plt.show() 
