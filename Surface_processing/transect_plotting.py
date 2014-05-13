import numpy as np
from matplotlib import pyplot as plt # for ploting
import matplotlib as mpl
site = 'Helheim'
output = r'/geog/data/sirius/epsilon/ggwillc/filter_transects/%s/%s_transect_plots' %(site, site)
transect_number = 1
filter_type = 'max'
method = 'Maximum'
window_size = 241
px = 0.5
kernel = window_size*px

data_path = r'/geog/data/sirius/epsilon/ggwillc/filter_transects/%s' %(site)

transect = "trn%i_%s_kernel_%i.txt" %(transect_number, filter_type, window_size)

transect_location = "%s/%s" %(data_path, transect)
print transect_location

# read in file
f_np = np.loadtxt(transect_location, skiprows=3, usecols = (0,1))

x = f_np [0:600,0]
y = f_np [0:600,1]


fig = plt.figure()
mpl.axes.set_default_color_cycle(['b'])
ax = fig.add_subplot(111) # add a subplot in position 111 (rows/cols/position) to figure called "ax"
ax.plot(x,y) 
ax.set_title("%s surface anomaly : Kernel %d m : Transect %i : %s " %(method, kernel, transect_number, site))
plt.xlabel('Metres along transect')
plt.ylabel('%s surface anomaly' %(method))
#name = "test" # set file name to save
#output = r'/home/staff/ggwillc/Desktop/'
#plt.savefig(output) 
plt.show() 


# line thickness

