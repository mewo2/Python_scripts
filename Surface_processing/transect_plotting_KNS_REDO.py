import numpy as np
import os
import sys
from matplotlib import pyplot as plt # for ploting
import matplotlib as mpl
import pylab

site = 'KNS'
output = r'/geog/data/sirius/epsilon/ggwillc/filter_transects/%s/%s_transect_plots' %(site, site)
transect_number = 3
method = 'Gaussian' # Maximum Median
filter_type = 'gaussian' # max median
window_size = 120 # 241
px = 0.5


def gaussian_settings(window_size):
	kernel = window_size
	return kernel

	
def normal_window_settings(window_size, px):
	kernel =  window_size*px
	return kernel

## ~~~~~~~~~~~~~~~~~~~

#Gaussian

kernel = gaussian_settings(window_size)
transect = "trn%i_%s_sigma_%ipx_REDO.txt" %(transect_number, filter_type, window_size)
name = "%s_surface_anomaly__Kernel_%dm__Transect_%i__%s_REDO.pdf" %(method, kernel, transect_number, site)
plot_title = "%s difference surface : Kernel %d m : Transect %i : %s" %(method, kernel, transect_number, site)

#Max
'''
kernel = normal_window_settings(window_size, px)
transect = "trn%i_%s_kernel_%i_REDO.txt" %(transect_number, filter_type, window_size)
name = "%s_surface_anomaly__Kernel_%dm__Transect_%i__%s_REDO.pdf" %(method, kernel, transect_number, site)
plot_title = "%s difference surface : Kernel %d m : Transect %i : %s" %(method, kernel, transect_number, site)
'''
#Median
'''
kernel = normal_window_settings(window_size, px)
transect = "trn%i_%s_kernel_%i_NEG_REDO.txt" %(transect_number, filter_type, window_size)
print transect
name = "%s_surface_anomaly__Kernel_%dm__Transect_%i__%s_REDO.pdf" %(method, kernel, transect_number, site)
plot_title = "%s negative surface anomaly : Kernel %d m : Transect %i : %s" %(method, kernel, transect_number, site)
'''

#Median pos neg
'''
kernel = normal_window_settings(window_size, px)
transect = "posneg_diff_trn%i_%s_%ipx_REDO.txt" %(transect_number, filter_type, window_size)
name = "%s_difference_surface__Kernel_%dm__Transect_%i__%s_REDO.pdf" %(method, kernel, transect_number, site)
plot_title = "%s difference surface : Kernel %d m : Transect %i : %s" %(method, kernel, transect_number, site)
'''
## ~~~~~~~~~~~~~~~~~~~

data_path = r'/geog/data/sirius/epsilon/ggwillc/filter_transects/%s' %(site)
transect_location = "%s/%s" %(data_path, transect)

if os.path.isfile(transect_location):
	print "file exists"
else:
	print "%s DOESN'T exist...\n" %(transect_location)
	sys.exit("File doesn't exist - check in code syntax")
	
# read in file
f_np = np.loadtxt(transect_location, skiprows=3, usecols = (0,1))

x = f_np [0:600,0]
y = f_np [0:600,1]
x_m = f_np [0:600,0]/2

fig = plt.figure(figsize=(14,5))
mpl.axes.set_default_color_cycle(['b'])
ax = fig.add_subplot(111) # add a subplot in position 111 (rows/cols/position) to figure called "ax"
ax.plot(x_m,y) 

ax.set_title(plot_title)

plt.xlabel('Distance along transect (m)')
plt.ylabel('%s surface anomaly (m)' %(method))

pylab.ylim([-40,25])

output_path_complete = "%s/%s" %(output,name)
plt.savefig(output_path_complete, format='pdf') 
#plt.show() 

print "Plot saved --- script complete"

# line thickness

