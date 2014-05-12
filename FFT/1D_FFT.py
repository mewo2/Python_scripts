from __future__ import division

import numpy as np
from matplotlib import pyplot as plt # for ploting

# read in file
f_np = np.loadtxt('trn2_gaussian_sigma_120px.txt', skiprows=3, usecols = (0,1))

x = f_np [:,0]
y = f_np [:,1]

fft_1d = np.absolute(np.fft.fft(y))
fft_1d_snip = fft_1d[1:100] # value 2 - 101 (value 1 is noise)
wavelengths = len(x)/np.arange(1,100)

fig = plt.figure()
ax = fig.add_subplot(111) # add a subplot in position 111 (rows/cols/position) to figure called "ax"

# plot in freq. space
ax.plot(fft_1d_snip) 

# plot wavelengths (nb/  in pixels)
#ax.plot(wavelengths,fft_1d_snip) 

ax.set_title("1D FFT")
#name = "test" # set file name to save
#output = r'/home/st relative to wavelengtaff/ggwillc/Desktop/'
#plt.savefig(output) 
plt.show() 
