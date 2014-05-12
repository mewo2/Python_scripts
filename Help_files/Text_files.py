f = open('trn2_gaussian_sigma_120px.txt', 'r')  #  open file as "read"
f.read()     # read out the file into memory
f.close()	 # close the read
data = f.read() # read into object called data
data.split('\n')[3:] # split the file at every line (\n) and ignore the first 3
lines


see also:

import numpy as np

#np.loadtxt()
f_np = np.loadtxt('trn2_gaussian_sigma_120px.txt', skiprows=3, usecols = (0,1)
f_np.shape
print f_np
