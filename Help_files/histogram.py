####
#### Open file/plot it/display it
#### @author Chris 29/11/13
####

import matplotlib as matplotlib
import numpy as NP
import matplotlib.pyplot as plt

print "here"

with open('test_numbers.csv') as f:
  v = NP.loadtxt(f, delimiter=",", dtype='float', comments="#", skiprows=1, usecols=None)
  
v_hist = NP.ravel(v)   # 'flatten' v
fig = plt.figure()
ax1 = fig.add_subplot(111)

n, bins, patches = ax1.hist(v_hist, bins=100, normed=1, facecolor='green')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.show()
  
