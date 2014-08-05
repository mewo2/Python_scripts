####
# Python: quick plots of arrays
####

from matplotlib import pyplot as plt # for ploting
import random
import numpy as np
from numpy import *

image_array_subsample = np.ndarray(shape = (20,20)) * random.random() # create a random array
plt.clf() # clear any previous plots
display = plt.imshow(image_array_subsample) # attach plot to a label
plt.colorbar(display) # attach a stretched color bar legend
plt.show()