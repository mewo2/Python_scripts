###############
Loops in PYTHON
###############

# Print a sequence of values
for i in [1,10,2,3]:
	print i

# Print a range of values
for i in range(10):
	print i

# Print every second number from a range of values	
for i in range(1,10,2):
	print i

# Print every third number from a range of values	
for i in range(1,10,3):
	print i

# Print a range of values using a decimal time step (0.5 in this case)
import numpy as np
for i in np.arange(1,10,0.5):
	print i
	
	