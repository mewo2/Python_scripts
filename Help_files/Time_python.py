########
### Time in python
########

#~~~~~~~~~~~~~
#Clock time:
#~~~~~~~~~~~~~

#GMT time:

from time import gmtime, strftime
strftime("%Y-%m-%d %H:%M:%S", gmtime())
# '2009-01-05 21:14:39'

#Local time:

from time import gmtime, strftime
strftime("%Y-%m-%d %H:%M:%S")
# '2009-01-05 23:14:39'

#~~~~~~~~~~~~~
#Script timer:
#~~~~~~~~~~~~~

import time as time 

# start timing
startTime = time.time()

## YOUR CODE HERE

endTime = time.time()
print "Script was running for: %f" %(endTime)

