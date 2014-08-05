###
# All about lists in python
###

from numpy import *

# Simple list functions
a = [1,2,3,4,5,6,7,8,9,10] # defines a 1D array
a.append(987) # appends 987 to array a
print a.count(987) # prints out the number of instances of 987
print a.index(987) # prints out the position of 987 in the array


# image image_array is a 2D array
# loop through image_array and for every cell, compare it to every other cell calculating variance and lag, inserting all values into lists "variance" and "lag"

variance = [] # define empty list
lag = [] # define empty list
var = 0.5 # define variable

image_array = array([(1,2,3),(4,5,6)])

 for i in range (len(image_array2)):
     for j in range (len(image_array2[i])):
		 variance.append(image_array2[i][j])
         lag_value = image_array2[i][j] - var
         lag.append(lag_value)

print variance
print lag

###
# loop through array and get values and positions of things in the array
###

 for i in range (len(image_array2)):
     for j in range (len(image_array2[i])):
         pos_now_row = i
         pos_now_col = j
		 value = image_array[i][j]
         print "%i,%i" %(pos_now_row,pos_now_col) # position
		 print "value = %f" %(value)

###
#Print last position in array
#Useful conditional
#[len(image_array2)-1][len(image_array2[i])-1] are used as python array lengths start at 1 but array positions start at 0
###		 
for i in range (len(image_array2)):
    for j in range (len(image_array2[i])):		
		value = image_array2[len(image_array2)-1][len(image_array2[0])-1]
		print "%f" %(value)

# loop through image_array and for every cell, compare it with the next cell to calculate variance and lag, inserting all values into lists "variance" and "lag"
# includes conditional to deal with boundaries so that when at last position in column, the next cell is the first column of the next row
# works for any EQUALLY SIZED array (can be tested by declaring a bigger array) - additional boundary condition needed if cols are not equal for all rows (perhaps instead of "if position is at j max, go to first col of next row" something like "if next cell does not exist, go to first col of next row" - SEE NEXT CODE SNIPPET)

variance = []
lag = []

#image_array2 = array([(10,12,38),(5,19,13),(2,41,23),(4,12,1)])
#image_array2 = array([(10,12,38),(5,19,13,26),(2,41,23),(4,12,1)])

for i in range (len(image_array2)):
    for j in range (len(image_array2[i])):		
		
		#print image_array2[i][j]
		
		rows = len(image_array2) # gets no. rows (count starts at 1)
		cols = len(image_array2[0]) # gets no. cols (count starts at 1)
		print "%d" %(rows)
		print "%d" %(cols)
		
		rows2 = len(image_array2)-1 # gets no. rows in terms of position (count starts at 0)
		cols2 = len(image_array2[0])-1 # gets no. cols in terms of position (count starts at 0)
		print "%d" %(rows2)
		print "%d" %(cols2)
		
		# prevent out of bounds exception
		if i == len(image_array2)-1 and j == len(image_array2[0])-1:
			print "End of array reached"
		else: 
			print "Within array"
						
			if j < len(image_array2[0])-1: # within the boundaries of j
				print "j < max"
				print "i = %i" %(i)
				point_value = image_array2[i][j]			
				next_point_value = image_array2[i][j+1]
				
				rows_now = i
				cols_now = j
				cols_next = j+1
				print "Position now: %i,%i" %(rows_now, cols_now)
				print "Next position co-ordinate: %i,%i" %(rows_now,cols_next)
				
				variance_value = point_value - next_point_value
				variance.append(variance_value)
			
			if j == len(image_array2[0])-1: # IF EXCEEDING J (COLS) EXTENT
				print "j == max"
				print "i = %i" %(i)
				point_value = image_array2[i][j] # POSITION STILL AT MAX POSSIBLE POINT
				a = i+1 # MOVE TO NEXT ROW
				b = 0 # GO TO 1ST COLUMN
				next_point_value = image_array2[a][b] #GET VALUE FROM NEXT ROW | FIRST COLUMN
				print "Position now: %i,%i" %(i, j)
				print "Next position co-ordinate: %i,%i" %(a,b)
				
				variance_value = point_value - next_point_value
				variance.append(variance_value)
			
# loop through image_array and for every cell, compare it with the next cell to calculate variance and lag, inserting all values into lists "variance" and "lag"
# includes conditional to deal with boundaries so that when at last position in column, the next cell is the first column of the next row
# works for any ANY EQUALLY OR UNEQUALLY SIZED array (can be tested by declaring a bigger array) - tests if next cell is within the extent and shifts to next line and 1st column if boundary error is recorded

variance = []
lag = []

#image_array2 = array([(10,12,38),(5,19,13),(2,41,23),(4,12,1)]) # equal dimensions
image_array2 = array([(10,12,38),(5,19,13,26),(2,41,23),(4,12,1)]) # unequal dimensions

for i in range (len(image_array2)):
    for j in range (len(image_array2[i])):		
		
		#print image_array2[i][j]
		
		rows = len(image_array2) # gets no. rows (count starts at 1)
		cols = len(image_array2[0]) # gets no. cols (count starts at 1)
		print "%d" %(rows)
		print "%d" %(cols)
		
		rows2 = len(image_array2)-1 # gets no. rows in terms of position (count starts at 0)
		cols2 = len(image_array2[0])-1 # gets no. cols in terms of position (count starts at 0)
		print "%d" %(rows2)
		print "%d" %(cols2)
		
		if i == len(image_array2)-1 and j == len(image_array2[0])-1: # CHECKS TO SEE IF AT END OF ARRAY
			print "End of array reached"
		else: 
			print "Within array"
			
			try: # TESTS IF STILL WITHIN COLUMN (J) LIMITS
				next_cell = image_array2[i][j+1]
				next_cell_test = 1
			except IndexError: # CATCHES THE OUTOFBOUNDS ERROR
				next_cell_test = 0
						
			if next_cell_test == 1:  # WITHIN J EXTENT 
				print "in bounds - next cell acquired"
				print "j < max"
				print "i = %i" %(i)
				point_value = image_array2[i][j]			
				next_point_value = image_array2[i][j+1]
				
				rows_now = i
				cols_now = j
				cols_next = j+1
				print "Position now: %i,%i" %(rows_now, cols_now)
				print "Next position co-ordinate: %i,%i" %(rows_now,cols_next)
				
				variance_value = point_value - next_point_value
				variance.append(variance_value)
			
			elif next_cell_test == 0: # EXCEEDING J EXTENT
				print "out of bounds - drop to next line"
				print "j == max"
				print "i = %i" %(i)
				point_value = image_array2[i][j] # POSITION STILL AT MAX POSSIBLE POINT
				a = i+1 # MOVE TO NEXT ROW
				b = 0 # GO TO 1ST COLUMN
				next_point_value = image_array2[a][b] #GET VALUE FROM NEXT ROW | FIRST COLUMN
				print "Position now: %i,%i" %(i, j)
				print "Next position co-ordinate: %i,%i" %(a,b)
				
				variance_value = point_value - next_point_value
				variance.append(variance_value)
			
# loop through image_array and for every cell, compare it with EVERY OTHER cell to calculate variance and lag, inserting all values into lists "variance" and "lag"
# works for array of ANY dimension
# variance list is converted to absolute values

variance = []
lag = []

image_array = array([(10.0,12.0,38.0),(5.0,19.0,13.0),(2.0,41.0,23.0),(4.0,12.0,1.0)]) # equal dimensions
#image_array = array([(10,12,38),(5,19,13,26),(2,41,23),(4,12,1)]) # unequal dimensions

for i in range (len(image_array)):
	for j in range (len(image_array[i])):		
		
		pos_now_row = i
		pos_now_col = j
		value = image_array[i][j]
		print "%i,%i" %(pos_now_row,pos_now_col) # position
		print "value = %f" %(value)

		rows = len(image_array) # gets no. rows (count starts at 1)
		cols = len(image_array[0]) # gets no. cols (count starts at 1)
		print "Number of rows (size starting at 1): %d" %(rows)
		print "Number of columns (size starting at 1): %d" %(cols)

		rows2 = len(image_array)-1 # gets no. rows in terms of position (count starts at 0)
		cols2 = len(image_array[0])-1 # gets no. cols in terms of position (count starts at 0)
		print "Number of rows (in terms of index starting at 0): %d" %(rows2)
		print "Number of cols (in terms of index starting at 0): %d" %(cols2)

		for i2 in range (len(image_array)):
			for j2 in range (len(image_array[i2])):
					
				pos_now_row_2 = i2
				pos_now_col_2 = j2
					
				var_value = value - image_array[i2][j2]
				var_value_ABSOLUTE = math.fabs(var_value)
				variance.append(var_value_ABSOLUTE)
				
				lag_val_pre = ((pos_now_row_2 - pos_now_row)*(pos_now_row_2 - pos_now_row))+((pos_now_col_2 - pos_now_col)*(pos_now_col_2 - pos_now_col))
				lag_val = math.sqrt(lag_val_pre)
				lag.append(lag_val)

####
# loop through the lists
####

for i in range (len(lag)):
    print lag[i]
	
	# loop through multiple lists of same length
for i in range (len(lag)):
    print lag[i]
    print variance[i]	
	
	# loop through any lists as long as same length (lag set as the length but variance read)
for i in range (len(lag)):
    print variance[i]
	
####
# combine to lists into a list of 2 columns
####

list1 = [0,1,2,3,4,5]
list2 = [10,43,987,-5,76,2]

list_cols = zip(list1, list2)

####
# Sorting lists (using sorting() with a key)
####

# The following sorts a 2 column list by the 1st column

print sorted(list_cols, key=lambda value: value[0])sorted(list_cols, key=lambda value: value[0])

# which is the same as this:

def sort_key(value):
    return value[0]

print sorted(list_cols, key=sort_key)

# lambda just provides an alternative syntax for function definition. The result is a function object, just like the one created by def. However, there are certain things that lambda functions can't do -- like defining new variables. They're good (depending on who you ask) for creating small one-use functions, such as this one.

# Once you understand that, then all you have to know is that key accepts a function, calls it on every value in the sequence passed to sorted, and sorts the values according to the order that their corresponding key values would take if they were sorted themselves.