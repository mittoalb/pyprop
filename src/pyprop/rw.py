import fabio
#Read an input data
def readata(filename):
	return fabio.open(filename)

#write an output data
def writedata(img,filename,ext='edf'):
	img.save(filename)
	
	
def next_greater_power_of_2(val): 
	"""
	Compute the next greater power of 2
	"""
	return 2**(val-1).bit_length()
