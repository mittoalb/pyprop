import fabio
from libtiff import TIFF
import glymur


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
	
	
def read_tif(name):
	"""
	Read a tif file
	"""
	tif = TIFF.open(name, mode='r')
	return tif.read_image()
	
	
def write_tif(name,data):
	"""
	Write a tif file
	"""
	tif = TIFF.open(name, mode='w')
	tif.write_image(data)
	
	
def read_jp2(name):
	"""
	Read a JPG2000 file
	"""
	jp2 = glymur.Jp2k(name)
	return jp2	
	
def write_jp2(name,data):
	"""
	Write a JPG2000 file
	"""
	jp2=glymur.Jp2k(name, data)			
