import fabio


#Read an input data
def readata(filename):
	return img.data, img.header = fabio.open(filename)

#write an output data
def writedata(img,filename,ext='edf'):
	img.save(filename)
