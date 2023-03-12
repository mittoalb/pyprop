

	
	




if __name__ == '__main__':

	"""
	Entrance matrix must be divisible by 4 without rest...
	"""

	from PyMca import EdfFile
	from math import sqrt
	import sys
	energy = 30. #in keV
	px_size = 3e-6 #in m
	distance = 2. #in m

	#fields1 = tester(px_size,energy)
	#fields1, d, b = testSphere(energy,300,46,px_size)
	#energy = 3.3 #in keV
	#px_size = 0.098e-6 #in m
	#distance = 0.000119 #in m

	pro = propagators(px_size,energy,distance,1)
	#fields1 = tester(px_size,energy)
	
	
	fields1 = simplewire()
	#write_edf("objectwire.edf",np.imag(fields1))
	#sys.exit()
	#Read speckle field data
	#fields1 = np.loadtxt('/data/id17/inhouse/MITTONE/SPECKLES/3DG/Projection_gold_1e8_3um.txt')
	x,y = fields1.shape
	#Reduce effect
	#fields1 = fields1*0.2 + np.ones((x,y),dtype=float)
	#fields1 /= (np.max(np.max(fields1))+0.0)	
	#Generate statistical noise

	#fields1 = tester(pixel_size)
	#x,y = fields1.shape


	

	#field1 = data.read()
	#fields1 = fields1[:-1,:-1]
	print fields1.shape
	
	pro.preparefresnel2D(fields1)
	#Single energy
	#fieldIntensity, fieldPhase = pro.fpro()
	
	#Polychromatic beam
	numberofevents = 1000
	path = "./"
	MD = 1e-4
	radix = "testwire_"
	savefileevery = 1000
	fieldIntensity = pro.polybeamMC(fields1,numberofevents,MD,path,radix,savefileevery)
	
	#blurring
	fieldIntensity = scipy.ndimage.filters.gaussian_filter(fieldIntensity, 1.0)	
	
	#with noise
	photons = 5000
	mu, sigma = photons, sqrt(photons)
	noise = np.random.normal(mu, sigma, x*y)
	noise = np.reshape(noise,(x,y))
	#fieldIntensity *= noise
	


	#test = np.asarray(fieldPhase)
	#print test.shape
	write_edf("testwire.edf",fieldIntensity) #path for saving propagated intensity
	#write_edf("/data/id17/inhouse/MITTONE/PBISIM/D90keVTheory.edf",d) #path for saving propagated intensity
	
	""" #NOT USED BY M.Eckermann#
	#Test polychromaticity - parameters
	numberofevents = 1000
	#Level of monochromaticity
	MD = 1e-4
	path = "./simpleobject_1_volume/"
	radix = "Propagated_PMMA_gold_1e8_3um_1em4_F11m_S0p6m_GaussNoise_Thick0p01cm_"
	savefileevery = 1000
	
	#Launch calculation
	distance = 0.6
	propagatedspeckles = pro.polybeamMC(fields1,numberofevents,MD,path,radix,savefileevery)
	
	fields1 = np.zeros((x,y))
	#Read speckle field data
	fields = np.loadtxt('/data/id17/inhouse/MITTONE/SPECKLES/3DG/Projection_sandpaper_1e8_3um.txt')
	fields1[0:125,0:125] = fields #pad speckle data for 500x500 px images
	fields1[0:125,125:250] = fields
	fields1[0:125,250:375] = fields
	fields1[0:125,375:500] = fields
	fields1[125:250,:] = fields1[0:125,:]
	fields1[250:375,:] = fields1[0:125,:]
	fields1[375:500,:] = fields1[0:125,:]
	#fields1 = tester(px_size)
	x,y = fields1.shape
	#Reduce effect
	#fields1 = fields1*0.2 + np.ones((x,y),dtype=float)
	#fields1 /= (np.max(np.max(fields1))+0.0)

	pro2 = propagators(px_size,energy,distance,1)
	
	propagatedspeckles /= (np.max(np.max(np.real(propagatedspeckles)))+0.0)	
	fields1 *= propagatedspeckles #Add fields
	
	pro2.preparefresnel2D(np.real(fields1))
	fieldIntensity, fieldPhase = pro2.fpro()
	done = pro2.polybeamMC(fields1,numberofevents,MD,path,radix,savefileevery)
	write_edf("./simpleobject_1_volume/Detected_intensity_simpleobject_1_ex.edf",np.asarray(fieldIntensity,dtype=np.float))
	write_edf("./simpleobject_1_volume/Detected_phase_simpleobject_1_ex.edf",np.asarray(fieldPhase,dtype=np.float))
	from matplotlib import pylab as plt
	from scipy.special import jv

	plt.figure(1)
	plt.imshow(fieldIntensity,cmap='Greys_r')
	plt.title("Fresnel Diffraction")
	plt.xlabel("X [um]")
	plt.ylabel("Intensity [a.u.]")
	plt.show()
	
	plt.figure(2)
	plt.plot(fieldIntensity[400,:])
	plt.title("Fresnel Diffraction")
	plt.xlabel("X [um]")
	plt.ylabel("Intensity [a.u.]")
	plt.show()
	"""

