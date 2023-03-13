import math
import numpy as np


def simplewire():
	#50 um thickness wire
	datafieldR = np.zeros(shape=(400,400,64),dtype=float)#Refraction
	datafieldA = np.zeros(shape=(400,400,64),dtype=float)#Absorption
	field1 = np.ones(shape=(400,400),dtype=float)#Absorption
	wireradius = 8 #in pixels of 3 um
	usamp = 8 #oversampling ratio
	#Fill the layers
	pixelsize = 3.0e-6 #in m
	wl = (12.398/30.0)*1e-10
	WN = (2*np.pi) / wl
	
	for i in range(0,wireradius*usamp):
		rr = int(math.sqrt(wireradius**2-(8-(i+0.0)/usamp)**2))	
		ffactor = math.sqrt(wireradius**2-(8-(i+0.0)/usamp)**2) - rr #Parzial pixel calculation
		datafieldR[:,200-rr:200+rr,i] = 2.94e-7 * (ffactor * (pixelsize/usamp))#PMMA AT 30 keV
		datafieldA[:,200-rr:200+rr,i] = 1.18e-10 * ((4.*np.pi)/wl) * ffactor * (pixelsize/usamp)#PMMA AT 30 keV
		
	datafieldR = np.sum(datafieldR,axis=2)*2
	datafieldA = np.sum(datafieldA,axis=2)*2
	#print datafieldR
	transm = np.exp(-datafieldA+1j*datafieldR*WN)
	
	return field1*transm
	
def creategrid():
	xps = 200
	yps = 200
	grid = np.zeros(shape=[xps,yps],dtype=float)
	
	for i in range(0,xps,5):
		for j in range(0,yps,5):
			grid[i:i+3,j:j+3] = 1
			grid[i+3:i+4,j+3:j+4]=0
	return grid



"""Simple Sphere"""
def testSphere(energy,size,R,pixelsize):
	#size = 200
	from math import sin, asin
	field = np.ones(shape=(size,size),dtype=complex)
	delta = np.zeros(shape=(size,size),dtype=float)
	beta = np.zeros(shape=(size,size),dtype=float)
	wl = (12.398/energy)*1e-10
	WN = (2*np.pi) / wl
	#print WN
	for i in range(0,size):
		for j in range(0,size):
			if sqrt((i-size/2)**2+(j-size/2)**2) < R: #Values for PMMA at 30 keV
				rfactor = abs(sqrt((i-size/2)**2+(j-size/2)**2))
				#print rfactor
				tk = 2.*sin(np.pi/2.-asin(rfactor/(0.0+R)))#*thickness
				if i==size/2 and j==size/2:
					tk = tmp
				tmp = tk
				#print tk
				delta[i,j] = 6.70778E-07 * tk * pixelsize#2.94029E-07
				beta[i,j] = (3.39441E-10 * (4.*np.pi)/wl)*tk*pixelsize #100 um thickness1.17987E-10

				
	transm = np.exp(-beta-1j*delta*WN)		
	field *= transm			
	return field, delta*WN, beta

"""
Create field for testing purposes
"""
def tester(pixel_size,energy):
	xnpoints = 1200
	ynpoints = 1200

	detector_sizeX = pixel_size * xnpoints
	detector_sizeY = pixel_size * ynpoints

	position_x = np.linspace(-detector_sizeX/2,detector_sizeX/2,xnpoints)
	position_y = np.linspace(-detector_sizeY/2,detector_sizeY/2,ynpoints)

	fields1 = np.ones(shape=[xnpoints,ynpoints],dtype=complex)
	wl = (12.398/energy)*1e-10
	img = fabio.open('/data/id17/inhouse/MITTONE/LINALGEBRA/simpleobject_3_volume/Projected_beta_Sphere_3300eV_CORRECT.edf'
	
	#edf=EdfFile.EdfFile('/data/id17/inhouse/MITTONE/LINALGEBRA/simpleobject_3_volume/Projected_beta_Sphere_3300eV_CORRECT.edf',access='r')
	LAC_=edf.GetData(0)
	LAC = np.pad(img.data,((400,400),(400,400)),'edge') #following three lines: for making image quadratic
	edf=EdfFile.EdfFile('/data/id17/inhouse/MITTONE/LINALGEBRA/simpleobject_3_volume/Projected_delta_Sphere_3300eV_CORRECT.edf',access='r')
	delta_=edf.GetData(0)
	delta = np.pad(delta_,((400,400),(400,400)),'edge')
	#for NRMG: do the following padding:
	#LAC = np.pad(LAC, ((1448,1449),(1448,1449)),"edge")
	#delta = np.pad(delta, ((1448,1449),(1448,1449)),"edge")
	WN = (2*np.pi) / wl
	#delta = delta*WN
	print wl, WN
	transm = np.exp(-LAC-1j*delta*WN)

	fields1 *= transm#2e-7 ->Absorption
	return fields1	
	
def loadMyModel(energy):
	wl = (12.398/30.0)*1e-10
	wl1 = (12.398/energy)*1e-10
	WN = (2*np.pi) / wl
	delta, h = read_edf('RescaledTESLARe.edf')
	beta, h = read_edf('RescaledGAUSS.edf')
	#Rescale for energy - WRONG FOR BETA!!!
	delta = delta * (wl1/wl)**2
	beta = beta * (wl1/wl)**4
	
	sy, sx = delta.shape
	field = np.ones(shape=(sy,sx),dtype=complex)
	#transm = np.exp(-beta-1j*delta*WN)
	transm = np.exp(-beta-1j*delta)#*WN)		
	field *= transm			
	return field
