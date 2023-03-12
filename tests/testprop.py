"""
author: A. Mittone
last modification 05/09/2017
PBI simulator

ADJUST:
- some settings in lines 289-291
- path where to save in line 327
- some settings in lines 248-264
"""
import numpy as np
import scipy.ndimage
import copy
import random
from utility import write_edf, read_edf
from math import sqrt
import sys
#import scipy
from scipy import signal
import PBIpropagator as pb

sys.path.insert(0, '/data/id17/map/ID17_Rec/GITrepo/dev_v1.4/other/')
import db
#from PBIpropagator import propagator


	
######################################################################################
#TEST FIELDS
######################################################################################
"""
Create field for testing purposes

def tester(pixel_size,energy):
	xnpoints = 1200
	ynpoints = 1200

	detector_sizeX = pixel_size * xnpoints
	detector_sizeY = pixel_size * ynpoints

	position_x = np.linspace(-detector_sizeX/2,detector_sizeX/2,xnpoints)
	position_y = np.linspace(-detector_sizeY/2,detector_sizeY/2,ynpoints)

	fields1 = np.ones(shape=[xnpoints,ynpoints],dtype=complex)
	wl = (12.398/energy)*1e-10
	edf=EdfFile.EdfFile('/data/id17/inhouse/MITTONE/LINALGEBRA/simpleobject_3_volume/Projected_beta_Sphere_3300eV_CORRECT.edf',access='r')
	LAC_=edf.GetData(0)
	LAC = np.pad(LAC_,((400,400),(400,400)),'edge') #following three lines: for making image quadratic
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
"""
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
	#PMMA 90 keV
	#delta = 3.28e-8
	#beta = 2.2e-11	
	#PMMA 30 keV
	#delta = 2.94e-7
	#beta = 1.18e-10
	#PMMA 70 keV
	#delta = 7.41e-8
	#beta = 3.77e-11
	#PMMA 70 keV
	#delta = 5.41423E-08
	#beta = 3.05746E-11			
	#PMMA 20 keV
	#delta = 6.70778E-07
	#beta = 3.39441E-10
		

				
	transm = np.exp(-beta-1j*delta*WN)		
	field *= transm			
	return field, delta*WN, beta

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
	
	
def simplewire(energy,pixelsize):
	import math
	#50 um thickness wire
	wireRadiusR = 50.0e-6 #in m
	
	wireradius = int(wireRadiusR/pixelsize) #in pixels of 3 um
	print "wire radius approximated to", wireradius * pixelsize * 1e6, "um"
	usamp = 8 #oversampling ratio

	datafieldR = np.zeros(shape=(400,400,wireradius*usamp),dtype=float)#Refraction
	datafieldA = np.zeros(shape=(400,400,wireradius*usamp),dtype=float)#Absorption
	field1 = np.ones(shape=(400,400),dtype=float)#Absorption

	#Fill the layers
	wl = (12.398/energy)*1e-10
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
	
	
	
def DoubleSlit(slitsize,slitdistance,pixelsize):
	field1 = np.ones(shape=(400,400),dtype=float)#Absorption
	datafieldA = np.zeros(shape=(400,400),dtype=float)	
	datafieldR = np.zeros(shape=(400,400),dtype=float)
	
	datafieldA += 100.
	datafieldR += 100.
	
	ys,xs = field1.shape
	
	slitdistance /= pixelsize
	
	print "Slits size approximated to:", int((slitsize/pixelsize)*2.)
	#Slit 1
	datafieldA[:,int(xs/2-(slitsize/pixelsize)-slitdistance/2.):int(xs/2+(slitsize/pixelsize)-slitdistance/2.)] = 0.0
	datafieldR[:,int(xs/2-(slitsize/pixelsize)-slitdistance/2.):int(xs/2+(slitsize/pixelsize)-slitdistance/2.)] = 0.0
	#Slit 2
	datafieldA[:,int(xs/2-(slitsize/pixelsize)+slitdistance/2.):int(xs/2+(slitsize/pixelsize)+slitdistance/2.)] = 0.0
	datafieldR[:,int(xs/2-(slitsize/pixelsize)+slitdistance/2.):int(xs/2+(slitsize/pixelsize)+slitdistance/2.)] = 0.0

	transm = np.exp(-datafieldA+1j*datafieldR)
	return field1*transm
######################################################################################


def blade(energy,thickness):#simulate tungsten blade
	#GET MATERIAL PARAMETERS
	thickness *= 10. #Convert in cm
	TEST = db.cdb()
	formula="WC"
	Energy = energy
	Density = 15.63

	d,b = TEST.getVal(formula,Energy,Density)
	print d, b

	field1 = np.ones(shape=(400,400),dtype=float)#Absorption
	datafieldA = np.zeros(shape=(400,400),dtype=float)	
	datafieldR = np.zeros(shape=(400,400),dtype=float)
	
	wl = (12.398/energy)*1e-10
	WN = (2*np.pi) / wl

	datafieldR[:,0:int(200)] = d * thickness
	datafieldA[:,0:int(200)] = b * ((4.*np.pi)/wl) * thickness 
	#AIR 
	TEST2 = db.cdb()
	formula="N4O"
	Energy = energy
	Density = 1.2e-3

	d,b = TEST2.getVal(formula,Energy,Density)

	datafieldR[:,int(200):] = d * thickness
	datafieldA[:,int(200):] = b * ((4.*np.pi)/wl) * thickness 

	print d, b

	ys,xs = field1.shape

	transm = np.exp(-datafieldA+1j*datafieldR)
	return field1*transm

def zpt():
	import sys
	sys.path.insert(0, '/data/id17/map/ID17_Rec/GITrepo/dev_v1.4/other/')
	from OpticElements import createZonePlate
	focus=2
	energy=4
	px=0.2e-6
	size=[0.2e-3,0.2e-3]
	n=16
	field1=createZonePlate(energy,focus,n,px,size)
	return field1

if __name__ == '__main__':

	"""
		Entrance matrix must be divisible by 4 without rest...
	"""

	#for i in range(1,400):
	
	energy = 30. #in keV
	px_size = 3.1e-6 #in m
	#nn = int(i)/10.
	thickness = 0.1e-3
	distance = 2.3 #in m
	source_distance = 140. #in m

	#fields1 = simplewire(energy,px_size)
	#fields1 = blade()#DoubleSlit(3e-6,10e-6,px_size)
	fields1 = simplewire(energy,px_size)
	#fields1 = blade(energy,thickness)#zpt()	
	#DO PADDING
	oxps, oyps = fields1.shape
	fields1 = np.pad(fields1,[oyps,oxps],'edge')	
	#x,y = fields1[0:1000,0:1000].shape
	#fields1 = fields1[0:1000,0:1000]
	#write_edf("ZPTplate.edf",fields1) #path for saving propagated intensity
	#Polychromatic beam
	numberofevents = 2
	MD = 1e-4

	pro = pb.propagators(px_size,energy,distance,1)
	pro.preparefresnel2D(fields1)
	pro.SRsource([20e-6,100e-6],source_distance) #prepare the SR source 20 x 100 um
	write_edf("originalwire.edf",np.real(fields1[oyps:-oyps,oxps:-oxps])) #path for saving propagated intensity
	fieldIntensity = pro.polybeamMC(fields1,numberofevents,MD)
	#fieldIntensity = pro.AddDetectorBlurring(20.,fieldIntensity)
	#fieldIntensity = pro.AddDetectorNoise(500,fieldIntensity)
	name = "PropagatedWIRE30keVNoNoise" + str('%4.4f' % 2.3) + "m.edf"
	fieldIntensity = fieldIntensity[oyps:-oyps,oxps:-oxps]
	write_edf(name,fieldIntensity) #path for saving propagated intensity
	
