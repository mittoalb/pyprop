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

from math import sqrt
import sys
#import scipy
from scipy import signal

from pyprop import rw
from pyprop import propagator as pb
from pyprop import models
from math import sqrt


#from PBIpropagator import propagator

if __name__ == '__main__':

	"""
		Entrance matrix must be divisible by 4 without rest...
	"""
	
	energy = 30. #in keV
	px_size = 1e-6 #in m
	#nn = int(i)/10.
	thickness = 0.1e-3
	distance = 0.5 #in m
	source_distance = 37. #in m
	source_size_V, source_size_H = 8.5e-6, 140e-6 #in um - SIGMAS
	counts_on_pixel = 50000
	delta = 6.70778e-07 #PMMA
	beta  = 3.39441e-10 #PMMA
	field_size = 400
	sphere_radius = 20.

	print(delta,beta)
	#SPHERE EXAMPLE - (energy,size,R,pixelsize,delta,beta)
	fields1, _, _ = models.testSphere(energy, field_size, sphere_radius, px_size, delta, beta)

	#DO PADDING
	oxps, oyps = fields1.shape
	fields1 = np.pad(fields1,[oyps,oxps],'edge')	
	#Polychromatic beam
	numberofevents = 200 #number of energies to sum
	MD = 1e-2	     #energy bandwidth

	pro = pb.propagators(px_size,energy,distance,1)
	pro.preparefresnel2D(fields1)
	pro.SRsource([source_size_V, source_size_H],source_distance) #prepare the SR source 20 x 100 um

	rw.write_tif("OriginalField.tif",np.real(fields1[oyps:-oyps,oxps:-oxps]))
	

	fieldIntensity = pro.polybeamMC(fields1,numberofevents,MD)
	fieldIntensity = pro.AddDetectorBlurring(1.4,fieldIntensity)
	fieldIntensity = pro.AddDetectorNoise(counts_on_pixel,fieldIntensity)
	name = "Out" + str('%4.4f' % 2.3) + "m.tif"
	fieldIntensity = fieldIntensity[oyps:-oyps,oxps:-oxps]
	
	rw.write_tif(name,fieldIntensity)
