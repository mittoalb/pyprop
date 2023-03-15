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

#sys.path.insert(0, '/data/id17/map/ID17_Rec/GITrepo/dev_v1.4/other/')
#import db
#from PBIpropagator import propagator




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
	fields1 = models.simplewire()#energy,px_size)
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

	rw.write_tif("originalwire.tif",np.real(fields1[oyps:-oyps,oxps:-oxps]))
	

	fieldIntensity = pro.polybeamMC(fields1,numberofevents,MD)
	#fieldIntensity = pro.AddDetectorBlurring(20.,fieldIntensity)
	#fieldIntensity = pro.AddDetectorNoise(500,fieldIntensity)
	name = "PropagatedWIRE30keVNoNoise" + str('%4.4f' % 2.3) + "m.tif"
	fieldIntensity = fieldIntensity[oyps:-oyps,oxps:-oxps]
	
	rw.write_tif(name,fieldIntensity)
