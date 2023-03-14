__version__ = "0.1"
__author__ = "Alberto Mittone"
__lastrev__ = "Last revision: 12/10/17"

import numpy as np

def ZonePlatesPar(energy,focus,n):
	"""
	Basic formulas for Fresnel zone plates
	energy: 	in keV
	focus: 		fix the focus distance
	n: 		number of zones
	
	Returns:
	rn: 	zone radius -> n=0 is excluded automatically
	Drn: 	zone area
	NA:	Numerical aperture
	Rres: 	Rayleigh resolution
	k:	Number of rings
	
	#Beware: only if NA << 1, the spherical aberration can be neglected
	Spherical aberrations are not considered here
	
	"""
	wl = 1e-10*(12.398/(energy))
	#radius nth zone:
	rn = np.zeros(n)
	Drn = np.zeros(n)
	#Ring counter
	k = 0
	for i in range(1,n,2):
		rn[k] = np.sqrt(i*wl*focus)
		Drn[k] = (wl*focus) / (2.*rn[k])
		k += 1
	NA = wl / (2.*Drn)
	Rres = (0.61*wl) / NA
	return rn[0:k], Drn[0:k], NA[0:k], Rres[0:k], k

def createZonePlate(energy,focus,n,px,size):
	"""
	Create a Zone Plate object
	Input:
	energy 	in keV, 
	focus 	in m, 
	n: 	number of zones
	px: 	pixel size in m
	size: 	size of the zone plates in m
	Returns:
	matrix: defined zone plate - approximated absorption/no absorption
	"""
	sy, sx = np.asarray((np.asarray(size)/px),dtype=int)
	rn, Drn, NA, Rres, k = ZonePlatesPar(energy,focus,n)
	rn = np.asarray(rn/px,dtype=int)
	Drn = np.asarray(Drn/px,dtype=int)
	
	#Get odd dimensions of the matrix
	if sy%2 == 0:
		sy += 1
	if sx%2 == 0:
		sx += 1
	#Create matrix field
	field = np.zeros([sy,sx],dtype=float)
	#Final field
	Ffield = np.zeros([sy,sx],dtype=float)
	
	#Create arrays field
	x = np.linspace(0,sx-1,sx)
	y = np.linspace(0,sy-1,sy)
	nt = np.ones(sx)
	#Center arrays
	x -= (sx-1)/2
	y -= (sy-1)/2
	#Define ring condition
	for i in range(0,k):
		print("rn:",rn[i], "Drn:", Drn[i])
	def Cond(x,y,nt,rn,Drn):
		tmpx = np.sqrt(np.outer(nt,x)**2 + np.outer(y,nt)**2) > rn - Drn/2
		tmpy = np.sqrt(np.outer(nt,x)**2 + np.outer(y,nt)**2) < rn + Drn/2
		return tmpx*tmpy
	#Loop over rings - get only odd numbers
	for i in range(0,k):
		field[Cond(x,y,nt,rn[i],Drn[i])] = 1
		Ffield += field
		#Reset to avoid multiple sums
		field[Cond(x,y,nt,rn[i],Drn[i])] = 0
	return Ffield
