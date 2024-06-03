"""
author: A. Mittone
modified by: M. Eckermann
PBI simulator

ADJUST:
- some settings in lines 289-291
- path where to save in line 327
- some settings in lines 248-264
"""
import numpy as np
import scipy.ndimage

from pyprop.rw import *

import copy
import random
import fabio
from pyprop import rw
from math import sqrt

class propagators:
	"""
		PBI propagation - Main class
	"""
	def __init__(self,pixelsize,energy,distance,over):
		self.px = pixelsize
		self.E = energy
		self.z = distance
		self.wl = (12.398/self.E)*1e-10
		self.oversampling = over
		self.peakenergy = energy

	def binning(self,data):
		"""
			Perform binning to restore image size after oversampling
		"""
		if self.oversampling > 1:
			sy = int(len(data[:,0])/self.oversampling)
			sx = int(len(data[0,:])/self.oversampling)
			binned = np.zeros(shape=(sy,sx),dtype=np.float)
			metabin = np.zeros(shape=(len(data[:,0]),sx),dtype=np.float)
			for i in range(0,sx):
				metabin[:,i] = np.sum(data[:,self.oversampling*i:self.oversampling*(i+1)],axis=1)
			for j in range(0,sy):
				binned[j,:] = np.sum(metabin[self.oversampling*j:self.oversampling*(j+1),:],axis=0)
			data = binned/(self.oversampling*self.oversampling)
		else:
			"Do nothing"
		return data

	def readdata(self,datareal,dataimg):
		"""
			Read absorption and phase information
		"""
		return datareal + 1j*dataimg

	def SRsource(self,size,distance):
		"""
			Define the Synchrotron radiation source
			size: size in m [Y,X]
			distance source-object: in m
			#sampling: number of points per pixel which with the source is modelized
		"""
		#Coherence length
		#wl = (12.398/self.E)*1e-10
		self.DDratio = self.z / distance
		self.SDdistance = distance
		
		#Define Gaussian source - adjust for the distance
		self.source = np.ones(shape=[self.yps,self.xps],dtype=float)
		size[0], size[1] = size[0]*self.DDratio, size[1]*self.DDratio 
		self.source = self.gauss2D(self.source,size)

		self.kernel_coh = self.source
		self.kernel_coh /= np.sum(self.kernel_coh) #Normalize the kernel
		rw.write_tif("WGL-source.tif",self.kernel_coh)
		

	def preparefresnel2D(self,dataimage):
		"""
			Prepare 2D fresnel diffraction fields via convolution by Fourier transform
		"""
		#Padding of the image
		oxps, oyps = dataimage.shape
		#dataimage = np.pad(dataimage,[oyps,oxps],'edge')
		#oxps, oyps = dataimage.shape
		"""
		interpolation ->
			1 nearest interpolation
			2 bilinear interpolation
			3 cubic interpolation
		"""
		if self.oversampling != 1: ####BUGGEDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD do not use
			#Separating component and oversampling
			dataimagereal = np.real(dataimage)
			dataimageimag = np.imag(dataimage)
			dataimagereal = scipy.ndimage.zoom(dataimagereal, self.oversampling, order=2)
			dataimageimag = scipy.ndimage.zoom(dataimageimag, self.oversampling, order=2)
		
			#recombined oversampled image
			dataimage = dataimagereal + 1j*dataimageimag
		
		xps, yps = dataimage.shape
		#xps = len(dataimage[0,:])
		#yps = len(dataimage[:,0])
		self.px /= self.oversampling
		Dxsize = self.px * xps
		Dysize = self.px * yps
		
		P_x = np.linspace(-Dxsize/2,Dxsize/2,xps)
		P_y = np.linspace(-Dysize/2,Dysize/2,yps)
		#print P_x
		#xoffset = P_x[0]
		#yoffset = P_y[0]

		#fft_sizex = xps
		#fft_sizey = yps
    
		fft_deltax = 1.0/xps/(P_x[1] - P_x[0])
		fft_deltay = 1.0/yps/(P_y[1] - P_y[0])

		if np.mod(xps,2) == 1:
			fft_offsetx = -fft_deltax*(float(xps-1)/2.0)#fft_offsetx = -fft_deltax*float(xps-1)/2.0
		else:
			fft_offsetx = -fft_deltax*float(xps)/2.0

		if np.mod(yps,2) == 1:
			fft_offsety = -fft_deltay*(float(yps-1)/2.0)#fft_offsety = -fft_deltay*float(yps-1)/2.0 
		else:
			fft_offsety = -fft_deltay*float(yps)/2.0

		#print fft_deltax, fft_offsetx
		F1 = np.fft.fft2(dataimage)
		#import scipy
		#from scipy import fftpack
		#wfou_fft = scipy.fftpack.fftshift(F1)
		#wfou_fft = np.fft.fftshift(F1)
		
		wfou_fft_x = np.fft.ifftshift(np.arange(start=fft_offsetx, stop = -fft_offsetx, step=fft_deltax, ))
		wfou_fft_y = np.fft.ifftshift(np.arange(start=fft_offsety, stop = -fft_offsety, step=fft_deltay, ))
		
		[U, V] = np.meshgrid(wfou_fft_y,wfou_fft_x)
		self.frq2 = U**2 + V**2
		
		self.wfou_fft = F1#wfou_fft
		self.xps = xps
		self.yps = yps

	
	def fpro(self):
		"""
			Only propagate and back-fourier transform for speed reason
		"""
		wfou_fft2 = copy.copy(self.wfou_fft)
		#Calculate effective propagation distance
		self.z1 = (self.z*self.SDdistance)/(self.z + self.SDdistance)
		wfou_fft2 *= np.exp(-1j * np.pi * self.wl * self.z1 * self.frq2)
	
		#back FT
		#wfou_fft2 = np.fft.ifftshift(wfou_fft2)
		propagated = np.fft.ifft2(wfou_fft2)

		fieldIntensity = np.abs(propagated)**2
		fieldPhase = np.arctan2(np.real(propagated), \
				       np.imag(propagated))
		
		#back to original size
		fieldIntensity = self.binning(fieldIntensity)
		fieldPhase = self.binning(fieldPhase)

		return fieldIntensity, fieldPhase
	

	def gauss(self,sigma,nob):
		"""
			Needs sigma and number of bins
		"""
		FWHM = 2.355 * sigma
		step = FWHM / nob
		self.x = np.linspace(-FWHM/2.,FWHM/2.,nob)
		val = np.zeros(nob)
		for i in range(0,nob):
			rr = -FWHM/2 + i*step
			val[i] = np.exp(-rr**2/(2*sigma**2))
		return val
		
	def gauss2D(self,field,sigma,center=None):
		'''
			Take an empty square input field, it fills it with a gaussian centered in the middle
		'''
		xf, yf = field.shape
		#Check the size
		if xf%2 == 0 and yf%2 == 0:
			field = np.zeros(shape=[yf+1,xf+1],dtype=float)
			xf, yf = field.shape
		elif xf%2 != 0 and yf%2 == 0:
			field = np.zeros(shape=[yf+1,xf],dtype=float)
			xf, yf = field.shape
		elif xf%2 == 0 and yf%2 != 0:
			field = np.zeros(shape=[yf,xf+1],dtype=float)
			xf, yf = field.shape			
			
		x = np.arange(0, xf, 1, float)
		y = x[:,np.newaxis]

		if center is None:
			x0 = y0 = xf // 2
		else:
			x0 = center[0]
			y0 = center[1]
			
		#Get asymmetric sigmas
		sigma_y, sigma_x = sigma
		sigma_y /= self.px
		sigma_x /= self.px
		return (1/(2*np.pi*sigma_x*sigma_y) * np.exp(-((x-x0)**2/(2*sigma_x**2) + (y-y0)**2/(2*sigma_y**2))))				


	def gaussMC(self,center,FHWM):
		"""
			Gives weight following a Gaussian distribution
		"""
		sigma = FHWM / 2.355
		x = random.gauss(center,sigma)
		y = np.exp(-(x-center)**2/(2*sigma**2))
		return x, y

	def beamdistribution(self,disttype,sigma,nob):
		"""
			Not used yes, to be implemented....maybe
		"""
		if disttype == "gauss":
			self.bdistr = self.gauss(sigma,nob)
		else:
			print(disttype, "Not supported.")

	def polybeam(self,field1,numberofwaves,sigma,path,radix,savefileevery):
		"""
			Calculation of polychromatic beam with gauss shape
		"""
		step = 2.355*sigma/numberofwaves
		weight = self.gauss(sigma,numberofwaves)
		weight /= sum(weight)
		#oxps, oyps = fields1.shape
		#fields1 = np.pad(fields1,[oyps,oxps],'edge')	
		
		self.preparefresnel2D(fields1)
		#Initialize 0 field of intensity
		self.fieldIntensity = np.zeros(shape=[len(fields1[0,:]),len(fields1[:,0])])
		self.E += self.x[0]
		for n in range(0,numberofwaves):
			self.E += step
			tmp, fieldPhase1 = self.fpro()
			#Update of the field
			self.fieldIntensity += tmp*weight[n]
			if np.mod(n,savefileevery) == 0:
				print("# Waves: ", n + 1, "Energy value: ", self.E, "Weight: ", weight[n])
				outname = path + radix + str('%4.4d' % n) + ".edf"
				rw.write_tif(outname,self.fieldIntensity)



	def polybeamMC(self,fields1,numberofevents,MD):
		"""
			Temporal coherence effect
			MonteCarlo sum of waves, it assumes a gaussian distribution of the beam
			field1: is the original geometry of the simulation
			numberofevents: MC number of events
			MD: Energy polychmaticity DeltaE/E
		"""
		#Check energy conservation
		#oxps, oyps = fields1.shape
		#fields1 = np.pad(fields1,[oyps,oxps],'edge')	
		#norm = np.sum(fields1)
		self.preparefresnel2D(fields1)
		#Initialize 0 field of intensity
		#Ny, Nx = fields1.shape
		self.fieldIntensity = np.zeros_like(fields1,dtype=float)
		self.WLFWHM = self.peakenergy * MD
		#Total Weigth
		self.TWe = 0.0
		
		for n in range(0,int(numberofevents)):
			E, weight = self.gaussMC(self.peakenergy,self.WLFWHM)
			self.wl = (12.398/E)*1e-10	
			tmp, fieldPhase1 = self.fpro()
			#Update of the field
			self.fieldIntensity += tmp*weight
			self.TWe += weight
			print("Summing up:", E, weight)

		#Normalize
		self.fieldIntensity /= self.TWe
		#Consider the effect of spatial coherence
		self.fieldIntensity = scipy.signal.fftconvolve(self.fieldIntensity,self.kernel_coh,mode='same')
		return self.fieldIntensity
		
	def AddDetectorNoise(self,photons,fieldIntensity):
		"""
			Generate Poisson noise and rescale the results on the number of counts on the detector
		"""	
		mu, sigma = photons, sqrt(photons)
		x, y = fieldIntensity.shape
		noise = np.random.normal(mu, sigma, x*y)
		noise = np.reshape(noise,(x,y))
		fieldIntensity *= noise
		return fieldIntensity
		
	def AddDetectorBlurring(self,sigma,fieldIntensity):
		"""
			Simulate the blurring related to the detection system
		"""	
		fieldIntensity = scipy.ndimage.filters.gaussian_filter(fieldIntensity, sigma)
		return fieldIntensity

