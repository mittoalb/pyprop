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
from PyMca import EdfFile
import copy
import random
#from utility import write_edf
#from utility import read_edf

import fabio
#import cupy as cp

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

	def bin_data(self, data):
		if self.oversampling > 1:
			sy = len(data[:,0]) // self.oversampling
			sx = len(data[0,:]) // self.oversampling
			binned = np.zeros(shape=(sy, sx), dtype=np.float)
			metabin = np.zeros(shape=(len(data[:,0]), sx), dtype=np.float)
			for i in range(sx):
			    metabin[:,i] = np.sum(data[:,self.oversampling*i:self.oversampling*(i+1)], axis=1)
			for j in range(sy):
			    binned[j,:] = np.sum(metabin[self.oversampling*j:self.oversampling*(j+1),:], axis=0)
			data = binned / (self.oversampling*self.oversampling)
		return data

	def read_data(self, data_real, data_imag):
		return data_real + 1j * data_imag

	def prepare_fresnel_2D(self, data_image):
		if self.oversampling != 1:
			# Separating component and oversampling
			data_image_real = np.real(data_image)
			data_image_imag = np.imag(data_image)
			data_image_real = scipy.ndimage.zoom(data_image_real, self.oversampling, order=2)
			data_image_imag = scipy.ndimage.zoom(data_image_imag, self.oversampling, order=2)
			# recombined oversampled image
			data_image = data_image_real + 1j * data_image_imag

		xps, yps = data_image.shape
		Dxsize = self.px * xps
		Dysize = self.px * yps

		P_x = np.linspace(-Dxsize/2, Dxsize/2, xps)
		P_y = np.linspace(-Dysize/2, Dysize/2, yps)
		xoffset = P_x[0]
		yoffset = P_y[0]

		fft_size_x = xps
		fft_size_y = yps

		fft_delta_x = 1.0 / (xps * (P_x[1] - P_x[0]))
		fft_delta_y = 1.0 / (yps * (P_y[1] - P_y[0]))

		if np.mod(xps, 2) == 1:
			fft_offset_x = -fft_delta_x * (float(xps - 1) / 2.0)
		else:
			fft_offset_x = -fft_delta_x * float(xps) / 2.0

		if np.mod(yps, 2) == 1:
			fft_offset_y = -fft_delta_y * (float(yps - 1) / 2.0)
		else:
			fft_offset_y = -fft_delta_y * float(yps) / 2.0

		F1 = np.fft.fft2(data_image)
		wfou_fft_x = np.fft.ifftshift(np.arange(start=fft_offset_x, stop=-fft_offset_x, step=fft_delta_x))
		wfou_fft_y = np.fft.ifftshift(np.arange(start=fft_offset_y, stop=-fft_offset_y, step=fft_delta_y))
		U, V = np.meshgrid(wfou_fft_y, wfou_fft_x)
		self.frq2 = U ** 2 + V ** 2
		self.wfou_fft = F1

	"""
		Only propagate and back-fourier transform for speed reason
	"""	
	def fpro(self):
		wfou_fft2 = copy.copy(self.wfou_fft)
		wfou_fft2 *= np.exp(-1j * np.pi * self.wl * self.z * self.frq2)
	
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
	
	"""
		Needs sigma and number of bins
	"""
	def gauss(self,sigma,nob):
		FWHM = 2.355 * sigma
		step = FWHM / nob
		self.x = np.linspace(-FWHM/2.,FWHM/2.,nob)
		val = np.zeros(nob)
		for i in range(0,nob):
			rr = -FWHM/2 + i*step
			val[i] = np.exp(-rr**2/(2*sigma**2))
		return val
		
	"""
		Gives weight
	"""
	def gaussMC(self,center,FHWM):
		sigma = FHWM / 2.355
		x = random.gauss(center,sigma)
		y = np.exp(-(x-center)**2/(2*sigma**2))
		return x, y
	"""
		Not used yes, to be implemented....maybe
	"""
	def beamdistribution(self,disttype,sigma,nob):
		if disttype == "gauss":
			self.bdistr = self.gauss(sigma,nob)
		else:
			print disttype, "Not supported."
	"""
		Calculation of polychromatic beam with gauss shape
	"""
	def polybeam(self,field1,numberofwaves,sigma,path,radix,savefileevery):
		step = 2.355*sigma/numberofwaves
		weight = self.gauss(sigma,numberofwaves)
		weight /= sum(weight)
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
				print "# Waves: ", n + 1, "Energy value: ", self.E, "Weight: ", weight[n]
				outname = path + radix + str('%4.4d' % n) + ".edf"
				saved=EdfFile.EdfFile(outname)
				saved.WriteImage({},self.fieldIntensity, Append=0, DataType='FloatValue')
	"""
		MonteCarlo sum of waves
	"""
	def polybeamMC(self,field1,numberofevents,MD,path,radix,savefileevery):
		#Check energy conservation
		norm = np.sum(fields1)
		self.preparefresnel2D(fields1)
		#Initialize 0 field of intensity
		self.fieldIntensity = np.zeros_like(field1,dtype=float)
		WLFWHM = self.peakenergy * MD
		print self.peakenergy, WLFWHM
		if numberofevents != 0:
			for n in range(0,int(numberofevents)):
				E, weight = self.gaussMC(self.peakenergy,WLFWHM)
				self.wl = (12.398/E)*1e-10	
				tmp, fieldPhase1 = self.fpro()
				#Update of the field
				self.fieldIntensity += tmp*weight
				if np.mod(n,savefileevery) == 0:
					print "# Events: ", n+1, "Energy value: ", E, "Wavelenght:", self.wl, "Weight: ", weight
					outname = path + radix + str('%4.4d' % n) + ".edf"					
					prop = np.sum(self.fieldIntensity)
					#self.fieldIntensity *= (norm / prop)
					#saved=EdfFile.EdfFile(outname)
					#saved.WriteImage({},self.fieldIntensity, Append=0, DataType='FloatValue')
					#write_edf(outname,np.real(self.fieldIntensity))
		else:
			n = 0
			while 1:
				E, weight = self.gaussMC(self.peakenergy,WLFWHM)
				self.wl = (12.398/E)*1e-10	
				tmp, fieldPhase1 = self.fpro()
				#Update of the field
				self.fieldIntensity += tmp*weight
				if np.mod(n,savefileevery) == 0:
					print "# Events: ", n+1, "Energy value: ", E, "Wavelenght:", self.wl, "Weight: ", weight
					outname = path + radix + str('%4.4d' % n) + ".edf"
					#saved=EdfFile.EdfFile(outname)
					#saved.WriteImage({},self.fieldIntensity, Append=0, DataType='FloatValue')
					#write_edf(outname,np.real(self.fieldIntensity))
				n += 1
		#Normalize
		self.fieldIntensity /= numberofevents	
		return self.fieldIntensity

