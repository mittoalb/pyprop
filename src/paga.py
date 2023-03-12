import PhaseRetrieval
from tomopy_phase import retrieve_phase
from PhaseRetrieval import phase_retrieval
from utility import read_edf, write_edf, next_greater_power_of_2
import numpy
import scipy
from scipy import signal

def fft_frequencies(Nx):
	fx = numpy.zeros(Nx,"f")
	fx[0:(Nx)/2+1] = numpy.arange((Nx)/2+1)
	fx[(Nx)/2+1:Nx] = -(Nx)/2+1 + numpy.arange(Nx-(Nx)/2-1)
	return fx /Nx	


def PaganinPhaseRetrieval(data,px,PL):
	Ny, Nx = data.shape
	f1=fft_frequencies(Nx)/px
	f2=fft_frequencies(Ny)/px
	FF = (f2*f2)+(f1*f1)[:,None]
	kernel = 1.0/((1.0+FF*PL*PL)*Nx*Ny)
	padsize = int(next_greater_power_of_2(Nx))
	Sx, Sy = kernel.shape
	padding = padsize - Nx
	#kernel = kernel[(Sy-KERNEL_size)/2:(Sy+KERNEL_size)/2,(Sx-KERNEL_size)/2:(Sx+KERNEL_size)/2]
	data = scipy.signal.fftconvolve(numpy.pad(data,padding,'edge'),kernel,mode='same')
	NNx, NNy = data.shape
	data = data#[(NNy-Ny)/2:(NNy+Ny)/2,(NNx-Nx)/2:(NNx+Nx)/2]
	return data

data, h  = read_edf('PropagatedWIRE30keVNoNoise2.3000m.edf')
px = 3.1e-4
dist = 230.
energy = 30.
alpha = 1.0e-6
padding = False
PL = 800.
#out = PaganinPhaseRetrieval(data,px,PL)
out = retrieve_phase(data, px, dist, energy, alpha, padding)
write_edf('PRSum30keVWireNoNoise.edf', out)
