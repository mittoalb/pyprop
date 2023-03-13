# -*- coding: utf-8 -*-
"""
Author: A. Mittone
Last revision: 21/10/16
Module for fft transform
"""
import numpy as np
import os
import pyfftw
import multiprocessing as mp

nthread = int(mp.cpu_count()-4)

def fftw(a):
	out = np.zeros_like(a)
	fft = pyfftw.builders.fft(a, s=None, axes=(-2, -1), overwrite_input=False, planner_effort='FFTW_MEASURE', threads=nthread, auto_align_input=False, auto_contiguous=False, avoid_copy=True)
	out = fft()
	return out

# --------------------------------------------------------------------


def ifftw(a):
	out = np.zeros_like(a)
	fft = pyfftw.builders.ifft(a, s=None, axes=(-2, -1), overwrite_input=False, planner_effort='FFTW_MEASURE', threads=nthread, auto_align_input=False, auto_contiguous=False, avoid_copy=True)
	out = fft()
	return out

# --------------------------------------------------------------------


def fftw2(a):
	out = np.zeros_like(a)
	fft = pyfftw.builders.fft2(a, s=None, axes=(-2, -1), overwrite_input=False, planner_effort='FFTW_MEASURE', threads=nthread, auto_align_input=False, auto_contiguous=False, avoid_copy=True)
	out = fft()
	return out

# --------------------------------------------------------------------

def ifftw2(a):
	out = np.zeros_like(a)
	fft = pyfftw.builders.ifft2(a, s=None, axes=(-2, -1), overwrite_input=False, planner_effort='FFTW_MEASURE', threads=nthread, auto_align_input=False, auto_contiguous=False, avoid_copy=True)
	out = fft()
	return out
