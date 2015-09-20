# -*- coding: utf-8 -*-

'''

Microscopy tools module
Subset of different small functions for electron microscopy

'''

import numpy as np
from scipy.fftpack import fftshift,fft2
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def fft(img):
    '''
    
    Calculate and plot FFT
    
    Parameters
    ----------
    img : array_like
        image data

    Returns
    -------
    fft_data : ndarray
        Complex array of calculated FFT

    '''
    
    fft_data = fftshift(fft2(img))    
    plt.imshow(np.log(np.absolute(fft_data)), cmap=cm.binary_r)