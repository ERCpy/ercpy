'''

Microscopy tools module
Subset of different small functions for electron microscopy

'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import utils
from matplotlib.patches import Rectangle
from scipy.fftpack import ifft2,fftshift,fft2,ifftshift
from scipy.signal import correlate2d, medfilt
from skimage.feature import register_translation
from scipy.misc import imresize


__all__ = ['fft']


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
    
    
def align_img(img_one, img_two, method = 'imgreg', show=False, roi=True, sb_filtering=False, filt_size= 200,
              binning=2, feducial=1, manualxy = None, **kwargs):
    '''
    Function to align images or holograms using X-correlation
    Parameters
    ----------
    img_one : ndarray
        The refrence image
    img_two : ndarray
        An image to align
    method : string
        Either 'imgreg' to use image registartion from skimage, 
        or 'xcorr' to use crosscorrelation in real space from scipy,
        or 'feducial' to use feducial markers,
        or 'manual' to displace the images manually
#        or 'gui_shift' to shift image intractively
    show : boolean or string
        Set True to plot the results, set 'diff' to show difference of the images
    roi : boolean
        Set true to do alignament on ROI instead of whole image
    sb_filtering : boolean
        Set True, to apply for holograms
    filt_size : int
        Size of the filter for main band filtering
        used only for alignment of holograms
    binning : int
        Binning for the images during alignment used in 'xcorr' method
    feducial : int
        Number of feducial markers for 'feducial' method
    manualxy : tuple of 2 int
        Coordiantes x,y for manual alignment
        
    Returns
    -------
    img_algn : ndarray
        Aligned image img_two
    (xdrift, ydrift): tuple
        drift correction coordinates

    Notes
    -----
    * Manual alignamnt method requiers cooridnates provided using 'manualxy' parameter
    * X-correlation method is slow! Use only small images and small ROI (use of binning is recomended)
    * Binning is implemented for "xcorr" method only
    
    See Also
    --------

    '''
    
    (ry,cx) = img_one.shape
#    img_one = img_one.astype(float)
#    img_two = img_two.astype(float)
    
    # --- IFFT main band:
    if sb_filtering:
        fft_img_one = fftshift(fft2(img_one))
        (xx,yy) = np.meshgrid(np.linspace(-ry/2, ry/2-1, ry), np.linspace(-cx/2, cx/2-1, cx))
        rr = np.sqrt(xx**2+yy**2)
        mask = np.zeros((ry,cx))
        mask[rr<filt_size] = 1
        img_one_m = np.absolute(ifft2(ifftshift(fft_img_one*mask)))
        # --- Processing second image
        fft_img_two = fftshift(fft2(img_two))
        img_two_m = np.absolute(ifft2(ifftshift(fft_img_two*mask)))
    else:
        img_one_m = img_one
        img_two_m = img_two
      
    # --- Use ROI if True
    if roi:
        # --- GUI based assignment of ROI
        f, ax = plt.subplots(1, 1)
        ax.imshow(img_one_m, cmap=cm.binary_r)
        rect = utils.RoiRect()
        if hasattr(f.canvas.manager, 'window'): f.canvas.manager.window.raise_()
        plt.waitforbuttonpress(100)
        plt.waitforbuttonpress(5)
        plt.close(f)
    else:
        rect = Rectangle((0,0), 1, 1,fc='none', ec='r')
        rect.x0 = 0
        rect.y0 = 0
        rect.y1 = ry-1
        rect.x1 = cx-1

    # --- Select allignment method
    if method is 'imgreg':
        
        img_one_roi = img_one_m[rect.y0:rect.y1, rect.x0:rect.x1]
        img_two_roi = img_two_m[rect.y0:rect.y1, rect.x0:rect.x1]
        
#        px_rescale_y = np.float(img_one_m.shape[0])/np.float(img_one_roi.shape[0])
#        px_rescale_x = np.float(img_one_m.shape[1])/np.float(img_one_roi.shape[1])
        upsample = 4
        # --- Upsampled image registration for ROI
        shift, error, diffphase = register_translation(img_one_roi, img_two_roi, upsample)
        
        ydrift = shift[0]
        xdrift = shift[1]
        print(shift)
        
#        # --- Accounting for change in pixel size   <- is not needed since ROI doesn't change px-size!!     
#        ydrift = ydrift*px_rescale_y
#        xdrift = xdrift*px_rescale_x    

    elif method is 'xcorr': # slow X-cor for small images only!
        # TODO: Check if the metod is working properly
        # --- Selecting ROI and X-correlating
        img_one_m = imresize(img_one_m, 1.0/binning)
        img_two_m = imresize(img_two_m, 1.0/binning)
        template = img_two_m[rect.y0/binning:rect.y1/binning, rect.x0/binning:rect.x1/binning]
        cc = correlate2d(img_one_m, template, boundary='symm', mode='same')
        imax = np.argmax(np.absolute(cc))
        ypeak, xpeak = np.unravel_index(imax, cc.shape) # coordinates of X-corr max
#        ydrift = (rect.y0/binning-ypeak-template.shape[0]-1)*binning
#        xdrift = (rect.x0/binning-xpeak-template.shape[1]-1)*binning
        ydrift = (ypeak-(rect.y1-rect.y0)/2)*binning
        xdrift = (xpeak-(rect.x1-rect.x0)/2)*binning
        
    elif method is 'feducial': # alignment using feducial markers
        # TODO: add multiple marker alignments
        img_one_roi = img_one_m[rect.y0:rect.y1, rect.x0:rect.x1]
        img_two_roi = img_two_m[rect.y0:rect.y1, rect.x0:rect.x1]
        
        f, ax = plt.subplots(1, 1)
        ax.imshow(img_one_roi, cmap=cm.binary_r)
        ax.set_title('Please set feducial marker position for 1st image')
        marker_one = utils.RoiPoint()
        if hasattr(f.canvas.manager, 'window'): f.canvas.manager.window.raise_()
        plt.waitforbuttonpress(100)
        plt.waitforbuttonpress(1)
        plt.close(f)
        
        f, ax = plt.subplots(1, 1)
        ax.imshow(img_two_roi, cmap=cm.binary_r)
        ax.set_title('Please set feducial marker position for 2nd image')
        marker_two = utils.RoiPoint()
        if hasattr(f.canvas.manager, 'window'): f.canvas.manager.window.raise_()
        plt.waitforbuttonpress(100)
        plt.waitforbuttonpress(1)
        plt.close(f)
        
        xdrift = marker_one.x0 - marker_two.x0
        ydrift = marker_one.y0 - marker_two.y0
        
    elif method is 'manual':
        if manualxy:
            xdrift = manualxy[0]
            ydrift = manualxy[1]
        else:
            raise ValueError('Method manual requires shifts provided in manualxy argument.')
            
#    elif method is 'gui_shift': #TODO: create 'gui_shift' method
        
    else:
        raise ValueError('Wrong method argument! Check doc.')

    print("xydrift = %d, %d" % (xdrift, ydrift) )
    
    img_algn = np.roll(img_two, np.int(ydrift), axis=0)
    img_algn = np.roll(img_algn, np.int(xdrift), axis=1)
    
    if isinstance(show, basestring):
        if show is 'diff':
            f, ax = plt.subplots(1,1)
            ax.imshow(img_algn - img_one, cmap=cm.binary_r)
            ax.set_title(('xydrift = ', str(xdrift), str(ydrift)))
    elif show:
        f, ax = plt.subplots(1,1)
        ax.imshow(img_algn, cmap=cm.binary_r)
    return (img_algn, (xdrift, ydrift))

    
def rm_duds(img, sigma=8.0, median_k=5):
    '''
    Removes dud pixels from images
    
    Parameters
    ----------
    img : ndarray
        The image
    sigma : float
    
    median_k : int
        Size of median kernel
        
    Returns
    -------
    img_nodud : ndarray
        Image with removed dud pixels (e.g. X-Rays spikes)

    Notes
    -----
    
    See Also
    --------

    '''
    # TODO: check if 1D for spectrum works as well
    img_mf = medfilt(img, median_k) # median filtered image
    diff_img = np.absolute(img-img_mf)
    mean_diff = sigma*np.sqrt(np.var(diff_img))
    duds = diff_img > mean_diff
    img[duds] = img_mf[duds]
    
    n_duds = np.sum(duds) # dud pixels
    print "The number of pixels changed = %d" % n_duds
    
    return (img, duds)