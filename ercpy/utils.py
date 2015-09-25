import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from skimage import draw
from scipy.fftpack import ifft2,fftshift,fft2,ifftshift
from scipy.signal import correlate2d, medfilt
import matplotlib.cm as cm
from skimage.feature import register_translation
from skimage.feature.register_translation import _upsampled_dft

def fib(n):    # write Fibonacci series up to n
    a, b = 0, 1
    while b < n:
        print b,
        a, b = b, a+b

def fib2(n): # return Fibonacci series up to n
    result = []
    a, b = 0, 1
    while b < n:
        result.append(b)
        a, b = b, a+b
    return result

class RoiRect(object):
    ''' Class for getting a mouse drawn rectangle
    Based on the example from:
    http://matplotlib.org/users/event_handling.html#draggable-rectangle-exercise
    Note that:
    
    * It makes only one roi
    
    '''
    def __init__(self):
        self.ax = plt.gca()
        self.rect = Rectangle((0,0), 1, 1,fc='none', ec='r')
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.ax.add_patch(self.rect)
        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        print 'press'
        self.x0 = event.xdata
        self.y0 = event.ydata
        self.rect.set_linestyle('dashed')
        self.set = False

    def on_release(self, event):
        print 'release'
        self.x1 = event.xdata
        self.y1 = event.ydata
        self.rect.set_width(self.x1 - self.x0)
        self.rect.set_height(self.y1 - self.y0)
        self.rect.set_xy((self.x0, self.y0))
        self.rect.set_linestyle('solid')
        self.ax.figure.canvas.draw()
        self.set = True
        self.ax.figure.canvas.mpl_disconnect(self.on_press)
        self.ax.figure.canvas.mpl_disconnect(self.on_release)
        self.ax.figure.canvas.mpl_disconnect(self.on_motion)
        
    def on_motion(self, event):
        # on motion will move the rect if the mouse
        if self.x0 is None: return
        if self.set: return
        # if event.inaxes != self.rect.axes: return
        self.x1 = event.xdata
        self.y1 = event.ydata
        self.rect.set_width(self.x1 - self.x0)
        self.rect.set_height(self.y1 - self.y0)
        self.rect.set_xy((self.x0, self.y0))
        self.ax.figure.canvas.draw()

#class roi_rect_new(object):
#    ''' Class for getting a mouse drawn rectangle
#    '''
#    def __init__(self):
#        self.ax = plt.gca()
#        self.rect = Rectangle((0,0), 1, 1, facecolor='None', edgecolor='green')
#        self.x0 = None
#        self.y0 = None
#        self.x1 = None
#        self.y1 = None
#        self.ax.add_patch(self.rect)
#        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
#        self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
#        self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
#    def on_press(self, event):
#        print 'press'
#        self.x0 = event.xdata
#        self.y0 = event.ydata    
#        self.x1 = event.xdata
#        self.y1 = event.ydata
#        self.rect.set_width(self.x1 - self.x0)
#        self.rect.set_height(self.y1 - self.y0)
#        self.rect.set_xy((self.x0, self.y0))
#        self.rect.set_linestyle('dashed')
#        self.ax.figure.canvas.draw()
#    def on_motion(self,event):
#        if self.on_press is True:
#            return
#        self.x1 = event.xdata
#        self.y1 = event.ydata
#        self.rect.set_width(self.x1 - self.x0)
#        self.rect.set_height(self.y1 - self.y0)
#        self.rect.set_xy((self.x0, self.y0))
#        self.rect.set_linestyle('dashed')
#        self.ax.figure.canvas.draw()
#    def on_release(self, event):
#        print 'release'
#        self.x1 = event.xdata
#        self.y1 = event.ydata
#        self.rect.set_width(self.x1 - self.x0)
#        self.rect.set_height(self.y1 - self.y0)
#        self.rect.set_xy((self.x0, self.y0))
#        self.rect.set_linestyle('solid')
#        self.ax.figure.canvas.draw()
#        print self.x0,self.x1,self.y0,self.y1
#        return [self.x0,self.x1,self.y0,self.y1]
        
def poly_to_mask(vertex_row_coords, vertex_col_coords, shape):
    '''
    Creates a poligon mask
    '''
    fill_row_coords, fill_col_coords = draw.polygon(vertex_row_coords, vertex_col_coords, shape)
    mask = np.zeros(shape, dtype=np.bool)
    mask[fill_row_coords, fill_col_coords] = True
    return mask

def wrap_to_pi(angle):
    """
    Wrap a given angle in radians to the range -pi to pi.
    
    @param angle : The angle to be wrapped
    @param type angle : float
    @return : Wrapped angle
    @rtype : float
    """
    return np.mod(angle+np.pi,2.0*np.pi)-np.pi
    
def wrap(angle):
    '''
    Wrap a given angle in radians to the range 0 to 2pi.
    
    @param angle : The angle to be wrapped
    @param type angle : float
    @return : Wrapped angle
    @rtype : float
    '''
    return angle % (2 * np.pi )
    
def align_img(img_one, img_two, method = 'imgreg', roi=True, sb_filtering=False, filt_size= 200, feducial=1, manualxy = None, **kwargs):
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
        or 'manual' to displace the images manually.
    roi : boolean
        Set true to do alignament on ROI instead of whole image
    sb_filtering : boolean
        Set True, to apply for holograms
    filt_size : int
        Size of the filter for main band filtering
        used only for alignment of holograms
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
    Manual alignamnt method requiers cooridnates provided using 'manualxy' parameter
    
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
        rect = RoiRect()
        f.canvas.manager.window.raise_()
        plt.waitforbuttonpress(5)
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
        
        px_rescale_y = np.float(img_one_m.shape[0])/np.float(img_one_roi.shape[0])
        px_rescale_x = np.float(img_one_m.shape[1])/np.float(img_one_roi.shape[1])
        upsample = np.min((px_rescale_x, px_rescale_y))
        # --- Upsampled image registration for ROI
        shift, error, diffphase = register_translation(img_one_roi, img_two_roi, upsample)
        
        ydrift = shift[0]
        xdrift = shift[1]
        print(shift)
        
        # --- Accounting for change in pixel size        
        ydrift = ydrift*px_rescale_y
        xdrift = xdrift*px_rescale_x    

    elif method is 'xcorr': # slow X-cor for small images only!
    
        # --- Selecting ROI and X-correlating
        template = img_one_m[rect.y0:rect.y1, rect.x0:rect.x1]
        cc = correlate2d(template,img_two_m)
        imax = np.argmax(np.absolute(cc))
        ypeak, xpeak = np.unravel_index(imax, cc.shape) # coordinates of X-corr max
        ydrift = rect.y0-ypeak-template.shape[0]-1
        xdrift = rect.x0-xpeak-template.shape[1]-1
    else:
        raise ValueError('Wrong method argument! Check doc.')

    print("xydrift = %d, %d" % (xdrift, ydrift) )
    
    img_algn = np.roll(img_two, np.int(ydrift), axis=0)
    img_algn = np.roll(img_algn, np.int(xdrift), axis=1)
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