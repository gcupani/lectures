import matplotlib.pyplot as plt
import numpy as np

figsize = (12,8)

def extrema(arr, mode='max'):
    """ @brief Find the location of maximum or minimum in an array
        @param arr Array
        @param mode Mode (`max` or `min`)
    """
    
    if mode not in ['max', 'min']:
        print("Mode must be either `max` or `min`. Using `max`.")
        mode = 'max'
    y, x = np.unravel_index(getattr(np, 'arg%s' % mode)(arr), arr.shape)
    return x, y


def plot_graph(x, y=None, new=True, mode='plot', xlabel=None, ylabel=None, **kwargs):
    """ @brief Plot a graph with proper formatting 
        @param img Image
        @param **kwargs Keyword arguments of imshow
    """

    if y is None:
        if len(x)==2:
            x, y = x
        else:
            x, y = range(len(x)), x
    if new: plt.figure(figsize=figsize)
    if mode not in ['plot', 'scatter', 'errorbar']:
        print("Only `plot`, `scatter` or `errorbar` available for graphs. Using `plot`.")
        mode = 'plot'
    getattr(plt, mode)(x, y, **kwargs)
    if xlabel!=None: plt.xlabel(xlabel)
    if ylabel!=None: plt.ylabel(ylabel)


def plot_img(img, new=True, **kwargs):
    """ @brief Plot an image with proper formatting 
        @param img Image
        @param **kwargs Keyword arguments of imshow
    """

    if new: plt.figure(figsize=figsize)
    plt.imshow(img, **kwargs)
    plt.colorbar()
    plt.xlabel('x pixel')
    plt.ylabel('y pixel')
    
    
def plot_img_medrange(img, hrange=200, **kwargs):
    """ @brief Plot a sequence of images with proper formatting 
        @param **kwargs Keyword arguments of imshow
    """

    med = np.median(img)
    vmin, vmax = med-hrange, med+hrange
    plot_img(img, vmin=vmin, vmax=vmax, **kwargs)
    
    
def reduce(img, mbias, mflat):
    """ @brief Reduce a science image
        @param img Science image
        @param mbias Master bias
        @param mflat Master flat
    """
    
    debiased = img-mbias
    ff = debiased/mflat
    return debiased, ff

class Star():
    
    def __init__(self, img, reg, rad=[12,24,36]):
        self.img = np.array(img, dtype=float)
        self.reg = reg
        
        self.sub_img = img[reg]
        shape = self.sub_img.shape
        self.xc, self.yc = extrema(self.sub_img)    
        self.rows, self.cols = np.ogrid[-self.yc:shape[0]-self.yc, -self.xc:shape[1]-self.xc]

    def mask(self, rad):
        if isinstance(rad, list) and len(rad)>1:
            self.targ_mask = self.rows**2+self.cols**2<rad[0]**2    
            self.bkg_mask = np.logical_and(self.rows**2+self.cols**2>rad[-2]**2,self.rows**2+self.cols**2<rad[-1]**2)    
        else:
            self.targ_mask = self.rows**2+self.cols**2<rad**2    
            self.bkg_mask = ~self.targ_mask
        
    def photometry(self, gain=0.6, ron=28.8, ron_npix=17500, rad=[12,24,36]):

        # Mask sub-image
        self.mask(rad)
        self.targ = np.copy(self.sub_img)
        self.bkg = np.copy(self.sub_img)
        self.targ[~self.targ_mask] = np.nan
        self.bkg[~self.bkg_mask] = np.nan
        npix = np.nansum(self.targ_mask)
    
        # Subtract background
        bkg_mean = np.nanmean(self.bkg)
        bkg_std = np.nanstd(self.bkg)
        self.bkg = bkg_mean*gain 
        self.bkg_noise = bkg_std*gain 
        self.targ_bkgsub = self.targ-bkg_mean
        bkg_subfact = np.nansum(self.targ_mask)/np.nansum(self.bkg_mask)
    
        # Compute RON subtraction. factor
        ron_subfact = np.nansum(self.targ_mask)/ron_npix
        
        # Compute flux and error
        self.flux = np.nansum(self.targ_bkgsub)*gain
        self.error = np.sqrt(self.flux + (1+ron_subfact)*ron**2*npix + (1+bkg_subfact)*self.bkg_noise**2*npix)
        self.snr = self.flux/self.error 
        
    def magnitude(self, ref):
        self.mag = ref.mag - 2.5*(np.log10(self.flux)-np.log10(ref.flux))    