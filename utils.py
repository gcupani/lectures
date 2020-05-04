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
    if mode not in ['plot', 'scatter']:
        print("Only `plot` and `scatter` available for graphs. Using `plot`.")
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
    