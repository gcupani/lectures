from astropy.io import fits
from IPython import display
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import scipy.fftpack
import scipy.signal as sgn

def argmax(arr):
    """Find the maximum point of a 1-d or 2-d array"""
    if len(arr.shape) == 1:
        return np.argmax(arr)
    else:
        return np.unravel_index(np.argmax(arr), arr.shape)

def flux_cal(sci_names, mbias, mflat, ref_xy, ref_width, p_range, signal, noise, mag_ref):
    signal_ref, noise_ref = phot_circ(sci_names, mbias, mflat, ref_xy, ref_width, p_range)
    return mag_ext(signal, noise, signal_ref, mag_ref)

def fft_filt(data, cut=1e-7):
    x = np.linspace(0,2*np.pi,len(data))
    w = scipy.fftpack.rfft(data)
    spectrum = w**2
    f = scipy.fftpack.rfftfreq(1000, x[1]-x[0])
    cutoff_idx = spectrum < (spectrum.max()*cut)
    w2 = w.copy()
    w2[cutoff_idx] = 0
    return scipy.fftpack.irfft(w2)

def get_data(name):
    """Get data from a FITS frame"""
    hdul = fits.open(name)
    return hdul[0].data

def get_header(name):
    """Get header from a FITS frame"""
    hdul = fits.open(name)
    return hdul[0].header

def mag_ext(signal, noise, signal_ref, mag_ref):
    mag = mag_ref - 2.5*(np.log10(signal)-np.log10(signal_ref))
    err_up = mag - mag_ref + 2.5*(np.log10(signal+noise)-np.log10(signal_ref))
    err_down = mag_ref - 2.5*(np.log10(signal-noise)-np.log10(signal_ref)) - mag
    return mag, err_up, err_down

def mask_circ(sci_red, am, rad=12):
    """Create a circular mask"""
    shape = sci_red.shape
    y, x = np.ogrid[-am[0]:shape[0]-am[0], -am[1]:shape[1]-am[1]]
    return x**2 + y**2 < rad**2

def phot_circ(sci_names, mbias, mflat, xy, width, p_range=[], gain=0.6,
              ron=20.0, filter_window=3, **kwargs):
    """Perform circular photometry"""
    signal = []
    noise = []
    targ_a = []
    bkg_a = []

    for i, n in enumerate(sci_names):
        plt.figure(figsize=(16,8))

        # Reduction
        sci = get_data(n)
        sci_red = (sci - mbias) / mflat

        # Subimage definition
        sub = np.s_[xy[0]:xy[0]+width, xy[1]:xy[1]+width]
        sci_sub = sci_red[sub]

        sci_filt = sgn.medfilt(sci_sub, [filter_window,filter_window])

        # Circular masking
        targ = np.copy(sci_sub)
        bkg = np.copy(sci_sub)
        mask_t = (targ, bkg)
        am = argmax(sci_filt)
        mask = mask_circ(sci_sub, am, **kwargs)
        targ[~mask] = np.nan
        bkg[mask] = np.nan

        if i in p_range:
            #targ_a.append(targ)
            #bkg_a.append(bkg)

            plot_medrange(sci_sub, new=False, argmax=am)
            plt.text(x=0, y=0, s=i, color='black', ha='left', va='top', fontsize='xx-large')
            display.clear_output(wait=True)
            display.display(plt.gcf())
        plt.close()

        # Background removal
        bkg_med = np.nanmedian(bkg)
        targ_sub = targ-bkg_med

        # Signal and noise computation
        s = np.nansum(targ_sub)*gain
        ron2 = ron**2*np.sum(mask)
        bkgn = np.nanstd(bkg)*gain
        totn = np.sqrt(s+2*bkgn**2+2*ron2)
        signal = np.append(signal, s)
        noise = np.append(noise, totn)

    return signal, noise


def plot(arr, title=None, argmax=None, new=True, vmin=None, vmax=None, scatter=False, selection=None, **kwargs):
    """Plot a 1-d or 2-d array"""
    if len(arr.shape) == 1:
        plot_1d(arr, title, argmax, new, vmin, vmax, scatter=scatter, **kwargs)
    elif len(arr.shape) == 2:
        plot_2d(arr, title, argmax, new, vmin, vmax, selection=selection, **kwargs)
    else:
        pass

def plot_1d(arr, title=None, argmax=None, new=True, vmin=None, vmax=None, xlabel='pixel', ylabel='ADU', scatter=False, **kwargs):
    """Plot a 1-d array, possibly with its maximum highlighted"""
    if new:
        plt.figure(figsize=(16,8))
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    if scatter:
        plt.scatter(range(len(arr)), arr, **kwargs)
    else:
        plt.plot(range(len(arr)), arr, **kwargs)
    if argmax is not None: plt.scatter(argmax, np.max(arr), color='r')
    plt.ylim(vmin, vmax)


def plot_2d(arr, title=None, argmax=None, new=True, vmin=None, vmax=None, xlabel='x pixel', ylabel='y pixel', selection=None,
            **kwargs):
    """Plot a 2-d array, possibly with its maximum highlighted"""
    if new:
        plt.figure(figsize=(16,8))
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    plt.imshow(arr, vmin=vmin, vmax=vmax, **kwargs)
    plt.colorbar()
    if selection is not None:
        rect = patches.Rectangle(selection['xy'], selection['width'], selection['height'], edgecolor='r', facecolor='none')
        plt.gca().add_patch(rect)
    if argmax is not None: plt.scatter(argmax[1], argmax[0], color='r')

def plot_medrange(arr, title=None, hwid=200, argmax=None, new=True, **kwargs):
    """Same as `plot`, but with a colorbar centered around the median value"""
    plot(arr, title, argmax, new, np.nanmedian(arr)-hwid, np.nanmedian(arr)+hwid, **kwargs)

def plot_sequence(seq, plot_func, argmax=None, **kwargs):
    """Plot a sequence of arrays in the same window"""
    for i, arr in enumerate(seq):
        plt.figure(figsize=(16,8))
        if argmax is None:
            plot_func(arr, new=False, **kwargs)
        else:
            plot_func(arr, new=False, argmax=argmax[i], **kwargs)
        plt.text(x=0, y=0, s=i, color='black', ha='left', va='top', fontsize='xx-large')
        display.clear_output(wait=True)
        display.display(plt.gcf())
        plt.close()

def plot_time_series(t, f, err=None, title=None, new=True, ylog=False, xlabel='time', ylabel='Flux (ph)', fmt='o', **kwargs):
    """Plot a time series"""
    if new:
        plt.figure(figsize=(16,8))
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
    if ylog:
        plt.gca().set_yscale('log')
    plt.errorbar(t, f, yerr=err, fmt=fmt, **kwargs)
