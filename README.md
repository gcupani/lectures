# Lectures at the University of Trieste

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gcupani/lectures/master)

## 2019 - Programming @ Astro Lab: exercises with Python
Each exercise can be run both as a Jupyter notebook (preferred) and as a Python script. You can launch the Jupyter notebooks also by cliking on "launch binder" above. 
The code is written in Python 3.5. Packages required: Astropy, IPython, Matplotlib, NumPy, SciPy. 
Be careful to edit the paths to fetch your own data. 

- **SNR_regimes**: Simple computation of how the signal-to-noise ratio varies with the target magnitude in different regimes (shot-noise limited, background-noise limited, and readout-noise limited observation).
- **image_handling**: Simple introduction to FITS image handling. Requires a test image (e.g. one frame of SZ Lyn).
- **SZ_Lyn_2015** (previously **image_reduction_analysis**): Extraction of the SZ Lyn light curve from observations taken at INAF-OATs in 2015. Requires a set of SZ Lyn images, a set of bias images and a set of flat images. Fully detailed and commented. **Beware**: this script loads all images into the RAM and may cause problems if images are too many.
- **SZ_Lyn_2019**: same as above, but for the observations taken in 2019. More sparsely commented. Uses functions in **functions.py**. This script load images into the RAM one by one and is much less memory intensive. 
