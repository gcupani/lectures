#!/usr/bin/env python
# coding: utf-8

# ### Quick start with FITS
# A simple image frame contains only a primary **header data unit (HDU)**. We will use **Astropy** to open it.
# > Edit the path to point to *your* frame!

# In[4]:


from astropy.io import fits
name = '/Users/guido/Dropbox/lectures/2019/SZ_Lyn_2015/Autosave Image -001.fit'
hdul = fits.open(name)
hdul.info()


# The HDU contains both **metadata** (header) and **data** (image).   
# The header is like a dictionary and its elements (**cards** or **keywords**) can be called by name:

# In[5]:


hdr = hdul[0].header
hdr


# In[6]:


hdr['EXPTIME']


# You can modify keywords and add comments to keywords:

# In[7]:


hdr['OBJECT'], hdr.comments['OBJECT']


# In[8]:


hdr['OBJECT'] = 'V* SZ Lyn'
hdr.comments['OBJECT'] = 'Variable Star of delta Sct type'
hdr['OBJECT'], hdr.comments['OBJECT']


# The image is an **array** and its pixels can be addressed by index:

# In[9]:


img = hdul[0].data
img


# In[10]:


img[0][0]


# You can of course save your changes (in a new frame, possibly):

# In[11]:


hdul.writeto('new_frame.fits', overwrite=True)
hdul_test = fits.open('new_frame.fits')
hdul_test[0].header['OBJECT']


# ### Learn how to speak array
# An image is (quite obviously) a **2-d matrix** of pixels. You can think of it as an **array of (sub-)arrays**. Each sub-array is a row of pixels. The array of all rows is the whole matrix. ```shape``` gives the size of the matrix (rows times columns).

# In[12]:


img.shape


# A sub-array can be called with a single index: 

# In[13]:


row0 = img[0]
row0


# Then you can call an element (or a slice) of the sub-array:

# In[14]:


row0[2]


# In[15]:


row0[::2]


# You can do the same with two indices:

# In[16]:


img[0][2], img[0,2]


# Notice that using two sets of brackets is *not* the same as using one. Each set of brackets is applied *after* the previous one, while operations within brackets are performed *all at once*.

# In[17]:


img[0][:], img[0,:]


# In[18]:


img[:][0], img[:,0]


# ### A plot is worth a thousand printouts
# We will use **Matplotlib** to visualize data.

# In[19]:


import matplotlib.pyplot as plt
plt.imshow(img)
plt.show()


# To extract a sub-image, you can slice the array:

# In[20]:


sub_img = img[510:530, 630:660]
plt.imshow(sub_img)
plt.show()


# You can use ```figure``` to create a new figure instead of overplotting it. ```show``` is used just at the end to display all figures. Compare:

# In[21]:


plt.imshow(img)
plt.imshow(sub_img)
...
plt.show()  # Only last image is shown


# In[22]:


plt.figure()  # First figure
plt.imshow(img)
plt.figure()  # Second figure
plt.imshow(sub_img)
...
plt.show()  # All images are shown in separate figures


# Based on the above example, let's make a simple function to plot a list of images (*bigger* images, with added colorbar and axis labels):

# In[23]:


def plot_list(imgs):
    for i in imgs:
        plt.figure(figsize=(16,8))
        plt.imshow(i)
        plt.colorbar()
        plt.xlabel('x pixel')
        plt.ylabel('y pixel')
    plt.show()
plot_list([img, sub_img])


# Beware that the input to this particular function *must* be a list, even if you want to plot a single image:

# In[24]:


#plot_list(sub_img)  # WRONG!
plot_list([sub_img])


# ```imshow``` works only with 2-d matrices (images). To plot simple arrays (like image rows or columns), use ```plot``` for a continuous plot or ```scatter``` for a scatter plot:

# In[25]:


row = sub_img[7,:]
col = sub_img[:,13]
row_range = range(sub_img.shape[1])
col_range = range(sub_img.shape[0])
plt.figure()
plt.plot(row_range, row)
plt.figure()
plt.scatter(col_range, col)
plt.show()


# Let's improve our plotting function to work with both rows/columns and full images. Note that the part that makes the plot (```plot```) has been splitted from the part that loops through all arrays (```plot_list```):

# In[26]:


def plot_list(arrs):
    for a in arrs:
        plt.figure(figsize=(16,8))
        plot(a)
    plt.show()

def plot(arr):
    if len(arr.shape) == 1:
        plt.plot(range(len(arr)), arr)
        plt.xlabel('pixel')
        plt.ylabel('ADU')
    elif len(arr.shape) == 2:
        plt.imshow(arr)
        plt.colorbar()
        plt.xlabel('x pixel')
        plt.ylabel('y pixel')
    else:
        pass


col = sub_img[:,13]
plot_list([sub_img, row, col])


# > **Your turn now**: write the code to extract a 11x11-pixel window centered in (7, 13) from ```sub_img``` and to plot in separate figures all the rows and columns of this new sub image.

# ### And now for some serious data crunching
# We will use **NumPy** to do statistics on arrays. A simple test case: find the global maximum of an array (remember from analysis: the *maximum point* is the index where the array reaches its *maximum value*).  

# In[27]:


import numpy as np
row_mv = np.max(row) 
row_mp = np.argmax(row)
row_mv, row_mp


# In[28]:


plt.plot(range(len(row)), row)
plt.scatter(row_mp, row_mv, color='r') 
plt.show()


# You can use ```where``` to check the result:

# In[29]:


np.where(row==row_mv)


# In the case of an image, ```argmax``` returns the maximum point of the flattened array (i.e. the array obtained by concatenating all rows):

# In[30]:


sub_img_mv = np.max(sub_img)
sub_img_mp = np.argmax(sub_img)
sub_img_mv, sub_img_mp


# To obtain a 2-tuple with the row and column of the maximum point, you have to ‘unravel’ the result onto the shape of the array: 

# In[31]:


(sub_img_mrow, sub_img_mcol) = np.unravel_index(np.argmax(sub_img), sub_img.shape)
sub_img_mrow, sub_img_mcol


# You obtain something similar with ```where```:

# In[32]:


np.where(sub_img==sub_img_mv)


# Let's again improve our plotting function to display also the maximum point (and maximum value, for single arrays). Note that we can redefine ```plot``` without changing ```plot_list```:

# In[33]:


def argmax(arr):
    if len(arr.shape) == 1:
        return np.argmax(arr)
    else:
        return np.unravel_index(np.argmax(arr), arr.shape)

def plot(arr):
    am = argmax(arr)
    if len(arr.shape) == 1:
        plt.plot(range(len(arr)), arr)
        plt.scatter(am, np.max(arr), color='r')  # In 1-d, we plot location and maximum value
        plt.xlabel('pixel')
        plt.ylabel('ADU')
    elif len(arr.shape) == 2:
        plt.imshow(arr)
        plt.colorbar()
        plt.scatter(am[1], am[0], color='r')  # In 2-d, we only plot location
        plt.xlabel('x pixel')
        plt.ylabel('y pixel')
    else:
        pass
    
plot_list([sub_img, row, col])


# > **Your turn now**: with help from the NumPy Reference (https://docs.scipy.org/doc/numpy-1.16.1/reference/), compute the arithmetic mean and the median of ```sub_img```. What can we infer from the difference between the two values? **Extra points**: compute the median of each rows of ```img``` separately and plot all the values as a function of the row index.
