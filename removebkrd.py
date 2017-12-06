# -*- coding: utf-8 -*-
"""
Created on Fri May 26 12:47:42 2017

@author: Adam
"""
# Program #
"""
This program takes an input of a .fits image and outputs a .fits image without the background radiaton 
This is done by creating a histogram of pixel values along each line (horizontally and vertically), creating
a histogram of the values and then fitting a gaussian to the lowest value peak. 

This yields an estimate for the background radiation along that line and the associated deviation. These 
values are created for the entire image, and then a linear regression is fitted to the average values, with 
the line deviations used to deduce error.

The horizontal and vertical regressions are combined to form a planar background, and this is removed from the 
image.
"""
# Libraries #

from astropy.io import fits
import numpy as np
from scipy import stats

# Global Variables #

bkrd_row = []
bkrd_col = []
dev_row = []
dev_col= []
bkrd = []

# Functions #

def reject_outliers(data, m = 3.):
    d = np.abs(data[1] - np.median(data[1]))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    for i in range(0,len(data[0])):
        if s[i] > m:
            data = np.delete(data,i,1)
    return data

# Main #

fitsname = input('Enter the file name of the .fits image: ')
fitsname += '.fits'

image = fits.open(fitsname, memmap=True)
image_array = image[1]

length = len(image_array[:,0])
width = len(image_array[0,:])

bkrd_row = np.shape(length,1)
bkrd_col = np.shape(width,1)
dev_row = np.shape(length,1)
dev_col= np.shape(width,1)

# Background in rows

for i in range(0,length):
    
    image_temp = np.empty([1,width])
    for k in range(0,width):
        image_temp[0,k] = k+1
    image_temp[1] = image_array[i,:]
    filt_image_temp = reject_outliers(image_temp)
    md = np.median(filt_image_temp)
    sigma = np.std(filt_image_temp)
    mu = np.mean(filt_image_temp)
    
    if mu < md:
        bkrd_row[i,1] = mu  
    else:
        bkrd_row[i,1] = md
    bkrd_row[i,0] = i
            
    dev_row[i,1] = sigma
    dev_row[i,0] = i
           
# Background in Columns           
           
for j in range(0,width):
    
    image_temp = np.empty([1,length])
    for k in range(0,length):
        image_temp[0,k] = k+1
    image_temp[1] = image_array[:,j]
    filt_image_temp = reject_outliers(image_temp)
    md = np.median(filt_image_temp)
    sigma = np.std(filt_image_temp)
    mu = np.mean(filt_image_temp)
    
    if mu < md:
        bkrd_col[j,1] = mu  
    else:
        bkrd_col[j,1] = md
    bkrd_row[j,0] = j
            
    dev_row[j,1] = sigma
    dev_row[j,0] = j
           
# Row Regression

row_Bslope, row_base, row_R2, row_pval, row_sigma = stats.linregress(bkrd_row[:,0],bkrd_row[:,1])

# Columnn Regression

col_Bslope, col_base, col_R2, col_pval, col_sigma = stats.linregress(bkrd_col[:,0],bkrd_col[:,1])

base = (row_base + col_base)/2
       
# Background array values

bkrd = np.shape(length,width)

for i in range(0,length):
    for j in range(0,width):
        bkrd[i,j] = base +row_Bslope*(1/2+i) + col_Bslope*(1/2+j)
            
# Removal of background from image

fixed_image = image_array - bkrd

minval = np.abs(min(fixed_image))

fixed_image += minval

image[1] = fixed_image
     
# Save image as background removed (Brem) .fits file

fitsname -= '.fits'
fitsname += 'BkrdRem.fits'

image.writeto(fitsname)