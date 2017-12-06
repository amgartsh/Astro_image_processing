# -*- coding: utf-8 -*-
"""
Created on Fri May 26 22:46:21 2017

@author: Adam """


# Program Details
"""
 This program is a first attempt at creating a FellWalker-style algorithm for finding the stellar sources of a .fits image. It works on a previously edited stationary
 image that has had its background removed (assumed planar), so that clumps can be easily found and analyzed by the iterative walk, and deemed significant or noise-based.
 
 The initial pass of all of the pixels will find points outside of the 1-sigma deviance of the background, and label these as potentially significant. All others will be
 deemed as noise-related. Next, the iterative walks will parse the potentially significant points into separate sources. 
 
 After all potentially significant points have been sorted, adjacent clumps that are not separable by variation will be amalgamated until we're left with only fully
 distinguishable sources. Then, the boundaries of each clump will be smoothed to remove some unneccessary noise-based differentials. In the last portion of the program,
 all unuseable clumps will be removed from the pool of good, useable sources. """

 
 # Libraries
 
from astroio.py import fits
import numpy as np
import math
 
 # Functions
 
def reject_outliers(data, m = 3.):
    md = np.median(data)
    for i in range(0,len(data[:,0])):
        for j in range(0,len(data[0,:])):
            d = np.abs(data[i,j] - md)
            s = d/md if md else 0.
            if s > m:
                data[i,j] = None
    return data  

 # Global Variables

clump_ass = []
sig_pix = []
peak_array = []
peak_pairs = []

 # Main
 
 # Import image
 
fitsname = input('Enter the file name of the .fits image: ')
fitsname += '.fits'

image = fits.open(fitsname, memmap=True)
image_array = image[1] 
length = len(image_array[:,0])
width = len(image_array[0,:])                                                         # CAA
clump_ass = np.shape(length,width)
                   
# Specification of Initial Pass Filter

bkrd = reject_outliers(image_array)
mu = np.mean(bkrd)
dev = np.std(bkrd)

# Test specs, input by tester

print('The mean value of the background is '), mu, ('and the standard deviation is '), dev, '\n'
ini_pass = input('Enter the initial pass filter value as it relates to background standard deviation: ')

# Specification of Minimum Initial Gradient Path and Pixel Averaging Length

min_grad = input('Enter the minumum gradient defining the edges of each source: ')
pix_mean = input('Enter the number of pixels over which to calculate the gradient: ')

# Specification of Minimum Significant Peak Value

min_peak = input('Enter the minimum intensity above background (in std. dev. of the background) defining a significant peak: ')

# Specification of Maximum Distance between Peaks to Merge

maxjump = input('Enter the maximum pixel distance over which to amalgamate peaks into a single source: ')

# Specification of Minimum Dip to Separate Sources

mindip = input('Enter the minimum decrease in intensity (std devs of larger peak) defining separate sources: ')

# Specification of mean wavelength of filter, aperture diameter (circular)

wavelength = input('Enter the mean wavelength of the bandwidth filter: ')
diameter = input('Enter the detector aperture diameter: ')

# Initial Pass of filter
# Possibly in clump => clump_ass = 0. Failure to pass filter => clump_ass = -1

pass_filt = mu + dev*ini_pass

for i in range(0,length):
    for j in range(0,width):
        if image_array[i,j] < pass_filt:
            clump_ass[i,j] = -1
        elif image_array[i,j] >= pass_filt:
            clump_ass[i,j] = 0
            sig_pix = np.append(sig_pix,[i,j])
                     
# Clump identification

while len(sig_pix) > 0:
    temp_path = []
    temp_path = np.append(temp_path, sig_pix[0,:])
    i = sig_pix[0,0]
    j = sig_pix[0,1]
    min_peak = min_peak * dev + mu
    c_a = 1
    peak = 0
    x = i
    y = j
    
    while peak == 0:
        maxval = image_array[x,y]
        nn = 1
        
        while (x == i) and (y == j):
            for k in range(-nn, nn):
                for l in range(-nn, nn):
                    if image_array[i + k, j + l] > maxval:
                        maxval = image_array[i + k, j + l]
                        x = i + k
                        y = j + l
            if nn < 4 and (x == i) and (y == j):
                nn += 1
            else:
                break
            
        if (x == i) and (y == j):
            peak_array = np.append(peak_array,[i,j,image_array[i,j]])
            peak += 1
            
            if clump_ass[i,j] == 0:
                if image_array[i,j] < min_peak:
                    for k in range(0,len(temp_path[:,0])):
                        clump_ass[temp_path[k,:]] = -1
                if image_array[i,j] >= min_peak:
                    clump_ass[i,j] = c_a
                    c_a += 1
                
            for k in range(0, len(temp_path[:,0])):
                x_iden = temp_path[k,0]
                y_iden = temp_path[k,1]
                clump_ass[x_iden, y_iden] = clump_ass[i,j]                       
                for l in range(0, len(sig_pix[:,0])):
                    if (x_iden == sig_pix[l,0]) and (y_iden == sig_pix[l,1]):
                        sig_pix = np.delete(sig_pix, sig_pix[l,:])
                        
            while len(temp_path[:,0]) > pix_mean:
                sum_grad = 0
                for k in range(0,pix_mean):
                    sum_grad += image_array[temp_path[k+1,:]] - image_array[temp_path[k,:]]
                mean_grad = sum_grad/pix_mean
                if mean_grad < min_grad:
                    clump_ass[temp_path[0,:]] = -1
                    temp_path = np.delete(temp_path, temp_path[0,:])
                if mean_grad >= min_grad:
                    break
                for k in range(0, length):
                    clump_ass[temp_path[k,:]] = -1
                             
            if len(temp_path[:,0]) < pix_mean:
                for k in range(0, len(temp_path[:,0])):
                    grad = image_array[temp_path[k+1,:]] - image_array[temp_path[k,:]]
                    if grad < min_grad:
                        clump_ass[temp_path[k,:]] = -2

        temp_path = np.append(temp_path,[x,y])
        i = x
        j = y
        
# Clump/Source Merging

for k in range(0, len(peak_array[:,0]) - 1):
    for l in range(1, len(peak_array[:,0])):
        while k < l:
            distance = math.sqrt((peak_array[k,0]-peak_array[l,0])**2 + (peak_array[k,1]-peak_array[l,1])**2)
            if distance < maxjump:
                pair = np.empty([0,3])
                pair[0,0] = peak_array[k,0]
                pair[0,1] = peak_array[k,1]
                pair[0,2] = peak_array[l,0]
                pair[0,3] = peak_array[l,1]
                peak_pairs = np.append(peak_pairs,pair)
                
# Find Resolution distance between Fellwalker peaks based on the Rayleigh Criterion of the detector.
# Input: Filter bandwidth, circular aperture diameter, image size, image-solid angle ratio, peak pair array.
# Operations: Calculate Resolution angle for each peak pair, compare to detector's Rayleigh
# Output: Array of every resolved peak pair.

# input:
# variable: filter bandwidth = wavelength [wavelength]
# variable: circular aperture diameter = diameter [m]

rayleigh = 1.22*wavelength/diameter
res_pp = np.empty()
iden_changes = np.empty()

for i in range(0,len(peak_pairs)):

    dist_px = sqrt((peak_pairs[i,1]-peak_pairs[i,0])**2+(peak_pairs[i,3]-peak_pairs[i,2])**2)
    ang_dist = dist_px*(image_angle/image_size)
    
    if ang_dist >= rayleigh:
        res_pp = np.append(res_pp,[peak_pairs[i,0],peak_pairs[i,1]])
    if ang_dist < rayleigh:
    
        # change higher peak number to lower peak number
        if image_array[peak_pairs[i,0],peak_pairs[i,1]] >= image_array[peak_pairs[i,2],peak_pairs[i,3]]:
            iden_change = clump_ass[peak_pairs[i,2],peak_pairs[i,3]]
            new_iden = clump_ass[peak_pairs[i,0],peak_pairs[i,1]]
            
        else:
            iden_change = clump_ass[peak_pairs[i,0],peak_pairs[i,1]]
            new_iden = clump_ass[peak_pairs[i,2],peak_pairs[i,3]]
            
        iden_changes = np.append(iden_changes,[iden_change,new_iden])
        
for m in range(0,len(iden_changes)):
    for n in range(0,len(iden_changes)):
        x=iden_changes[m,0]
        y=iden_changes[n,1]
        if x==y:
            iden_changes[m,0] = iden_changes[n,0]
            iden_changes = np.delete(iden_changes,n,0     
            
for j in range(0,length):
    for k in range(0,width):
        for i in range(0,len(iden_changes[:,0])):
            if (clump_ass[j,k] == iden_changes[i,0]):
                clump_ass[j,k] = iden_changes[i,1]

## Clump Smoothing

pix_changes = 1
while pix_changes > 0:
    pix_changes = 0
    for j in range(0,length):
        for k in range(0,width):
            n_peaks = np.amax(clump_ass)
            hist_peak = np.empty()
            m =2
            adj_points = np.empty()
            for i in range(0,n_peaks):
                hist_peak=np.append(hist_peak,[i,0])
            for i in range(-m,m):
                for n in range(-m,m)
                dist = sqrt((j+i)**2+(k+n)**2)
                points= [i,n,clump_ass[j+i,k+n]]
                adj_points = np.append(adj_points,points)
            for c in range(0,len(adj_points)):
                if adj_points[c,0]=0 and adj_points[c,1]=0:
                    hist_peak[adj_points[c,2],1] += 2*sqrt(2)
                else:
                    weight = 1/(sqrt(adj_points[c,0]**2+adj_points[c,1]**2))
                    hist_peak[adj_points[c,2],1] += weight
            base = 0
            new_val = 0
            for d in range(0,len(hist_peak)):
                if hist_peak[d,1]>base:
                    new_val = hist_peak[d,0]
            if new_val != clump_ass[j,k]:
                pix_changes += 1
                clump_ass[j,k] = new_val

















