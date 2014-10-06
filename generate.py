#!/usr/bin/env python
# -*- coding: cp1252 -*-
"""
This example demonstrates how to create the 17 segment model for the left
ventricle recommended by the American Heart Association (AHA).
"""
##set working directory. must keep import csv file in the same 


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv
import itertools
from itertools import chain
import math


def bullseye_plot(ax, data, vlim=None, segBold=[]):
    """
    Bullseye representation for the left ventricle.

    Parameters
    ----------
    ax : axes
    data : list of int and float
        The intensity values for each of the 17 segments
    vlim : [min, max] or None
        Optional argument to set the Intensity range
    segBold: list of int
        A list with the segments to highlight


    Notes
    -----
    This function create the 17 segment model for the left ventricle according
    to the American Heart Association (AHA) [1]_

    References
    ----------
    .. [1] M. D. Cerqueira, N. J. Weissman, V. Dilsizian, A. K. Jacobs,
        S. Kaul, W. K. Laskey, D. J. Pennell, J. A. Rumberger, T. Ryan,
        and M. S. Verani, “Standardized myocardial segmentation and nomenclature
        for tomographic imaging of the heart,” Circulation, vol. 105, no. 4,
        pp. 539–542, 2002.
    """
	
    linewidth = 2
    data = np.array(data).ravel()

    if vlim is None:
        vlim = [data.min(), data.max()]

    theta = np.linspace(0, 2*np.pi, 768)
    r = np.linspace(0.2, 1, 4)

    # Create the bound for the segment 17
    for i in range(r.shape[0]):
        ax.plot(theta, np.repeat(r[i], theta.shape), '-k', lw=linewidth)

    # Create the bounds for the segments  1-12
    for i in range(6):
        theta_i = i*60*np.pi/180
        ax.plot([theta_i, theta_i], [r[1], 1], '-k', lw=linewidth)

    # Create the bounds for the segments 13-16
    for i in range(4):
        theta_i = i*90*np.pi/180 - 45*np.pi/180
        ax.plot([theta_i, theta_i], [r[0], r[1]], '-k', lw=linewidth)

    # Fill the segments 1-6
    r0 = r[2:4]
    r0 = np.repeat(r0[:, np.newaxis], 128, axis=1).T
    for i in range(6):
        # First segment start at 60 degrees
        theta0 = theta[i*128:i*128+128] + 60*np.pi/180
        theta0 = np.repeat(theta0[:, np.newaxis], 2, axis=1)
        z = np.ones((128, 2))*data[i]
        ax.pcolormesh(theta0, r0, z, vmin=vlim[0], vmax=vlim[1], cmap=plt.cm.Blues_r)
        if i+1 in segBold:
            ax.plot(theta0, r0, '-k', lw=linewidth+2)
            ax.plot(theta0[0], [r[2], r[3]], '-k', lw=linewidth+1)
            ax.plot(theta0[-1], [r[2], r[3]], '-k', lw=linewidth+1)

    # Fill the segments 7-12
    r0 = r[1:3]
    r0 = np.repeat(r0[:, np.newaxis], 128, axis=1).T
    for i in range(6):
        # First segment start at 60 degrees
        theta0 = theta[i*128:i*128+128] + 60*np.pi/180
        theta0 = np.repeat(theta0[:, np.newaxis], 2, axis=1)
        z = np.ones((128, 2))*data[i+6]
        ax.pcolormesh(theta0, r0, z, vmin=vlim[0], vmax=vlim[1], cmap=plt.cm.Blues_r)
        if i+7 in segBold:
            ax.plot(theta0, r0, '-k', lw=linewidth+2)
            ax.plot(theta0[0], [r[1], r[2]], '-k', lw=linewidth+1)
            ax.plot(theta0[-1], [r[1], r[2]], '-k', lw=linewidth+1)

    # Fill the segments 13-16
    r0 = r[0:2]
    r0 = np.repeat(r0[:, np.newaxis], 192, axis=1).T
    for i in range(4):
        # First segment start at 45 degrees
        theta0 = theta[i*192:i*192+192] + 45*np.pi/180
        theta0 = np.repeat(theta0[:, np.newaxis], 2, axis=1)
        z = np.ones((192, 2))*data[i+12]
        ax.pcolormesh(theta0, r0, z, vmin=vlim[0], vmax=vlim[1], cmap=plt.cm.Blues_r)
        if i+13 in segBold:
            ax.plot(theta0, r0, '-k', lw=linewidth+2)
            ax.plot(theta0[0], [r[0], r[1]], '-k', lw=linewidth+1)
            ax.plot(theta0[-1], [r[0], r[1]], '-k', lw=linewidth+1)

    #Fill the segments 17
    if data.size == 17:
        r0 = np.array([0, r[0]])
        r0 = np.repeat(r0[:, np.newaxis], theta.size, axis=1).T
        theta0 = np.repeat(theta[:, np.newaxis], 2, axis=1)
        z = np.ones((theta.size, 2))*data[16]
        ax.pcolormesh(theta0, r0, z, vmin=vlim[0], vmax=vlim[1], cmap=plt.cm.Blues_r)
        if 17 in segBold:
            ax.plot(theta0, r0, '-k', lw=linewidth+2)

    ax.set_ylim([0, 1])
    ax.set_yticklabels([])
    ax.set_xticklabels([])

def arr_av (arr_of_arrs) :
    """
    Assumes you are kind enough to pass in arrays of the same length. Returns
    an array of the averages of each value.
    """
    vals = []
    if(len(arr_of_arrs) == 0 ) : return vals
    for x in range(len(arr_of_arrs[0])) :
        vals.append(average(arr_of_arrs, x))
    return vals

def average(arr_of_arrs, pos) :
    """
    Assumes you are kind enough to pass in arrays of the same length. Returns
    the average of the values at position <pos> for each array in arr_of_arrs
    """
    if(len(arr_of_arrs) == 0 or pos > len(arr_of_arrs[0])) : return -1
    counter = 0
    for x in range(len(arr_of_arrs)) :        
        counter += arr_of_arrs[x][pos]
    return counter / len(arr_of_arrs);

arr_of_groups = []
with open('ss_rawdata.csv') as ssr :
    rdr = csv.reader(ssr)
    counter = 0
    group = -1
    group_arr = []
    for row in rdr :        
        # Don't include the first row
        if(counter > 0) :
            sub_arr = []
            rowcount = 0            
            for char in row :
                # Ignore first column. Second column defines groups - create new 2D array if new group
                if(rowcount == 1) :                    
                    if(group != -1 and int(char) != group) :
                        # Switch to next group
                        arr_of_groups.append(list(group_arr))
                        group_arr = []                                            
                    group = int(char)                    
                if(rowcount > 1) :
                    if(char.strip() != '') :
                        # Convert to int
                        sub_arr.append(int(char))
                    else :
                        # Translate blank as 0
                        sub_arr.append(0)
                rowcount += 1                
            group_arr.append(sub_arr)            
        counter += 1
    # Append the last group
    arr_of_groups.append(group_arr)



vlim = [1,12]#[data.min(), data.max()]

for x in range(len(arr_of_groups)) :
    arr = arr_of_groups[x]

    num_cols = int(math.sqrt(len(arr) + 1))
    
    # We need to create an x by x grid of graphs. I'm using n columns,
    # where n is calculated to be as square a grid as possible. To make
    # sure I have enough columns to fit everything in
    # (ie num rows of data / n, plus one more row if num rows of data / n has a remainder
    overrun = 1 if (len(arr) + 1) % num_cols > 0 else 0
    
    fig, ax = plt.subplots(figsize=(12, 8), nrows=int((len(arr)+1)/num_cols) + overrun, ncols=num_cols,
                               subplot_kw=dict(projection='polar'))

    #Convert 2D array of axes into a list
    chain = list(itertools.chain.from_iterable(ax))


    fig.canvas.set_window_title('Left Ventricle Bulls Eyes (AHA) Group ' + str(x + 1))

    for y in range (len(arr)):

        data = np.array(arr[y])#np.array(range(17)) + 1
        
        bullseye_plot(chain[y], data, vlim=vlim)
        chain[y].set_title('Patient ' + str(y+1))
        

    # Add average
    data = arr_av(arr)

    bullseye_plot(chain[len(arr)], data, vlim=vlim)
    chain[len(arr)].set_title("Average")

    #Add legend
    cm = plt.cm.Blues #plt.cm.jet

    #define the bins and normalize
    cNorm = mpl.colors.Normalize(vmin=vlim[0], vmax=vlim[1])

    ticks = [vlim[0], 0, vlim[1]]
    chain[len(arr)] = fig.add_axes([0.1, 0.01, 0.1, 0.025])
        
    cb = mpl.colorbar.ColorbarBase(chain[len(arr)], cmap=cm, norm=cNorm, ticks=ticks,
                                       orientation='horizontal')
plt.show()

