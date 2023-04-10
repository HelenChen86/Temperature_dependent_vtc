import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import pandas as pd
import shutil
import scipy.optimize as op

from pylab import *
from scipy.interpolate import interp1d





#reads the alldata file (as is from rattrap) and returns a temp list with the data in file 
def read_alldata_file(filename):
    temp_list = []
    try:
        file = open(filename, "r")
        start_reading = False
        for line in file:
            line = line.split('\t') # tab separated columns
            if start_reading:
                temp_list.append(np.array(line))
            if "Energy_(eV)" in line:
                # keeping track of data versus headers and comments
                Title = line
                Title[-1] = line[-1].format().replace('\n','')
                start_reading = True
    finally:
        file.close()
    return temp_list

#writing fitted data to file
def write_to_file(energy, fit, residual, name, delimiter='\t'): 
    name = name.strip('_norm_shifted.txt')
    # for csv file, change delimiter to ','
    if delimiter is ',':
        ext = '.csv'
    else:
        ext = '.txt'
    
    #check if there is an average counts directory and create one if not
    shifted_data_directory_name = 'Fitted Data'
    if not os.path.exists(f'{shifted_data_directory_name}/'):
        os.makedirs(f'{shifted_data_directory_name}/')
    
    # move to average cnts directory
    os.chdir(f'{shifted_data_directory_name}/')
    
    # Write the file
    with open(f'{name}_fit{ext}', 'w') as file:
        # Headers first
        file.write(f'Energy_(eV){delimiter}fit{delimiter}residual\n')
        # Then the data
        for i in range(len(energy)):
            file.write(f'{energy[i]}{delimiter}{fit[i]}{delimiter}{residual[i]}\n')
    # return to parent directory
    os.chdir('../')

    
#interpolate the end spectra to the energy grid in data
def interpolate(energy_oct, octa0, energy_tet, tetra0, energy0, data0):
    f_octa = interp1d(energy_oct, octa0, kind='cubic')
    f_tetra = interp1d(energy_tet, tetra0, kind='cubic')

    energy1 = energy0[30:-30]
    data1 = data0[30:-30]

    octa1 = f_octa(energy1)
    tetra1 = f_tetra(energy1)
    
    #convert the lists to arrays  
    energy1 = np.array(energy1, dtype=float)
    octa1 = np.array(octa1, dtype=float)
    tetra1 = np.array(tetra1, dtype=float)
    data1 = np.array(data1, dtype=float)

    return(energy1, octa1, tetra1, data1)


#specify the start and end of fitting range in energy
def get_index(energy1, start, end):
    start_fit0 = np.where(energy1 >= start)   #start = 9640eV
    end_fit0 = np.where(energy1 >= end)  #end = 9660eV

    start_fit1 = start_fit0[0]
    end_fit1 = end_fit0[0]

    start_fit2 = start_fit1[0]
    end_fit2 = end_fit1[0]


    start_fit = start_fit2.item()
    end_fit = end_fit2.item()
    
    return(start_fit, end_fit)


def model(coeff, octa, tetra):
    fit = (coeff * octa) + ((1 - coeff) * tetra)
    return fit



#define loss 
# must have paramter as first input!
def loss(coeff, octa, tetra, data):
    return np.sum(sqrt((data - model(coeff, octa, tetra))**2))

font = matplotlib.font_manager.FontProperties(family='Times New Roman',
                                               weight='normal',
                                               style='normal', size=24)

font_ticklabel = matplotlib.font_manager.FontProperties(family='Times New Roman',
                                               weight='normal',
                                               style='normal', size=28)

def plot(energy1, data1, octa1, tetra1, energy, fit, residual, name, namefile):
    namefile = namefile.rstrip('norm_shifted.txt')
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    
    plt.plot(energy1, data1, markersize=0.1, label = namefile)
    plt.plot(energy1, octa1, markersize=0.1, label = 'Octahedral')
    plt.plot(energy1, tetra1, markersize=0.1, label = 'Tetrahedral')
    plt.plot(energy, fit, markersize=0.1, color = 'r', label='Fit')
    plt.plot(energy, residual, markersize=0.1, label = 'Residual')

    # Move left y-axis and bottim x-axis to centre, passing through (0,0)
    #ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('zero')
    #ax.spines['bottom'].set_position('axes', -0.002)

    #specify axes range
    ##plt.axis([9638, 9668, -0.0002, 0.002])
    
    #new
    ##plt.plot(Energy, Spectra, markersize=0.1, label=name,)

    # nice formatting :D
    #plt.title(title, fontname='Times New Roman', fontsize=20)
    plt.xlabel('Energy (eV)', fontname='Times New Roman', fontsize=30)
    plt.ylabel('Intensity (a.u.)', fontname='Times New Roman', fontsize=30)
   
    leg = plt.legend(prop=font, loc='upper left', borderaxespad=0, labelspacing=0, handlelength=1.0, handletextpad=0.5)
    leg.get_frame().set_alpha(0)
    
    ax = plt.gca()
    ax.tick_params(axis='both', direction ="in", length=5, width=2, pad=10, labelsize=24, labelrotation=0)
   
    for label in ax.get_xticklabels():
        label.set_fontproperties(font_ticklabel)

    for label in ax.get_yticklabels():
        label.set_fontproperties(font_ticklabel)
        
    plt.yticks(np.arange(0, 0.0024, step=0.0004), np.arange(6))
    
    plt.axis([9638.5, 9661.5, -0.0002, 0.0024])
    
    figure = plt.gcf()

    figure.set_size_inches(8, 6)
    plt.savefig(f'{namefile}.png', bbox_inches='tight', dpi=300)
    ##plt.show()
    
    #plt.savefig(f'{name}.png', bbox_inches='tight')
