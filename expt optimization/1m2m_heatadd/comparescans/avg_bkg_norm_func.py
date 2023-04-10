import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import pandas as pd
import shutil
import scipy.optimize as op


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


#writes data to .txt file in the folder 'Averaged & Normalized' - (creates one if not found)
#args: array-energy, array-data, name
#output: write to .txt file
def write_to_file(Energy, avg_cnts, name, delimiter='\t'): 
    
    # for csv file, change delimiter to ','
    if delimiter is ',':
        ext = '.csv'
    else:
        ext = '.txt'
    
    #check if there is specified directory and create one if not
    average_cnts_directory_name = 'Averaged & Normalized'
    if not os.path.exists(f'{average_cnts_directory_name}/'):
        os.makedirs(f'{average_cnts_directory_name}/')
    
    # move to average cnts directory
    os.chdir(f'{average_cnts_directory_name}/')
    
    # Write the file
    with open(f'{name}{ext}', 'w') as file:
        # Headers first
        file.write(f'Energy_(eV){delimiter}Live_Counts\n')
        # Then the data
        for i in range(len(Energy)):
            file.write(f'{Energy[i]}{delimiter}{avg_cnts[i]}\n')
    
    # return to parent directory
    os.chdir('../')
    
#plot the spectra with style
def plot(Energy, Spectra, title, save=False):

    fig, ax = plt.subplots(1, figsize=(9, 6))

    ax.plot(Energy, Spectra, "darkcyan", markersize=0.1)  # label=folder_list[i])

    # nice formatting :D
    plt.title(title, fontsize=20)
    plt.xlabel('Energy (eV)', fontsize=16)
    plt.ylabel('arb. units', fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax.legend(fontsize=16)

    if save:
        plt.savefig(f'{directory}/{title}.png', bbox_inches='tight')

        
#find the index of a given value in an array
#returns the index for the value
def index_from_array(array, val):
    dif = (array - val)**2
    return np.argmin(dif)


#find the index for ROI in an array
#args: array- Energy, float - start_energy and end_energy
#returns the indices for start and end of ROI
def get_E_region(Energy, start_energy, end_energy):
   # get index from energy values
    start_index = index_from_array(Energy, start_energy)
    end_index = index_from_array(Energy, end_energy)
    return start_index, end_index


#loss for constant background subtraction
#args: arrays
def loss(b, y):
    return np.sum((y - b)**2/100)


#Integral Normalizes the data in the given range
#args: array- avg_cnts, array-Energy, tuple? (indices)- E_range
#returns normalized data as array
def Normalize(avg_cnts, Energy, E_range):
    
    # integral normalizes
    I = 0
    start, end = E_range
    energy = Energy[start:end]
    
    # find intergal in specified energy range
    for i in range(end - start - 1):
        dE = energy[i+1]-energy[i]
        I += dE*avg_cnts[i]
    
    # divide by sum
    return avg_cnts/I


