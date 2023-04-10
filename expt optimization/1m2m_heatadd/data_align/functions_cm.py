import numpy as np
import matplotlib.pyplot as plt
import os
from pylab import *
from math import pi 


#convert energy(eV) to angle(deg) using Bragg's law
#agrs: (list-energy, float-dspacing)
#output: array-angle

def convert_to_angle(energy, d): # 'd' in nm
    #energy = energy.tolist()
    angle = []
    for i in energy:
        angle_a = (1239.8425/(2*d*i))
        angle_b = math.asin(angle_a)*(180/pi)
        angle.append(angle_b)
    angle = np.array(angle, dtype=float)
    return(angle)  


#convert angle(deg) to energy(eV) using Bragg's law
#agrs: (list-angle, float-dspacing)
#output: array-energy

def convert_to_energy(angle, d): # 'd' in nm
    energy = (1239.8425/(2*d))*(1/sin(angle*pi/180))
    return(energy)


#Read in reference spectrum from file
#returns list of values in 1st column (energy) as 'energy_ref' 
#and 2nd column (normalized counts per live) as 'data_ref'

def read_file():
    energy_ref = []
    data_ref = []
    #read the reference spectra from a file and append to respective lists
    for i in open('file_ref_shifted.txt', 'r'):
        segs = i.split()
        energy_ref.append(float(segs[0]))
        data_ref.append(float(segs[1]))
    return(energy_ref, data_ref)


#Define the reading function for data files
#reads the data file and returns a temp list of energy and data values
#data files should have only two columns, 1st column header should read 'energy_(eV)'

def read_data_file(filename):
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


#check for (and make if not found) directory called 'Shifted Data' and 
#write each shifted data to a new .txt file as: originalfilename_shifted.txt

def write_to_file(Angle, Data, name, delimiter='\t'): 
    name = name.strip('.txt')
    # for csv file, change delimiter to ','
    if delimiter is ',':
        ext = '.csv'
    else:
        ext = '.txt'
    
    #check if there is a shifted data directory and create one if not
    shifted_data_directory_name = 'Shifted Data'
    if not os.path.exists(f'{shifted_data_directory_name}/'):
        os.makedirs(f'{shifted_data_directory_name}/')
    
    # move to shifted data directory
    os.chdir(f'{shifted_data_directory_name}/')
    
    # Write the file
    with open(f'{name}_shifted{ext}', 'w') as file:
        # Headers first
        file.write(f'Energy_(eV){delimiter}counts\n')
        # Then the data
        for i in range(len(Angle)):
            file.write(f'{Angle[i]}{delimiter}{Data[i]}\n')
    # return to parent directory
    os.chdir('../')

    
#shift the data (spectra) to match the reference spectrum wrt the center of mass of kbeta peak
#takes in 2 lists (angle, data), 2 arrays (angle_ref_new and data_ref_new) and name of the sample (file)
#output:shows a plot of the shifted data and writes the shifted data in separate .txt files 

def spectra_shift(angle, data, angle_ref_new, data_ref_new,  name):
    #change to arrrays, the other two are already arrays
    angle_new = np.array(angle, dtype=float)
    data_new = np.array(data, dtype=float)

    #shift the data using center of mass to allign the peaks
    #the cm of the peaks
    def center_of_mass(x, y):
        cm = np.sum(x*y)/np.sum(y)
        return(cm)
    
    #specify the ROI of the peak and calculate center of mass
    cm_data = center_of_mass(angle_new[0:200], data_new[0:200])
    cm_ref = center_of_mass(angle_ref_new[0:200], data_ref_new[0:200])
    
    #find the shift needed to allign the cm
    shift = cm_data - cm_ref

    #shift the angle axis
    angle_shifted = angle_new - shift
    data_shifted = data_new
   
    #convert the shifted angle to energy
    energy_shifted = convert_to_energy(angle_shifted, 0.0653269)  # d = 0.0653269 for Ge 555
 
    #show plot of the shifted data
    plt.plot(energy_shifted, data_shifted)
    plt.axis([9640, 9670, 0, 0.002])
     
    # note: this will overwrite any files with the same name
    # name specified from original directory/file name from whence the batches came from
    write = True
    if write:
        write_to_file(energy_shifted, data_shifted, name)

    