import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from math import pi

#convert energy(eV) to angle(deg)
def convert_to_angle(energy, d = 0.0653269): # 'd' in nm  - - d for Ge 555 is 0.0653269 nm
    angle = []
    for i in energy:
        angle_a = (1239.8425/(2*d*i))
        angle_b = math.asin(angle_a)*(180/pi)
        angle.append(angle_b)
    angle = np.array(angle, dtype=float)
    return(angle)  


#convert angle(deg) to energy(eV)
def convert_to_energy(angle, d = 0.0653269): # 'd' in nm - - d for Ge 555 is 0.0653269 nm
    energy = (1239.8425/(2*d))*(1/sin(angle*pi/180))
    return(energy)


#write the final shifted data to a new file
def write_to_file(Angle, Data, name, delimiter='\t'): 
    # for csv file, change delimiter to ','
    if delimiter is ',':
        ext = '.csv'
    else:
        ext = '.txt'
    
    # Write the file
    with open(f'{name}_shifted{ext}', 'w') as file:
        # Headers first
        #file.write(f'energy_(eV){delimiter}counts\n')
        # Then the data
        for i in range(len(Angle)):
            file.write(f'{Angle[i]}{delimiter}{Data[i]}\n')

            
def convert_shift(peak_lit, peak_my):  # args are both peak positions (single floats) in eV
    angle_lit = convert_to_angle([peak_lit])
    angle_my = convert_to_angle([peak_my])
    shift = angle_lit - angle_my
    return(shift)


def spectra_shift(angle, shift):  #both args are ndarray
    angle_new = angle + shift
    return(angle_new)


def read_file(file):
    energy = []
    data = []

#read spectra from a file and append to respective lists
    for i in open(f'{file}', 'r'):
        segs = i.split()
        energy.append(float(segs[0]))
        data.append(float(segs[1]))
    return(energy, data)


