#!/usr/bin/env python

# execute: python3 JDFTX_2_DEEPMD_v2.py filename.jdftxout 
# Rajan Khadka 10/5/2023
# make sure to have numpy, dpdata, os, sys, ase packages in python environment
# collects the raw data, shuffles them 
# Use deepmd-kit/data/raw/raw_to_set.sh to divide the collected data into multiple sets

# Importing the modules
import numpy as np
import dpdata
import os
import sys
from ase.data import chemical_symbols, atomic_numbers

# cat md.* > compiled.jdftxout # To compile all the AIMD runs
# we will use compiled.jdftxout to create the training set


atom_types = ['Co', 'Si', 'O']  # Add the types of elements that you use 
                                # in your project even if certain datas has only few elements

# Conversion from jdftx to deepmd or lammps
ene_conv = 27.2114              # hartree to eV conversion
len_conv = 0.529177             # bohr to Angstrom conversion
force_conv = ene_conv/len_conv  # E = fd, so f = E/d 

filename=sys.argv[1]            # providing filename as first argument

print("Note:")
print('-'*100)
print(">Uses *.jdftxout file which is jdftxoutput of single or multiple runs combined")
print(">Only combine jdftxout for same system eg. CoSi")
print(">Only combine jdftxout with same number of atoms")
print(">Use deepmd-kit/data/raw/raw_to_set.sh to divide the collected data into multiple sets")
print(">Rajan Khadka, 10/4/2023")
print('-'*100)

with open(filename) as file:
    lines = file.readlines()
    total_lines = len(lines)
    Etot = []                    # Extracting the total energy
    Force = []                   # Extracting the forces fx, fy, fz on ions
    Stress = []                  # Extracting the stress tensor in Cartesian coordinates
    Lattice = []                 # Extractin the lattice vectors
    Coord = []                   # Extracting ionic position x, y, z

    for i, line in enumerate(lines):
        # no. of atoms:
        if "Initialized" in line:
            natoms = int(line.strip().split()[4])
            nspecies = int(line.strip().split()[1])
            
    for i, line in enumerate(lines):
        
        # Energy information:
        if "Etot =    " in line:
            Etot.append(float(lines[i].strip().split()[2]))
        
        # Lattice vector information:
        if "# Lattice vectors:" in line:
            x = list(map(float, lines[i+2].strip().replace('[', '').replace(']', '').split()))
            y = list(map(float, lines[i+3].strip().replace('[', '').replace(']', '').split()))
            z = list(map(float, lines[i+4].strip().replace('[', '').replace(']', '').split()))
            Lattice.append(x + y + z)
        
        if "unit cell volume =" in line:
                volume = float(lines[i].strip().split()[4])
                
        # Stress tensor information:
        if "# Stress tensor in Cartesian coordinates" in line:
            x = list(map(float, lines[i+1].strip().replace('[', '').replace(']', '').split()))
            y = list(map(float, lines[i+2].strip().replace('[', '').replace(']', '').split()))
            z = list(map(float, lines[i+3].strip().replace('[', '').replace(']', '').split()))
            xx = [i*volume for i in x]
            yy = [i*volume for i in y]
            zz = [i*volume for i in z]
            Stress.append(xx + yy + zz)
            
        # Coordinate information:
        if "Ionic positions in cartesian coordinates" in line:
            clines = lines[i+1:i+1+natoms]
            
            # [2:5] because that is the location of x, y, z coordinates:
            coord = np.array([list(map(float, line.strip().split()[2:5])) for line in clines]).reshape(-1)
            Coord.append(coord)
            
            # atom-type:
            allatoms = [line.strip().split()[1] for line in clines]

        # Force information:
        if "Forces in Cartesian coordinates" in line:
            flines = lines[i+1:i+1+natoms]
            
            # [2:5] because that is the location of x, y, z coordinates:
            force = np.array([list(map(float, line.strip().split()[2:5])) for line in flines]).reshape(-1)
            Force.append(force)

def automapping():
    # Identifying unique atom types in the simulation:
    atom_types = list(set(allatoms)); 

    # Assigning integer numbers to unique elements
    mapping = {atom: i for i, atom in enumerate(atom_types)}; 

    # Replacing elements with their mapping
    atm_mapped = [mapping[atom] for atom in allatoms]

    # Writing type.raw:
    with open("type.raw", 'w') as f:
        for value in atm_mapped:
            f.write(f"{value}\n")

    # Writing type_map.raw:
    with open("type_map.raw", 'w') as f:
        for i in atom_types:
            f.write(f"{i}\n")
    return 

def manualmapping():
    atomid = [atomic_numbers[i] for i in allatoms]

    # Writing type_map.raw:
    with open("type.raw", 'w') as f:
        for value in atomid:
            f.write(f"{value}\n")

    # Writing type.raw:
    with open("type_map.raw", 'w') as f:
        for i in atom_types:
            f.write(f"{i}\n")
    return

# Choose the type of mapping
automapping()             # Automatically identifies atom types and assigns the integer
# manualmapping()           # Uses the provided atom types, and assigns atomic number as integers


print(">Saving the text data file\n")
print(">Converting the .raw to .npy\n")

permuted_indices = np.random.permutation(len(Etot))

# Energy data:
Etot = np.array(Etot)*ene_conv
np.savetxt("unshuffled_energy.raw", Etot)
np.savetxt("energy.raw", Etot[permuted_indices])
# np.save("energy.npy", Etot)

# Lattice data:
Lattice = np.array(Lattice)*len_conv
np.savetxt("unshuffled_box.raw", Lattice)
np.savetxt("box.raw", Lattice[permuted_indices])
# np.save("box.npy", Lattice)

# position data:
Coord = np.array(Coord)*len_conv
np.savetxt("unshuffled_coord.raw", Coord)
np.savetxt("coord.raw", Coord[permuted_indices])
# np.save("coord.npy", Coord)

# Force data:
Force = np.array(Force)*force_conv
np.savetxt("unshuffled_force.raw", Force)
np.savetxt("force.raw", Force[permuted_indices])
# np.save("force.npy", Force)

# Stress data:
Stress = np.array(Stress)*ene_conv
np.savetxt("unshuffled_virial.raw", Stress)
np.savetxt("virial.raw", Stress[permuted_indices])
# np.save("virial.npy", Stress)

print("# Summary of jdftx to DeepMD\n")
print('-'*50)
print("No. of Elements   :", nspecies)
print("Elements involved :", list(set(allatoms)))
print("All Elements      :", atom_types)
# print("Element to int map:", mapping)
print("No. of atoms      :", natoms)
print("No. of frames     :", Force.shape[0])
print("force.raw         :", Force.shape)
print("coord.raw         :", Coord.shape)
print("energy.raw        :", Etot.shape)
print("box.raw           :", Lattice.shape)
print("virial.raw        :", Stress.shape)
print('-'*50)
    
with open("conversionInfo.dat", 'w') as f:
    f.write("# Summary of jdftx to DeepMD\n")
    f.write('-'*50 + '\n')
    f.write("No. of Elements   :" + f'{nspecies}'+ '\n')
    f.write("Elements involved :" + f'{atom_types}'+ '\n')
#     f.write("Element to int map:" + f'{mapping}'+ '\n')
    f.write("No. of atoms      :" + f'{natoms}'+ '\n')
    f.write("No. of frames     :" + f'{Force.shape[0]}'+ '\n')
    f.write("force.raw         :" + f'{Force.shape}'+ '\n')
    f.write("coord.raw         :" + f'{Coord.shape}'+ '\n')
    f.write("energy.raw        :" + f'{Etot.shape}'+ '\n')
    f.write("box.raw           :" + f'{Lattice.shape}'+ '\n')
    f.write("virial.raw        :" + f'{Stress.shape}'+ '\n')
    f.write('-'*50)

# if not os.path.exists("set.000"):
#     os.mkdir("set.000")
# os.system('mv *.npy set.000')

if not os.path.exists("unshuffled_raw_data"):
    os.mkdir("unshuffled_raw_data")

os.system('mv -f unshuffled_*.raw unshuffled_raw_data')

if not os.path.exists("shuffled_raw_data"):
    os.mkdir("shuffled_raw_data")
    
os.system('mv -f virial.raw force.raw energy.raw coord.raw box.raw shuffled_raw_data')

if not os.path.exists("training_data"):
    os.mkdir("training_data")
os.system('mv -f type.raw type_map.raw shuffled_raw_data unshuffled_raw_data training_data')

if not os.path.exists("validation_data"):
    os.mkdir("validation_data")

