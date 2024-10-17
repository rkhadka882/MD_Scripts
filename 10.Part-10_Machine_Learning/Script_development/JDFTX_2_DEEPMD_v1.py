#!/usr/bin/env python

# Importing the modules
import numpy as np
import dpdata
import os
import sys

# cat md.* > compiled.jdftxout # To compile all the AIMD runs
# we will use compiled.jdftxout to create the training set


def ls():
    """ list out the file in a directory"""
    return os.system('ls')

# Conversion from jdftx to deepmd or lammps
ene_conv = 27.2114              # hartree to eV conversion
len_conv = 0.529177             # bohr to Angstrom conversion
force_conv = ene_conv/len_conv  # E = fd, so f = E/d 

filename= sys.argv[1]           # providing filename as first argument

print("Note:")
print('-'*100)
print(">Uses *.jdftxout file which is jdftxoutput of single or multiple runs combined")
print(">Only combine jdftxout for same system eg. CoSi")
print(">Only combine jdftxout with same number of atoms")
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
        cd s
            
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
        
        
        # Stress tensor information:
        if "# Stress tensor in Cartesian coordinates" in line:
            xx = list(map(float, lines[i+1].strip().replace('[', '').replace(']', '').split()))
            yy = list(map(float, lines[i+2].strip().replace('[', '').replace(']', '').split()))
            zz = list(map(float, lines[i+3].strip().replace('[', '').replace(']', '').split()))
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

# Identifying unique atom types in the simulation:
atom_types = list(set(allatoms)); 

# Assigning integer numbers to unique elements
mapping = {atom: i for i, atom in enumerate(atom_types)}; 

# Replacing elements with their mapping
atm_mapped = [mapping[atom] for atom in allatoms]

# Writing type_map.raw:
with open("type_map.raw", 'w') as f:
    for value in atm_mapped:
        f.write(f"{value}\n")

# Writing type.raw:
with open("type.raw", 'w') as f:
    for i in atom_types:
        f.write(f"{i}\n")

print(">Saving the text data file\n")
print(">Converting the .raw to .npy\n")

# Energy data:
Etot = np.array(Etot)*ene_conv
np.savetxt("energy.raw", Etot)
np.save("energy.npy", Etot)

# Lattice data:
Lattice = np.array(Lattice)*len_conv
np.savetxt("box.raw", Lattice)
np.save("box.npy", Lattice)

# position data:
Coord = np.array(Coord)*len_conv
np.savetxt("coord.raw", Coord)
np.save("coord.npy", Coord)

# Force data:
Force = np.array(Force)*force_conv
np.savetxt("force.raw", Force)
np.save("force.npy", Force)

# Stress data:
Stress = np.array(Stress)*force_conv/len_conv**2
np.savetxt("virial.raw", Stress)
np.save("virial.npy", Stress)

print("# Summary of jdftx to DeepMD")
print('-'*50)
print("No. of Elements   :", nspecies)
print("Elements involved :", atom_types)
print("Element to int map:", mapping)
print("No. of atoms      :", natoms)
print("No. of frames     :", Force.shape[0])
print("force.raw         :", Force.shape)
print("coord.raw         :", Coord.shape)
print("energy.raw        :", Etot.shape)
print("box.raw           :", Lattice.shape)
print("virial.raw        :", Stress.shape)
print('-'*50)


if not os.path.exists("set.000"):
    os.mkdir("set.000")
os.system('mv *.npy set.000')

if not os.path.exists("set.raw"):
    os.mkdir("set.raw")
    
os.system('mv virial.raw force.raw energy.raw coord.raw box.raw set.raw')

if not os.path.exists("training_data"):
    os.mkdir("training_data")
os.system('mv type.raw type_map.raw set.000 set.raw training_data')

if not os.path.exists("validation_data"):
    os.mkdir("validation_data")


