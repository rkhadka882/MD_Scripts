#!/usr/bin/env python
# coding: utf-8

#Importing the modules
import numpy as np
import dpdata
import os

print("Note:\n\
---------------------------------------------------------------------------------------------\n\
> No. of atoms should be same in both NPT and NVT data\n\
> Combine the NPT and NVT AIMD runs separately for a given system, use combine_AIMD.sh\n\
> To execute: python JDFTX_2_DEEPMD_v4.py \n\
> Output: type.raw type_map.raw force.raw, energy.raw, coord.raw, box.raw which are shuffled\n\
> Use /deepmd-kit/data/raw/raw_2_set.sh to create multiple sets of data\n\
> This script assumes the following: \n\
>> NPT and NVT runs have same no. of atoms and same atom types only \n\
>> jdftx ionic position are dumped alphabetically for all ionic steps not random \n\
---------------------------------------------------------------------------------------------\n\
\n\
")

# Ask the user for a folder name
folder = input("Enter the folder to put all data eg. CoSi_AIMD: ")

#file_AIMD = input("Enter the combined NPT filename eg. CoSi_AIMD.jdftxout: ")
file_AIMD = f'./{folder}/AIMD.jdftxout'

frame_sep = input("Enter the frame separation in AIMD data eg. 10: ")
frame_sep = int(frame_sep)

#frame_sep = 5 			# Separation between the AIMD data that is taken for training 

# Check if directory exists. If not, create it.
if not os.path.exists(folder):
    os.makedirs(folder)
    print(f"Folder '{folder}' has been created.")
else:
    print(f"Folder '{folder}' already exists.")

atom_types = ['Co', 'Si', 'O']  # CHANGE HERE, Add the types of elements that you use 
                                # in your project even if certain datas has only few elements

# Conversion from jdftx to deepmd or lammps
ene_conv = 27.2114              # hartree to eV conversion
len_conv = 0.529177             # bohr to Angstrom conversion
force_conv = ene_conv/len_conv  # E = fd, so f = E/d 

# Single atom energies:
#E_Co = -148.9280175 * ene_conv # hartree to eV
#E_Si = -3.600491552 * ene_conv # hartree to eV
#E_O = -15.72003172 * ene_conv # hartree to eV

E_Co = -148.5065685 * ene_conv # hartree to eV
E_Si = -3.80475254 * ene_conv # hartree to eV
E_O = -15.8884064 * ene_conv # hartree to eV

# Function definitions:

def get_atoms_info(line):
    natoms = int(line.strip().split()[4])
    nspecies = int(line.strip().split()[1])
    return natoms, nspecies

def get_energy(line):
    return float(line.strip().split()[2])

def get_lattice_vectors(lines, i):
    x = list(map(float, lines[i].strip().replace('[', '').replace(']', '').split()))
    y = list(map(float, lines[i+1].strip().replace('[', '').replace(']', '').split()))
    z = list(map(float, lines[i+2].strip().replace('[', '').replace(']', '').split()))

    box = np.zeros((3,3))
    box[0, :] = np.array(x)
    box[1, :] = np.array(y)
    box[2, :] = np.array(z)
    box_corr = box.T
    x_corr = box_corr[0, :].tolist()
    y_corr = box_corr[1, :].tolist()
    z_corr = box_corr[2, :].tolist()

    return x_corr + y_corr + z_corr

def get_volume(line):
    return float(line.strip().split()[4])

def get_stress_tensor(lines, i, volume):
    x = list(map(float, lines[i].strip().replace('[', '').replace(']', '').split()))
    y = list(map(float, lines[i+1].strip().replace('[', '').replace(']', '').split()))
    z = list(map(float, lines[i+2].strip().replace('[', '').replace(']', '').split()))
    xx = [i*volume for i in x]
    yy = [i*volume for i in y]
    zz = [i*volume for i in z]
    return xx + yy + zz

def get_coordinates(lines, natoms, i):
    clines = lines[i:i+natoms]
    coord = np.array([list(map(float, line.strip().split()[2:5])) for line in clines]).reshape(-1)
    allatoms = [line.strip().split()[1] for line in clines]
    return coord, allatoms

def get_forces(lines, natoms, i):
    flines = lines[i:i+natoms]
    force = np.array([list(map(float, line.strip().split()[2:5])) for line in flines]).reshape(-1)
    return force

def automapping(allatoms):
    """
    Use the atom_types in the input and create type.raw
    and type_map.raw files
    """
    # Identifying unique atom types in the simulation:
    #atom_types = list(set(allatoms)); 

    # Assigning integer numbers to unique elements
    mapping = {atom: i for i, atom in enumerate(atom_types)}; 

    # Replacing elements with their mapping
    atm_mapped = [mapping[atom] for atom in allatoms]

    # Writing type_map.raw:
    with open("type.raw", 'w') as f:
        for value in atm_mapped:
            f.write(f"{value}\n")

    # Writing type.raw:
    with open("type_map.raw", 'w') as f:
        for i in atom_types:
            f.write(f"{i}\n")
    return atm_mapped

def all_convert(filename, frame_sep):
    """
    Collects the Force, enegy, box, virial datas from JDFTXoutput file
    """
    with open(filename) as file:
        lines = file.readlines()

        Etot, Force, Stress, Lattice, Coord, dump_index = [], [], [], [], [], []

        # Getting the dump indexing:
        for i, line in enumerate(lines):
            if "Dumping 'md.ionpos' ... done" in line:
                dump_index.append(i)

        start_ind = dump_index[0:-1]
        end_ind = dump_index[1:]
        total_frame = len(start_ind)
        cframe = 0

        for i in range(total_frame):                          # Looping over all frames
            frame = lines[start_ind[i]:end_ind[i]]            # Get the frame
            
            for linee in frame:
                if "Converged" in linee:
                    cframe = cframe + 1
                    good_frame = frame

                    for j, nline in enumerate(good_frame):
                        if "Initialized" in nline:
                            natoms, nspecies = get_atoms_info(nline)
                        if "F =    " in nline:
                            Etot.append(get_energy(nline))
                        if "# Lattice vectors:" in nline:
                            Lattice.append(get_lattice_vectors(good_frame, j+2))
                        if "unit cell volume =" in nline:
                            volume = get_volume(nline)
            #             if "# Stress tensor in Cartesian coordinates" in line:
            #                 Stress.append(get_stress_tensor(lines, i+1, volume))
                        if "Ionic positions in cartesian coordinates" in nline:
                            coord, allatoms = get_coordinates(good_frame, natoms, j+1)
                            Coord.append(coord)
                        if "Forces in Cartesian coordinates" in nline:
                            Force.append(get_forces(good_frame, natoms, j+1))
    uframe = total_frame-cframe
    print("Converged frames:", cframe, "Unconverged frames:", uframe)
    
    # Converting the list to array
    Etot = np.array(Etot)
    Lattice = np.array(Lattice)
    Coord = np.array(Coord)
    Force = np.array(Force)
    
    # Sampling the AIMD data:
    AIMD_frame_no = len(Etot) // frame_sep
    Frame_index = np.arange(AIMD_frame_no)*frame_sep

    Etot = Etot[Frame_index]
    Lattice = Lattice[Frame_index]
    Coord = Coord[Frame_index]
    Force = Force[Frame_index]

    # Shuffling the data:
    print(">Force, Energy, Coordinate, Lattice converted, shuffled, and saved")
    permuted_indices = np.random.permutation(len(Etot))
    Etot = np.array(Etot)[permuted_indices]*ene_conv
    Lattice = np.array(Lattice)[permuted_indices]*len_conv
    Coord = np.array(Coord)[permuted_indices]*len_conv
    Force = np.array(Force)[permuted_indices]*force_conv
#     Stress = np.array(Stress)[permuted_indices]*ene_conv
    
    # Energy correction:
    Etot = Etot - Energy_atom(allatoms)

    # Saving the data:
    np.savetxt("energy.raw", Etot)
    np.savetxt("box.raw", Lattice)
    np.savetxt("coord.raw", Coord)
    np.savetxt("force.raw", Force)
#     np.savetxt("virial.raw", Stress)
    
    return Etot, Force, Coord, Lattice, allatoms, natoms, nspecies, cframe, uframe

def Energy_atom(allatoms):
    """ Calculate the energy of atoms to be subtracted from bulk """
    # Assigning integer numbers to unique elements
    mapping = {atom: i for i, atom in enumerate(atom_types)}; 
    # Replacing elements with their mapping
    atm_mapped = [mapping[atom] for atom in allatoms]
    nCo = atm_mapped.count(0)
    nSi = atm_mapped.count(1)
    nO = atm_mapped.count(2)
    print('Co:', nCo, 'Si:', nSi,'O:', nO)
    E_atom = nCo*E_Co + nSi*E_Si + nO*E_O
    return E_atom

# AIMD Conversion:
print("!!!Converting the AIMD data!!!")
print('-'*50)

# Force, energy, virial, coord:
Etot_AIMD, Force_AIMD, Coord_AIMD, Lattice_AIMD, allatoms_AIMD, natoms, nspecies, Cframe, Uframe = all_convert(file_AIMD, frame_sep)

# Energy correction:
#Etot_AIMD = Etot_AIMD - Energy_atom(allatoms_AIMD)

# type, type_map
atom_map_AIMD = automapping(allatoms_AIMD)

#---------------------------------------------------------------------------------------------

# Pringin the summary information:
print("# Summary of jdftx to DeepMD\n")
print('-'*50)
print("Converged frames  :", Cframe)
print("Unconverged frames:", Uframe)
print("No. of Elements   :", nspecies)
print("Elements involved :", list(set(allatoms_AIMD)))
print("All Elements      :", atom_types)
# print("Element to int map:", mapping)
print("No. of atoms      :", natoms)
print("No. of frames     :", Force_AIMD.shape[0])
print("force.raw         :", Force_AIMD.shape)
print("coord.raw         :", Coord_AIMD.shape)
print("energy.raw        :", Etot_AIMD.shape)
print("box.raw           :", Lattice_AIMD.shape)
# print("virial.raw        :", Stress_AIMD.shape)
print('-'*50)
    
with open("conversionInfo.dat", 'w') as f:
    f.write("# Summary of jdftx to DeepMD\n")
    f.write('-'*50 + '\n')
    f.write("Converged frames  :" + f'{Cframe}' + '\n')
    f.write("Unconverged frames:" + f'{Uframe}' + '\n')
    f.write("No. of Elements   :" + f'{nspecies}'+ '\n')
    f.write("Elements involved :" + f'{list(set(allatoms_AIMD))}'+ '\n')
    f.write("All Elements      :" + f'{atom_types}'+ '\n')
    f.write("No. of atoms      :" + f'{natoms}'+ '\n')
    f.write("No. of frames     :" + f'{Force_AIMD.shape[0]}'+ '\n')
    f.write("force.raw         :" + f'{Force_AIMD.shape}'+ '\n')
    f.write("coord.raw         :" + f'{Coord_AIMD.shape}'+ '\n')
    f.write("energy.raw        :" + f'{Etot_AIMD.shape}'+ '\n')
    f.write("box.raw           :" + f'{Lattice_AIMD.shape}'+ '\n')
#     f.write("virial.raw        :" + f'{Stress_AIMD.shape}'+ '\n')
    f.write('-'*50)
    
# Cleaning up 
os.system(f'mv -f *.raw conversionInfo.dat {folder}')

# Creating the training and validation data sets:
if os.path.exists(f'{folder}/training_data'):
    os.system(f'rm -r {folder}/training_data')
    
if os.path.exists(f'{folder}/validation_data'):
    os.system(f'rm -r {folder}/validation_data')

if not os.path.exists(f'{folder}/training_data'):
    os.makedirs(f'{folder}/training_data')
    print(f"'training_data' has been created.")

if not os.path.exists(f'{folder}/validation_data'):
    os.makedirs(f'{folder}/validation_data')
    print(f"'validation_data' has been created.")

# Copying the type.raw and type_map.raw to training_data and validation_data
os.system(f'cp {folder}/type* {folder}/training_data')
os.system(f'cp {folder}/type* {folder}/validation_data')


n_sets = int(Etot_AIMD.shape[0] // 4.5)

os.system(f'cd {folder}; /gpfs/u/home/TMGN/TMGNrjkh/scratch/DeePMD_calculation/deepmd-kit/data/raw/raw_to_set.sh {n_sets}')

os.system(f'mv {folder}/set.000 {folder}/validation_data')
os.system(f'mv {folder}/set.* {folder}/training_data')
