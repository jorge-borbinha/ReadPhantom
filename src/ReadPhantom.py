#!/usr/bin/env python3
# @author:  Jorge Cebola Borbinha
# @github:  jorge-borbinha
# @website: jorge-borbinha.github.io
#
# Copyright (C) 2025 Jorge Cebola Borbinha
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License v3 (AGPLv3) as 
# published by the Free Software Foundation.
# For more information, see https://www.gnu.org/licenses/agpl-3.0.html
#
# This software is for educational and research purposes only. The user is 
# solely responsible for ensuring data privacy and must de-identify all patient
# data before use. The authors provide this software "as is" without any
# warranty. By using this software, you agree to these terms.
#

import numpy as np
import pandas as pd
import struct
import os

def read_progress(z):
    """Prints a progress message based on the current slice number. Serves as 
    visual indicator of the program's progress during file reading/writing. It 
    prints a message to the console when the number of processed slices reaches
    specific milestones.
    
     Parameters
     ----------
     z : int
         Current no. of slices that have been processed.
    
     Returns
     -------
     None
         Function only prints to the console.
     """
     
    if z in [50, 100, 200, 300, 350, 400]:
        print(f" > Number of slices read: {z}")
        


def create_vox_file (vox_file, arr_material, arr_density, n_vox_x, n_vox_y, n_vox_z,
                       vox_res_x, vox_res_y, vox_res_z):
    """Creates a .vox voxel phantom file in the format required by the PENELOPE/
    penEasy simulation framework. It first writes a fixed and then appends the 
    material and density data for each voxel. The voxel data is written 
    efficiently using NumPy's `savetxt` function.

    Parameters
    ----------
    vox_file : str
        The name of the output .vox file (e.g., 'phantom.vox').
    arr_material : numpy.ndarray
        A 1D array of integers representing the material ID for each voxel.
    arr_density : numpy.ndarray
        A 1D array of floats representing the density (in g/cm**3) for each voxel.
    n_vox_x : int
        The number of voxels along the x-axis.
    n_vox_y : int
        The number of voxels along the y-axis.
    n_vox_z : int
        The number of voxels along the z-axis.
    vox_res_x : float
        The physical size of a voxel along the x-axis /cm.
    vox_res_y : float
        The physical size of a voxel along the y-axis /cm.
    vox_res_z : float
        The physical size of a voxel along the z-axis /cm.

    Returns
    -------
    None
        The function writes directly to a file.
    """
    
    # Create the phantom.vox file
    print('Creating the VOX file...\n')
    with open(vox_file, 'w') as f:
        # VOX file header
        f.write('[SECTION VOXELS HEADER v.2008-04-13]\n')
        f.write(f' {n_vox_x:4d}{n_vox_y:4d}{n_vox_z:4d}\n')
        f.write(f' {vox_res_x:7.5f} {vox_res_y:7.5f} {vox_res_z:7.5f}\n')
        f.write(' 1\n')
        f.write(' 2\n')
        f.write(' 0\n')
        f.write('[END OF VXH SECTION]\n')

    # Combine arr_meterial and arr_density to use np.savetxt(), more efficient
    combined_array_mat_den = np.column_stack((arr_material, arr_density))
    with open('phantom.vox','a') as f1:
        np.savetxt(f1, combined_array_mat_den, fmt = '%3d %7.4f')
    
    print(f'The {vox_file} file was created.\n')

    def create_ct_den_mat_file(file_name, arr_density, arr_material, arr_organ, n_vox_x, n_vox_y, n_vox_z, 
                                 len_x, len_y, len_z, order_indices, loop_order):
        """
        Writes CT voxel data (density, material, and organ ID) to a file in GNUPLOT 
        format, which includes a header and array data. The resulting file is 
        compatible with GNUPLOT scripts available with the PENELOPE/penEasy 
        framework. The function is optimized for efficiency by building the entire 
        data section in memory before writing to the file in a single bulk 
        operation, which is significantly faster than writing line-by-line in
        nested loops.

        Parameters
        ----------
        file_name : str
            The full path and filename for the output file (e.g., 'ct-den-matXY.dat').
        arr_density : numpy.ndarray
            A 1D array containing the density values for each voxel.
        arr_material : numpy.ndarray
            A 1D array containing the material IDs for each voxel.
        arr_organ : numpy.ndarray
            A 1D array containing the organ IDs for each voxel.
        n_vox_x : int
            The number of voxels in the X-dimension.
        n_vox_y : int
            The number of voxels in the Y-dimension.
        n_vox_z : int
            The number of voxels in the Z-dimension.
        len_x : float
            The total length of the phantom in the X-dimension / cm.
        len_y : float
            The total length of the phantom in the Y-dimension / cm.
        len_z : float
            The total length of the phantom in the Z-dimension /cm.
        order_indices : callable
            Lambda function that takes the loop counters and returns the
            correctly ordered (x, y, z) tuple corresponding to the original flattened 
            array indexing order. Allows for the same nested loop to be used for 
            all ct-den-mat files.
        loop_order : tuple
            Tuple containing the dimensions, in the correct order, to loop over
            for generating the file (e.g., (n_vox_z, n_vox_y, n_vox_x) for an XY file).

        Returns
        -------
        None
            The function writes directly to a file.
        """
        
        with open(file_name, 'w') as f:
            # Write headers of ct-den-mat file
            f.write('#  CT structure (GNUPLOT format).\n')
            f.write(f"#  CT enclosure limits:  XL = {0.0:.6e} cm,  XU = {len_x:.6e} cm\n")
            f.write(f"#                       YL = {0.0:.6e} cm,  YU = {len_y:.6e} cm\n")
            f.write(f"#                       ZL = {0.0:.6e} cm,  ZU = {len_z:.6e} cm\n")
            f.write(f"#  Numbers of voxels:    NVX = {n_vox_x}, NVY = {n_vox_y}, NVZ = {n_vox_z}\n")
            f.write('#\n')
            f.write('#\n')
            f.write('#  columns 1 to 3: bin indices IX, IY and IZ\n')
            f.write('#  4th column: density (g/cm**3).\n')
            f.write('#  5th column: material. 6th column: organ ID\n')
            f.write('#  CT structure (GNUPLOT format).\n')

            # Write lines of ct-den-mat files in nested loops to "lines" list
            # Loop order is the list that allows for the same nested loop to be used for all ct-den-mat files
            lines = []
            for i_loop in range(loop_order[0]):
                for j_loop in range(loop_order[1]):
                    for k_loop in range(loop_order[2]):
                        # Calculate counter, order_indices lambda function is what gives us usual (x, y, z) iteration
                        # in order to access/perform calculations on arr_density and other arrays
                        x, y, z = order_indices(i_loop, j_loop, k_loop)
                        counter = x + n_vox_x * (y + n_vox_y * z)

                        lines.append(f" {x+1:3d} {y+1:3d} {z+1:3d} {arr_density[counter]:7.5e} {arr_material[counter]:4d} {arr_organ[counter]:4d}\n")
                    
                    lines.append(' \n')
                
                lines.append(' \n')

                read_progress(i_loop)
            
            # Write all lines at once. 
            f.writelines(lines)
        
        print(f"\nFile {file_name} created.\n")
        
        print("--- FILE VERIFICATION CT-DEN-MAT ---")
        print("\n--- First 10 lines of the ct-den-mat file: ---\n")
        with open(file_name, 'r') as f:
            for i in range(20):  # 7 header lines + 10 data lines
                print(f.readline().strip())
                
        print("\n--- First 10 data lines that do not have two zeros: ---\n")
        with open(file_name, 'r') as f:
            # Skip the 7 header lines
            for _ in range(10):
                next(f)
            # Read and process data lines, skipping the zero-value ones
            count = 0
            for line in f:
                # Split the line into material and density values
                values = line.strip().split()
                if not values:
                    continue
                density = int(values[3])
                material = float(values[4])
                # Check if both values are not zero.
                if material != 0 or density != 0.0:
                    print(line.strip())
                    count += 1
                    if count >= 10:
                        break



def create_ct_den_mat_file(file_name, arr_density, arr_material, arr_organ, n_vox_x, n_vox_y, n_vox_z, 
                         len_x, len_y, len_z, order_indices, loop_order):
    """
    Writes CT voxel data (density, material, and organ ID) to a file in GNUPLOT 
    format, which includes a header and array data. The resulting file is 
    compatible with GNUPLOT scripts available with the PENELOPE/penEasy 
    framework. The function is optimized for efficiency by building the entire 
    data section in memory before writing to the file in a single bulk 
    operation, which is significantly faster than writing line-by-line in
    nested loops.
    
    Parameters
    ----------
    file_name : str
        The full path and filename for the output file (e.g., 'ct-den-matXY.dat').
    arr_density : numpy.ndarray
        A 1D array containing the density values for each voxel.
    arr_material : numpy.ndarray
        A 1D array containing the material IDs for each voxel.
    arr_organ : numpy.ndarray
        A 1D array containing the organ IDs for each voxel.
    n_vox_x : int
        The number of voxels in the X-dimension.
    n_vox_y : int
        The number of voxels in the Y-dimension.
    n_vox_z : int
        The number of voxels in the Z-dimension.
    len_x : float
        The total length of the phantom in the X-dimension / cm.
    len_y : float
        The total length of the phantom in the Y-dimension / cm.
    len_z : float
        The total length of the phantom in the Z-dimension /cm.
    order_indices : callable
        Lambda function that takes the loop counters and returns the
        correctly ordered (x, y, z) tuple corresponding to the original flattened 
        array indexing order. Allows for the same nested loop to be used for 
        all ct-den-mat files.
    loop_order : tuple
        Tuple containing the dimensions, in the correct order, to loop over
        for generating the file (e.g., (n_vox_z, n_vox_y, n_vox_x) for an XY file).
    
    Returns
    -------
    None
        The function writes directly to a file.
    """
    
    with open(file_name, 'w') as f:
        # Write headers of ct-den-mat file
        f.write('#  CT structure (GNUPLOT format).\n')
        f.write(f"#  CT enclosure limits:  XL = {0.0:.6e} cm,  XU = {len_x:.6e} cm\n")
        f.write(f"#                       YL = {0.0:.6e} cm,  YU = {len_y:.6e} cm\n")
        f.write(f"#                       ZL = {0.0:.6e} cm,  ZU = {len_z:.6e} cm\n")
        f.write(f"#  Numbers of voxels:    NVX = {n_vox_x}, NVY = {n_vox_y}, NVZ = {n_vox_z}\n")
        f.write('#\n')
        f.write('#\n')
        f.write('#  columns 1 to 3: bin indices IX, IY and IZ\n')
        f.write('#  4th column: density (g/cm**3).\n')
        f.write('#  5th column: material. 6th column: organ ID\n')
        f.write('#  CT structure (GNUPLOT format).\n')
    
        # Write lines of ct-den-mat files in nested loops to "lines" list
        # Loop order is the list that allows for the same nested loop to be used for all ct-den-mat files
        lines = []
        for i_loop in range(loop_order[0]):
            for j_loop in range(loop_order[1]):
                for k_loop in range(loop_order[2]):
                    # Calculate counter, order_indices lambda function is what gives us usual (x, y, z) iteration
                    # in order to access/perform calculations on arr_density and other arrays
                    x, y, z = order_indices(i_loop, j_loop, k_loop)
                    counter = x + n_vox_x * (y + n_vox_y * z)
    
                    lines.append(f" {x+1:3d} {y+1:3d} {z+1:3d} {arr_density[counter]:7.5e} {arr_material[counter]:4d} {arr_organ[counter]:4d}\n")
                
                lines.append(' \n')
            
            lines.append(' \n')
    
            read_progress(i_loop)
        
        # Write all lines at once. 
        f.writelines(lines)
    
    print(f"\nFile {file_name} created.\n")
    
    # print("--- FILE VERIFICATION CT-DEN-MAT ---")
    # print("\n--- First 10 lines of the ct-den-mat file: ---\n")
    # with open(file_name, 'r') as f:
    #     for i in range(20):  # 7 header lines + 10 data lines
    #         print(f.readline().strip())
            
    # print("\n--- First 10 data lines that do not have two zeros: ---\n")
    # with open(file_name, 'r') as f:
    #     # Skip the 7 header lines
    #     for _ in range(11):
    #         next(f)
    #     # Read and process data lines, skipping the zero-value ones
    #     count = 0
    #     for line in f:
    #         # Split the line into material and density values
    #         values = line.strip().split()
    #         if not values:
    #             continue
    #         density = float(values[3])
    #         material = int(values[4])
    #         # Check if both values are not zero.
    #         if material != 0 or density != 0.0:
    #             print(line.strip())
    #             count += 1
    #             if count >= 10:
    #                 break



# Main program
def readPhantom():
    """
    This program reads a voxel phantom file, processes the data, and
    generates files for visualization (.vox) and simulation (ct-den-mat).
    
    """
    
    # Load phantom file name
    print('Please insert the name of the phantom file:')
    phantom_file = input('> ').strip()

    # Check if phantom file is in binary or ASCII
    while True:
        try:
            print('\nIs your file in binary or in ASCII? (type 0 for binary or 1 for ASCII)')
            file_type = int(input('> '))
            if file_type in [0, 1]:
                break
            else:
                print('The number you typed is not valid. Type 0 for bin or 1 for ASCII.')
        except ValueError:
            print('The number you typed is not valid. Type 0 for bin or 1 for ASCII.')

    # Verify if phantom file exists, try opening phantom_file.
    if file_type == 0:  # Binary file
        try:
            with open(phantom_file, 'rb') as f:
                pass
        except FileNotFoundError:
            print(f"Error: Binary file '{phantom_file}' not found.")
            return

    elif file_type == 1:  # ASCII file
        try:
            with open(phantom_file, 'rb') as f:
                pass
        except FileNotFoundError:
            print(f"Error: ASCII file '{phantom_file}' not found.")
            return


    # Ascertain addictional info about the phantom
    print('\nType in the number of voxels of the phantom in x,y,z.')
    n_vox_x, n_vox_y, n_vox_z = map(int, input('> ').split())

    print('\nType in the voxel resolution in x,y,z /cm.')
    vox_res_x, vox_res_y, vox_res_z = map(float, input('> ').split())

    print('\nType in the number of different materials in the phantom.')
    material_num = int(input('> '))

    print('\nType in the number of organ IDs - i.e. number of lines in organlist.dat.')
    organ_num = int(input('> '))
    
    # n_vox_x = 299
    # n_vox_y = 137
    # n_vox_z = 348
    # vox_res_x = 0.21
    # vox_res_y = 0.21
    # vox_res_z=0.8
    # material_num=53
    # organ_num=141
    
    # Calculate phantom statistics, present them to user and ask user for permission to continue
    n_vox_tot = n_vox_x * n_vox_y * n_vox_z
    len_x = vox_res_x * n_vox_x
    len_y = vox_res_y * n_vox_y
    len_z = vox_res_z * n_vox_z
    
    # Check number of materials
    if material_num <= 1:
        raise ValueError('The number of materials must be larger than 1.')
        return

    print('\n>>> Characteristics of your phantom file:')
    print(f' > Number of voxels in x,y,z: {n_vox_x} {n_vox_y} {n_vox_z}')
    print(f' > Voxel resolution in x,y,z /cm: {vox_res_x:.5f} {vox_res_y:.5f} {vox_res_z:.5f}')
    print(f' > Phantom size (approximate value) in x,y,z /cm: {len_x:.5f} {len_y:.5f} {len_z:.5f}')
    print(f' > Total number of voxels: {n_vox_tot}')
    print('\nPlease check if these values are correct. Do you wish to continue? (y/n)')
    
    while True:
        choice_continue = input().strip().lower()
        if choice_continue in ['y', 'n']:
            break
        print(' Please type y or n.')

    if choice_continue == 'n':
        raise SystemExit(' The script will stop at your request.')
        return
     
    
    # Read the phantom and organlist files into databases
    print('\nReading the phantom and organlist files...')
    
    # Load organlist file name
    print("\nWhat is the name of the organlist file?")
    organlist_file = input('> ').strip()
    
    print('\nShould any rows be skipped when reading the organlist file (includes headers)? If yes, how many? Default is 4.')
    organlist_skip_rows = int(input('> ').strip())
    if not organlist_skip_rows:
        organlist_skip_rows=4
    
    print('\nWhat are the headers of the columns in the organlist file? Write the names separated by commas.\
It is strongly advised that the names contain no spaces.\
Default names are: "Organ_ID", "Organ", "Material_ID", "Density").\
Warning: Number of column headers provided must be equal to number of columns in the organlist file.\n')
    organlist_df_headers = input('> ').strip().split(',')
    # Check if the list contains at least one empty string, or if the original input was empty
    # Split the input
    organlist_df_headers = [header.strip() for header in organlist_df_headers]
    if not organlist_df_headers or any(not header for header in organlist_df_headers):
        organlist_df_headers = ['Organ_ID', 'Organ', 'Material_ID', 'Density']
        print(f'The names for the columns in the organlist file are: {organlist_df_headers}')
    else:
        pass
    

    # Read the organlist file using pandas_fixed width file
    # The `organlist.dat` file is assumed to be space-separated
    try:
        organlist_df = pd.read_fwf(organlist_file, skiprows=organlist_skip_rows,
                                   names=organlist_df_headers)
    except FileNotFoundError:
        print("Error: Organlist file provided not found or is empty.")
        return
    
    
    # Initialize numpy arrays for marterial_ID, organ_ID and density
    arr_material = np.zeros(n_vox_tot, dtype=np.int16)
    arr_density = np.zeros(n_vox_tot, dtype=np.float32)
    arr_organ = np.zeros(n_vox_tot, dtype=np.int16)
    
    
    # Read the phantom data based on file type (bin or ASCII)
    print(f"\nLoading {phantom_file}...\n")

    if file_type == 0:  # Binary file
        # Read the entire binary file and interpret it as a stream of single bytes (characters)
        # Then convert each byte to its ASCII integer value
        try:
            with open(phantom_file, 'rb') as f:
                phantom_binary_data = f.read()
                arr_organ = np.array([struct.unpack('B', byte)[0] for byte in phantom_binary_data], dtype=np.int32)
        except FileNotFoundError:
            print(f"Error: Binary file '{phantom_file}' not found.")
            return

    elif file_type == 1:  # ASCII file
        # Cannot use numpy loadtxt because the number of columns(i.e. numbers) per line is inconsistent.
        # Read all lines as a string, split the string using the spaces
        try:
            with open(phantom_file, 'rb') as f:
                ph_file_str = f.read()
                arr_organ = np.array(ph_file_str.split(), dtype=np.int16)
        except FileNotFoundError:
            print(f"Error: ASCII file '{phantom_file}' not found.")
            return

    # Check Data quality, e.g. if the number of voxels read matches the expected total
    if len(arr_organ) != n_vox_tot:
        raise ValueError(f"Number of voxels read ({len(arr_organ)}) does not match expected total ({n_vox_tot}).\
                         Check the phantom characteristics.")
    
    # Check for negative densities
    if np.any(arr_density < 0):
        print('There are negative densities in the organlist file.\n')
        print('For the software to function properly, this cannot happen.\n')
        print('Each density will be chaged to its absolute value.\n\n')
        arr_density = np.abs(arr_density)
                        

    print(f"Finished loading the phantom file. {n_vox_tot} voxels were read.\n")
    
    
    # Start processing the data for loading to other files
    print('Processing the phantom data...\n')

    print('In the organlist file, what are the names of the columns that correspond to the Organ ID, Material ID and Density?\
These are the same names as provided before. Write the names in order, separated by commas (case sensitive).\
Default is "Organ_ID","Material_ID", "Density".')
    column_headers = input('> ').strip().split(',')
    # Check if the list contains at least one empty string, or if the original input was empty
    # Split the input
    column_headers = [header.strip() for header in column_headers]
    if not column_headers or any(not header for header in column_headers):
        column_headers = ['Organ_ID', 'Material_ID', 'Density']
        print(f'The names of the columns are {column_headers}.\n')
    else:
        pass
    
    print('\nWhat is the name of the .vox phantom file you want to create? (default is "phantom.vox")')
    vox_file = input('> ').strip()
    if not vox_file:
        vox_file='phantom.vox' 
        print('The name of the .vox phantom file is "phantom.vox" \n')
    
    print('\nWhat is the name of the ct-den-matXY, XZ and YZ visualization files you want to create?\
Write the names in order, separated by commas (case sensitive). (default is "ct-den-matXY.dat", XZ and YZ)')
    ct_den_mat_files = input('> ').strip().split(',')
    ct_den_mat_files = [header.strip() for header in ct_den_mat_files]
    if not ct_den_mat_files or any(not header for header in ct_den_mat_files):
        ct_den_mat_files=['ct-den-matXY.dat','ct-den-matXZ.dat','ct-den-matYZ.dat']
        print(f'The name of the ct-den-mat files is {ct_den_mat_files}.\n')
    else:
        pass
    
    
    # Create arrays equivalent to arr_organ, arr_material and arr_density, with the material ID and density values
    for row in organlist_df.itertuples(index=False):
        organ_indices = np.where(arr_organ == getattr(row, column_headers[0]))
        arr_material[organ_indices] = getattr(row, column_headers[1])
        arr_density[organ_indices] = getattr(row, column_headers[2])

    #Calculate and report max values found in phantom file.
    max_material = np.max(arr_material)
    max_density = np.max(arr_density)
    max_organ = np.max(arr_organ)

    print(f'  > Number of materials: {material_num}')
    print(f'  > Maximum value of material ID: {max_material}')
    print(f'  > Maximum value of organ (tag) ID: {max_organ}')
    print(f'  > Maximum density: {max_density:.5e}\n')
    print('  The phantom data was processed.\n')
   
    
    # Create .VOX file
    create_vox_file (vox_file, arr_material, arr_density, n_vox_x, n_vox_y, n_vox_z, 
                       vox_res_x, vox_res_y, vox_res_z)
        
    
    # Create ct-den-mat.dat files
    print('Creating the ct-den-mat.dat files...\n')
    
    # Create ct-den-matXY.dat
    order_indices_xy = lambda z, y, x: (x, y, z)
    loop_order_xy = (n_vox_z, n_vox_y, n_vox_x)
    create_ct_den_mat_file(ct_den_mat_files[0], arr_density, arr_material, arr_organ, n_vox_x, n_vox_y, n_vox_z, 
                           len_x, len_y, len_z, order_indices_xy, loop_order_xy)
    
    # Create ct-den-matXZ.dat
    order_indices_xz = lambda y, x, z: (x, y, z)
    loop_order_xz = (n_vox_y, n_vox_x, n_vox_z)
    create_ct_den_mat_file(ct_den_mat_files[1], arr_density, arr_material, arr_organ, n_vox_x, n_vox_y, n_vox_z, 
                           len_x, len_y, len_z, order_indices_xz, loop_order_xz)
    
    # Create ct-den-matYZ.dat
    order_indices_yz = lambda x, z, y: (x, y, z)
    loop_order_yz = (n_vox_x, n_vox_z, n_vox_y)
    create_ct_den_mat_file(ct_den_mat_files[2], arr_density, arr_material, arr_organ, n_vox_x, n_vox_y, n_vox_z, 
                           len_x, len_y, len_z, order_indices_yz, loop_order_yz)

    print(' All Files for visualization and simulation with this phantom have been created.\n')


if __name__ == '__main__':
    readPhantom()
    
    
    
