# ReadPhantom: Voxel phantom data processing for PENELOPE simulations

readPhantom is a Python script to read, process, and convert voxel phantom data into files suitable for simulation and visualization with the PENELOPE/penEasy Monte Carlo framework [1,2]. :atom_symbol:

The `readPhantom.py script` is a powerful tool designed to streamline the workflow for researchers and students working with voxel phantoms. It automates the conversion of raw phantom data into structured .vox files required by PENELOPE/penEasy [3,4].

The program is optimized for efficiency, especially when handling large datasets with millions of voxels. It uses Python's numpy and pandas libraries for fast data processing, avoiding slow, line-by-line processing.

## Contextualization
Voxel phantoms are models of the human anatomy consisting of several voxels (i.e. 3D pixels) grouped together to model the human body, assembled from Computed Tomography or Magnetic Resonance Imaging. Monte Carlo (MC) simulations are numerical computations that use random sampling to obtain numerical results, being a common tool in methematics, economics, physics, finance, etc. When coupled with voxel phantoms, MC simulations are widely employed to model radiation with biologic matter and the human body [1, 2]. 

Nonetheless, on the one hand voxel phantom files are almost exclusively available as raw data. On the other hand, different MC softwares only allow voxel phantom loading in specific formats. Therefore, when working with a MC software, you very probably will have to create your own voxel phantom file.

If you want to learn more about voxel phantoms, MC simulations and the origin/purpose of the script, please follow the references.

__Note:__ This program is a modern Python re-implementation of a similar tool originally written in Fortran in 2017, also available in this repository. The ReadMe for the Fortran version is in Annex 1 of Reference [1].

:arrow_right: __References:__

[1] Borbinha J. Organ Dose Estimates in Thorax CT: Voxel Phantom Organ Matching With Individual Patient Anatomy. MSc Thesis, 2017. Available from: [http://hdl.handle.net/10362/](https://run.unl.pt/bitstream/10362/29982/1/Borbinha_2017.pdf) 

[2] Borbinha J, Di Maria S, Madeira P, Belchior A, Baptista M, Vaz P. Increasing organ dose accuracy through voxel phantom organ matching with individual patient anatomy. Radiat Phys Chem. 2019 Jun;159:35–46. doi: [10.1016/j.radphyschem.2019.02.014](https://www.doi.org/10.1016/j.radphyschem.2019.02.014)

[3] FNEA (2019), PENELOPE 2018: A code system for Monte Carlo simulation of electron and photon transport: Workshop Proceedings, Barcelona, Spain, 28 January – 1 February 2019, OECD Publishing, Paris, doi: [doi.org/10.1787/32da5043-en](https://doi.org/10.1787/32da5043-en) 

[4] Sempau J., Badal A., Brualla L. A PENELOPE-based system for the automated Monte Carlo simulation of clinacs and voxelized geometries-application to far-from-axis fields. Med Phys. 2011;38(11):5887–5895. doi:[10.1118/1.3643029](https://doi.org/10.1118/1.3643029)

## Description

### Input files:

- __Phantom file:__ Contains organ IDs for each voxel, ordered continuously along x, y, and z axes.
  - Can be either in either ASCII or binary format.
  - May come in many formats: '.dat','.txt','.csv', etc.
  - In older phantom files, there may be a header present, that needs to be removed, for the script to work properly.

- __Organlist file:__ Lists at least Organ ID, Material ID and Density in columns, not necessarily in this order.
  - May also come in many formats: '.dat', '.xlsx','.txt','.csv', etc.
  - The organlist file can be created by the user or come with the phantom as a package.
  - The columns (i.e. at least Organ ID, Material ID and Density) should all be a fixed width, as in the example image.
An example of organlist file is in the next image. I recommend you define the Organ ID for air outside body as 0, the organlists included with phantom files usually don't include this Organ ID. As for the Material ID for air outside phantom, it would depend on your computational dosimetry objectives, but the general recommendation is Material ID (air outside phantom) != Material ID (air inside phantom).

![Simple example of organlist created by user](images/organlist_ubuntu.png)

### Output files:

- __.vox file:__ This is the voxel phantom format required for MC simulation and read directly by PENELOPE.
  - The default is 'phantom.vox', but the script prompts you to choose another name.
  - To avoid the mixing of several files from different phantoms, I recommend you rename the 'phantom.vox' file with indicative names.
  - This file has a 7 line header

Example header for International Commission on Radiation Protection (ICRP) adult female reference phantom:

```
[SECTION VOXELS HEADER V.2008_04_13]
299 137 348                                                             (number of voxels in x,y,z) 
0.1775 0.1775 0.484                                                     (voxel resolution in x,y,z /cm)  
1                                                                       (1st column for material)  
2                                                                       (2nd column for density)  
0                                                                       (no blank line after xy cycle)  
[END OF VXH SECTION] 
```

  - It’s  possible  to  check  if  the  file  was  correctly  created  by  checking  its  number  of  lines. 
The ICRP adult female reference phantom file should have 299 * 137 * 348 + 7 lines, corresponding to the total number of voxels in 
the phantom (299 * 137 * 348) plus the number of lines of header (7). 

- __ct-den-mat files (.dat):__ These  are  visualization  files,  which  allow  visualization  of  the  phantom  in  the  x,  y  and  z  planes.
  - The  ct_den_mat.dat  files  list,  for  each  voxel,  the  x,  y  and  z  indices  (column  1,  2  and  3), as well as the density, Material ID and organ ID  (columns  4, 5 and 6).  Using  gnuplot  scripts (command_line  driven  graphic  utility) provided in the PENELOPE/penEasy framework,  it’s  possible  to  visualize  the  phantom  slice  per  slice.  The  corresponding  gnuplot  files  are:  visualizeVoxelsDensityXY/XZ/YZ.gpl.  To  visualize  the  phantom:

```
gnuplot
>load!“visualizeVoxelsDensityXY.gpl”
```

## Prerequisites

To run this script, you'll need Python 3.x, as well as the numpy and pandas libraries. These are included with the anaconda distribution. You can install them with pip: 

```
pip install numpy pandas
```

## Usage

To get started, place the script and your input files in the same directory. Then, simply run the script from your terminal:

```
python readPhantom.py
```

The script will then guide you through a command-line interface, consisting of a series of interactive prompts, to gather all the necessary phantom information and file names.

## Configuration

The script is highly user-friendly and configurable through the command-line interface. Most parameters, like file names and phantom dimensions, are provided during runtime.

__Key Prompts   Default Values:__

- __Input File Type:__ 0 for binary or 1 for ASCII.

- __Voxel Dimensions:__ The number of voxels in x, y, and z.

- __Voxel Resolution:__ The physical size of each voxel in x, y, and z, in centimeters.

- __Organlist File Details:__ You'll be asked for the name of your organlist.dat file, the number of header rows to skip, and the column names.

- __Output File Names:__ You can specify custom names for the output files or accept the defaults:

- __VOX file:__ Defaults to phantom.vox.

- __GNUPLOT files:__ Defaults to ct-den-matXY.dat, ct-den-matXZ.dat, and ct-den-matYZ.dat.

If you just press Enter on a prompt with a default value, the script will automatically use that default.



