# ReadPhantom: Voxel phantom data processing for PENELOPE simulations

A Python script to read, process, and convert voxel phantom data into files suitable for simulation and visualization with the PENELOPE/penEasy framework. :atom_symbol: ⚛️

The `ReadPhantom.py` script is a powerful tool designed to streamline the workflow for researchers and students working with voxel phantoms. It automates the conversion of raw phantom data—which can be either in ASCII or binary format—into structured output files required by PENELOPE/penEasy. This program is a modern Python re-implementation of a similar tool originally written in Fortran.

The script handles two primary input files:

- A phantom file containing the organ ID for each voxel, ordered continuously along x, y, and z axes.

- An organlist file that provides the corresponding Material ID and density for each Organ ID.

The script then produces two main types of output files:

- A .vox file, which is the voxel phantom format read directly by PENELOPE.

- ct-den-mat files (.dat), which are formatted for visualization using GNUPLOT scripts that come with the PENELOPE and penEasy distributions.

The program is optimized for efficiency, especially when handling large datasets with millions of voxels. It uses Python's numpy and pandas libraries for fast data processing, avoiding slow, line-by-line processing.

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

__Key Prompts & Default Values:__

- __Input File Type:__ 0 for binary or 1 for ASCII.

- __Voxel Dimensions:__ The number of voxels in x, y, and z.

- __Voxel Resolution:__ The physical size of each voxel in x, y, and z, in centimeters.

- __Organlist File Details:__ You'll be asked for the name of your organlist.dat file, the number of header rows to skip, and the column names.

- __Output File Names:__ You can specify custom names for the output files or accept the defaults:

- __VOX file:__ Defaults to phantom.vox.

- __GNUPLOT files:__ Defaults to ct-den-matXY.dat, ct-den-matXZ.dat, and ct-den-matYZ.dat.

If you just press Enter on a prompt with a default value, the script will automatically use that default.



