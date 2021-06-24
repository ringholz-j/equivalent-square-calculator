# Equivalent Square Calculator

## Description
In the field of radiation therapy and dosimetry the correction factor of specific radiation field detectors are in most cases determined for square fields only. So to retrieve these factors also for different formed fields such as rectangular fields, the equivalent square fields for these fields have to be known. Because for rectangular fields the pure geometric relation such as geometric mean and Sterling equation are still in use today, we developed a physical method based on the pencil beam model according to Ahnesjö et al. [[1]](#1) and Nyholm et al. [[2]](#2),[[3]](#3). We implement this method as Python library and GUI.
We implement both definitions for equivalent square: equal axis dose and equal TPR<sub>20,10</sub>.<br>
The assignments between WFF (with flattening filter) square fields and WFF rectangular fields, between FFF (flattening filter free) square fields and FFF rectangular fields and between WFF square fields and FFF rectangular fields are possible.<br>
The grafical user interface (GUI) includes a calculater for the equivalent square, a difference plotter for the comon geometric mean and Sterling equation, and a table generator for equivalent squares with an export function in excel-format.

## Installation
The package is provided as pure Python Version for both the library and the grafical user interface, and also as standalone GUI version for the operating systems Windows, Linux and Mac OSX.

### Python
The file `equivalent_square_lib.py` includes the library and `equivalent_square_gui.py` the grafical user interface. A Python version of minimum 3.6.9 is required. 

#### Required Python packages
* `decimal`
* `matplotlib`
* `numpy`
* `pandas`
* `pandastable`
* `scipy`
* `tkinter`


### Windows
The standalone version of the GUI can be installed using the installation file `EquivalentSquareCalculator-1.0.0-amd64.msi`.
#### Requirements
* Tested for Windows 10
* disc space: 300MB
### Linux
The standalone version of the GUI can be installed unpacking the ZIP-file `exe.linux-x86_64-3.6.zip` in your desired location. After that open the file `EquivalentSquareCalculator.desktop` with an arbitrary text editor and and insert in the following row the path of the `EquivalentSquareCalculator` file

  ```sh
  Exec=/YOURPATH/EquivalentSquareCalculator
  ```
If you move the file `EquivalentSquareCalculator.desktop` to the application folder:

  ```sh
  /home/YOURUSERNAME/.local/share/applications/
  ```
you find the EquivalentSquareCalculator in your launchpad.

#### Requirements
* Tested for: Ubuntu 18
* disc space: 270MB

### OSX
The standalone version of the GUI can be installed unpacking the ZIP-file `exe.macosx-10.9-x86_64-3.9.zip` and move file `EquivalentSquareCalculator.app`to your application folder. Now the Equivalent Square Calculator should appear in your launchpad.

#### Requirements
* Tested for: OSX 10.15.7 Catalina
* disc space: 370MB

## Usage

### Python library
The first part of the library contains the functions according to Ahnesjö et al. [[1]](#1) and Nyholm et al. [[2]](#2),[[3]](#3) describing the pencil beam kernel. The following functions (`Pz_tri_part_b` ... `Dz_rect_tpr`) integrating the kernel over the considered rectangular field dimensions and calculating the required axis dose of the fields. The functions `newton_f` ... `newton` define the newton method to find the square field with the best dose agreement to the considered rectangular field.<br><br>


### GUI


## References
<a id="1">[1]</a> 
A. Ahnesjö, M. Saxner, and A. Trepp, A Pencil Beam Model for Photon Dose Calculation, Medical Physics 19, 263–273 (1992).

<a id="2">[2]</a> 
T. Nyholm, J. Olofsson, A. Ahnesjo ̈, and M. Karlsson, Photon Pencil Kernel Param- eterisation Based on Beam Quality Index, Radiotherapy and Oncology 78, 347–351 (2006).

<a id="3">[3]</a> 
T. Nyholm, J. Olofsson, A. Ahnesj ̈o, and M. Karlsson, Corrigendum to “Photon Pencil Kernel Parameterisation Based on Beam Quality Index” [Radiother. Oncol. 78 (2006) 347–351], Radiotherapy and Oncology 98, 286 (2011).
