# Equivalent Square Calculator


<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li><a href="#description">Description</a></li>
    <li><a href="#installation">Installation</a></li>
    <ul>
      <li><a href="#python">Python</a></li>
      <li><a href="#windows">Windows</a></li>
      <li><a href="#python">Linux</a></li>
      <li><a href="#osx">OSX</a></li>
    </ul>
    <li><a href="#usage">Usage</a></li>
    <ul>
      <li><a href="#python-library">Python Library</a></li>
      <li><a href="#gui">GUI</a></li>
    </ul>
    <li><a href="#references">References</a></li>
  </ol>
</details>


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
The standalone version of the GUI can be installed downloading the newiest installation file <a href='https://github.com/ringholz-j/equivalent-square-calculator/raw/main/install/windows/EquivalentSquareCalculator-1.0.0-amd64.msi'>`EquivalentSquareCalculator-1.0.0-amd64.msi`</a>.
#### Requirements
* Tested for Windows 10
* disc space: 300MB
### Linux
The standalone version of the GUI can be installed downloading and unpacking the newest release ZIP-file <a href='https://github.com/ringholz-j/equivalent-square-calculator/raw/main/install/linux/exe.linux-x86_64-3.6.zip'>`exe.linux-x86_64-3.6.zip`</a> in your desired location. After that open the file `EquivalentSquareCalculator.desktop` with an arbitrary text editor and and insert in the following row the path of the `EquivalentSquareCalculator` file

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
The standalone version of the GUI can be installed unpacking the ZIP-file <a href='https://github.com/ringholz-j/equivalent-square-calculator/raw/main/install/osx/EquivalentSquareCalculator.zip'>`EquivalentSquareCalculator.zip`</a> and move file `EquivalentSquareCalculator.app`to your application folder. Now the Equivalent Square Calculator should appear in your launchpad.

#### Requirements
* Tested for: OSX 10.15.7 Catalina
* disc space: 160MB

## Usage

### Python library
The first part of the library contains the functions according to Ahnesjö et al. [[1]](#1) and Nyholm et al. [[2]](#2),[[3]](#3) describing the pencil beam kernel. The following functions (`Pz_tri_part_b` ... `Dz_rect_tpr`) integrating the kernel over the considered rectangular field dimensions and calculating the required axis dose of the fields. The functions `newton_f` ... `newton` define the newton method to find the square field with the best dose agreement to the considered rectangular field.<br><br>

#### class EquivalentSquare:
Class for assigning WFF square fields to WFF rectangular fields using the equal axis dose definition.

##### Input variables:
* `dx`: x-dimension of the considered rectangular field in cm (required)
* `dy`: y-dimension of the considered rectangular field in cm (required)
* `z:` depth in water in cm; default z=10cm (optional)
* `tpr2010`: quality index TPR<sub>20,10</sub> of used LINAC; default tpr2010=671 (optional)
* `epsilon`: maximal deviation in integrated axis dose (cm<sup>2</sup> g<sup>-1</sup>)) between rectangular and square field; default epsilon=0.000001 (optional)
* `max_iter`: maximal amount of iterations for the newton optimation; default max_iter=100 (optional)
* `no_kernel`: if false also the integrated kernal dose of the square field is saved.

##### calculated output class variables:
* `equi_sq`: calculated equivalent square round to 0.01cm
* `equi_sq_raw`: unround calculated equivalent square
* `geometric_mean`: Geometric Mean
* `geo_dif`: difference between equivalent square and geometric mean
* `geo_dif_rel`: relative difference between equivalent square and geometric mean
* `sterling`: Output from Sterling equation
* `sterling_dif`: difference between equivalent square and output from Sterling equation
* `sterling_dif_rel`: relative difference between equivalent square and output from Sterling equation


#### class EquivalentSquareFFF:
Class for assigning FFF square fields to FFF rectangular fields using the equal axis dose definition.

##### Input variables:
* `dx`: x-dimension of the considered rectangular field in cm (required)
* `dy`: y-dimension of the considered rectangular field in cm (required)
* `z:` depth in water in cm; default z=10cm (optional)
* `tpr2010`: quality index TPR<sub>20,10</sub> of used LINAC; default tpr2010=671 (optional)
* `epsilon`: maximal deviation in integrated axis dose (cm<sup>2</sup> g<sup>-1</sup>)) between rectangular and square field; default epsilon=0.000001 (optional)
* `max_iter`: maximal amount of iterations for the newton optimation; default max_iter=100 (optional)
* `no_kernel`: if false also the integrated kernal dose of the square field is saved.
* `energy`: Selection of the implemented FFF-profiles, used as weighting function in the calculation; obtained for ELEKTA Versa HD LINACS at a depth of 1cm; possibilities: ['6mv','10mv']; default energy='6mv' (optional)

##### calculated output class variables:
* `equi_sq`: calculated equivalent square round to 0.01cm
* `equi_sq_raw`: unround calculated equivalent square
* `geometric_mean`: Geometric Mean
* `geo_dif`: difference between equivalent square and geometric mean
* `geo_dif_rel`: relative difference between equivalent square and geometric mean
* `sterling`: Output from Sterling equation
* `sterling_dif`: difference between equivalent square and output from Sterling equation
* `sterling_dif_rel`: relative difference between equivalent square and output from Sterling equation


#### class EquivalentSquareFFFWFF:
Class for assigning WFF square fields to FFF rectangular fields using the equal axis dose definition.

##### Input variables:
* `dx`: x-dimension of the considered rectangular field in cm (required)
* `dy`: y-dimension of the considered rectangular field in cm (required)
* `z:` depth in water in cm; default z=10cm (optional)
* `tpr2010`: quality index TPR<sub>20,10</sub> of used LINAC; default tpr2010=671 (optional)
* `epsilon`: maximal deviation in integrated axis dose (cm<sup>2</sup> g<sup>-1</sup>)) between rectangular and square field; default epsilon=0.000001 (optional)
* `max_iter`: maximal amount of iterations for the newton optimation; default max_iter=100 (optional)
* `no_kernel`: if false also the integrated kernal dose of the square field is saved.
* `energy`: Selection of the implemented FFF-profiles, used as weighting function in the calculation; obtained for ELEKTA Versa HD LINACS at a depth of 1cm; possibilities: ['6mv','10mv']; default energy='6mv' (optional)

##### calculated output class variables:
* `equi_sq`: calculated equivalent square round to 0.01cm
* `equi_sq_raw`: unround calculated equivalent square
* `geometric_mean`: Geometric Mean
* `geo_dif`: difference between equivalent square and geometric mean
* `geo_dif_rel`: relative difference between equivalent square and geometric mean
* `sterling`: Output from Sterling equation
* `sterling_dif`: difference between equivalent square and output from Sterling equation
* `sterling_dif_rel`: relative difference between equivalent square and output from Sterling equation


#### class EquivalentSquareTPR:
Class for assigning WFF square fields to WFF rectangular fields using the equal TPR<sub>20,10</sub> definition.

##### Input variables:
* `dx`: x-dimension of the considered rectangular field in cm (required)
* `dy`: y-dimension of the considered rectangular field in cm (required)
* `tpr2010`: quality index TPR<sub>20,10</sub> of used LINAC; default tpr2010=671 (optional)
* `epsilon`: maximal deviation in TPR<sub>20,10</sub> between rectangular and square field; default epsilon=0.000001 (optional)
* `max_iter`: maximal amount of iterations for the newton optimation; default max_iter=100 (optional)

##### calculated output class variables:
* `equi_sq`: calculated equivalent square round to 0.01cm
* `equi_sq_raw`: unround calculated equivalent square
* `geometric_mean`: Geometric Mean
* `geo_dif`: difference between equivalent square and geometric mean
* `geo_dif_rel`: relative difference between equivalent square and geometric mean
* `sterling`: Output from Sterling equation
* `sterling_dif`: difference between equivalent square and output from Sterling equation
* `sterling_dif_rel`: relative difference between equivalent square and output from Sterling equation
* `tpr`: TPR<sub>20,10</sub> of rectangular field
* `tpr_sq`: TPR<sub>20,10</sub> of the rectangular field
* `tpr_dif`: TPR<sub>20,10</sub> difference between rectangular and square field 


#### class EquivalentSquareFFFTPR:
Class for assigning FFF square fields to FFF rectangular fields using the equal TPR<sub>20,10</sub> definition.

##### Input variables:
* `dx`: x-dimension of the considered rectangular field in cm (required)
* `dy`: y-dimension of the considered rectangular field in cm (required)
* `tpr2010`: quality index TPR<sub>20,10</sub> of used LINAC; default tpr2010=671 (optional)
* `epsilon`: maximal deviation in TPR<sub>20,10</sub> between rectangular and square field; default epsilon=0.000001 (optional)
* `max_iter`: maximal amount of iterations for the newton optimation; default max_iter=100 (optional)
* `energy`: Selection of the implemented FFF-profiles, used as weighting function in the calculation; obtained for ELEKTA Versa HD LINACS at a depth of 1cm; possibilities: ['6mv','10mv']; default energy='6mv' (optional)

##### calculated output class variables:
* `equi_sq`: calculated equivalent square round to 0.01cm
* `equi_sq_raw`: unround calculated equivalent square
* `geometric_mean`: Geometric Mean
* `geo_dif`: difference between equivalent square and geometric mean
* `geo_dif_rel`: relative difference between equivalent square and geometric mean
* `sterling`: Output from Sterling equation
* `sterling_dif`: difference between equivalent square and output from Sterling equation
* `sterling_dif_rel`: relative difference between equivalent square and output from Sterling equation
* `tpr`: TPR<sub>20,10</sub> of rectangular field
* `tpr_sq`: TPR<sub>20,10</sub> of the rectangular field
* `tpr_dif`: TPR<sub>20,10</sub> difference between rectangular and square field


#### class EquivalentSquareFFFWFFTPR:
Class for assigning FFF square fields to FFF rectangular fields using the equal TPR<sub>20,10</sub> definition.

##### Input variables:
* `dx`: x-dimension of the considered rectangular field in cm (required)
* `dy`: y-dimension of the considered rectangular field in cm (required)
* `tpr2010`: quality index TPR<sub>20,10</sub> of used LINAC; default tpr2010=671 (optional)
* `epsilon`: maximal deviation in TPR<sub>20,10</sub> between rectangular and square field; default epsilon=0.000001 (optional)
* `max_iter`: maximal amount of iterations for the newton optimation; default max_iter=100 (optional)
* `energy`: Selection of the implemented FFF-profiles, used as weighting function in the calculation; obtained for ELEKTA Versa HD LINACS at a depth of 1cm; possibilities: ['6mv','10mv']; default energy='6mv' (optional)

##### calculated output class variables:
* `equi_sq`: calculated equivalent square round to 0.01cm
* `equi_sq_raw`: unround calculated equivalent square
* `geometric_mean`: Geometric Mean
* `geo_dif`: difference between equivalent square and geometric mean
* `geo_dif_rel`: relative difference between equivalent square and geometric mean
* `sterling`: Output from Sterling equation
* `sterling_dif`: difference between equivalent square and output from Sterling equation
* `sterling_dif_rel`: relative difference between equivalent square and output from Sterling equation
* `tpr`: TPR<sub>20,10</sub> of rectangular field
* `tpr_sq`: TPR<sub>20,10</sub> of the rectangular field
* `tpr_dif`: TPR<sub>20,10</sub> difference between rectangular and square field 


### GUI
#### Calculator Tab
Gui for calculating equivalent squares for both definitions equal axis dose and equal TPR<sub>20,10</sub>. For comparison also the output for geometric mean and the Sterling equation is calculated.

##### Input parameters
* x: x-dimension of the considered rectangular field
* y: y-dimension of the considered rectangular field
* depth: depth in water in cm
* TPR2010: the quality index TPR<sub>20,10</sub> of the used LINAC
* mode: assignment mode between square and rectangular fields; options: [WFF-WFF,FFF-FFF,FFF,WFF]
* FFF-profile: Only relevant for modes including FFF-fields; selection of the implemented FFF-profiles, used as weighting function in the calculation; obtained for ELEKTA Versa HD LINACS at a depth of 1cm; options: [6mv,10mv]


#### Differences Plots Tab
In this subprogram the differences between our equivalent square calculation method and the comon geometric mean respectively Sterling equation can be plotted as color map. Optional the plot can also be done using an one colored limed plot. The plots can be saved as png-file.

##### Input parameters
###### Required
* Smin: minimal fieldsize in cm
* Smax: maximal field size in cm
* increment: increment between the plotted pixel in each direction in cm
* depth: depth in water in cm
* TPR2010: the quality index TPR<sub>20,10</sub> of the used LINAC
* method: selection of the equivalent square definition; options: [Axis dose, TPR2010]
* mode: assignment mode between square and rectangular fields; options: [WFF-WFF,FFF-FFF,FFF,WFF]
* FFF-profile: Only relevant for modes including FFF-fields; selection of the implemented FFF-profiles, used as weighting function in the calculation; obtained for ELEKTA Versa HD LINACS at a depth of 1cm; options: [6mv,10mv]

###### Optional
* representation: 
  * relative: if true the ploted differences are relative
  * absolute values: if true the plotted differences are absolute values (positive)
* color map: selection of color map; all comon `matplotlib`color maps can be typed in
  * min color: if true the minimal value for the 'minimal' color of the colorbar can be typed in; false the smalest value of the array is choosen
  * max color: if true the maximal value for the 'maximal' calor of the colorbar can be typed in; false the largest value of the array is choosen
* Limit Plot:
  * show: if true a limit plot is shown instead of the normal one; this means all values biger or equal to the limit are shown in red
  * limit: the limit for the limit plot can be set

#### Create Tables Tab
A table for equivalent squares can be calculated and shown. Afater the calculation the table can be saved using the excel format.

##### Input parameters
* Smin: minimal fieldsize in cm
* Smax: maximal field size in cm
* increment: increment between the plotted pixel in each direction in cm
* depth: depth in water in cm
* TPR2010: the quality index TPR<sub>20,10</sub> of the used LINAC
* method: selection of the equivalent square definition; options: [Axis dose, TPR2010]
* mode: assignment mode between square and rectangular fields; options: [WFF-WFF,FFF-FFF,FFF,WFF]
* FFF-profile: Only relevant for modes including FFF-fields; selection of the implemented FFF-profiles, used as weighting function in the calculation; obtained for ELEKTA Versa HD LINACS at a depth of 1cm; options: [6mv,10mv]


## References
<a id="1">[1]</a> 
A. Ahnesjö, M. Saxner, and A. Trepp, A Pencil Beam Model for Photon Dose Calculation, Medical Physics 19, 263–273 (1992).

<a id="2">[2]</a> 
T. Nyholm, J. Olofsson, A. Ahnesjo ̈, and M. Karlsson, Photon Pencil Kernel Param- eterisation Based on Beam Quality Index, Radiotherapy and Oncology 78, 347–351 (2006).

<a id="3">[3]</a> 
T. Nyholm, J. Olofsson, A. Ahnesj ̈o, and M. Karlsson, Corrigendum to “Photon Pencil Kernel Parameterisation Based on Beam Quality Index” [Radiother. Oncol. 78 (2006) 347–351], Radiotherapy and Oncology 98, 286 (2011).
