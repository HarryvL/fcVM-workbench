# fcVM-workbench
Finite element collapse analysis based on the von Mises plasticity model for use with [FreeCAD](https://freecad.org)

<img src="https://github.com/HarryvL/fcVM-workbench/blob/main/pictures/Embankment_with_Ditch_Example_Load_Displacement.png" height="250"/> <img src="https://github.com/HarryvL/fcVM/blob/main/pictures/Embankment_with_Ditch_Example_Failure_Mechanism.png" height="220" raw=true/> 

### Background
fcVM is a finite element workbench and solver for performing collapse analysis of structures and soil bodies. It is based on the theory of elastoplasticity and gives insight in ductility and reserve strength beyond first yield. The theory underpinning fcVM can be found in the document included in this repository.

<img src="https://github.com/HarryvL/fcVM-workbench/blob/main/pictures/Plate_with_hole_Example_Load_Displacement.png" height="250"/><img src="https://github.com/HarryvL/fcVM-workbench/blob/main/pictures/Plate_with_hole_Example_Failure_Mechanism.png" height="200"/>

### Prerequisites
* FreeCAD >= v0.19

### Installation
1. Create a fcVM-workbench directory in the FreeCAD Mod directory, the location of which can be found by typing App.getUserAppDataDir() in the FreeCAD python console.
1. Copy the entire fcVM-workbench repository into this new fcVM-workbench directory.
1. Restart FreeCAD
1. If all went well, the fcVM workbench should now be an option in the FreeCAD workbench dropdown menu

```bash
App.getUserAppDataDir()
├── Mod
│   ├── fcVM-workbench
│   │   ├── control files
│   │   │   ├── Block_Disp_Control_Example.inp
│   │   │   └── etc.
│   │   ├── documentation
│   │   │   └── Numerical Analysis of Soil Structure Interaction HvL.PDF
│   │   ├── freeCAD files
│   │   │   ├── Block_Disp_Control_Example.FCStd
│   │   │   └── etc.
│   │   ├── icons
│   │   │   └── fcFEM.svg
│   │   ├── output files
│   │   │   ├── Block_Disp_Control_Example.out
│   │   │   └── etc.
│   │   ├── pictures
│   │   │   ├── Embankment_with_Ditch_Example_Failure_Mechanism.png
│   │   │   └── etc.
│   │   ├── source code
│   │   │   ├── .gitignore
│   │   │   ├── fcVM.FCMacro
│   │   │   └── fcVM.py
│   │   ├── user_interface
│   │   │   └── fcVM.ui
│   │   ├── .gitignore
│   │   │
│   │   ├── Init.py
│   │   │
│   │   ├── InitGui.py
│   │   │
│   │   ├── README.md
│   │   │
│   │   └── dummy.py
│   │
```

<img src="https://github.com/HarryvL/fcVM-workbench/blob/main/pictures/Pit_Example_Load_Displacement.png" height="250"/><img src="https://github.com/HarryvL/fcVM-workbench/blob/main/pictures/Pit_Example_Failure_Mechanism.png" height="250"/>

### Dependencies
fcVM imports (from) the following packages: numpy, scipy (version 1.11.3), numba, matplotlib, scikit-sparse, pyvista and meshio. On Linux and Mac-OS systems these can be installed with the usual package managers (e.g. conda or mamba). On Windows this requires more effort:
1. download Miniforge3: [Miniforge3](https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe)
1. run the installer: Miniforge3-Windows-x86_64.exe
1. find and run Miniforge on your system - this opens a Miniforge Prompt: (base) C:\"path">
1. create a new virtual environment: (base) C:\"path"> mamba create --name fcVM (or any other name of your choice)
1. change into the new environment: (base) C:\"path"> mamba activate fcVM (or the other name you chose)
1. instal freeCAD and dependencies: (fcVM) C:\"path"> mamba install freecad scipy=1.11.3 numba matplotlib scikit-sparse pyvista meshio (with spaces and no commas)
1. check with python if the dependencies can be imported: (fcVM) C:\"path"> python
1. \>>> import scipy.sparse
1. \>>> import sksparse.cholmod
1. if no problems, quit python and start freecad: (fcVM) C:\"path"> freecad
1. if you encounter a "black screen" then follow the advice [here](https://forum.freecad.org/viewtopic.php?t=36087&start=40#p669458), i.e. rename all opengl32sw.dll files on your system to opengl32.dll

### Documentation
Please refer to source code for in-line comments and to the FreeCAD forum (TBC)

<img src="https://github.com/HarryvL/fcVM-workbench/blob/main/pictures/Tubes_Example_Load_Displacement.png" height="250"/><img src="https://github.com/HarryvL/fcVM-workbench/blob/main/pictures/Tubes_Example_Failure_Mechanism.png" height="200"/>

### Current Limitations
fcVM only works with 2nd order volumetric meshes generated with GMSH or Netgen (10-noded tetrahedrals). 

### Licence information

Copyright (c) 2024 - Harry van Langen <hvlanalysis@gmail.com>  


This program is free software; you can redistribute it and/or modify  
it under the terms of the GNU Lesser General Public License (LGPL)    
as published by the Free Software Foundation; either version 2 of     
the License, or (at your option) any later version.                   
for detail see the LICENCE text file.                                 
                                                                         
This program is distributed in the hope that it will be useful,       
but WITHOUT ANY WARRANTY; without even the implied warranty of        
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
GNU Library General Public License for more details.                  
                                                                         
You should have received a copy of the GNU Library General Public     
License along with this program; if not, write to the Free Software   
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  
USA                                                                   
