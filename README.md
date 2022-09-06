# DMFTtools
This is a collection of fortran modules and procedures to support quantum many-body, specifically Dynamical Mean-Field THeory, calculations. 

There are large portions of useful software missing. Anyone is welcome to contribute or to test the software. 


### Dependencies

* gfortran > 5.0 **OR** ifort  > 13.0
* cmake > 3.0.0    
* scifor  ( https://github.com/QcmPlab/SciFortran )   
* MPI ( https://github.com/open-mpi/ompi )  [optional, recommended]

The `scifor` library must be installed and loaded in the system. Info are contained in the linked repository. 


## BUILD

Installation is available using CMake. 
Clone the repo:

`git clone https://github.com/QcmPlab/DMFTtools dmft_tools`

### Make
From the repository directory (`cd dmft_tools`) make a standard out-of-source CMake-Make compilation:

`mkdir build`  
`cd build`  
`cmake ..`      (*)  
`make`     
`make install`   


### Ninja

If `ninja` ( https://ninja-build.org ) be available in your system, you can use it to build and install the library. 

From the repository directory (`cd dmft_tools`) make a standard out-of-source CMake-Ninja compilation:

`mkdir build`  
`cd build`  
`cmake -GNinja ..`      (*)  
`ninja`     
`ninja install`   


(*) *In some cases CMake fails to find the MPI Fortran compiler, even if it is effectively installed and loaded into the system. An easy fix is to setup and export the `FC=mpif90` environment variable before invoking `cmake`.* 


## INSTALL

Installation is completed after the build step using either: 

`make install`  

or   

`ninja install`  
 
To complete the installation you should follow the instructions printed on the screen, which suggest different ways to load the library in your system.


The `CMake` compilation can be controlled using the following additional variables, default values between `< >`:   

* `-DPREFIX=prefix directory <~/opt/scifor>` 

* `-DUSE_MPI=<yes>/no`  

* `-DVERBOSE=yes/<no> `  

* `-DBUILD_TYPE=<RELEASE>/TESTING/DEBUG`  


## UNINSTALL

`Cmake` does not officially provide uninstall procedure in the generate Makefile. Yet, it is possible to construct one such uninstall mode in a simple way. SCIFOR provides a way to uninstall the files generated inside any out-of-source compilation calling: 
`make uninstall`  


### CONTACT

For any information contact the author as:  
adriano DOT amaricci @ gmail DOT com

--

***LICENSE***  
Copyright (C) Adriano Amaricci, Lorenzo Crippa, Giacomo Mazza, Gabriele Bellomia, Samuele Giuli, Massimo Capone

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL) as published by
the Free Software Foundation, either version 3 of the License, or any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU LGPL for more details.

You should have received a copy of the GNU LGPL along with this program.  If not, see <http://www.gnu.org/licenses/>.


