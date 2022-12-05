# DMFTtools  
This is a collection of fortran modules and routines to support quantum many-body calculations, with a strong focus on Dynamical Mean-Field Theory.

There are still many useful features missing. Anyone is welcome to contribute or to test the software. 


### Dependencies

* [GNU Fortran (`gfortran`)](https://gcc.gnu.org/fortran/) > 5.0 **OR** [Intel Fortran Compiler Classic (`ifort`)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html)  > 13.0
* [CMake](https://cmake.org/) ≥ 3.0 [> 3.16 for ninja support] 
* [Make](https://www.gnu.org/software/make/) **OR** [Ninja](https://ninja-build.org/) ≥ 1.10 
* [SciFortran](https://github.com/QcmPlab/SciFortran), our scientific fortran library
* [MPI](https://github.com/open-mpi/ompi)  [optional, recommended]
* [ScaLAPACK](https://github.com/Reference-ScaLAPACK/scalapack)  [optional, recommended] 

If any of the required libraries is not available in your system, or the version requirement is not satisfied, please install/upgrade them. Except for SciFortran (for which you should follow the installation instructions reported in the linked README), we generally advice for pre-packaged versions, as provided by either [`apt`](https://en.wikipedia.org/wiki/APT_(software)), [`pip`](https://pypi.org/project/pip/), [`brew`](https://formulae.brew.sh/), [`conda`](https://docs.conda.io/en/latest/) or [`spack`](https://spack.io/). The latter may provide the best options for HPC environments (trivial installation without sudo, easy management of interdependent multiple versions, automatic generation of environment modules, etc.), but choose freely according to your needs.

## BUILD

Our build system relies on CMake. Clone the repo:

`git clone https://github.com/QcmPlab/DMFTtools dmft_tools`

<details>
<summary> Using <tt>make</tt> (click to expand) </summary>
Default CMake workflow, with widest version support (CMake > 3.0).

```
mkdir build 
cd build  
cmake .. 
make
```      

</details>

<details>
<summary> Using <tt>ninja</tt> (click to expand)</summary>

If a fortran-capable[^3] version of `ninja` ( https://ninja-build.org ) is available in your system (and CMake can[^4] take advantage of it), you can use it to build the library at lightning, multi-threaded, speed. 

```
mkdir build    
cd build  
cmake -GNinja ..  
ninja
```       

</details>

[^3]: Ninja did not support fortran before version 1.10, although Kitware has long mantained a fortran-capable fork, which might be obtained easily as a [Spack package](https://packages.spack.io/package.html?name=ninja-fortran). Nevertheless we note that as of fall 2022 `pip install ninja --user` [ships Ninja v1.10.2](https://pypi.org/project/ninja/), hence obtaining a suitable official Ninja release should be trivial.

[^4]: This depends on your CMake version. Comparing [this](https://cmake.org/cmake/help/v3.16/generator/Ninja.html#fortran-support) to [this](https://cmake.org/cmake/help/v3.17/generator/Ninja.html#fortran-support) would suggest that CMake started supporting Ninja's fortran features only after v3.17 but we have verified that at least v3.16.3 (current version shipped by `apt` on Ubuntu 20.04 LTS) does indeed work. For more information you can take a look to a [related SciFortran issue](https://github.com/QcmPlab/SciFortran/issues/16). 

## INSTALL

System-wide installation is completed after the build step using either: 

```
make install
```  

or   

```
ninja install
```  
 
To actually link the library we provide some alternatives: 

* A generated [pkg-config](https://github.com/freedesktop/pkg-config) file to, installed to `~/.pkgconfig.d/dmft_tools.pc`  
* A generated [environment module](https://github.com/cea-hpc/modules), installed to `~/.modules.d/dmft_tools/<PLAT>/<VERSION>`  
* Two generated `bash` scripts at `<PREFIX>/bin/`, to be sourced for permanent loading (user or global).

which you can choose among by following the instructions printed on screen.

## UNINSTALL

CMake does not officially provide uninstall procedures in the generated Make/Ninja files. Hence DMFTtools supplies a homebrew method to remove the generated files by calling (from the relevant build folder): `make uninstall` / `ninja uninstall`. 


### CONTACT

If you encounter bugs or difficulties, please [file an issue](https://github.com/QcmPlab/DMFTtools/issues/new/choose). For any other communication, please reach out to:    
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



