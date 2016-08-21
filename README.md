# libbsdf
## Overview
libbsdf is a C++ library for BSDF, BRDF, and BTDF.
libbsdf provides basic functions: file readers, interpolation, data structure and so on.
[Eigen 3][1] is used for arrays, vectors, and matrices.
libbsdf is licensed under the terms of the Mozilla Public License, version 2.0.
See the LICENSE file.

## Make it build on OSX 10.10 with MacPorts
* ports to install:
  * eigen3
* ``cmake -DEIGEN3_INCLUDE_DIR=/opt/local/include/eigen3 -DCMAKE_CXX_FLAGS=-std=c++11 ..``

## Install it on Windows using MSYS2
* install MSYS2
* ``pacman -S mingw-w64-(i686|x86_64)-libbsdf`` (choose correct architecture)
* Be aware that this is work in progress, as it is not yet as performant as it ought to be!

## Future Plans
* Adds FindLibbsdf.cmake
* Prepares sample codes

[1]: http://eigen.tuxfamily.org/index.php?title=Main_Page "Eigen"
