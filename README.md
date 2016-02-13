# libbsdf
## Overview
libbsdf is a C++ library for BSDF, BRDF, and BTDF.
libbsdf provides basic functions: file readers, interpolation, data structure and so on.
[Eigen 3][1] is used for arrays, vectors, and matrices.
libbsdf is licensed under the terms of the Mozilla Public License, version 2.0.
See the LICENSE file.

## Make it build on OSX 10.10 with MacPorts
* ports to install:
** eigen3
** qt4-mac
** OpenSceneGraph (this needs an upgrade for the current version disables Qt support)
* ``cmake -DEIGEN3_INCLUDE_DIR=/opt/local/include/eigen3 -DCMAKE_CXX_FLAGS=-std=c++11 ..``


## Future Plans
* Adds spline or rational interpolation
* Adds virtual constructors to create a copy of lb::Bsdf, lb::Brdf, or lb::Btdf
* Adds FindLibbsdf.cmake
* Prepares sample codes

[1]: http://eigen.tuxfamily.org/index.php?title=Main_Page "Eigen"
