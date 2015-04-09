# libbsdf
## Overview
libbsdf is a C++ library for BSDF, BRDF, and BTDF.
libbsdf provides basic functions: file readers, interpolation, data structure and so on.
[Eigen 3][1] is used for arrays, vectors, and matrices.
libbsdf is licensed under the terms of the Mozilla Public License, version 2.0.
See the LICENSE file.

## Future Plans
* Adds spline or rational interpolation
* Adds virtual constructors to create a copy of lb::Bsdf, lb::Brdf, or lb::Btdf
* Adds FindLibbsdf.cmake
* Prepares sample codes

[1]: http://eigen.tuxfamily.org/index.php?title=Main_Page "Eigen"
