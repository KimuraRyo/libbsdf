# libbsdf

[![Linux](https://github.com/KimuraRyo/libbsdf/actions/workflows/linux.yml/badge.svg)](https://github.com/KimuraRyo/libbsdf/actions/workflows/linux.yml)

## Overview

libbsdf is a C++ library for BSDF, BRDF, and BTDF.
libbsdf provides basic functions: file IO, interpolation, data structure and so on.
[Eigen 3][Eigen] is used for arrays, vectors, and matrices.
libbsdf is licensed under the terms of the Mozilla Public License, version 2.0.
See the LICENSE file.

## Building libbsdf

C++ compiler that supports C++17 is required.
CMake is used as the build system.
The search path for [Eigen 3][Eigen] is set through CMake variables: `Eigen3_DIR`.
To enable the fitting function to reflectance models, enable LIBBSDF_ENABLE_FITTING and set the Ceres Solver path to the Ceres_DIR variable.

### Data structure

- Scatter and specular data
  - BRDF
  - BTDF
  - Specular reflectance
  - Specular transmittance
- Color
  - Monochrome
  - RGB
  - XYZ
  - Spectrum
- Parameterization
  - Spherical coordinate system
  - Specular coordinate system
  - Half-difference coordinate system

### File reader

| Format                                            | Extension  | Measured Data                      |
| ------------------------------------------------- | ---------- | ---------------------------------- |
| [Surface Scattering Distribution Data][Spec_SSDD] | .ssdd      |                                    |
| Integra Diffuse Distribution                      | .ddr, .ddt |                                    |
| Integra Specular Distribution                     | .sdr, .sdt |                                    |
| Zemax BSDF                                        | .bsdf      | [RPC Photonics][Data_RPC]          |
| LightTools BSDF                                   | .bsdf      |                                    |
| ASTM E1392-96(2002)                               | .astm      | [Cornell University][Data_Cornell] |
| MERL BRDF                                         | .binary    | [MERL][Data_MERL]                  |

### File writer

| Format                                            | Extension  |
| ------------------------------------------------- | ---------- |
| [Surface Scattering Distribution Data][Spec_SSDD] | .ssdd      |
| Integra Diffuse Distribution                      | .ddr, .ddt |

[Eigen]: http://eigen.tuxfamily.org/index.php?title=Main_Page "Eigen"
[Data_Cornell]: http://www.graphics.cornell.edu/online/measurements/reflectance/
[Data_MERL]: http://www.merl.com/brdf
[Data_RPC]: http://www.rpcphotonics.com/bsdf-data-optical-diffusers/
[Spec_SSDD]: https://github.com/KimuraRyo/libbsdf/blob/master/doc/SsddFileFormatSpecification.md
