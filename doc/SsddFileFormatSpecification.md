# Surface scattering distribution data (SSDD) file format

Version: 0.2
Revised: 2023-08-14

The SSDD file format is designed to store tabular data of measured or generated bidirectional scattering distribution function (BSDF). A bidirectional reflectance distribution function (BRDF), bidirectional transmittance distribution function (BTDF), specular reflectance, and specular transmittance can be contained in an SSDD file.

The extension of file is **.ssdd**.

## Purpose

The goal of SSDD file format is to save scattering properties in a file with the effective conversion of data, colors, parameters, and so on. Existing file formats depend on a specific software or measurement device. They supports a small number of property types. As a result, the file saving and conversion of arbitrary BSDF are difficult. The SSDD file solves this problem by supporting multiple types of scattering data, color models, and parameterizations of incoming and outgoing directions.

**Following properties are supported:**

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

## Format

The file format consists of a file header and data blocks. The file header is stored before data blocks. A data block contains the meta-data and tabular data of a BRDF, BTDF, specular reflectance, or specular transmittance. A file header and at least one data block exist in a file.

**Overall SSDD file structure:**

```
[file header]
[meta-data 1]
[tabular data 1]
...
[meta-data 4]
[tabular data 4]
```

Each data block has different data type. For example, two BRDFs are not stored at once. The file header and meta-data are in ASCII form. Tabular data is in ASCII or binary form.

Comments in ASCII sections start at **#** followed by a space character and end at a new line (**\n**). The string of each name property is the sequence of UTF-8 encoded characters and terminated by a new line (**\n**). The format of date is [ISO 8601][ISO_8601] (YYYY-MM-DD).

### File header

File information is stored in a file header as ASCII text. The file header must start with **VERSION \<*number*\>**, which is the entry for the file version. Others are in random order and optional.

| Name      | Type      | Value                                     | Description |
|-|-|-|-|
| VERSION   | string    |\<*major*\>.\<*minor*\>                    | Version of SSDD file. The string consists of the major version number followed by a dot and the minor version number. |
| SOFTWARE  | string    |  <nobr>e.g. *BSDFProcessor-1.3.0*</nobr>  | Name of software used to create this file. |
| API       | string    | e.g. *libbsdf-1.0.0*                      | Name of API used to create this file. |
| DATE      | string    | e.g. *2020-01-16*                         | Date of file creation. |

### Data block

A BRDF, BTDF, specular reflectance, or specular transmittance is stored in a data block. The data structure consists of meta-data and tabular data. Meta-data has properties of data and the format type of tabular data. Properties include a scattering type, color model, parameterization type, etc. Tabular data is expressed as the list of colors. An SSDD file can contain up to four data blocks.

Specular reflectance and transmittance are supported in addition to BRDF and BTDF as date types, because BRDF/BTDF of specular surface is the Dirac delta function. In other words, the value is infinite at the specular direction in infinitesimal projected solid angle. This causes difficulty for the tabular representation and measurement of specular BRDF/BTDF. Specular reflectance and transmittance for each incoming direction are practical workarounds for this problem.

#### Meta-data

Properties of data are stored in ASCII form. The following is the list of entries.

| Name              | Type      | Value                                                                     | Description |
|-|-|-|-|
| DATA_TYPE         | string    | *brdf*<br> *btdf*<br> *specular_reflectance*<br> *specular_transmittance* | Type of data. This entry is located in the starting position of meta-data. |
| COLOR_MODEL       | string    | *monochrome*<br> *rgb*<br> *xyz*<br> *spectrum*                           | Color model of data. ***monochrome*** is a scalar model. ***rgb*** and ***xyz*** have tristimulus values. The color space of ***rgb*** is not defined in the current SSDD version. ***xyz*** is the CIE XYZ color space. ***spectrum*** has multiple channels with the **WAVELENGTH_LIST** entry. |
| WAVELENGTH_LIST   | float[]   | e.g. *400 500 600 700*                                                    | List of wavelengths expressed in nanometers (nm). This is required if **COLOR_MODEL** is ***spectrum***. Numbers are in ascending order. |
| PARAM_TYPE        | string    | *spherical_coordinate_system*<br> *specular_coordinate_system*<br> *half_difference_coordinate_system* | BRDF/BTDF parameterization type of incoming and outgoing directions. Parameter lists are contained in **PARAM0_LIST**, **PARAM1_LIST**, **PARAM2_LIST**, **PARAM3_LIST**, and **PARAM4_LIST**. Details are explained in [Parameterization types](#parameterization-types). |
| REDUCTION_TYPE    | string[]  | *bilateral_symmetry*<br> *reciprocity*                                    | Reduction types of **PARAM3_LIST**. ***bilateral_symmetry*** uses a symmetry along the incident plane. ***reciprocity*** uses the Helmholtz reciprocity principle and is only valid for ***half_difference_coordinate_system***. The range is reduced to [0, 180]. If both types are specified, the range is [0, 90]. |
| PARAM0_LIST       | double[]   | e.g. *0 15 30 45 60 75 90*                                                | List of parameters. The definition of parameter depends on **PARAM_TYPE**, see [Parameterization types](#parameterization-types) for details. |
| PARAM1_LIST       | double[]   | e.g. *0 90 180 270 360*                                                   | List of parameters. The definition of parameter depends on **PARAM_TYPE**, see [Parameterization types](#parameterization-types) for details. |
| PARAM2_LIST       | double[]   | e.g. *0 0.5 1 1.5 2 3 4 5 7.5 10 12.5 15 20 30 45 60 75 90*               | List of parameters. The definition of parameter depends on **PARAM_TYPE**, see [Parameterization types](#parameterization-types) for details. |
| PARAM3_LIST       | double[]   | e.g. *0 45 90 135 180 225 270 315 360*                                    | List of parameters. The definition of parameter depends on **PARAM_TYPE**, see [Parameterization types](#parameterization-types) for details. |
| PARAM4_LIST       | double[]   | e.g. *0 -5.06 -10.52 -16.87 -24.73 -34.9 -48.18*                          | List of optional parameters.  This is valid if **PARAM_TYPE** is ***specular_coordinate_system***, see [Parameterization types](#parameterization-types) for details. |
| NAME              | string    | e.g. *metallic paint*                                                     | Name of data. |
| SOURCE_TYPE       | string    | *measured*<br> *generated*<br> *edited*                                   | Way of creating data. Types are predefined for scattering data measured by a device, generated artificially, and edited for specific purposes. |
| DEVICE            | string    | e.g. *GCMS-4*                                                             | Name of the measurement device. |
| CREATION_DATE     | string    | e.g. *2020-01-16*                                                         | Date of data creation. |
| MEASUREMENT_DATE  | string    | e.g. *2019-12-26*                                                         | Date of measurement. |
| DATA              | string    | *ascii*<br> *binary*                                                      | Format type of the following tabular data. This entry is located in the ending position of meta-data. |

Meta-data must start with **DATA_TYPE \<*type*\>** and end with **DATA \<*mode*\>**.
**DATA_TYPE**, **COLOR_MODEL**, **PARAM0_LIST**, **DATA** are mandatory entries.
**PARAM_TYPE** must be used for BRDF/BTDF. **WAVELENGTH_LIST**, **PARAM1_LIST**, **PARAM2_LIST**, and **PARAM3_LIST** are required for specific settings. Others are in random order and optional.

**Entries must be set in the following order:**

```
VERSION
(Optional entries of file header in random order)
DATA_TYPE
COLOR_MODEL
WAVELENGTH_LIST
PARAM_TYPE
REDUCTION_TYPE
PARAM0_LIST
PARAM1_LIST
PARAM2_LIST
PARAM3_LIST
PARAM4_LIST
(Optional entries of meta-data in random order)
DATA
```

#### Parameterization types

Incoming and outgoing directions of BRDF/BTDF are expressed as three or four parameters. If ***brdf***/***btdf*** is specified in the **DATA_TYPE \<*type*\>** entry, a BRDF/BTDF parameterization type (**PARAM_TYPE**) and parameter lists (**PARAM0_LIST**, **PARAM1_LIST**, **PARAM2_LIST**, and **PARAM3_LIST**) are used to construct data structure. **PARAM_TYPE** is not used for specular reflectance/transmittance (***specular_reflectance***/***specular_transmittance***) because incoming directions in a spherical coordinate system are sufficient as parameters.

Followings are currently supported parameterization types and parameter lists. The angle unit of parameter is in degrees. Numbers are in ascending order.

***spherical_coordinate_system*** consists of two spherical coordinates in each of the incoming and outgoing directions.

| Name          | Definition                                | Range                 | Description |
|-|-|-|-|
| PARAM0_LIST   | Incoming polar angles                     | [0, 90]               ||
| PARAM1_LIST   | <nobr>Incoming azimuthal angles</nobr>    | <nobr>[0, 360]</nobr> | This property is required for anisotropic data. |
| PARAM2_LIST   | Outgoing polar angles                     | [0, 90]               ||
| PARAM3_LIST   | Outgoing azimuthal angles                 | [0, 360]              ||

***specular_coordinate_system*** has the incoming direction of spherical coordinates and the reparameterized outgoing direction. The latter is spherical coordinates rotated toward the specular direction.

| Name          | Definition                                | Range                   | Description |
|-|-|-|-|
| PARAM0_LIST   | Incoming polar angles                     | [0, 90]                 ||
| PARAM1_LIST   | <nobr>Incoming azimuthal angles</nobr>    | [0, 360]                | This property is required for anisotropic data. |
| PARAM2_LIST   | Specular polar angles                     | [0, 180]                | Angles between outgoing and specular directions. |
| PARAM3_LIST   | Specular azimuthal angles                 | [0, 360]                | Rotation angles of outgoing directions around specular directions. |
| PARAM4_LIST   | Offset angles                             | <nobr>[-90, 90]</nobr>  | Offsets from specular directions. The size of the list must be equal to the number of incoming polar angles (**PARAM0_LIST**). This parameter list is efficient to represent the BTDF with refraction. |

[***half_difference_coordinate_system***][half_diff_cs] is based on the halfway vector between incoming and outgoing directions. Half polar and azimuthal angles are defined using halfway and normal directions. Difference polar and azimuthal angles are parameters between incoming and halfway directions.

| Name          | Definition                                | Range                 | Description |
|-|-|-|-|
| PARAM0_LIST   | Half polar angles                         | [0, 90]               | Angles between halfway directions and the normal direction. |
| PARAM1_LIST   | Half azimuthal angles                     | <nobr>[0, 360]</nobr> | Rotation angles of halfway directions around the normal direction. This property is required for anisotropic data. |
| PARAM2_LIST   | Difference polar angles                   | [0, 90]               | Angles between incoming and halfway directions. |
| PARAM3_LIST   | <nobr>Difference azimuthal angles</nobr>  | [0, 360]              | Rotation angles of incoming directions around halfway directions. |

Specular reflectance (***specular_reflectance***) and transmittance (***specular_transmittance***) are defined in each incoming direction. The parameterization type is a spherical coordinate system.

| Name          | Definition                                | Range     | Description |
|-|-|-|-|
| PARAM0_LIST   | Incoming polar angles                     | [0, 90]   ||
| PARAM1_LIST   | <nobr>Incoming azimuthal angles</nobr>    | [0, 360]  | This property is required for anisotropic data. |
| PARAM2_LIST   | Unused                                    |||
| PARAM3_LIST   | Unused                                    |||

#### Tabular data

Data is stored below the **DATA \<*mode*\>** entry. Values are in ascii or binary form.
The size of data depends on **COLOR_MODEL** and **PARAM_TYPE**. The number of stored values is derived by multiplying the size of color channels, **PARAM0_LIST**, **PARAM1_LIST**, **PARAM2_LIST**, and **PARAM3_LIST**. Data consists of the list of colors. A color is a monochromatic value, RGB, XYZ, or spectrum for each pair of incoming and outgoing directions. In ascii form, a color is stored in a line. Values in a line are separated by one or more spaces or tabs and terminated by a new line (**\n**). The binary format of floating-point number is IEEE 754 in little-endian ordering.

The index of a value is got by the following function:

```C
index = param0_index
      + param0_list_size * param1_index
      + param0_list_size * param1_list_size * param2_index
      + param0_list_size * param1_list_size * param2_list_size * param3_index
```

## Example

Lambertian BRDF with a sample point.

```
VERSION 0.2

DATA_TYPE brdf
COLOR_MODEL monochrome
PARAM_TYPE spherical_coordinate_system
PARAM0_LIST 0
PARAM2_LIST 0
PARAM3_LIST 0
DATA ascii
0.3183
```

Material containing BRDF, BTDF, specular reflectance, and specular transmittance.

```
VERSION 0.2

# Bluish diffuse BRDF
DATA_TYPE brdf
COLOR_MODEL rgb
PARAM_TYPE half_difference_coordinate_system
REDUCTION_TYPE bilateral_symmetry reciprocity
PARAM0_LIST 0 90
PARAM2_LIST 0 90
PARAM3_LIST 0 90
DATA ascii
# PARAM3: 0
# PARAM2: 0
0 0.0524275 0.159779
0 0.0524275 0.159779
# PARAM2: 90
0 0.0524275 0.159779
0 0.0524275 0.159779
# PARAM3: 90
# PARAM2: 0
0 0.0524275 0.159779
0 0.0524275 0.159779
# PARAM2: 90
0 0.0524275 0.159779
0 0.0524275 0.159779

# Orange diffuse BTDF in the CIE XYZ color space
DATA_TYPE btdf
COLOR_MODEL xyz
PARAM_TYPE specular_coordinate_system
REDUCTION_TYPE bilateral_symmetry
PARAM0_LIST 0 90
PARAM2_LIST 0 180
PARAM3_LIST 0 180
DATA ascii
# PARAM3: 0
# PARAM2: 0
0.106 0.082 0.011
0.106 0.082 0.011
# PARAM2: 180
0.106 0.082 0.011
0.106 0.082 0.011
# PARAM3: 180
# PARAM2: 0
0.106 0.082 0.011
0.106 0.082 0.011
# PARAM2: 180
0.106 0.082 0.011
0.106 0.082 0.011

# Simple specular reflectance
DATA_TYPE specular_reflectance
COLOR_MODEL spectrum
WAVELENGTH_LIST 400 450 500 550 600 650 700
PARAM0_LIST 0 90
DATA ascii
0.05 0.05 0.05 0.05 0.05 0.05 0.05
0.05 0.05 0.05 0.05 0.05 0.05 0.05

# Simple specular transmittance
DATA_TYPE specular_transmittance
COLOR_MODEL monochrome
PARAM0_LIST 0 90
DATA ascii
0.05
0.05
```

[ISO_8601]: https://www.iso.org/iso-8601-date-and-time-format.html "ISO 8601"
[half_diff_cs]: https://www.cs.princeton.edu/~smr/papers/brdf_change_of_variables/ "Rusinkiewicz, S. 1998. A New Change of Variables for Efficient BRDF Representation."
