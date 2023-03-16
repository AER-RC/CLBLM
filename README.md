# CLBLM

---
**Contents**

1. [Introduction](#intro)
2. [Cloning the Latest Release](#cloning)
3. [General LNFL/CLBLM File Information](#general)
    1. [Platforms on which CLBLM can be run](#platforms)
    2. [Issues relating to unformatted files on UNIX and LINUX systems](#unformatted)
	  3. [LNFL/CLBLM Naming Convention](#nomenclature)
4. [Instructions and Tips for Running LNFL](#runlnfl)
	1. [Input files for LNFL](#lnflin)
	2. [Output files for LNFL](#lnflout)
	3. [Sequence for running LNFL](#lnflseq)
5. [Instructions and Tips for Compiling and Running CLBLM](#runclblm)
	1. [Required input files for CLBLM](#clblmin)
	2. [Layer numbering scheme](#laynum)
	3. [Output files for CLBLM](#clblmout)
	4. [Sequence for running CLBLM](#clblmseq)
6. [Tests](#tests)
7. [Frequently Asked Questions](#faq)

# Introduction <a name="intro"></a>

CLBLM (Community Line-By-Line radiative transfer Model) is an accurate and efficient line-by-line radiative transfer model derived from the Line-By-Line Radiation Transfer Model (LBLRTM). CLBLM is a modernized version of LBLRTM, written in Fortran 90 and utilizing netCDF for input and output files.

The [HITRAN database](http://cfa-www.harvard.edu/hitran) provides the basis for the line parameters used in CLBLM. These line parameters, as well as additional line parameters from other sources, are extracted for use in CLBLM by a line file creation program called LNFL. A line parameter database built from HITRAN and suitable for use with LNFL can be downloaded with the [AER Line File retrieval code](https://github.com/AER-RC/AER_Line_File) or directory from the [Zenodo repository](https://zenodo.org/record/4019178).

CLBLM uses the line parameters and [MT_CKD continuum](https://github.com/AER-RC/MT_CKD) in its calculations. The models and data are thus linked. For the latest release, the relationships are:

| CLBLM Release | MT_CKD Release | Line File |
| :---: | :---: | :---: |
| [v1](https://github.com/AER-RC/CLBLM/releases/tag/v1) | [4.1.1](https://github.com/AER-RC/MT_CKD/releases/tag/4.1.1) | [v3.8.1](https://zenodo.org/record/4019178/files/aer_v_3.8.1.tar.gz?download=1) |

If any build or run issues occur, please [create an issue](https://github.com/AER-RC/LBLRTM/issues) or contact the [AER-RC Group](https://github.com/AER-RC).

For more, please see the [Wiki page](https://github.com/AER-RC/CLBLM/wiki/)

# Cloning the Latest Release <a name="cloning"></a>

Assuming the output directory should be `CLBLM` and that the user has [created an SSH key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent) (we use RSA):

```
git clone --recursive git@github.com:AER-RC/CLBLM.git
```

Alternatively, users that have not setup an SSH key but have [created a personal access token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token#creating-a-fine-grained-personal-access-token) can use the HTTPS protocol:

```
git clone --recursive https://github.com/AER-RC/CLBLM.git
```

`--recursive` is important, because this repository is linked with our [common FORTRAN modules repository](https://github.com/AER-RC/aer_rt_utils) that are required in the model builds. The [cross section database](https://github.com/AER-RC/cross-sections) is also added as a submodule (it is not required for all model runs). If this keyword is forgotten, one can do:

```
git submodule init
git submodule update
```

in the `CLBLM` directory.

The current release is CLBLM v1, and it is recommended that this be the version that users clone and checkout (rather than the `master` branch). To do this, one needs to simply checkout the `v1` tag:

```
git checkout tags/v1
```

Instead of cloning, users can also download a CLBLM [tarball](https://github.com/AER-RC/CLBLM/archive/v1.zip) and unpack it:

```
tar xvf clblm_v1.tar.gz
```

# General LNFL/CLBLM File Information <a name="general"></a>

## Platforms on which CLBLM can be run <a name="platforms"></a>

It is recommended that LNFL and CLBLM be compiled in Fortran 90. LBLRTM has previously been run on DEC alpha, Cray, MS-DOS, and HP platforms.

## Issues relating to unformatted files on UNIX and LINUX systems <a name="unformatted"></a>

Unformatted files are often not compatible between systems due to differences in the way the bytes are written to the files (big-endian versus little-endian).  Note that the `byteswap` option available with most compilers will not work with most CLBLM unformatted output files because of the mixing of real and integer data within records.

## LNFL/CLBLM Naming Convention <a name="nomenclature"></a>

Specific information on the input/output files from LNFL and CLBLM is located in their respective instruction manuals, `lnfl_instructions` and `clblm_instructions`, and the examples provided in the code `tar` files.  

# Instructions and Tips for Running LNFL <a name="runlnfl"></a>

LNFL is used to generate a unformatted file (`TAPE3`) of all the line parameters required by LBLRTM.

## Input files for LNFL <a name="lnflin"></a>

1. `TAPE1`: The line parameter database in ASCII format (downloaded with the [AER Line File repository](https://github.com/AER-RC/AER_Line_File) or from [Zenodo](https://zenodo.org/record/4019178)).

2. `TAPE5`: LNFL input file. 
The `TAPE5` input file is read as formatted FORTRAN. As a consequence of the formatted read, any blank space will be read as "zero".  Thus, one may leave blanks for most of the parameters and within the code they will default to an acceptable value. Real numbers format input as either `E` or `F` format, with the entire number within the range specified in the input instructions.  Integers are read in with the `I` format and must be specified exactly in the integer format.

## Output files for LNFL <a name="lnflout"></a>

1. `TAPE3`: Unformatted LNFL output file containing the line parameters for LBLRTM.
2. `TAPE6`: Informational output file.
3. `TAPE7`: Optional output file containing ASCII version of the parameters contained in `TAPE3`.

## Sequence for running LNFL <a name="lnflseq"></a>

*	[Clone the latest LNFL code](https://github.com/AER-RC/LNFL#cloning-)  and download the latest line parameter database with the [AER Line File repository](https://github.com/AER-RC/AER_Line_File) or from [Zenodo](https://zenodo.org/record/4019178).
*	Compile LNFL using the makefiles found in the LNFL tar file.  Note: one needs to compile in the `build` directory.
*	Link the line parameter database to `TAPE1` in the LNFL working directory.
*	Remove `TAPE3` file from the LNFL working directory.
*	Edit necessary parameters in the `TAPE5` input file.  Note that the beginning and ending wavenumber (_v<sub>1</sub>_, _v<sub>2</sub>_) in `TAPE5` must extend at least 25 cm<sup>-1</sup> beyond each end of the desired spectral range for the LBLRTM calculations.
*	Run the LNFL executable.

# Instructions and Tips for Compiling and Running CLBLM <a name="runclblm"></a>

CLBLM is used to generate line-by-line upwelling and downwelling transmittances and radiances.

## Required input files for CLBLM <a name="clblmin"></a>
1. `scenes.nc`: NetCDF file containing profile information (ex. temperature, pressure, molecular amounts, etc.) required to run CLBLM.
2. `clblm_config.json`: JSON control file with model run parameters (ex. requested output, spectral inverval, instrument functions, etc.) , required to run CLBLM.

CLBLM does not have a limit for the spectral interval, unlike LBLRTM (runs must not exceed 2000 cm<sup>-1</sup>). 

Other input files are required if you are using the solar source function, cross sections, surface emissivity, etc. See the CLBLM instruction manual and examples provided.

## Layer numbering scheme <a name="laynum"></a>

The CLBLM convention is that layer 1 is at the highest pressure level (lowest altitude).  The layer information for a given run may be found in the scenes.nc input file.

## Output files for CLBLM <a name="clblmout"></a>
The CLBLM output is a netCDF file containing the desired output (ex. transmittances/radiances).

The file specifies the dimensions of the output as the number of points in the model run (_numPoints_) and for some runs, the number of vertical levels (_numLayers_). The global attributes include the (_spectralDataType_) such as "monochromatic" or "convolved", the starting and ending points of the model run (_v<sub>1</sub>_ and _v<sub>2</sub>_), and the spectral spacing of the points (_dv_). 


## Sequence for running CLBLM <a name="clblmseq"></a>
* [Clone the latest CLBLM code](https://github.com/AER-RC/CLBLM#cloning-the-latest-release-) and download the latest line parameter database with the [AER Line File repository](https://github.com/AER-RC/AER_Line_File) or from [Zenodo](https://zenodo.org/record/4019178).
* Compile CLBLM, scene_writer, and build_solar with the makefile in the CLBLM tar file.
```
make all
``` 
* Link the line parameter database (`TAPE3` from LNFL) to CLBLM/clblm_data/spectroscopy/TAPE3.
* Use the sceenwriter executable to build the profile input from an existing LBLRTM TAPE5
```
scene_writer TAPE5_user_defined_upwelling  user_archive/scene_files/scenes.nc
```
* For a solar radiance run, set the solar configuration options (ex. start/end wavenumber, solar irradiance data file) in solar_config.json and use the build_solar executable to build the solar source function file (SOLAR.RAD.nc).
```
build_solar solar_config.json
```
* Edit the JSON control file (clblm_config.json) to set all the configuration options (similar to LBLRTM TAPE5 parameters) for the model run. 
* Run the CLBLM executable.
```
clblm
```

# Tests <a name="tests"></a>

A [run example package](https://github.com/AER-RC/CLBLM/releases/tag/v1/clblm_v1.examples.tar) is provided separately from the code repository. It can be used to validate building and running of the model for select atmospheric specifications and model configurations. See `README.setup` in top level of the package for further direction.


# Frequently Asked Questions <a name="faq"></a>

1. **What is the difference between a line-by-line calculation and a band-model calculation?**

Absorption/emission spectra are comprised of a complicated array of spectral lines.  The HITRAN 2008 Database (Version 13.0) contains over 2,713,000 lines for 39 different molecules. In order to resolve these individual lines, a nominal spectral sampling rate of less than the mean line half width must be utilized.  Such highly resolved radiative transfer calculations are called line-by-line (LBL) calculations. The computational time associated with calculating broadband fluxes from LBL calculations is formidable.  A band model aims to simplify radiative transfer calculations by using approximations to represent the line-by-line characteristics of a particular spectral interval.  Band models are appropriate for situations where the desired spectral resolution is much smaller than the Lorentz and Doppler widths of the spectral lines. Such approximations are also of use in general circulation models.

2. **What are the standard units used in CLBLM calculations?**

| Output Variable | Units |
| :---: | :---: |
| Wavenumber | cm<sup>-1</sup> |
| Radiance | W cm<sup>-2</sup> sr<sup>-1</sup> /  cm-1 |
| Brightness Temperature | K |
| Analytic Jacobians (dR/dx) | <ul><li>molecules: W cm<sup>-2</sup> sr<sup>-1</sup> / cm<sup>-1</sup> / log(VMR)</li><li>temperature: W cm<sup>-2</sup> sr<sup>-1</sup> / cm-1  / K |

3. **Radiance Derivatives (Jacobians)**

CLBLM features a straighforward setup for Jacobians, in which Jacobians can be output directly from a single CLBLM run. In the clblm_config.json input file, set clblm_out to jacobians and specify which jacobians you want to run. See the CLBLM example package and instructions for more details.
```
"clblm-out":               {"convolved jacobians":"clblm_out/jacobian.nc","jacobian-list": ["T","H2O","Tskin","emis"]}
```
	
4. **Is it possible to scale the profile of one or more species?**

Yes. Sceen_writer is able to handle the scaling provided in the LBLRTM TAPE5.

5. **Does CLBLM include heavy molecule parameters (cross-sectional species)?**

Heavy molecules (such as CCL4, F11, and others) can be included in CLBLM calculations. Sceen_writer handles in the inclusion of cross-section molecules included in the LBLRTM TAPE5. An additional file (`FSCDXS`) and directory (`xs`) are required for these calculations and can be obtained from the [cross sections repository](https://github.com/AER-RC/cross-sections) (which is cloned with LBLRTM if the directions in the [Cloning][#cloning] section are followed) or from the CLBLM example tar file (available in the [CLBLM v1 Release](https://github.com/AER-RC/CLBLM/releases/download/v1/clblm_v1.examples.tar)).

6. **Format of external surface emissivity/reflectivity files**

Sea surface spectral emissivity and reflectivity files are provided with the example (available in the [CLBLM v1 Release](https://github.com/AER-RC/CLBLM/releases/download/v1/clblm_v1.examples.tar)). The files must have the file names of `EMISSIVITY` and `REFLECTIVITY`. The format is as follows:

| Parameter | Format | Description |
| :--- | :---: | :--- |
| `V1EMIS` | `E10.3` | Initial emissivity/reflectivity frequency value [cm<sup>-1</sup>] |
| `V2EMIS` | `E10.3` | Finial emissivity/reflectivity frequency value [cm<sup>-1</sup>] |
| `DVEMIS` | `E10.3` | Frequency Increment [cm<sup>-1</sup>] |
| `NLIMEM` | `I5` | Number of spectral emissivity/reflectivity points in the file |
| `ZEMIS` | `E15.7` | Emissivity at each spectral point |

**NOTE**: It is assumed that the spectral emissivity/reflectivity points are equally spaced and there is a maximum number of points (see LBLRTM instructions).

7. **Absorption due to clouds/aerosols**

Absorption due to clouds and aerosols is not available with the initial release of CLBLM.

8. **Solar Radiance**

The CLBLM package contains two solar irradiance datasets; an average solar irradiance file with no temporal variability and a multicomponent solar irradiance file that allows the user to account for variability in the 11-year solar cycle, facular brighening and sunspot darkening. The solar irradiance data is derived from the NRLSSI2 model and covers a wavenumber range of 100-860000. Parameters for solar radiance calculations are set in the clblm_config.json file. For example,
```
"solar-irradiance":           {"option": 2, "cycle-frac" :0.382576, "facula-var": 1.0, "spot-var": 1.0}
```
Solar radiance runs with CLBLM require a solar source function file named `SOLAR.RAD.nc`. The user will generate `SOLAR.RAD.nc` using the `build_solar` executable to extract solar irradiance data from one of the above datasets. The parameters for `build_solar` are specified in the solar_config.json. For example,
```
"inputs":                     [{"path": "/project/p2326/broot/build_comb_solar_rad_multi_comp_50000plus.nc",
                                "start-wavenumber":22000.0,
                                "end-wavenumber":24000.0}],

```
See the CLBLM example package and instructions for more details.

9. **Line coupling/mixing**

Line coupling parameters are utilized in LBLRTM for O<sub>2</sub>, CO<sub>2</sub> and CH<sub>4</sub>. The line coupling parameters are provided in the AER line parameter database (available in the [AER Line File repository](https://github.com/AER-RC/AER_Line_File) or on [Zenodo](https://zenodo.org/record/4019178)) and are written to the line parameter input file (`TAPE3`) by LNFL.

10. **What is the appropriate reference for CLBLM calculations in journal articles and presentations?**

11. **How do you calculate fluxes?**

CLBLM features built-in flux calculations, such that radiative fluxes and heating rates can be output directly from a single CLBLM run. In the clblm_config.json input file, add the flux-flags group and specify the flux parameters, such as the spectral interval for the flux output and the number of quadrature angles to use for the flux calculation. See the CLBLM example package and instructions for more details.
```
"flux-flags":                  {"flux_flag":true, "dv_flux":10.0,"nang":3}
```

