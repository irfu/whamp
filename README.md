WHAMP
=====
[![ci-build](https://github.com/irfu/whamp/actions/workflows/ci-build.yml/badge.svg)](https://github.com/irfu/whamp/actions/workflows/ci-build.yml)
[![DOI](https://zenodo.org/badge/12994595.svg)](https://zenodo.org/doi/10.5281/zenodo.11639728)

WHAMP - Waves in Homogeneous Anisotropic Magnetized Plasma.
Fortran code, originally written by Kjell Rönnmark, that calculates
general wave dispersion relation in plasmas.

Package includes also a matlab code to run WHAMP from inside MATLAB
and to visualize distribution functions.


# WHAMP - Wave in Homogeneous, Anisotropic, Multi-component Plasma

This is a modified version of the WHAMP dispersion relation solver that allows model input from files and provides additional output parameters. The original WHAMP code can be downloaded from [http://www.tp.umu.se/forskning/space/WHAMP/](http://www.tp.umu.se/forskning/space/WHAMP/).

## Compilation

To build the code, simply run:

```bash
make
```

This will create the following executables:
- `whamp` - Main program
- `whamp_engine_test` - Test program for the engine
- `libwhamp.a` - Static library

## Usage

Run WHAMP with a plasma model file:

```bash
./whamp -file <model_file>
```

### Example

```bash
./whamp -file ../Models/Ex3
```

## Input Model File Format

The model file contains plasma parameters in the following order (with optional comments starting with `#`):

```
# Density values (m^-3) for 10 species
n(1)  n(2)  ...  n(10)        

# Temperature values (eV) for 10 species  
t(1)  t(2)  ...  t(10)        

# Loss cone parameter D for 10 species
d(1)  d(2)  ...  d(10)        

# Temperature anisotropy A for 10 species
a(1)  a(2)  ...  a(10)        

# Loss cone parameter B for 10 species
b(1)  b(2)  ...  b(10)        

# Atomic mass for 10 species (0=electron, 1=proton, 16=oxygen, etc.)
ass(1)  ass(2)  ...  ass(10)  

# Drift velocity for 10 species (v_drift/v_thermal)
vd(1)  vd(2)  ...  vd(10)     

# Gyrofrequency (kHz)
fce                           

# PZL parameter (1=log scale, 0=linear scale)
pzl                           
```

### Example Model File (Ex3)

```
# WHAMP plasma model file - Example 3
# Density values (m^-3) for 10 species
1.4e6 1.4e6 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
# Temperature values (eV) for 10 species  
.002  .002 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
# Loss cone parameter D for 10 species
1.   1.   1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 
# Temperature anisotropy A for 10 species
1.   1.   1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 
# Loss cone parameter B for 10 species
0.0  0.0  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
# Atomic mass for 10 species (0=electron, 1=proton, etc.)
0.  1.0  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
# Drift velocity for 10 species
0.   .0  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
# Gyrofrequency (kHz)
3.3
# PZL parameter
0
```

## Example Output

```bash
./whamp -file ../Models/Ex3
# TOTAL PLASMA FREQ.:   10.62634 KHZ; SPEC 1 GYRO FREQ.:   3.30000 KHZ; SPEC 1 PLASMA FREQ.:   10.62345 KHZ; 
# SPECIES 1 V_TH/C:   0.002798;  SPECIES 1 BETA:   0.009007 (ratio of thermal velocity to gyrofrequency)
# e-   DN= 1.40000E+06  T=  0.00200  D=1.00  A=1.00  B=0.00 VD= 0.00
# H+   DN= 1.40000E+06  T=  0.00200  D=1.00  A=1.00  B=0.00 VD= 0.00
#INPUT: 
  
p0z.01f.1
#OUTPUT: 

pzf
    0.0000000      0.0100000    5.4631239E-01 -1.01E-12    

#INPUT: 
  
o
#OUTPUT: 

pzf/ey
    0.0000000      0.0100000    5.4631239E-01 -1.01E-12    
 EX= 0.7071  0.0000  EY=-0.0000  0.7071  EZ=-0.0000  0.0000  enden= 0.504E+02 enfl_p=      NaN enfl_z= -.206E-05  enden= 0.196E-01 enfl_p=      NaN enfl_z= 0.349E-12   

#INPUT: 

```

## Interactive Commands

After loading a model, WHAMP enters interactive mode where you can:

- Enter wave parameters like `p0z.01f.1` (parallel wavenumber, perpendicular wavenumber, frequency)
- Use output commands like `pzf/e` to get detailed wave field information

## Available Output Parameters

The output parameters that can be specified include:

| Parameter | Description |
|-----------|-------------|
| `e` | Electric field components (ex, ey, ez) |
| `b` | Magnetic field components (bx, by, bz) |
| `f` | Frequency (real, imaginary) |
| `g` | Group velocity |
| `h` | \|E\|/\|B\| ratio [mV/nT] |
| `l` | \|Bp\|/\|Bz\| ratio |
| `m` | Im[bx]/Re[by] |
| `n` | Ellipticity bx/by |
| `o` | General ellipticity |
| `p` | Perpendicular wavenumber |
| `s` | Spatial growth (sp, sz) |
| `u` | Energy ratio between total wave energy and electric field energy |
| `v` | Poynting flux (in μW/m² for ⟨E²⟩=0.5(mV/m)²) |
| `z` | Parallel wavenumber |
| `x` | Phase of bz against bx |
| `y` | Energy density and flux of each component |

## Changes from Original WHAMP

This version includes:
- Support for input from model files with comments
- Additional output parameters including plasma beta
- Enhanced plasma parameter display
- Improved input file parsing with comment support

## References

- Original WHAMP: [http://www.tp.umu.se/forskning/space/WHAMP/](http://www.tp.umu.se/forskning/space/WHAMP/)
- See `whamp_quick_reference.pdf` for detailed


Some articles that use WHAMP:
https://ui.adsabs.harvard.edu/public-libraries/Iw3OYJ1iQa2WSeG8SfHxyg



