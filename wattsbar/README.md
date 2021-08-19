
Utility Program to create VERA HDF file for
Watts Bar Reactor using user data.

Watts Bar is a standard Westinghouse 4-Loop PWR.

All of the geometry parameters are hard-wired in the "Mod_geom" module.
These parameters can be modified if you want to apply to a different reactor.


## Data Format

Pin data (pin power, pin exposures, etc.) are arranged in 4D arrays

Dimensions:
```
    pin(maxassm,kdfuel,npin,npin)

    npin   - number of pins across the assembly (17)
    kdfuel - number of axial levels (non-uniform) (49)
    maxasm - number of assemblies in qtr-core (56)
```

Order of data:
```
    pin(na,k,i,j)
```


## Assemblies

There are 56 assemblies in qtr-core symmetry.

The assemblies are numbered in the following order:
```
     1  2  3  4  5  6  7  8
     9 10 11 12 13 14 15 16
    17 18 19 20 21 22 23 24
    25 26 27 28 29 30 31 32
    33 34 35 36 37 38 39  0
    40 41 42 43 44 45 46  0
    47 48 49 50 51 52  0  0
    53 54 55 56  0  0  0  0
```

## Axial levels

There are 49 non-uniform axial levels.
1 is the bottom plane with fuel.
49 is the top plane with fuel.

(there are usually more planes in the VERA models, but only data in planes
with fuel is written to the data file)

## Pin data

Pin data is on an (17,17) grid with ordering (i,j).

The I-direction goes from left-to-right
```
   1, 2, 3, 4, ... 17
```

The j-direction goes from top to bottom
```
  1
  2
  3
  ...
  17
```

# User Data

User data is specified in a text file
There is one line per pin and axial level.
Each line has the following format:

```
  na, k, i, j, pow(na,k,i,j) /
```
Refer to the sample file "power.inp" for an example.


# Usage

```
> makewb.exe power.inp power.h5
```

This will create a VERA HDF file which can be read by VERAView.




