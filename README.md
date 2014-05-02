VERAout-tools
=============

HDF Tools to post-process VERA output files

VERA stands for the "Virtual Environment for Reactor Analysis" and is 
a product of the Consortium for Advanced Simulation of LWR's (CASL).
More information about CASL can be found at [http://www.casl.gov]

The latest version of this repository can be found at 
[https://github.com/palmtag/VERAout-tools]

This repository contains the following source files:

* Mod_hdftools.f90 - Wrapper routines for performing basic HDF functions
* hunittest.f90 - Simple program to perform unit tests
* makehdf.f90 - Simple program to create a sample VERAout file for testing
* mpactread.f90 - Simple program to read VERAout file created by MPACT


#### Hint: How to install HDF5 Library on your computer

This section contains information on installing HDF5 on your computer system.

You must have a fortran compiler installed on your system.

Information on the latest HDF5 source can be found at:
   http://www.hdfgroup.org/HDF5/release/obtain5.html

```
All:

  HVER=hdf5-1.8.12    # latest source on 2014/04/29
  wget http://www.hdfgroup.org/ftp/HDF5/current/src/$HVER.tar.gz
  tar vxfz  $HVER.tar.gz
  cd $HVER

Linux:

  # export F9X=ifort     # add this if using Intel Fortran (not recommended)
  ./configure --enable-fortran --enable-cxx --prefix=/opt/hdf5
  make
  make check
  sudo make install
  ls -Fl /opt/hdf5

Cygwin:

  ./configure --enable-fortran --enable-cxx --prefix=/usr/local/lib
  make
  make check
  make install
```

