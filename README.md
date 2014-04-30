VERAout-tools
=============

HDF Tools to post-process VERA output files



#### Example: How to install HDF5 Library on your computer

Download latest HDF5 source:
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
  dir /opt

Cygwin:

  ./configure --enable-fortran --enable-cxx --prefix=/usr/local/lib
  make
  make check
  make install
```

