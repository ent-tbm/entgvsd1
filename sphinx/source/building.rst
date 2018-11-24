.. _building

Building EntGVSD
================

Install Prerequisites
---------------------

EntGVSD requires the following prerequisites:

* GNU Fortran / C / C++ version 4.9 or later.

* `NetCDF-C <https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html>`_, version 4 or later

* `NetCDF-Fortran <https://www.unidata.ucar.edu/software/netcdf/docs/building_netcdf_fortran.html>`_, version 4 or later

* (For now) `Python <https://www.python.org>`_ version 3.5 or later

There are a variety of ways these prerequisites can be installed.
Some suggestions are provided below; assuming you have already
required the required Fortran and C compilers.

Spack
`````

The easiest way to install these prerequisites is with `Spack
<https://spack.io>`_ (Linux or MacOS).  If you need to install Spack,
do:

.. code-block:: bash

   git clone https://github.com/spack/spack.git
   . spack/share/spack/setup-env.sh

Now you can build the packages:

.. code-block:: bash

   spack install netcdf-fortran
   spack install python@3:
   spack install py-netcdf^python@3:
   spack install cmake
   spack load netcdf-fortran
   spack load netcdf
   spack load python@3:
   spack load cmake

If needed, locations can be identified as follows:

.. code-block:: bash

   NETCDF4_FORTRAN_ROOT=`spack location -i netcdf-fortran`
   NETCDF4_C_ROOT=`spack location -i netcdf`

MacPorts
````````

If you're on a Macintosh, MacPorts may be used as an alternative to
Spack.  If you need to install MacPorts, follow instructions for
`user-level MacPorts <https://github.com/citibeth/usermacports>`_.
Then:

.. code-block:: bash

   port install python37
   port select -- set python3 python37
   port install py37-netcdf4
   port install netcdf-fortran
   port install cmake

Locations can be identified as follows:

   NETCDF4_FORTRAN_ROOT=`which port`/../..
   NETCDF4_C_ROOT=$NETCDF4_FORTRAN_ROOT


Run CMake
---------

CMake works by a two-step process:

 1. Run ``cmake``, which creates a ``Makefile``.

 2. Run ``make install`` to build and install the project.

In Step 1, CMake identifies the location of your project's
dependencies --- in this case, various NetCDF4 libraries.  CMake
locates dependencies in a variety of places: standard system
locations, paths specified by environment variables, etc.  In the best
of worlds, CMake is "automagically" able to find and choose the
verions of the dependencies that you wanted.  This is accomplished as
follows:

.. code-block:: bash

   mkdir build; cd build
   cmake .. -DCMAKE_INSTALL_PREFIX:PATH=`pwd`

If this step succeeds, CMake generates a ``Makefile``.  However, it is
still important to verify that CMake found the dependencies you were
hoping it would find.  Display these paths with:

.. code-block:: bash

   $ egrep 'CMAKE_Fortran_COMPILER:|NETCDF4' CMakeCache.txt 
   CMAKE_Fortran_COMPILER:FILEPATH=/home2/rpfische/spack-tools/opt/spack/linux-x86_64/gcc-4.8.5/gcc-4.9.3-jfebnuuusdch34j7pvdnvlxffe2rmoe4/bin/gfortran
   NETCDF4_C_BINARY_DIR:PATH=/home2/rpfische/spack6/opt/spack/linux-centos7-x86_64/gcc-4.9.3/netcdf-4.4.0-7hecfhzw4sj7wj2h5ioxmiv7dxvpcjeh/bin
   NETCDF4_C_INCLUDE_DIR:PATH=/home2/rpfische/spack6/opt/spack/linux-centos7-x86_64/gcc-4.9.3/netcdf-4.4.0-7hecfhzw4sj7wj2h5ioxmiv7dxvpcjeh/include
   NETCDF4_C_LIBRARY:FILEPATH=/home2/rpfische/spack6/opt/spack/linux-centos7-x86_64/gcc-4.9.3/netcdf-4.4.0-7hecfhzw4sj7wj2h5ioxmiv7dxvpcjeh/lib/libnetcdf.so
   NETCDF4_FORTRAN_INCLUDE_DIR:PATH=/home2/rpfische/spack6/opt/spack/linux-centos7-x86_64/gcc-4.9.3/netcdf-fortran-4.4.4-p2cmykx3iwkc2tqa6reuih75t4iysbuc/include
   NETCDF4_FORTRAN_LIBRARY:FILEPATH=/home2/rpfische/spack6/opt/spack/linux-centos7-x86_64/gcc-4.9.3/netcdf-fortran-4.4.4-p2cmykx3iwkc2tqa6reuih75t4iysbuc/lib/libnetcdff.so

If CMake did *not* find the dependencies you were hoping it would
find, this needs to be addressed by telling CMake where to find them
(you must fill in the ``...`` below:

.. code-block:: bash

   cmake .. -DCMAKE_INSTALL_PREFIX:PATH=`pwd`
   FC=`which gfortran` cmake .. -DCMAKE_INSTALL_PREFIX:PATH=. -DNETCDF4_C_ROOT=<...> -DNETCDF4_FORTRAN_ROOT=<...>

CMake is a standard and widely-used system, with abundant on-line
documentation and help, and every CMake-based package works about the
same.  See `here <https://cmake.org>`_ for more information.

Build And Install
-----------------

Once the ``Makefile`` has been generated, it is time to build and
install the EntGVSD supportin code.  This is done by:

.. code-block:: bash

   make
   make install

The final version of the support code will be installed in your
``build/bin``, ``build/lib`` and ``build/include`` directories.  The
script ``build/bin/entgvsd``, which compiles and immediately runs any
of the Fortran-based EntGVSD programs, is the only entry point to this
installation required by the user.

The installation can now be tested by running the EntGVSD test program:

.. code-block:: bash

   $ bin/entgvsd ../bnu/test.f90
   Test program succeeded!


Install xent
------------

The version of the ``bin/entgvsd`` script required is linked to the
version of the Fortran programs being run in the ``bnu`` directory.
If multiple versions of the EntGVSD source are checked out, it is
important to run the correct version of ``entgvsd`` for a given
program in ``bnu``.  This onerous task can be automated using the
``xent`` script, also installed in ``build/bin``.  This can be tested
as follows:

.. code-block:: bash

   $ bin/xent ../bnu/test.f90
   Test program succeeded!

A single ``xent`` script works for any version of EntGVSD; it works by
identifying the correct version of ``entgvsd`` to run based on the
location of the Fortran program being run.  One can therefore copy the
``xent`` program (once) into a directory in the user's ``$PATH``, and
then use it for all EntGVSD programs in any source checkout.
