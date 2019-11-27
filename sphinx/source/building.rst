.. _building

Building EntGVSD
================

Install Prerequisites
---------------------

EntGVSD requires the following prerequisites:

* GNU Fortran / C / C++ version 4.9 or later.

* `NetCDF-C <https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html>`_, version 4 or later

* `NetCDF-Fortran <https://www.unidata.ucar.edu/software/netcdf/docs/building_netcdf_fortran.html>`_, version 4 or later

* (Optional) `Python <https://www.python.org>`_ version 3.5 or later

There are a variety of ways these prerequisites can be installed.
Some suggestions are provided below; assuming you have already
required the required Fortran and C compilers.

Spack
`````

The easiest way to install these prerequisites is with `Spack
<https://spack.io>`_ (Linux or MacOS), using a pre-configured Spack Environment.


.. code-block:: bash

   cd ~
   git clone git@github.com:citibeth/spack.git -b efischer/giss2
   . spack/share/spack/setup-env.sh
   spack -e entgvsd-gibbs concretize -f
   spack -e entgvsd-gibbs install
   spack -e entgvsd-gibbs env loads

These commands build all the prerequisites needed for working with EntGVSD, and install them as environment modules.  Those environment modules may be loaded with:

.. code-block:: bash

   source ~/spack/var/spack/environments/entgvsd-gibbs/loads-x

.. note::

   1. Spack can be cloned into any location.  From here on, we will
      assume WLOG it has been installed in ``~/spack``.

   1. If you have trouble installing prerequisites with Spack, *please*
      and ask questions on the `Spack discussion
      group<https://groups.google.com/forum/#!forum/spack>`.  This will
      get you better help, faster, than contacting the EntGVSD authors
      directly.


Download EntGVSD Source
-----------------------

Download the EntGVSD source from the Simplex git server.  One must be
connected to the NASA network for this to work (replace
``<ndcusername>`` with your NDC username):

.. code-block:: bash

   cd ~/git
   git clone <ndcusername>@simplex.giss.nasa.gov:/giss/gitrepo/entgvsd1.git -b develop

.. note::

   1. EntGVSD can be cloned into any location.  From here on, we will
      assume WLOG it has been installed in ``~/git/entgvsd``.

   1. To gain access to Simplex, contact Igor Aleinov
      *igor.aleinov@nasa.gov*.


Build EntGVSD Library
---------------------

EntGVSD is structured as a set of Fortran programs that can be run
from the command line, to perform each step in the process.  These
programs, in the ``src/`` directory, can be thought of as "scripts"
because each one is self-contained in a single source file; and as
with a scripting language, they can be edited and run immediately,
without explicit compilation.

Supporting the scripts is the EntGVSD library, which is built and
"installed" within the EntGVSD directory structure.  This is built as
follows:

.. code-block:: bash

   cd ~/git/entgvsd1
   mkdir build
   cd build
   FC=`which gfortran` cmake .. -DCMAKE_INSTALL_PREFIX:PATH=`pwd`
   make install

.. admonition:: OPTIONAL:

   The ``xent`` may be used to conviently launch Fortran scritpts.  It
   should be added to your ``.bashrc`` file as follows:

   .. code-block:: bash

      export PATH=$PATH:~/git/entgvsd1/build/bin

   Alternately, you can just copy it to an existing directory in your
   ``PATH`` (eg ``~/sh``):

   .. code-block:: bash

      cp ~/git/entgvsd1/build/bin/xent ~/sh


   .. note::

      If you clone EntGVSD more than once, you still only need one
      copy of ``xent``, they are all the same.




bin/xent ../bnu/test.f90
cd ../bnu
./run_all.sh


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

   FC=`which gfortran` cmake .. -DCMAKE_INSTALL_PREFIX:PATH=`pwd` -DNETCDF4_C_ROOT=<...> -DNETCDF4_FORTRAN_ROOT=<...>

CMake is a standard and widely-used system, with abundant on-line
documentation and help, and every CMake-based package works about the
same.  See `here <https://cmake.org>`_ for more information.

.. note::

   The following command worked when using MacPorts and GCC 7:

   .. code-block:: bash

      FC=gfortran-mp-7 cmake .. -DCMAKE_INSTALL_PREFIX:PATH=`pwd`


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
