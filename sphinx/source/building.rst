.. _building

Build EntGVSD
=============

Install Prerequisites
---------------------

EntGVSD requires the following prerequisites:

* GNU Fortran / C / C++ version 4.9 or later.  Also works with Intel Compiler version 18.

* `NetCDF-C <https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html>`_, version 4 or later

* `NetCDF-Fortran <https://www.unidata.ucar.edu/software/netcdf/docs/building_netcdf_fortran.html>`_, version 4 or later

* `CMake <https://cmake.org>`_, version 3.0 or later

* `cURL <https://curl.haxx.se>`_ command line utility

* (Optional) `Python <https://www.python.org>`_ version 3.5 or later.  It should have at least the packages `netCDF4-python <https://unidata.github.io/netcdf4-python/netCDF4/index.html>`_ installed.  This is only needed to convert final files to the GISS ModelE format.

* (Optional) `R <https://www.r-project.org>`_.  It should have at least the packages *sp*, *maps*, *sdmtools*, *plotrix*, *spam*, *fields*, *maptools*, *rworldmap* and *netcdf*.  See `CRAN <https://cran.r-project.org>`_ for more details.  R is used to plot the results in an easy-to-check format.

There are a variety of ways these prerequisites can be installed.  For
example, they may be installed by hand.  Some suggestions are provided
below; assuming you have already required the required Fortran and C
compilers.

Spack
`````

One way to install these computing prerequisites, along with useful
utilities, is through `Spack <https://spack.io>`_ (Linux or MacOS).
See :ref:`spackenv` for further details.

Download EntGVSD Source
-----------------------

Download the EntGVSD source from GitHub:

.. code-block:: bash

   cd ~/git
   git clone git@github.com:ent-tbm/entgvsd1.git -b develop

.. note::

   To use GitHub, you need to generate a public/private SSH keypair on
   each machine on which you will use git.  If you have not done so
   already, follow these `instructions
   <https://help.github.com/en/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account>`_
   to generate a keypair and install it to your *GitHub* account.



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
   cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$(pwd)
   make install

.. admonition:: OPTIONAL:

   The ``xent`` may be used to conviently launch Fortran scripts.  It
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


Fetch Input Data and Create the Makefile
----------------------------------------

The main EntGVSD process is structured as a series of Fortran scripts,
to be run in order, starting with the capital letter `B`.  For example:

.. code-block:: bash

   B01_bnu_laimax.F90
   B02_lc_modis_entpftrevcrop.F90
   B03_regrid_snowice.F90
   B04_veg_height.F90
   B05_carrer_mean.F90
   B06_albmodis_gridfill.F90
   B07_soil_albedo.F90
   B08_lc_laimax.F90
   B09_lc_lai_doy.F90
   B10_lc_lai_monthly.F90
   B11_reclass_annual.F90
   B12_reclass_doy.F90
   B13_reclass_monthly.F90
   B14_regrid.F90
   B15_regrid_controls.F90
   B16_trim.F90
   B17_checksum.F90
   B18_modele.F90


The EntGVSD creates a Makefile to run these in sequence.  To download
all necessary input data and create the Makefile, run the *mkgen*
script.

.. code-block:: bash

   cd ~/git/entgvsd1/src
   ./mkgen


.. note::

   1. Downloading input files can take a while; and can also get
      stuck, depending on the condition of the network and NCCS.

   1. The input data files and their subdirectory structures used to
      produce the Ent GVSD, are mirrored at the `NCCS Data Portal
      <https://portal.nccs.nasa.gov/datashare/GISS/Ent_TBM/EntGVSD/inputs/>'_.

   1.  The input files are not automatically downloaded with a git
       clone of the code, due to their size.  These are pre-processed
       data files that are read by the ``B*.F90`` fortran programs that
       generate the Ent GVSD. The *mkgen* script downloads the input
       files to their correct directories in your EntGVSD clone and
       also avoids repeating if previously downloaded.

   1. Input files are stored in compressed (gzip) form on the
      dataportal, and are uncompressed immediately after downloading.
      Uncompressed files can be markedly larger than their compressed
      form, sometimes up to 50-100X.

   1. *mkgen* may take a long time, due to downloading the files.
      If it is stopped in the middle, simply restart it agian.

   1. In addition to downloading datafiles, the *mkgen* script
      generates dependency files in the ``mkfiles/`` directory, which
      indicate the input and ouput files of each EntGVSD script.
      These are not used for the ``Makefile``; however they offer a
      definitive reference of what files each step uses and produces.

Run EntGVSD
============

Once EntGVSD has been built, the Fortran programs can be run with:

.. code-block:: bash

   cd ~/git/entgvsd1/src
   make

This will run the steps, in order, and is expected to take a few days.
If you alter a Fortran script in the `src/` directory, recompilation is not necessary.  However, if you alter any code in the `slib/` directory, be sure to recompile by repeating the *make install* command:

.. code-block:: bash

   cd ~/git/entgvsd1/build
   make install

In order to force rerun of a step ; say, step ``B01_bnu_laimax``, do:

.. code-block:: bash

   cd ~/git/entgvsd1/src
   rm ../outputs/B01_bnu_laimax.txt
   make

.. note::

   This will rerun the desired step, plus all subsequent steps, which
   are assumed to depend on all previous steps.

To run a single program at a time, such as B11_reclass.F90:

.. code-block:: bash

   cd ~/git/entgvsd1/src
   ../build/bin/entgvsd B11_reclass.F90


Input / Output Records
----------------------

Each step of EntGVSD, when it runs, writes out a file ending in
``.mk``, which details the input and ouptut files used by that
program.  These ``.mk`` files are written twice:

1. When `mkgen` is run, they are written in the `mkfiles/` directory.

1. When the programs are run for real, they are written again, in the
   `outputs/` directory.

Looking in these ``.mk`` files is useful to give a definitive answer
on what files each program opens.


Modifying Parameters
====================

User-editable parameters are in the file ``slib/ent_params.f90``.
Once parameter(s) in this file are changed, the following steps must
take place to make sure they take effect:

.. code-block:: bash

   cd ~/git/entgvsd1/build
   make install
   

.. note::

   1. The ``ent_params.f90`` file is NOT checked into git.  It is a
      user configuration file.

   1. To revert to default values as stored in git, do:

      .. code-block:: bash

         cd ~/git/entgvsd1/slib
         rm ent_params.f90
         cd ../build
         FC=$(which gfortran) cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$(pwd)

   1. The parameters ``LAI_YEAR`` and ``sLAI_YEAR`` must match.  One
      is a string, one is an integer.

   1. Changing the ``LAI_YEAR`` parameter will cause ``2004`` to be
      replaced by a different year, everywhere it is appropriate in
      input filenames, output filenames, metadata and folders ---
      except for ``B20_plots.R``, where the year must be changed manually.

Rerun EntGVSD
=============

If EntGVSD has already run and you wish to re-run it with a "clean"
slate, the following steps are will do so:

.. code-block:: bash

   cd ~/git/entgvsd1
   rm -rf outputs build
   mkdir build
   cd build
   FC=$(which gfortran) cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$(pwd)
   make install
   cd ../src
   ./mkgen
   make

.. note::

   As long as the downloaded data files in the `inputs/` directory are
   not deleted, this procedure will not need to re-download them.

   
Modifying Parameters
====================

User-editable parameters are in the file ``slib/ent_params.f90``.
Once parameter(s) in this file are changed, the following steps must
take place to make sure they take effect:

.. code-block:: bash

   cd ~/git/entgvsd1/build
   make install
   

.. note::

   1. The ``ent_params.f90`` file is NOT checked into git.  It is a
      user configuration file.

   1. The parameters ``LAI_YEAR`` and ``sLAI_YEAR`` must match.  One
      is a string, one is an integer.

   1. Changing the ``LAI_YEAR`` parameter will cause ``2004`` to be
      replaced by a different year, everywhere it is appropriate in
      input filenames, output filenames, metadata and folders ---
      except for ``B20_plots.R``, where the year must be changed manually.

Pre-Processsed Raw Data Files
============================

Code to pre-process original source data files (many of which serve as
input to EntGVSD) are in the ``data/`` directory, created and
downloaded by the *mkgen* script.  These codes have been run
previously and their output pre-processed files are provided; but
unlike the scripts in ``src/``, the codes do not come with a curated
build system.  They are provided as-is, for reference.

Accompanying the code are a number of data files from the original data sources.  
They may be downloaded by running the ``entdata'' script in each subdirectory of ``data/``.  For example:

.. code-block:: bash

   cd ~/git/entgvsd1/data/climstats
   ./entdata

The contents of the data directory are described here.

.. note::

   **TODO**: Add link to new page named data.rst to describe the data
    directory

