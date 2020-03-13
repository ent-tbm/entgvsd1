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

An option to install these computing prerequisites is through `Spack
<https://spack.io>`_ (Linux or MacOS). A pre-configured Spack Environment for the EntGVSD v1 can be set up through these installation instructions:


.. code-block:: bash

   cd ~ #cd to your home directory.
   bash #Use the bash shell.
   git clone https://github.com/nkianggiss/spack.git -b nasasgiss/entgvsd1 #Clone the Spack branch developed for the Ent GVSD v1.
   source ~/spack/share/spack/setup-env.sh         #Set up the environment variables in the user's shell.

The needed software packages are listed in files named "spack.yaml" here:

.. code-block:: bash

~/spack/var/spack/environments/environments/<environment name>/spack.yaml

where <environment name> is a subdirectory for a particular known computing environment.  The spack.yaml files for each <environment name> differ only in an optional include section for customizable yaml files that specify  preferred versions of the software packages. Customized yaml files are put in this directory, which contains example configurations:

.. code-block:: bash

~/spack/var/spack/environments/environments/config/

Otherwise, the default is to download the most recent versions.

Select your <environment name> to use, and continue setting up your environment:

.. code-block:: bash

   spack -e <entgvsd-environment> concretize -f    #Discern what software packages are needed and their dependencies.
   spack -e <entgvsd-environment> install          #Download and build the software packages.
   spack -e <entgvsd-environment> env loads        #Generate a file "loads". This is a list of modules to load.

The "loads" file may have duplicates that slow downloading them.  Commands to eliminate duplicates in the loads files are:

.. code-block:: bash

   cd ~/spack/var/spack/environments/entgvsd-gibbs
   sort loads | uniq >loads2
   cp loads2 loads

The above commands build all the prerequisites needed for working with EntGVSD, and install them as environment modules.  To load the environment modules, the user should make a version of the below environment file, "loads-x".  This environment file performs the functions that might be in the user’s .profile or other environment files, including:

1.  Purges modules, loads modules
2.  Sets environment paths such as optionally python path, library path (avoiding the need to edit a .profile or .bashrc file).
3.  Sets FC environment variable to select which fortran to use.

After putting your loads-x into your spack environment directory, then load the environment:

.. code-block:: bash

   source ~/spack/var/spack/environments/<environment name>/loads-x

.. note::

   1. Spack can be cloned into any location.  From here on, we will
      assume without loss of generality it has been installed in ``~/spack``.

   1. If you have trouble installing prerequisites with Spack, *please*
      and ask questions on the `Spack discussion
      group<https://groups.google.com/forum/#!forum/spack>`_.  This will
      get you better help, faster, than contacting the EntGVSD authors
      directly.


Download EntGVSD Source
-----------------------

Download the EntGVSD source from the Simplex git server.  

If you are outside the NASA network, download a snapshot of the code from:
* link TBA on NASA-approved git site.


If you are inside the NASA network, (replace
``<ndcusername>`` with your NDC username):

.. code-block:: bash

   cd ~/git
   git clone <ndcusername>@simplex.giss.nasa.gov:/giss/gitrepo/entgvsd1.git -b develop

.. note::

   1. EntGVSD can be cloned into any location.  From here on, we will
      assume without loss of generality it has been installed in ``~/git/entgvsd``.

   2. To gain access to Simplex, contact Igor Aleinov
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
   FC=$(which gfortran) cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$(pwd)
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


Fetch Input Data and Create the Makefile
----------------------------------------

The main EntGVSD process is structured as a series of Fortran scripts,
to be run in order, starting with the capital letter `B`.  For example:
| B01_bnu_laimax.F90
| B02_lc_modis_entpftrevcrop.F90
|  ...


The EntGVSD creates a Makefile to run these in sequence.  To download all necessary input data and
create the Makefile, run the ``mkgen`` script.

.. code-block:: bash

   cd ~/git/entgvsd1/src
   ./mkgen

Downloading input files can take a while; and can also get stuck, depending on the condition of 
the network and NCCS.

.. note::

   1. The input data files and their subdirectory structures used to produce the Ent GVSD, are mirrored at 
      the 'NCCS Data Portal. 
      <https://portal.nccs.nasa.gov/datashare/GISS/Ent_TBM/EntGVSD/inputs/>'_

   2.  The input files are not automatically downloaded with a git clone of the code, due to their size.  
       These are pre-processed data files that are read by the B*.F90 fortran programs that generate the 
       Ent GVSD. The ``mkgen`` script downloads the input files to their correct directories in your 
       EntGVSD clone and also avoids repeating if previously downloaded. 

   3. Input files are stored in compressed form on the dataportal
      (gzip format), and are uncompressed immediately after
      downloading.  Uncompressed files can be markedly larger than
      their compressed form, sometimes up to 50-100X.

   4. ``mkgen`` may take a long time, due to downloading the files.
      If it is stopped in the middle, simple restart it agian.

   5. In addition to downloading datafiles, the ``mkgen`` script
      generates dependency files in the ``mkfiles/`` directory, which
      indicate the input and ouput files of each EntGVSD script.
      These are not used for the ``Makefile``.

Run EntGVSD
============

Once EntGVSD has been built, the fortran programs can be run, with simply:

.. code-block:: bash

   cd ~/git/entgvsd1/src
   make

This will run the steps, in order, and is expected to take a few days.
In order to force rerun of a step ; say, step ``B01_bnu_laimax``, do:

.. code-block:: bash

   cd ~/git/entgvsd1/src
   rm ../outputs/B01_bnu_laimax.txt
   make

.. note::

   This will rerun the desired step, plus all subsequent steps (which
   are assumed to depend on all previous steps).

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

Code to pre-process original source data files (many of which serve as input to EntGVSD)
are in the ``data/`` directory, created and downloaded by the ``mkgen`` script.  These codes 
have been run previously and their output pre-processed files are provided; but unlike the 
scripts in ``src/``, the codes do not come with a
curated build system.  They are provided as-is, for reference.

Accompanying the code are a number of data files from the original data sources.  
They may be downloaded by running the ``entdata'' script in each subdirectory of ``data/``.  For example:

.. code-block:: bash

   cd ~/git/entgvsd1/data/climstats
   ./entdata

The contents of the data directory are described here.
##Add link to new page named data.rst to describe the data directory ##

