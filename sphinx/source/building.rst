.. _building

Build EntGVSD
=============

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
   git clone https://github.com/spack/spack.git -b efischer/giss2
   . spack/share/spack/setup-env.sh
   spack -e entgvsd-gibbs concretize -f
   spack -e entgvsd-gibbs install
   spack -e entgvsd-gibbs env loads

These commands build all the prerequisites needed for working with EntGVSD, and install them as environment modules.  Those environment modules may be loaded with:

.. code-block:: bash

   source ~/spack/var/spack/environments/entgvsd-gibbs/loads-x

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
| B02_lc_laimax_modis_entpftrevcrop.F90
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
In order to force rerun of a step; say, step ``B01_bnu_laimax``, do:

.. code-block:: bash

   cd ~/git/entgvsd1/src
   rm ../outputs/B01_bnu_laimax.txt
   make

This will rerun the desired step, plus all subsequent steps (which are
assumed to depend on all previous steps).


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

