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
      assume WLOG it has been installed in ``~/spack``.

   1. If you have trouble installing prerequisites with Spack, *please*
      and ask questions on the `Spack discussion
      group<https://groups.google.com/forum/#!forum/spack>`_.  This will
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


Create the Makefile
-------------------

The main EntGVSD process is structured as a series of Fortran scripts,
to be run in order, starting with the capital letter `B`.  For example:
| B01_bnu_laimax.F90
| B02_lc_laimax_modis_entpftrevcrop.F90
|  ...

Before these scripts can be run, EntGVSD creates a Makefile to run
them in sequence.  This is done using the ``mkgen`` script.

.. code-block:: bash

   cd ~/git/entgvsd1/src
   ./mkgen

This step also downloads all required input files for EntGVSD from the
NCCS Dataportal.  Downloading input files can take a while; and can
also get stuck, depending on the condition of the network and NCCS.

.. note::

   1. The input files are stored in the EntGVSD inputs directory of
      the `NCCS Dataportal at
      <https://portal.nccs.nasa.gov/datashare/GISS/Ent_TBM/EntGVSD/inputs>`_

   1. Input files are stored in compressed form on the dataportal
      (gzip format), and are uncompressed immediately after
      downloading.  Uncompressed files can be markedly larger than
      their compressed form, sometimes up to 50-100X.

   1. ``mkgen`` may take a long time, due to downloading the files.
      If it is stopped in the middle, simple restart it agian.

   1. In addition to downloading datafiles, the ``mkgen`` script
      generates dependency files in the ``mkfiles/`` directory, which
      indicate the input and ouput files of each EntGVSD script.
      These are not used for the ``Makefile``.

Run EntGVSD
============

Once EntGVSD has been built, it can be run, with simply:

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


Data Files
==========

Code to produce data files (many of which serve as input to EntGVSD)
are in the ``data/`` directory.  These codes have been run in the
past; but unlike the scripts in ``src/``, they do not come with a
curated build system.  They are provided as-is, for reference.

Accompanying the code are a number of data files received by our
group.  They may be downloaded by running the ``entdata'' script in
each subdirectory of ``data/``.  For example:

.. code-block:: bash

   cd ~/git/entgvsd1/data/climstats
   ./entdata

