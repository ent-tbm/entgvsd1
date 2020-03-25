.. _spackenv

Installing the EntGVSD1 Spack Environment
=========================================

These instructions outline how to install the prerequisites for EntGVSD1 on
*discover* (or any computer), starting from scratch.  For the most
part, command lines that can be cut-n-pasted into *discover* are used.
This tutorial assumes the use of the *bash* command shell.



Path Setup
----------

We assume the following paths.  The software may be installed in any
place; and the user may set these paths however they like.  These
environment variables do NOT need to be set to use this software,
they are merely set for the purposes of this tutorial.  However, if
you set these variables as appropriate in your shell now, you can
cut-n-paste from later parts of the tutorial.

.. code-block:: bash

   # Location of shared Spack instance containing installed
   # environment
   SPACK=~eafisch2/spack

   # Use this Environment spec within $SPACK
   # The environment needs to be tailored to your HPC system.
   SPENV=entgvsd-discover12

   # Your username on the simplex git server
   SIMPLEX_USER=$USER

Build the Spack Environment
---------------------------

This builds the prerequisite software.  It only needs to be done once
by any user on the system; and shared by everyone else.

Clone Spack
```````````

Download the for of Spack containing tested recipes to build software
needed for ModelE-PISM coupling:

.. code-block:: bash

   cd $(dirname $SPACK)
   git clone git@github.com:citibeth/spack.git -b efischer/giss2 $(basename $SPACK))
   # Don't worry about errors on this command
   source $SPACK/var/spack/environments/$SPENV/loads-x

Spack Compiler Configuration
````````````````````````````

Spack compiler definitions are in ``~/..spack/linux/compilers.yaml``.
Follow directions on `Compiler Configuration in the Spack Manual
<https://spack.readthedocs.io/en/latest/getting_started.html#compiler-configuration>`_.
Here is a sample used on *discover*.

.. code-block:: none
   # SLES12
   compilers:
   - compiler:
       environment: {}
       extra_rpaths: [/usr/local/intel/2020/compilers_and_libraries_2020.0.166/linux/compiler/lib/intel64_lin,/usr/local/intel/2020/compilers_and_libraries_2020.0.166/linux/mpi/intel64/lib]
       flags:
          cflags: -gcc-name=/usr/local/other/gcc/9.2.0/bin/gcc
          cxxflags: -gxx-name=/usr/local/other/gcc/9.2.0/bin/g++
       modules: []
       operating_system: suse_linux11
       paths:
         # https://software.intel.com/en-us/node/522750 
         cc: /usr/local/intel/2020/compilers_and_libraries_2020.0.166/linux/bin/intel64/icc
         cxx: /usr/local/intel/2020/compilers_and_libraries_2020.0.166/linux/bin/intel64/icpc
         f77: /usr/local/intel/2020/compilers_and_libraries_2020.0.166/linux/bin/intel64/ifort
         fc: /usr/local/intel/2020/compilers_and_libraries_2020.0.166/linux/bin/intel64/ifort
       spec: intel@20.0.166
       target: x86_64

.. note::

   The Intel compiler must be paird with an appropraite version of
   GCC.  See, for example, this `Discussion on Spack
   <https://github.com/spack/spack/issues/8356>`_.

Once you have configured your compiler, try building something simple.

.. code-block:: 

    $SPACK/bin/spack --no-env install zlib

If this works --- congratulations, you are ready to go with Spack!  If not... please contact the Spack community for help.  Here's what you will see if it worked:

.. code-block:: none

   do_install <spack.pkg.builtin.zlib.Zlib object at 0x7ffff185d190> None install_deps=True
   Installing zlib dependencies
   deps = []
   ==> Installing zlib
   ==> Searching for binary cache of zlib
   ==> Warning: No Spack mirrors are currently configured
   ==> No binary for zlib found: installing from source
   ==> Fetching http://zlib.net/fossils/zlib-1.2.11.tar.gz
   ######################################################################## 100.0%
   ==> Staging archive: /gpfsm/dnb53/eafisch2/spack/var/spack/stage/zlib-1.2.11-lk267u47ez67rkzl7z5gnrdqvhca2n46/zlib-1.2.11.tar.gz
   ==> Created stage in /gpfsm/dnb53/eafisch2/spack/var/spack/stage/zlib-1.2.11-lk267u47ez67rkzl7z5gnrdqvhca2n46
   ==> No patches needed for zlib
   ==> Building zlib [Package]
   ==> Executing phase: 'install'
   ==> Successfully installed zlib
     Fetch: 0.76s.  Build: 4.40s.  Total: 5.15s.
   [+] /gpfsm/dnb53/eafisch2/spack/opt/spack/linux-sles12-x86_64/intel-20.0.166/zlib-1.2.11-lk267u47ez67rkzl7z5gnrdqvhca2n46


Use Spack to Build Environement
```````````````````````````````

.. code-block:: bash

   $SPACK/bin/spack -e $SPENV concretize -f
   $SPACK/bin/spack -e $SPENV install
   cd $SPACK/var/spack/environments/$SPENV
   $SPACK/bin/spack -e $SPENV env loads -r
   sort loads | uniq >loads2
   cp loads2 loads

.. note::

   The Spack environment entgvsd-discover12 (``$SPENV``) is meant to
   work on the *NCCS Discover* supercomputer, SLES12 version.  If this
   is being built on another system, then that environent should be
   copied, modified as appropriate for that system, checked in and
   submitted as a pull request.  Further details are out of the scope
   of this document; see `Spack Environments
   <https://spack.readthedocs.io/en/latest/environments.html>`_:

   .. code-block:: bash

      cd $SPACK/var/spack/environments
      cp -r twoway-discover twoway-mything
      nano twoway-mything/spack.yaml

Make Sure Spack is World Readable
`````````````````````````````````

When you are done building the prerequisite software, it is polite to
make it world readable for everyone, so others can use it too:

.. code-block:: bash

   chmod -R a+r $SPACK


Load the Spack Environment
``````````````````````````

The EntGVSD Spack environment is now built; and now needs to be loaded
upon start of a shell session.  That is done with the command:

.. code-block:: bash

    source $SPACK/var/spack/environments/$SPENV
