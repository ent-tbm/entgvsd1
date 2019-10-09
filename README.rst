

Files: ``/``
========

``CMakeLists.txt``, ``cmake/``
  The CMake build system, generates Makefile(s).  

``bnu/``
  Contains all main fortran programs and output directory lc_lai_ent

``scripts/``
  OLD PRELIM. Elizabeth’s unzipping and organizing of Carlo’s code.		

``sphinx/``
  Restructured Text documentation files.  See `Sphinx Document Generator <https://www.sphinx-doc.org/en/master/>`_

``inputs_manual/``
  Pre-process carrer/ and sw_albedo/ files for use.  The result is input files for the programs.

``slib/``
  Fortran modules and subroutines shared by the programs in ``bnu/``

Files: ``bnu/``
---------------

``run_all.sh``
  Main script to run the entire process.  Stops if any step encounters an error.

``makefile``
  Not used YET , run_all.sh generates the .mk files that list the
  input/output files.  Idea was that then this makefile could be run
  to run those .mk files for changing input/output but hasn’t been
  tested, yet.  MAIN.sh - OLD stuff


Files: OBSOLETE
---------------

``entgvsd_root.txt``
  Not used. OLD.
``README.md``
  github front page in “markdown” format.
``templates/``
  netcdf .cdl template files.  Now are generated in code, so not needed any more.
``bnu/lc_lai_ent/*.mk``
  .mk files are not really Makefiles, but lists the inputs and output files for each program.
