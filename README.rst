******************************************
Ent Global Vegetation Structure Dataset v1
******************************************

User Manual
===========

See the `full user manual <https://entgvsd.readthedocs.io>`_ for more details on using this project.

Introduction
============
The *Ent Global Vegetation Structure Dataset (Ent GVSD) v1* provides satellite- and climate data-derived global vegetation structure for boundary conditions and evaluation data sets for dynamic global vegetation models (DGVMs), tailored to the Ent Terrestrial Biosphere Model (Ent TBM) (Kiang et al 2006, Kim et al 2015, Ni-Meister et al 2010, Yang et al 2010), and formulated to be of general use for for demographic dynamic global vegetation models (dDGVMs).  The data fields provided are gridded maps of: 
1. Subgrid fractions of land cover.  Land cover types can be bare soil and patches of vegetation communities composed of identical individual plants, which can be one of the Ent TBM's 12-13 plant functional types (PFTs).
1. Height of each patch community.
1. Maximum annual leaf area index (LAI) of each community.
1. Temporal values of LAI of each community.
1. Spectral soil albedo for advanced land albedo calculations.  

The vegetation communities are comprised of the Ent TBM's 12-13 plant functional types (PFTs) described at :ref:`EntPFTs`  

These boundary conditions when used together can be used to derive plant population density given allometric relations between height and maximum leaf area per plant for each PFT, as rationalized according to pipe model theory {Shinozaki, 1964 #708}. Then temporally varying LAI can be prescribed for simulations of vegetation seasonal dynamics.  Canopy radiative transfer models that account for soil albedo underlying the canopy additionally can utilize this vegetation structure with the soil albedo data set.  

Versions of the Ent GVSD differ in spatial and temporal resolutions, levels of community heterogeneity, and input data sources.
:ref:`v1.0`:  Demonstration version derived at 0.5 x 0.5 degrees.  LAI: Moderate Resolution Imaging Spectroradiometer (MODIS) MCD12Q1 V005 L3 data product, year 2004.  Soil albedo: grey soil albedos for bare soil from Matthews (1983). In use in the GISS GCM ModelE2.
:ref:`v1.1`:  Primary version derived at the native resolutions of satellite observations. LAI: Beijing Normal University (BNU) (Yuan et al 2011) post-processing of the MODIS LAI product.  Soil albedo: MODIS-derived spectral soil albedo from Carrer et al. (2014).
:ref:`v1.2`:  Higher temporal resolution LAI for multi-harvest crops.  LAI: LAI3g (Zhu et al. 2013) post-processing of the MODIS LAI product.

A *NASA Technical Publication (####)* describes in detail the data sources and procedures to create the Ent GVSD Versions 1.0, 1.1, and 1.2 (v1.X).

The codes released here and the documented source data enable the users to regenerate these versions of the Ent GVSD.

Requirements
============

Main programs are written in Fortran95, and **ELIZABETH WILL WRITE THIS SECTION**.  Script in R are provided for automated plotting of maps of all steps in generation of the Ent GVSD v1.1

HOW-TO get the code
=====================
get_code


Ent GVSD v1.1 Files: ``/``
=========================

``CMakeLists.txt``, ``cmake/``
  The CMake build system, generates Makefile(s).  

``bnu/``
  Contains all main fortran programs and output directory lc_lai_ent

``slib/``
  Fortran modules and subroutines shared by the programs in ``bnu/``

``scripts/``
  OLD PRELIM. Elizabeth’s unzipping and organizing of Carlo’s code.  **REMOVE**

``sphinx/``
  Restructured Text documentation files.  See `Sphinx Document Generator <https://www.sphinx-doc.org/en/master/>`_

``inputs/``
  Pre-processed data files for input to the main programs. This directory is generated via a script to link to the data files downloaded separately.  


Files: ``bnu/``
---------------

``run_all.sh``
  Main script to run the entire process.  Stops if any step encounters an error.

``makefile``
  Not used YET , run_all.sh generates the .mk files that list the
  input/output files.  Idea was that then this makefile could be run
  to run those .mk files for changing input/output but hasn’t been
  tested, yet.  MAIN.sh - OLD stuff



