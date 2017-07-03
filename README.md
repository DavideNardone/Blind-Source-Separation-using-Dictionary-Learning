# Blind-Source-Separation-using-Dictionary-Learning

Blind source separation (BSS) is a technique for estimating individual source components from their mixtures at multiple sensors.
It is called blind because we don't use any other information besides the mixtures.

## Introduction

The sparse decomposition of images and audio signals found great use in the field of: Compression, Noise removal and also in the Sources separation.
This implies the decomposition of signals in the form of linear combinations with some elements of a redundant dictionary.
The dictionary may be either a fixed dictionary (Fourier, Wavelet, etc) or may be learned from a set of samples.
The algorithms based on learning the dictionary can be applied to a broad class of signals and have a better compression performance than methods based on fixed dictionary.

Here we present a Compressed Sensing (CS) approach with an adaptive dictionary for solving a **Determined Blind Source Separation (DBSS)**.
The proposed method has been developed by reformulating a DBSS as Sparse Coding (SC) problem.
The algorithm consist of few steps:
  1. Mixing matrix estimation;
  2. Sparse source separation and Source reconstruction;

A sparse mixture of the original source signals has been used for the estimation of the mixing matrix which has been used for the reconstruction of the source signals.
A 'block signal representation' is used for representing the mixture, in order to greatly improve the computation efficiency of the 'mixing matrix estimation' and the 'signal recovery' processes,
without any particular loss of accuracy in the separation.

For more information about BSS, please refer to: http://www.fon.hum.uva.nl/praat/manual/blind_source_separation.html

## Requisites and dependency

The software is mostly written in MATLAB with consistent C/C++ dependency. In order to make it works you must have installed:

  * A working gcc/g++ compiler;
  * A MATLAB version (tested with MATLAB 2016b)
  
This DBSS implementation depends in particular by **_SPArse Modeling Software_** (SPAM)[1] which is an optimization toolbox for solving various sparse estimation problems.
For more information, please refer to: http://spams-devel.gforge.inria.fr/

The performaces of the proposed algorithm have been evaluated by using BSSEVAL[2] library.

## Usage

In order to run this project, go into the project folder: `/sources` and run the file `main.m`. The latter, in addition to the main DBSS method also integrates two well known
BSS method 1) Fast-ICA and ERICA. We've added them into the project just to compare each other performance.

The `main.m` file contains some commented code-lines which have some different purpose such as: saving outputs or displaying different plots.


## Dataset

The dataset used for testing the algorithms are those for the Sixth Community-Based Signal Separation Evaluation Campaign,
SiSEC 2015. They are placed into the project folder: `/resources`. In the latter folder, other kind of datasets have been provided.

# Authors

Davide Nardone, University of Naples Parthenope, Science and Techonlogies Departement, Msc Applied Computer Science 
https://www.linkedin.com/in/davide-nardone-127428102

# Contacts

For any kind of problem, questions, ideas or suggestions, please don't esitate to contact me at: 
- **davide.nardone@studenti.uniparthenope.it**

# References

[1] SPArse Modeling Software (SPAM): http://spams-devel.gforge.inria.fr/
[2] BSSEVAL: http://bass-db.gforge.inria.fr/bss_eval/

