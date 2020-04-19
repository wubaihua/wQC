# wQC
A simple Fortran Quantum Chemistry Program

Last-update: 2020-4-13

Author: Baihua Wu (wubaihua@pku.edu.cn)



## Introduction

wQC is a simple quantum chemistry program writing by fortran. Currently it is still in the development stage. The author will continue to develop it.



## Compile

Be sure you have installed Inter Fortran Compiler (ifort). 

You should install Libcint( https://github.com/sunqm/libcint ) according to the method provided by the website.

Just type `make`.



## Features currently implemented

* Restricted Hartree-Fock for closed-shell system
* Population Analysis
* MP2
* DIIS for SCF



## TODO

* Unrestricted HF
* DFT
* CCSD
* CAS Method
* TDHF and TDDFT



## Acknowledgments

* wQC using the Libcint (J. Compet. Chem. 2015, 36, 1664–1671) to calculate single and double electronic integral and build the overlap, kinetic and potential matrix. Some code in cint.f90 refers to St Maxwell's website (https://gensoukyo.me/use-libcint/). Author thanks for his contribution. Author also wants to thank Dr. Qiming Sun, the author of Libcint, for some advices about using Libcint.
* Some code in def.f90 refers to the source code of Multiwfn 3.6 (http://sobereva.com/multiwfn/). Multiwfn is a quantum chemistry wavefunction analysis program developed by Dr. Tian Lu. Author  thanks for his contribution.
* Author also wants to thank Prof. F. Chen from USTB, author’s Enlightenment Tutor in Quantum Chemistry, for his careful guidance to author for so long.  Author also wants to thank Prof. H. Jiang and Prof. J.Liu from PKU, for their excellent course in PKU.

