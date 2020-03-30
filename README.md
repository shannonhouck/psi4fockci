# Psi4FockCI

<img src="https://travis-ci.com/shannonhouck/psi4fockci.svg?branch=master">

Code Author: Shannon E. Houck

Updated: Nov. 21, 2019

This code is a Python package that uses Psi4’s DETCI module to run Fock-space CI (RAS-nSF-IP/EA) calculations. 
The method handles spin and spatial degeneracies in molecular systems by solving the orbitals of a reference state 
that can be well-represented by a single determinant, and then using non-particle-conserving and non-spin-conserving 
operators to obtain the desired state.

<a href="https://shannonhouck.github.io/psi4fockci/build/index.html">Read the docs here!</a>
