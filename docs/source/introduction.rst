This code is a Python package that uses Psi4's DETCI module to run 
Fock-space CI (RAS-nSF-IP/EA) calculations. Details about this method can 
be found here. [TODO: LINK TO PAPER]

Installation
============

Clone the program from the GitHub repository::

    $ git clone https://github.com/shannonhouck/psi4_spinflip_wfn.git

Then navigate into the directory and use pip to install::

    $ cd psi4_spinflip_wfn
    $ pip install -e .

You can import this as a Python package and use it however you want! 
If you have pytest installed, you can use it to test your installation::

    $ cd psi4fockci
    $ pytest

Running CAS-nSF-IP/EA
=====================

The program is primarily run through the ``sf_cas`` function call. The 
function should be run as follows:

    1. Initialize a Psi4.core.Molecule object with the correct 
       reference charge/multiplicity.
    2. Determine the target charge/multiplicity.
    3. Run ``sf_cas``.

In order to run a 1SF-CAS/STO-3G 2-SF 1-IP calculation on N2+, one 
could set an input file up as follows::

    import psi4
    import spinflip
    from spinflip import sf_cas

    # setting up molecule
    n2 = psi4.core.Molecule.create_molecule_from_string("""
    0 7
    N 0 0 0
    N 0 0 2.5
    symmetry c1
    """)

    # set target charge and multiplicity
    charge = 1
    multiplicity = 2

    # set up additional options
    options = {"basis": "sto-3g"}

    # run SF-CAS
    e = sf_cas( charge, multiplicity, n2, conf_space="" , add_opts=options)

This file can subsequently be run from the command line::

    $ python example.dat

All relevant information is written to standard output.

Passing Keywords to Psi4
========================

Additional keywords can be passed to Psi4 using the ``add_opts`` keyword. 
These options should be put in the dictionary form usually taken by Psi4. 
For example, if I wanted to change the number of CI roots, I could specify 
it as follows::

    options = {"basis": "sto-3g", "num_roots": 10}
    e = sf_cas( charge, multiplicity, n2, conf_space="" , add_opts=options)

Adding Excitations
==================

Excitations are important, particularly for the nSF-IP/EA schemes. 
(Hole excitations are recommended for IP-type and particle excitations 
are recommended for EA-type; see the paper for details.) 
Excitations outside of the CAS space can be requested by setting the 
``conf_space`` keyword appropriately. The following keywords are valid:

    * ``""`` CAS-nSF-IP/EA (no additional excitations)
    * ``"h"`` RAS(h)-nSF-IP/EA (hole excitations)
    * ``"p"`` RAS(p)-nSF-IP/EA (particle excitations)
    * ``"S"`` RAS(S)-nSF-IP/EA (singles)
    * ``"SD"`` RAS(SD)-nSF-IP/EA (singles and doubles)
    * ``"SDT"`` RAS(SDT)-nSF-IP/EA (singles, doubles, and triples)


