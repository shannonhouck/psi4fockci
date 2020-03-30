This code is a Python package that uses Psi4's DETCI module to run 
Fock-space CI (RAS-nSF-IP/EA) calculations. The method handles spin and 
spatial degeneracies in molecular systems by solving the orbitals of a 
reference state that can be well-represented by a single determinant, and 
then using non-particle-conserving and non-spin-conserving operators to 
obtain the desired state. 
`A more detailed description of the method itself can be found here.
<The RAS-nSF-IP/EA Method>`_
`Further details about this method, including examples and analysis, 
can be found in this paper.
<https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b01268>`_

Installation
============

Clone the program from the GitHub repository::

    $ git clone https://github.com/shannonhouck/psi4fockci.git

Then navigate into the directory and use pip to install::

    $ cd psi4fockci
    $ pip install -e .

You can import this as a Python package and use it however you want! 
If you have pytest installed, you can use it to test your installation::

    $ cd psi4fockci
    $ pytest

Running CAS-nSF-IP/EA
=====================

The plugin can be run directly through Psi4's energy() call, as with 
any Psi4 plugin. The number of spin-flips and IP/EA to perform are 
determined automatically based on the given charge and multiplicity 
of the target state. In order to run a CAS-1SF-IP/STO-3G calculation, 
for example, one could set an input file up in the following way::

    molecule {
    0 7
    N 0.0 0.0 0.0
    N 0.0 0.0 1.3
    symmetry c1
    }

    set {
      basis cc-pVDZ
    }

    energy('psi4fockci', new_charge=1, new_multiplicity=1)

The input file can then be fun from the command line::

    $ psi4 example.dat

The program can also be run through the ``run_psi4fockci`` function call.
See the documentation of that function for information about the various 
options and keywords.

Passing Keywords to Psi4
========================

If running with Psi4, keywords for various modules can be set as normal 
in the input file::

    set detci {
      ci_maxiter 500
      num_roots 7
    }

Alternately, keywords can be passed to Psi4 using the ``add_opts`` keyword. 
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

    * ``""`` CAS-nSF-IP/EA (default, no additional excitations)
    * ``"h"`` RAS(h)-nSF-IP/EA (hole excitations)
    * ``"p"`` RAS(p)-nSF-IP/EA (particle excitations)
    * ``"S"`` RAS(S)-nSF-IP/EA (singles)
    * ``"SD"`` RAS(SD)-nSF-IP/EA (singles and doubles)
    * ``"SDT"`` RAS(SDT)-nSF-IP/EA (singles, doubles, and triples)


