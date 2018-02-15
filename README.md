# psi4_spinflip_wfn
Code Author: Shannon Houck

Updated: Jan 26, 2018

Currently, the molecule and alpha/beta specifications are set in the input file.

To run an SF-CAS/STO-3G 2-SF 1-IP calculation on N2+, for example,
 you would edit the example.dat file like so:

```
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
```

This file can subsequently be run from the command line:

```
$ psi4 example.dat
```

To run a SF-CAS(1x) or SF-CAS(S), set the conf_space variable; so, one
would write `sf_cas( charge, multiplicity, n2, conf_space="1x" )` 
in place of the sf_cas function call above. Additional things 
(ex. the CI wavefunction) can be returned using various keywords 
as outlined in the spinflip.py file comments.

If you have pytest installed, the tests in the tests/ directory can be run with:

```
$ cd tests/
$ pytest
```

