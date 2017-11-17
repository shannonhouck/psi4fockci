# psi4_spinflip_wfn
Code Author: Shannon Houck

Updated: Nov 14, 2017

Currently, the molecule and alpha/beta specifications are set in the input file.

To run an SF-CAS/STO-3G calculation, for example, on the septet state of N2,
 with -2a/+2b, you would edit the input.dat file like so:

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
charge = 0
multiplicity = 3

# set up additional options
options = {"basis": "sto-3g"}

# run SF-CAS
e = sf_cas( charge, multiplicity, n2, conf_space="" , add_opts=options)
```

To run a SF-CAS(1x), you set the conf_space variable; so, one
would write `sf_cas( charge, multiplicity, n2, conf_space="1x" )` 
in place of the sf_cas function call above.

The output from the calculation can be found in output.dat.
