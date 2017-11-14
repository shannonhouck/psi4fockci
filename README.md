# psi4_spinflip_wfn
Code Author: Shannon Houck

Updated: Nov 14, 2017

Currently, the molecule and alpha/beta specifications are set in the input file.

To run a calculation, for example, on the septet state of N2, with -2a/+2b, 
you would edit the input.dat file like so:

```
import spinflip
from spinflip import run_sf

# setting up molecule
n2 = psi4.core.Molecule.create_molecule_from_string("""
0 7
N 0 0 0
N 0 0 2.5
symmetry c1
""")
da = -2 # change in alpha electron count
db = 2 # change in beta electron count
# running the spin-flip calculation
run_sf( da, db, n2 )
```

The output from the calculation can be found in output.dat.
