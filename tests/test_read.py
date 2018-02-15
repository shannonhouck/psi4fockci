import sys, os
import psi4
sys.path.insert(1, '../')
import spinflip
from spinflip import sf_cas

# threshold for value equality
threshold = 1e-7
# setting up molecule
n2 = psi4.core.Molecule.create_molecule_from_string("""
0 7
N 0 0 0
N 0 0 2.5
symmetry c1
""")

# Test: SF-CAS/STO-3G with N2 (0,7 to -1,2)
def test_1():
  options = {"basis": "sto-3g"}
  e_1 = sf_cas( -1, 2, n2, add_opts=options, write_rohf_wfn="test_rohf.npz" )
  e_2 = sf_cas( -1, 2, n2, add_opts=options, read_rohf_wfn="test_rohf.npz" )
  assert (e_1 - e_2) < threshold


