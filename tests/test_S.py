import sys, os
import psi4
sys.path.insert(1, '../')
import spinflip
from spinflip import sf_cas

# threshold for value equality
threshold = 1e-7
# setting up molecule
o2 = psi4.core.Molecule.create_molecule_from_string("""
0 3
O
O 1 1.2
symmetry c1
""")

# Test: SF-xCIS/CC-PVDZ with O2 (0,3 to 0,1)
def test_1():
  options = {"basis": "cc-pvdz"}
  e = sf_cas( 0, 1, o2, conf_space="xcis", add_opts=options )
  expected = -149.604321051649
  assert (e - expected) < threshold

