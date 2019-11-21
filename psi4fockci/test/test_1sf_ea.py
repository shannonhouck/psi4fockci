import psi4
import psi4fockci
from psi4fockci import sf_cas

# threshold for value equality
threshold = 1e-7
# setting up molecule
n2 = psi4.core.Molecule.from_string("""
0 7
N 0 0 0
N 0 0 2.5
symmetry c1
""")

# Test: CAS-1SF-EA/CC-PVDZ with N2
def test_1():
  options = {"basis": "cc-pvdz"}
  e = sf_cas( -1, 2, n2, add_opts=options, localize=True )
  expected = -108.600832070267
  assert (e - expected) < threshold

