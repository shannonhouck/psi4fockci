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

# Test: RAS(p)-3SF/cc-pvdz with N2 (0,7 to 0,1)
def test_1():
  options = {"basis": "cc-pvdz", 'num_roots': 2, 'diis_start': 20}
  e = sf_cas( 0, 1, n2, conf_space="p", add_opts=options, localize=True )
  expected = -108.773240257969
  assert (e - expected) < threshold

