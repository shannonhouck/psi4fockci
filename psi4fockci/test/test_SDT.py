import psi4
import psi4fockci
from psi4fockci import sf_cas

# threshold for value equality
threshold = 1e-7
# setting up molecule
o2 = psi4.core.Molecule.from_string("""
0 3
O
O 1 1.2
symmetry c1
""")

# Test: RAS(SDT)-1SF/6-31G with O2 (0,3 to 0,1)
def test_1():
  options = {"basis": "6-31G", "num_roots": 4}
  e = sf_cas( 0, 1, o2, conf_space="SDT", add_opts=options, localize=True )
  expected = -149.730699517793994
  assert (e - expected) < threshold

