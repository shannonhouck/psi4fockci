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

# Test: SF-CAS/STO-3G with N2 (0,7 to 0,1)
def test_1():
  options = {"basis": "sto-3g"}
  e = sf_cas( 0, 1, n2, add_opts=options, localize=True )
  expected = -107.439176904454
  assert (e - expected) < threshold

# Test: SF-CAS/STO-3G with N2 (0,7 to 0,3)
def test_2():
  options = {"basis": "sto-3g"}
  e = sf_cas( 0, 3, n2, add_opts=options, localize=True )
  expected = -107.437970126831
  assert (e - expected) < threshold

# Test: SF-CAS/STO-3G with N2 (0,7 to -1,2)
def test_3():
  options = {"basis": "sto-3g"}
  e = sf_cas( -1, 2, n2, add_opts=options, localize=True )
  expected = -107.112752304674
  assert (e - expected) < threshold

# Test: SF-CAS/CC-PVDZ with N2 (0,7 to 0,1)
def test_4():
  options = {"basis": "cc-pvdz"}
  e = sf_cas( 0, 1, n2 )
  expected = -108.776024394853295
  assert (e - expected) < threshold

