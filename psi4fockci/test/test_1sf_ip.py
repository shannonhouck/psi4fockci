import psi4
import psi4fockci
from psi4fockci import run_psi4fockci

# threshold for value equality
threshold = 1e-7
# setting up molecule
n2 = psi4.core.Molecule.from_string("""
0 7
N 0 0 0
N 0 0 2.5
symmetry c1
""")

# Test: CAS-2SF-IP/CC-PVDZ with N2+
def test_1():
    psi4.core.clean()
    psi4.core.clean_options() 
    psi4.core.clean_variables() 
    options = {"basis": "cc-pvdz"}
    wfn = run_psi4fockci('psi4fockci', n2, new_charge=1, new_multiplicity=2, 
        add_opts=options)
    e = psi4.core.get_variable("CI ROOT 0 TOTAL ENERGY")
    expected = -108.250639579451
    assert (e - expected) < threshold


