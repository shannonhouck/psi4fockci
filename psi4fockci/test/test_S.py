import psi4
import psi4fockci
from psi4fockci import run_psi4fockci

# threshold for value equality
threshold = 1e-7
# setting up molecule
o2 = psi4.core.Molecule.from_string("""
0 3
O
O 1 1.2
symmetry c1
""")

# Test: RAS(S)-1SF/CC-PVDZ with O3
def test_1():
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.core.clean_variables()
    options = {"basis": "cc-pvdz"}
    wfn = run_psi4fockci('psi4fockci', o2, new_charge=0, new_multiplicity=1, 
        conf_space="S", add_opts=options)
    e = psi4.core.get_variable("CI ROOT 0 TOTAL ENERGY")
    expected = -149.506943097607547
    assert (e - expected) < threshold

