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

# Test: SF-CAS/STO-3G with N2 (0,7 to 0,1)
def test_1():
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.core.clean_variables()
    options = {"basis": "sto-3g"}
    wfn = run_psi4fockci('psi4fockci', n2, new_charge=0, new_multiplicity=1, 
        add_opts=options )
    e = psi4.core.variable("CI ROOT 0 TOTAL ENERGY")
    expected = -107.439176904454
    assert (e - expected) < threshold

# Test: SF-CAS/STO-3G with N2 (0,7 to 0,3)
def test_2():
    psi4.core.clean()
    psi4.core.clean_options() 
    psi4.core.clean_variables() 
    options = {"basis": "sto-3g"}
    wfn = run_psi4fockci('psi4fockci', n2, new_charge=0, new_multiplicity=3, 
        add_opts=options )
    e = psi4.core.variable("CI ROOT 0 TOTAL ENERGY")
    expected = -107.437970126831
    assert (e - expected) < threshold

# Test: SF-CAS/CC-PVDZ with N2 (0,7 to 0,1)
def test_3():
    psi4.core.clean()
    psi4.core.clean_options() 
    psi4.core.clean_variables() 
    options = {"basis": "cc-pvdz"}
    wfn = run_psi4fockci('psi4fockci', n2, new_charge=0, new_multiplicity=1,
        add_opts=options )
    e = psi4.core.variable("CI ROOT 0 TOTAL ENERGY")
    expected = -108.776024394853295
    assert (e - expected) < threshold

