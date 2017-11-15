import sys, os
sys.path.insert(1, '/usr/local/psi4/lib')
import psi4
from psi4 import *

# A method to run a spin-flip calculation
# Params:
#  del_alpha - change in number of alpha electrons
#  del_beta - change in number of beta electrons
#  mol - molecule to run the calculation on
#  opts - the options to pass into Psi4 (optional)
def sf_cas( del_alpha, del_beta, mol, opts={'scf_type': 'pk', 'basis': 'cc-pvdz', 'reference': 'rohf', 'guess': 'sad', 'diis_start': 20, 'e_convergence': 1e-12, 'd_convergence':1e-12} ):
  print("SF-CAS CALCULATION")
  # run rohf calculation on reference state
  print("RUNNING REFERENCE...\tCHARGE %i\tMULT %i" %(mol.molecular_charge(), mol.multiplicity()))
  psi4.set_options(opts)
  e_rohf, wfn_rohf = energy('scf', molecule=mol, return_wfn=True, options=opts)
  print("SCF (%i %i): %6.12f" %(mol.molecular_charge(), mol.multiplicity(), e_rohf))
  #
  # change number of alpha/beta electrons
  print("DOING SPIN-FLIP: %+i ALPHA, %+i BETA" % (del_alpha, del_beta))
  # sanity check!!
  if((del_alpha > wfn_rohf.soccpi()[0]) or (del_beta > wfn_rohf.soccpi()[0])):
    print("The number of alphas flipped (%i) cannot be greater than the number of singly occupied orbitals (%i)!" % (del_alpha, wfn_rohf.soccpi()[0]) )
    print("Exiting...")
    exit()
  # del_charge is the difference between total electrons
  del_charge = -(del_alpha + del_beta)
  mol.set_molecular_charge(mol.molecular_charge() + del_charge)
  # multiplicity is the difference between total alpha and beta
  na = wfn_rohf.nalpha() + del_alpha
  nb = wfn_rohf.nbeta() + del_beta
  mol.set_multiplicity(max(na, nb) - min(na, nb) + 1)
  #
  # set up reference wfn to pass into detci
  # run RHF calculation to initialize soccpi, doccpi, nalpha, nbeta, etc.
  psi4.set_options(opts)
  e_rohf_new, wfn_rohf_new = energy('scf', molecule=mol, return_wfn=True, options=opts)
  # fill wfn with values from reference calculation
  wfn_rohf_new.Ca().copy(wfn_rohf.Ca())
  wfn_rohf_new.Cb().copy(wfn_rohf.Cb())
  wfn_rohf_new.H().copy(wfn_rohf.H())
  #
  # set active space and docc space
  opts.update({'frozen_docc': [wfn_rohf.doccpi()[0]]})
  opts.update({'ras1': [0]})
  opts.update({'ras2': [wfn_rohf.soccpi()[0]]})
  opts.update({'ras3': [0]})
  opts.update({'ras4': [0]})
  #
  # run cas
  print("RUNNING CAS...\t\tCHARGE %i\tMULT %i" %(mol.molecular_charge(), mol.multiplicity()))
  psi4.set_options(opts)
  e_cas, wfn_cas = energy('detci', ref_wfn=wfn_rohf_new, return_wfn=True, molecule=mol)
  print("CAS (%i %i): %6.12f" %(mol.molecular_charge(), mol.multiplicity(), e_cas))

# Run a Spin-Flip RAS calculation with the correct excitation level.
def sf_ras( del_alpha, del_beta, mol, ex_level, opts={'scf_type': 'pk', 'basis': 'cc-pvdz', 'reference': 'rohf', 'guess': 'sad', 'diis_start': 20, 'e_convergence': 1e-12, 'd_convergence':1e-12} ):
  print("SF-RAS(%ix) CALCULATION" % ex_level)
  # run rohf calculation on reference state
  print("RUNNING REFERENCE...\tCHARGE %i\tMULT %i" %(mol.molecular_charge(), mol.multiplicity()))
  psi4.set_options(opts)
  e_rohf, wfn_rohf = energy('scf', molecule=mol, return_wfn=True, options=opts)
  print("SCF (%i %i): %6.12f" %(mol.molecular_charge(), mol.multiplicity(), e_rohf))
  #
  # change number of alpha/beta electrons
  print("DOING SPIN-FLIP: %+i ALPHA, %+i BETA" % (del_alpha, del_beta))
  # sanity check!!
  if((del_alpha > wfn_rohf.soccpi()[0]) or (del_beta > wfn_rohf.soccpi()[0])):
    print("The number of alphas flipped (%i) cannot be greater than the number of singly occupied orbitals (%i)!" % (del_alpha, wfn_rohf.soccpi()[0]) )
    print("Exiting...")
    exit()
  # del_charge is the difference between total electrons
  del_charge = -(del_alpha + del_beta)
  mol.set_molecular_charge(mol.molecular_charge() + del_charge)
  # multiplicity is the difference between total alpha and beta
  na = wfn_rohf.nalpha() + del_alpha
  nb = wfn_rohf.nbeta() + del_beta
  mol.set_multiplicity(max(na, nb) - min(na, nb) + 1)
  #
  # set up reference wfn to pass into detci
  # run RHF calculation to initialize soccpi, doccpi, nalpha, nbeta, etc.
  psi4.set_options(opts)
  e_rohf_new, wfn_rohf_new = energy('scf', molecule=mol, return_wfn=True, options=opts)
  # fill wfn with values from reference calculation
  wfn_rohf_new.Ca().copy(wfn_rohf.Ca())
  wfn_rohf_new.Cb().copy(wfn_rohf.Cb())
  wfn_rohf_new.H().copy(wfn_rohf.H())
  #
  # set active space and docc space
  opts.update({'frozen_docc': [0]})
  opts.update({'ex_level': ex_level})
  opts.update({'ras1': [wfn_rohf.doccpi()[0]]})
  opts.update({'ras2': [wfn_rohf.soccpi()[0]]})
  opts.update({'ras3': [wfn_rohf.frzvpi()[0]]})
  opts.update({'ras4': [0]})
  #
  # run cas
  print("RUNNING RAS(%i)...\t\tCHARGE %i\tMULT %i" %(ex_level, mol.molecular_charge(), mol.multiplicity()))
  psi4.set_options(opts)
  e_cas, wfn_cas = energy('detci', ref_wfn=wfn_rohf_new, return_wfn=True, molecule=mol)
  print("RAS (%i %i): %6.12f" %(mol.molecular_charge(), mol.multiplicity(), e_cas))

