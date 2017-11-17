import sys, os
sys.path.insert(1, '/usr/local/psi4/lib')
import psi4
from psi4 import *

# A method to run a spin-flip calculation
# Params:
#  new_charge - the target charge of the molecule
#  new_multiplicity - the target multiplicity of the molecule
#  mol - molecule to run the calculation on
#  conf_space - the configuration space to use ("" or "1x" currently supported)
#  opts - additional options to pass into Psi4 (optional)
def sf_cas( new_charge, new_multiplicity, ref_mol, conf_space="", add_opts={} ):
  # update options to include any additional opts from the user
  psi4.core.clean()
  opts = {'scf_type': 'pk', 'basis': 'cc-pvdz', 'reference': 'rohf', 'guess': 'sad', 'diis_start': 20, 'e_convergence': 1e-12, 'd_convergence':1e-12, 'mixed': False}
  opts.update(add_opts)
  if(conf_space == ""):
    print("SF-CAS CALCULATION")
  else:
    print("SF-CAS(%s) CALCULATION" % conf_space)
  # run rohf calculation on reference state
  print("RUNNING REFERENCE...\tCHARGE %i\tMULT %i" %(ref_mol.molecular_charge(), ref_mol.multiplicity()))
  psi4.set_options(opts)
  e_rohf, wfn_rohf = energy('scf', molecule=ref_mol, return_wfn=True, options=opts)
  print("SCF (%i %i): %6.12f" %(ref_mol.molecular_charge(), ref_mol.multiplicity(), e_rohf))
  #
  # change charge and multiplicity to new target values
  print("DOING SPIN-FLIP: CHARGE %i, MULTIPLICITY %i" % (new_charge, new_multiplicity))
  # copy molecule so original molecule passed in is unchanged
  mol = ref_mol.clone()
  mol.set_molecular_charge(new_charge)
  mol.set_multiplicity(new_multiplicity)
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
  # set active space and docc space based on configuration space input
  if(conf_space == ""):
    opts.update({'frozen_docc': [wfn_rohf.doccpi()[0]]})
    opts.update({'ras1': [0]})
    opts.update({'ras2': [wfn_rohf.soccpi()[0]]})
    opts.update({'ras3': [0]})
    opts.update({'ras4': [0]})
  elif(conf_space == "1x"):
    opts.update({'frozen_docc': [0]})
    opts.update({'ex_level': 1})
    opts.update({'ras1': [wfn_rohf.doccpi()[0]]})
    opts.update({'ras2': [wfn_rohf.soccpi()[0]]})
    opts.update({'ras3': [wfn_rohf.nmo() - wfn_rohf.soccpi()[0] - wfn_rohf.doccpi()[0]]})
    opts.update({'ras4': [0]})
  elif(conf_space == "2x"):
    opts.update({'frozen_docc': [0]})
    opts.update({'ex_level': 2})
    opts.update({'ex_allow': [0, 1]})
    opts.update({'ras1': [wfn_rohf.doccpi()[0]]})
    opts.update({'ras2': [wfn_rohf.soccpi()[0]]})
    opts.update({'ras3': [wfn_rohf.nmo() - wfn_rohf.soccpi()[0] - wfn_rohf.doccpi()[0]]})
    opts.update({'ras4': [0]})
  else:
    print("Configuration space $conf_space not supported. Exiting...")
    exit()
  #
  # run cas
  print("RUNNING CAS...\t\tCHARGE %i\tMULT %i" %(mol.molecular_charge(), mol.multiplicity()))
  psi4.set_options(opts)
  e_cas, wfn_cas = energy('detci', ref_wfn=wfn_rohf_new, return_wfn=True, molecule=mol)
  print("CAS (%i %i): %6.12f" %(mol.molecular_charge(), mol.multiplicity(), e_cas))
  return e_cas

