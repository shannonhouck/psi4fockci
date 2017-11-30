import sys, os
sys.path.insert(1, '/usr/local/psi4/lib')
import psi4
from psi4 import *

def sf_cas( new_charge, new_multiplicity, ref_mol, conf_space="", add_opts={} ):
    """
    A method to run a spin-flip electron addition calculation.

    Parameters
    ----------
    new_charge : int
        Target charge of the molecule.
    new_multiplicity : int
        Target multiplicity of the molecule.
    ref_mol : psi4.core.Molecule
        Molecule to run the calculation on.
    conf_space : str ("")
        Configuration space to use for CAS after spin-flip. 
        Options "", "S", and "1x" are currently supported.
    add_opts : dict ({})
        Additional options to pass into Psi4.

    Returns
    -------
    e_ci : double
        The SF-CAS([conf_space]) energy.

    """

    # printing initial information about the calculation
    print("SF-CAS(%s) CALCULATION" % conf_space)
    # update options to include any additional opts from the user
    opts = {'scf_type': 'pk',
            'basis': 'cc-pvdz',
            'reference': 'rohf',
            'guess': 'sad',
            'diis_start': 20,
            'e_convergence': 1e-12,
            'd_convergence': 1e-12,
            'mixed': False}
    opts.update(add_opts) # add additional options from user

    # run rohf calculation on reference state
    print("RUNNING REFERENCE...\tCHARGE %i\tMULT %i" %(ref_mol.molecular_charge(), ref_mol.multiplicity()))
    psi4.core.clean() # cleanup (in case Psi4 has been run before)
    psi4.set_options(opts)
    e_rohf, wfn_rohf = energy('scf', molecule=ref_mol, return_wfn=True, options=opts)
    psi4.core.print_variables()
    print("SCF (%i %i): %6.12f" %(ref_mol.molecular_charge(), ref_mol.multiplicity(), e_rohf))

    # change charge and multiplicity to new target values
    print("DOING SPIN-FLIP: CHARGE %i, MULTIPLICITY %i" % (new_charge, new_multiplicity))

    # copy molecule so original molecule passed in is unchanged
    mol = ref_mol.clone()
    mol.set_molecular_charge(new_charge)
    mol.set_multiplicity(new_multiplicity)

    # set up reference wfn to pass into detci
    #wfn_rohf_new = psi4.core.Wavefunction.build(mol, psi4.core.BasisSet.build(mol))
    wfn_rohf_new = psi4.core.ROHF.build(mol, psi4.core.BasisSet.build(mol))
    #psi4.core.ROHF.initialize(psi4.core.ROHF.build(mol, psi4.core.BasisSet.build(mol)))
    wfn_rohf_new.set_array('Ca', wfn_rohf.Ca())
    wfn_rohf_new.set_array('Cb', wfn_rohf.Cb())
    wfn_rohf_new.set_array('H', wfn_rohf.H())
    # calculate soccpi, doccpi
    n_soccpi = mol.multiplicity()-1
    wfn_rohf_new.soccpi()[0] = n_soccpi
    print(wfn_rohf_new.soccpi()[0])
    
    del_electrons = ref_mol.molecular_charge() - mol.molecular_charge()
    n_total = wfn_rohf.nalpha() + wfn_rohf.nbeta() + del_electrons
    wfn_rohf_new.doccpi()[0] = n_total - n_soccpi
    print(wfn_rohf_new.doccpi()[0])

    # set active space and docc space based on configuration space input
    if(conf_space == ""):
      opts.update({'frozen_docc': [wfn_rohf.doccpi()[0]]})
      opts.update({'ras1': [0]})
      opts.update({'ras2': [wfn_rohf.soccpi()[0]]})
      opts.update({'ras3': [0]})
      opts.update({'ras4': [0]})
    elif(conf_space == "S" or conf_space == "xcis"):
      opts.update({'frozen_docc': [0]})
      opts.update({'ex_level': 1})
      opts.update({'ras1': [wfn_rohf.doccpi()[0]]})
      opts.update({'ras2': [wfn_rohf.soccpi()[0]]})
      opts.update({'ras3': [wfn_rohf.nmo() - wfn_rohf.soccpi()[0] - wfn_rohf.doccpi()[0]]})
      opts.update({'ras4': [0]})
    elif(conf_space == "1x"):
      opts.update({'frozen_docc': [0]})
      opts.update({'ex_level': 0})
      opts.update({'val_ex_level': 1})
      opts.update({'ras3_max': 1})
      opts.update({'ras1': [wfn_rohf.doccpi()[0]]})
      opts.update({'ras2': [wfn_rohf.soccpi()[0]]})
      opts.update({'ras3': [wfn_rohf.nmo() - wfn_rohf.soccpi()[0] - wfn_rohf.doccpi()[0]]})
      opts.update({'ras4': [0]})
    #elif(conf_space == "2x"):
    #  opts.update({'frozen_docc': [0]})
    #  opts.update({'ex_level': 2})
    #  opts.update({'ex_allow': [0, 1]})
    #  opts.update({'val_ex_level': 1})
    #  opts.update({'ras1': [wfn_rohf.doccpi()[0]]})
    #  opts.update({'ras2': [wfn_rohf.soccpi()[0]]})
    #  opts.update({'ras3': [wfn_rohf.nmo() - wfn_rohf.soccpi()[0] - wfn_rohf.doccpi()[0]]})
    #  opts.update({'ras4': [0]})
    else:
      print("Configuration space %s not supported. Exiting..." % conf_space)
      exit()
    print(opts)

    # run cas
    print("RUNNING CAS...\t\tCHARGE %i\tMULT %i" %(mol.molecular_charge(), mol.multiplicity()))
    psi4.set_options(opts)
    e_cas, wfn_cas = energy('detci', ref_wfn=wfn_rohf_new, return_wfn=True, molecule=mol)
    print("CAS (%i %i): %6.12f" %(mol.molecular_charge(), mol.multiplicity(), e_cas))
    psi4.core.print_variables()
    psi4.core.clean_options() # more cleanup
    return e_cas

