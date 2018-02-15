import sys, os
sys.path.insert(1, '/usr/local/psi4/lib')
import psi4
from psi4 import *

def sf_cas( new_charge, new_multiplicity, ref_mol, conf_space="", add_opts={}, ref_rohf_wfn="NONE", return_ci_wfn=False, return_rohf_wfn=False, return_rohf_e=False ):
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
            'maxiter': 1000,
            'ci_maxiter': 50,
            'mixed': False}
    opts.update(add_opts) # add additional options from user

    # run rohf calculation on reference state
    print("RUNNING REFERENCE...\tCHARGE %i\tMULT %i" %(ref_mol.molecular_charge(), ref_mol.multiplicity()))
    psi4.core.clean() # cleanup (in case Psi4 has been run before)
    # storing MOs (in case we need to restart the job)
    psi4.set_options(opts)
    mol = ref_mol.clone() # clone molecule so original isn't modified
    if(ref_rohf_wfn == "NONE"):
        e_rohf, wfn_rohf = energy('scf', molecule=mol, return_wfn=True, options=opts)
    else:
        e_rohf, wfn_rohf = energy('scf', molecule=mol, ref_wfn=ref_rohf_wfn, return_wfn=True, options=opts)
    print("SCF (%i %i): %6.12f" %(mol.molecular_charge(), mol.multiplicity(), e_rohf))

    # change charge and multiplicity to new target values
    print("DOING SPIN-FLIP: CHARGE %i, MULTIPLICITY %i" % (new_charge, new_multiplicity))

    # copy molecule so original molecule passed in is unchanged
    mol.set_molecular_charge(new_charge)
    mol.set_multiplicity(new_multiplicity)

    # set up reference wfn to pass into detci
    # save orbital values from reference calculation
    doccpi = wfn_rohf.doccpi()[0]
    soccpi = wfn_rohf.soccpi()[0]
    nmo = wfn_rohf.nmo()

    # calculate soccpi and doccpi
    new_soccpi = mol.multiplicity() - 1
    del_electrons = ref_mol.molecular_charge() - mol.molecular_charge()
    n_total = wfn_rohf.nalpha() + wfn_rohf.nbeta() + del_electrons
    # set orbital occupations
    wfn_rohf.force_soccpi(psi4.core.Dimension([new_soccpi]))
    wfn_rohf.force_doccpi(psi4.core.Dimension([(int)((n_total - new_soccpi)/2)]))

    # set active space and docc space based on configuration space input
    if(conf_space == ""):
      opts.update({'frozen_docc': [doccpi]})
      opts.update({'ras1': [0]})
      opts.update({'ras2': [soccpi]})
      opts.update({'ras3': [0]})
      opts.update({'ras4': [0]})
    elif(conf_space == "S" or conf_space == "xcis"):
      opts.update({'frozen_docc': [0]})
      opts.update({'ex_level': 1})
      opts.update({'ras1': [doccpi]})
      opts.update({'ras2': [soccpi]})
      opts.update({'ras3': [nmo - soccpi - doccpi]})
      opts.update({'ras4': [0]})
    elif(conf_space == "1x"):
      opts.update({'frozen_docc': [0]})
      opts.update({'ex_level': 0})
      opts.update({'val_ex_level': 1})
      opts.update({'ras3_max': 1})
      opts.update({'ras1': [doccpi]})
      opts.update({'ras2': [soccpi]})
      opts.update({'ras3': [nmo - soccpi - doccpi]})
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

    # run cas
    print("RUNNING CAS...\t\tCHARGE %i\tMULT %i" %(mol.molecular_charge(), mol.multiplicity()))
    psi4.set_options(opts)
    e_cas, wfn_cas = energy('detci', ref_wfn=wfn_rohf, return_wfn=True, molecule=mol)
    print("CAS (%i %i): %6.12f" %(mol.molecular_charge(), mol.multiplicity(), e_cas))
    psi4.core.print_variables()
    psi4.core.clean_options() # more cleanup
    # return required output
    if((not return_ci_wfn) and (not return_rohf_wfn) and (not return_rohf_e)):
        return e_cas
    else:
        out = (e_cas,)
        if(return_ci_wfn):
            out = out + (wfn_cas,)
        if(return_rohf_wfn):
            out = out + (wfn_rohf,)
        if(return_rohf_e):
            out = out + (e_rohf,)
        return out

