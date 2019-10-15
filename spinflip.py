import sys, os
import glob
import shutil
sys.path.insert(1, '/usr/local/psi4/lib')
import psi4
from psi4 import *
import numpy as np
import numpy.linalg as LIN

# Used to compute various matrices...
# Includes S^(1/2) and S^(-1/2).
import numpy as np
from numpy import linalg as LIN

def sf_cas( new_charge, new_multiplicity, ref_mol, conf_space="", add_opts={}, return_ci_wfn=False, 
            return_rohf_wfn=False, return_rohf_e=False, read_rohf_wfn="", write_rohf_wfn="", 
            write_ci_vects=True, localize=False, frozen_docc=0, frozen_uocc=0):
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
    return_ci_wfn : bool (False)
        Whether to return the CI wavefunction object.
    return_rohf_wfn : bool (False)
        Whether to return the ROHF wavefunction object.
    return_rohf_e : bool (False)
        Whether to return the ROHF energy.
    read_rohf_wfn : str ("")
        Name of file (.npz) to read ROHF wavefunction info from.
        By default, no wavefunction is read in.
    write_rohf_wfn : str ("")
        Name of file (.npz) to read ROHF wavefunction info from.
        By default, no wavefunction is written.
        Note that you MUST run Psi4 with the -m (messy) flag for this to work!
    localize : bool (False)
        Whether to perform BOYS localization on the RAS 2 space before computing.
        Can help with visualization and analysis of orbitals.

    Returns
    -------
    e_ci : double
        The SF-CAS([conf_space]) energy.
    return_ci_wfn : psi4.core.Wavefunction
        (optional) The SF-CAS([conf_space]) wavefunction.
    return_rohf_e : double
        (optional) The ROHF energy.
    return_rohf_wfn : psi4.core.Wavefunction
        (optional) The ROHF wavefunction.

    """

    # printing initial information about the calculation
    print("SF-CAS(%s) CALCULATION" % conf_space)
    # default options for Psi4
    opts = {'scf_type': 'pk',
            'basis': 'cc-pvdz',
            'reference': 'rohf',
            'guess': 'sad',
            'maxiter': 1000,
            'ci_maxiter': 50,
            'mixed': False}
    opts.update(add_opts) # add additional options from user

    # run ROHF calculation on reference state or read it in
    psi4.core.clean() # cleanup in case Psi4 has been run before
    psi4.set_options(opts)
    mol = ref_mol.clone() # clone molecule so original isn't modified
    # read in ROHF guess wavefunction if provided
    if(read_rohf_wfn != ""):
        # set up options and run
        opts.update({'guess': 'read'})
        psi4.set_options(opts)
        print("RUNNING ROHF FROM REFERENCE...\tCHARGE %i\tMULT %i" %(ref_mol.molecular_charge(), ref_mol.multiplicity()))
        e_rohf, wfn_rohf = energy('scf', molecule=mol, return_wfn=True, options=opts, restart_file=read_rohf_wfn)
    # else, run ROHF on reference state
    else:
        print("RUNNING REFERENCE...\tCHARGE %i\tMULT %i" %(ref_mol.molecular_charge(), ref_mol.multiplicity()))
        e_rohf, wfn_rohf = energy('scf', molecule=mol, return_wfn=True, options=opts)
        print("SCF (%i %i): %6.12f" %(mol.molecular_charge(), mol.multiplicity(), e_rohf))

    # change charge and multiplicity to new target values
    print("DOING SPIN-FLIP: CHARGE %i, MULTIPLICITY %i" % (new_charge, new_multiplicity))

    # saving npz file of wavefunction
    if(write_rohf_wfn != ""):
        shutil.copy(glob.glob('./*.180.npz')[0], write_rohf_wfn)

    # update molecular charge and multiplicity
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

    # if we need to localize...
    if(localize):
        C = psi4.core.Matrix.to_array(wfn_rohf.Ca(), copy=True)
        ras1_C = C[:, :doccpi]
        ras2_C = C[:, doccpi:doccpi+soccpi]
        ras3_C = C[:, doccpi+soccpi:]
        loc = psi4.core.Localizer.build('BOYS', wfn_rohf.basisset(), psi4.core.Matrix.from_array(ras2_C))
        loc.localize()
        ras2_localized = psi4.core.Matrix.to_array(loc.L, copy=True)
        localized_orbs = np.column_stack((ras1_C, ras2_localized, ras3_C))
        new_Ca = psi4.core.Matrix.from_array(localized_orbs, name="Ca")
        new_Cb = psi4.core.Matrix.from_array(localized_orbs, name="Cb")
        wfn_rohf.Ca().copy(new_Ca)
        wfn_rohf.Cb().copy(new_Cb)

    # set active space and docc space based on configuration space input
    # Regular CAS configuration space
    # includes only active space configurations
    if(conf_space == ""):
      opts.update({'frozen_docc': [doccpi]})
      opts.update({'ras1': [0]})
      opts.update({'ras2': [soccpi]})
      opts.update({'ras3': [0]})
      opts.update({'ras4': [0]})
    # just (h) excitations
    elif(conf_space == "h"):
      opts.update({'ex_level': 0})
      opts.update({'val_ex_level': 1})
      opts.update({'ras3_max': 0})
      opts.update({'frozen_docc': [frozen_docc]})
      opts.update({'ras1': [doccpi - frozen_docc]})
      opts.update({'ras2': [soccpi]})
      opts.update({'ras3': [0]})
      opts.update({'ras4': [0]})
    # just (p) excitations
    elif(conf_space == "p"):
      opts.update({'ex_level': 0})
      opts.update({'val_ex_level': 0})
      opts.update({'ras3_max': 1})
      opts.update({'frozen_docc': [doccpi]})
      opts.update({'ras1': [0]})
      opts.update({'ras2': [soccpi]})
      opts.update({'ras3': [nmo - soccpi - doccpi - frozen_uocc]})
      opts.update({'frozen_uocc': [frozen_uocc]})
      opts.update({'ras4': [0]})
    # 1x configuration space
    # includes (h, p) excitations
    elif(conf_space == "1x"):
      opts.update({'frozen_docc': [0]})
      opts.update({'ex_level': 0})
      opts.update({'val_ex_level': 1})
      opts.update({'ras3_max': 1})
      opts.update({'frozen_docc': [frozen_docc]})
      opts.update({'ras1': [doccpi - frozen_docc]})
      opts.update({'ras2': [soccpi]})
      opts.update({'ras3': [nmo - soccpi - doccpi - frozen_uocc]})
      opts.update({'frozen_uocc': [frozen_uocc]})
      opts.update({'ras4': [0]})
    # S configuration space
    # includes (h, p, hp) excitations
    elif(conf_space == "S" or conf_space == "xcis"):
      opts.update({'frozen_docc': [0]})
      opts.update({'ex_level': 1})
      opts.update({'frozen_docc': [frozen_docc]})
      opts.update({'ras1': [doccpi - frozen_docc]})
      opts.update({'ras2': [soccpi]})
      opts.update({'ras3': [nmo - soccpi - doccpi - frozen_uocc]})
      opts.update({'frozen_uocc': [frozen_uocc]})
      opts.update({'ras4': [0]})
    elif(conf_space == "SD"):
      opts.update({'frozen_docc': [0]})
      opts.update({'ex_level': 2})
      opts.update({'frozen_docc': [frozen_docc]})
      opts.update({'ras1': [doccpi - frozen_docc]})
      opts.update({'ras2': [soccpi]})
      opts.update({'ras3': [nmo - soccpi - doccpi - frozen_uocc]})
      opts.update({'frozen_uocc': [frozen_uocc]})
      opts.update({'ras4': [0]})
    elif(conf_space == "SDT"):
      opts.update({'frozen_docc': [0]})
      opts.update({'ex_level': 3})
      opts.update({'frozen_docc': [frozen_docc]})
      opts.update({'ras1': [doccpi - frozen_docc]})
      opts.update({'ras2': [soccpi]})
      opts.update({'ras3': [nmo - soccpi - doccpi - frozen_uocc]})
      opts.update({'frozen_uocc': [frozen_uocc]})
      opts.update({'ras4': [0]})
    # Other configuration spaces aren't supported yet
    else:
      print("Configuration space %s not supported. Exiting..." % conf_space)
      exit()

    # run cas
    print("RUNNING CAS...\t\tCHARGE %i\tMULT %i" %(mol.molecular_charge(), mol.multiplicity()))
    psi4.set_options(opts)
    e_cas, wfn_cas = energy('detci', ref_wfn=wfn_rohf, return_wfn=True, molecule=mol)
    print("CAS (%i %i): %6.12f" %(mol.molecular_charge(), mol.multiplicity(), e_cas))

    # obtain eigenvectors if needed
    # partly based on Daniel Smith's answer on Psi4 forums
    if(write_ci_vects):
        wfn_cas_2 = psi4.core.CIWavefunction(wfn_rohf)
        n_roots = add_opts['NUM_ROOTS']
        C = np.zeros((wfn_cas_2.ndet(), n_roots))
        print(C.shape)
        for i in range(n_roots):
            dvec = wfn_cas_2.new_civector(i+1, 53, True, True)
            dvec.set_nvec(i+1)
            dvec.init_io_files(True)
            dvec.read(i,0)
            C[:, i] = np.array(dvec)
        np.savetxt('ci_vect.txt', C)

    psi4.core.print_variables() # printing Psi4 variables
    psi4.core.clean_options() # more cleanup

    # return output specified by the user
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

