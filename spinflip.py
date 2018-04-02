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

# Computes the orthogonalized square root of a given matrix.
# Params:
#  numpy A - Matrix to be square-rooted
#  int threshold_S - Threshold for S values
# Returns:
#  numpy matrix A^(1/2)
def matrix_sqrt( A, threshold_S ):
  # do an SVD of matrix A
  A_vals, A_vects = LIN.eig(A)
  # set up output matrix and populate with A_s values
  A_sqrt = np.zeros((len(A_vals), len(A_vals)))
  for i in range(0, len(A_vals)):
    if(A_vals[i] > threshold_S):
      A_sqrt[i,i] = np.sqrt(A_vals[i])
  # form the final output matrix and return
  return np.dot(A_vects, np.dot(A_sqrt, LIN.inv(A_vects)))

# Computes the orthogonalized square root inverse of a given matrix.
# Params:
#  numpy A - Matrix to be square-inverted
#  int threshold_S - Threshold for S values
# Returns:
#  numpy matrix A^(-1/2)
def matrix_inv_sqrt( A, threshold_S ):
  # do an SVD of matrix A
  A_vals, A_vects = LIN.eig(A)
  # set up output matrix and populate with A_s values
  A_sqrt = np.zeros((len(A_vals), len(A_vals)))
  for i in range(0, len(A_vals)):
    if(A_vals[i] > threshold_S):
      A_sqrt[i,i] = 1/np.sqrt(A_vals[i])
  # form the final output matrix and return
  return np.dot(A_vects, np.dot(A_sqrt, LIN.inv(A_vects)))


def sf_cas( new_charge, new_multiplicity, ref_mol, conf_space="", add_opts={}, return_ci_wfn=False, return_rohf_wfn=False, return_rohf_e=False, read_rohf_wfn="", write_rohf_wfn="", rotate_orbitals=False):
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

    # set active space and docc space based on configuration space input
    # Regular CAS configuration space
    # includes only active space configurations
    if(conf_space == ""):
      opts.update({'frozen_docc': [doccpi]})
      opts.update({'ras1': [0]})
      opts.update({'ras2': [soccpi]})
      opts.update({'ras3': [0]})
      opts.update({'ras4': [0]})
    # 1x configuration space
    # includes (h, p) excitations
    elif(conf_space == "1x"):
      opts.update({'frozen_docc': [0]})
      opts.update({'ex_level': 0})
      opts.update({'val_ex_level': 1})
      opts.update({'ras3_max': 1})
      opts.update({'ras1': [doccpi]})
      opts.update({'ras2': [soccpi]})
      opts.update({'ras3': [nmo - soccpi - doccpi]})
      opts.update({'ras4': [0]})
    # S configuration space
    # includes (h, p, hp) excitations
    elif(conf_space == "S" or conf_space == "xcis"):
      opts.update({'frozen_docc': [0]})
      opts.update({'ex_level': 1})
      opts.update({'ras1': [doccpi]})
      opts.update({'ras2': [soccpi]})
      opts.update({'ras3': [nmo - soccpi - doccpi]})
      opts.update({'ras4': [0]})
    # Other configuration spaces aren't supported yet
    else:
      print("Configuration space %s not supported. Exiting..." % conf_space)
      exit()

    # rotate orbitals for convergence
    # note that this ONLY works for ROHF orbitals right now!!
    # modify for UHF if needed later on
    # also probably need to add 1x and S cases (to include RAS1/3??)
    if(rotate_orbitals):
        if(conf_space == ""):
            # get S for orthogonalization
            S = psi4.core.Matrix.to_array(wfn_rohf.S(), copy=True)
            S_inv_sqrt = matrix_inv_sqrt(S, 1e-12)
            S_sqrt = matrix_sqrt(S, 1e-12)
            # get orthogonzlied Fock matrix
            Fa_np = psi4.core.Matrix.to_array(wfn_rohf.Fa(), copy=True)
            Fa_np = np.dot(S_inv_sqrt.T, np.dot(Fa_np, S_inv_sqrt))
            # get RAS 2 from orthogonalized orbitals
            orbs_uncanon = psi4.core.Matrix.to_array(wfn_rohf.Ca(), copy=True)
            orbs_full = np.dot(S_sqrt, orbs_uncanon)
            orbs = orbs_full[:, doccpi:doccpi+soccpi]
            # now, canonicalize orbitals (rotate RAS2 space)
            Fa_np = np.dot(orbs.T, np.dot(Fa_np, orbs))
            diag_f, vect_f = LIN.eigh(Fa_np)
            Ca_rotated = np.dot(orbs, vect_f)
            Cb_rotated = np.dot(orbs, vect_f)
            # concatenate
            Ca_full = np.column_stack((orbs_full[:, :doccpi], Ca_rotated, orbs_full[:, soccpi+doccpi:]))
            Cb_full = np.column_stack((orbs_full[:, :doccpi], Cb_rotated, orbs_full[:, soccpi+doccpi:]))
            # put back into unorthogonal form and set
            wfn_rohf.Ca().copy(psi4.core.Matrix.from_array(np.dot(S_inv_sqrt, Ca_full), name="Ca (Alpha)"))
            wfn_rohf.Cb().copy(psi4.core.Matrix.from_array(np.dot(S_inv_sqrt, Cb_full), name="Cb (Beta)"))
        else:
            # get S for orthogonalization
            S = psi4.core.Matrix.to_array(wfn_rohf.S(), copy=True)
            S_inv_sqrt = matrix_inv_sqrt(S, 1e-12)
            S_sqrt = matrix_sqrt(S, 1e-12)
            # get orthogonzlied Fock matrix
            Fa_np = psi4.core.Matrix.to_array(wfn_rohf.Fa(), copy=True)
            Fa_np = np.dot(S_inv_sqrt.T, np.dot(Fa_np, S_inv_sqrt))
            # get RAS 2 from orthogonalized orbitals
            orbs_uncanon = psi4.core.Matrix.to_array(wfn_rohf.Ca(), copy=True)
            orbs_full = np.dot(S_sqrt, orbs_uncanon)
            ras1 = orbs_full[:, :doccpi]
            ras2 = orbs_full[:, doccpi:doccpi+soccpi]
            ras3 = orbs_full[:, doccpi+soccpi:]
            # now, canonicalize orbitals (rotate RAS2 space)
            Fa_1 = np.dot(ras1.T, np.dot(Fa_np, ras1))
            Fa_2 = np.dot(ras2.T, np.dot(Fa_np, ras2))
            Fa_3 = np.dot(ras3.T, np.dot(Fa_np, ras3))
            diag_f1, vect_f1 = LIN.eigh(Fa_1)
            diag_f2, vect_f2 = LIN.eigh(Fa_2)
            diag_f3, vect_f3 = LIN.eigh(Fa_3)
            Ca_ras1 = np.dot(ras1, vect_f1)
            Ca_ras2 = np.dot(ras2, vect_f2)
            Ca_ras3 = np.dot(ras3, vect_f3)
            # concatenate
            Ca_full = np.column_stack((Ca_ras1, Ca_ras2, Ca_ras3))
            Cb_full = np.column_stack((Ca_ras1, Ca_ras2, Ca_ras3))
            # put back into unorthogonal form and set
            wfn_rohf.Ca().copy(psi4.core.Matrix.from_array(np.dot(S_inv_sqrt, Ca_full), name="Ca (Alpha)"))
            wfn_rohf.Cb().copy(psi4.core.Matrix.from_array(np.dot(S_inv_sqrt, Cb_full), name="Cb (Beta)"))
        

    # run cas
    print("RUNNING CAS...\t\tCHARGE %i\tMULT %i" %(mol.molecular_charge(), mol.multiplicity()))
    psi4.set_options(opts)
    e_cas, wfn_cas = energy('detci', ref_wfn=wfn_rohf, return_wfn=True, molecule=mol)
    print("CAS (%i %i): %6.12f" %(mol.molecular_charge(), mol.multiplicity(), e_cas))
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

