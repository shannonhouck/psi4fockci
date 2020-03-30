import psi4
import psi4.driver.p4util as p4util
import numpy as np
import numpy.linalg as LIN

def run_psi4fockci(name, molecule, **kwargs):
    """
    A method to run a RAS-nSF-IP/EA calculation.

    This runs a RAS-nSF-IP/EA calculation using Psi4's DETCI module. The 
    number of spin-flips and IP/EA is determined based on setting the 
    ``new_charge`` and ``new_multiplicity`` of the desired target state.
    Additional excitations are included by setting the ``conf_space`` 
    keyword; excitations through the CISDT level are currently supported.

    Parameters
    ----------
    name : str
        Name of method (for Psi4 interfacing)
    molecule : psi4.core.Molecule
        Molecule to run the calculation on.
    new_charge : int
        Target charge of the molecule.
    new_multiplicity : int
        Target multiplicity of the molecule.
    conf_space : str ("")
        Option for including additional excitations outside of the CAS space. 
        Defaults to CAS-nSF-IP/EA. Valid options include:
            * ``""`` CAS-nSF-IP/EA
            * ``"h"`` RAS(h)-nSF-IP/EA
            * ``"p"`` RAS(p)-nSF-IP/EA
            * ``"1x"`` RAS(h,p)-nSF-IP/EA
            * ``"S"`` RAS(S)-nSF-IP/EA
            * ``"SD"`` RAS(SD)-nSF-IP/EA
            * ``"SDT"`` RAS(SDT)-nSF-IP/EA
    add_opts : dict ({})
        Additional options to pass into Psi4.
    return_ci_wfn : bool (False)
        Whether to return the CI wavefunction object.
    return_rohf_wfn : bool (False)
        Whether to return the ROHF wavefunction object.
    return_rohf_e : bool (False)
        Whether to return the ROHF energy.
    read_rohf_wfn : bool (False)
        Whether to read a Psi4 ROHF wavefunction.
    rohf_wfn_in : psi4.core.Wavefunction
        The Psi4 ROHF reference wavefunction (pre-computed).
    write_rohf_wfn : str ("")
        Name of file (.npz) to write
    localize : bool (False)
        Perform BOYS localization on the RAS 2 space before DETCI call?
        Can help with visualization and analysis of orbitals.
    frozen_docc : int (0)
        Number of frozen core orbitals.
    frozen_vir : int (0)
        Number of frozen virtual orbitals.

    Returns
    -------
    return_ci_wfn : psi4.core.Wavefunction
        The SF-CAS([conf_space]) wavefunction.
    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    optstash = p4util.OptionsState(
        ['SCF', 'SCF_TYPE'],
        ['BASIS'],
        ['SCF', 'MAXITER'],
        ['DETCI', 'CI_MAXITER']
    )

    new_charge = kwargs.get('new_charge')
    new_multiplicity = kwargs.get('new_multiplicity')
    ref_mol = molecule

    if(not 'new_charge' in kwargs):
        print("ERROR: Please designate a target charge!")
        exit()
    if(not 'new_multiplicity' in kwargs):
        print("ERROR: Please designate a target multiplicity!")
        exit()

    conf_space = kwargs.get('conf_space', "")
    add_opts = kwargs.get('add_opts', {})
    read_rohf_wfn = kwargs.get('read_rohf_wfn', False)
    wfn_rohf_in = kwargs.get('wfn_rohf_in', None)
    write_rohf_wfn = kwargs.get('write_rohf_wfn', "")
    write_ci_vects = kwargs.get('write_ci_vects', False)
    localize = kwargs.get('localize', False)
    frozen_docc = kwargs.get('frozen_docc', 0)
    frozen_uocc = kwargs.get('frozen_vir', 0)

    print("Starting Psi4FockCI...\n")
    # default options for Psi4
    opts = {'scf_type': 'pk',
            'reference': 'rohf',
            'mixed': False,
            'maxiter': 1000,
            'ci_maxiter': 50 }
    opts.update(add_opts) # add additional options from user

    # run ROHF calculation on reference state or read it in
    psi4.core.clean()
    psi4.set_options(opts)
    mol = ref_mol.clone() # clone molecule so original isn't modified
    # read in ROHF guess wavefunction if provided
    if(read_rohf_wfn):
        # set up options and run
        psi4.set_options(opts)
        print("Using ROHF from reference...")
        wfn_rohf = wfn_rohf_in
        e_rohf = wfn_rohf.energy()
    # else, run ROHF on reference state
    else:
        print("Running reference...") 
        # TODO: Change to scf_helper call
        e_rohf, wfn_rohf = psi4.energy('scf', molecule=mol, return_wfn=True, 
                                  options=opts)
        print("SCF (%i %i): %6.12f" 
              %(mol.molecular_charge(), mol.multiplicity(), e_rohf))

    # saving npz file of wavefunction if needed
    if(write_rohf_wfn != ""):
        wfn_rohf.to_file(write_rohf_wfn)

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
    wfn_rohf.force_doccpi(psi4.core.Dimension([(int)((n_total-new_soccpi)/2)]))

    # if we need to localize...
    if(localize):
        C = psi4.core.Matrix.to_array(wfn_rohf.Ca(), copy=True)
        ras1_C = C[:, :doccpi]
        ras2_C = C[:, doccpi:doccpi+soccpi]
        ras3_C = C[:, doccpi+soccpi:]
        loc = psi4.core.Localizer.build('BOYS', wfn_rohf.basisset(), 
                                        psi4.core.Matrix.from_array(ras2_C))
        loc.localize()
        ras2_localized = psi4.core.Matrix.to_array(loc.L, copy=True)
        localized_orbs = np.column_stack((ras1_C, ras2_localized, ras3_C))
        new_Ca = psi4.core.Matrix.from_array(localized_orbs, name="Ca")
        new_Cb = psi4.core.Matrix.from_array(localized_orbs, name="Cb")
        wfn_rohf.Ca().copy(new_Ca)
        wfn_rohf.Cb().copy(new_Cb)

    # change charge and multiplicity to new target values
    n_sf = (ref_mol.multiplicity() - abs(del_electrons) - new_multiplicity)/2
    print("\nRunning RAS-SF-IP/EA...")
    print("  New Charge/Mult: (%i %i)" %(new_charge, new_multiplicity))
    print("  Spin-Flips: %i" % n_sf)
    print("  Electron Count: %i" % del_electrons)

    # set active space and docc space based on configuration space input
    # Regular CAS configuration space
    # includes only active space configurations
    if(conf_space == "" or conf_space == "CAS"):
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
    elif(conf_space == "1x" or conf_space == "h,p"):
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
    elif(conf_space == "s"):
        opts.update({'frozen_docc': [0]})
        opts.update({'ex_level': 1})
        opts.update({'frozen_docc': [frozen_docc]})
        opts.update({'ras1': [doccpi - frozen_docc]})
        opts.update({'ras2': [soccpi]})
        opts.update({'ras3': [nmo - soccpi - doccpi - frozen_uocc]})
        opts.update({'frozen_uocc': [frozen_uocc]})
        opts.update({'ras4': [0]})
    elif(conf_space == "sd"):
        opts.update({'frozen_docc': [0]})
        opts.update({'ex_level': 2})
        opts.update({'frozen_docc': [frozen_docc]})
        opts.update({'ras1': [doccpi - frozen_docc]})
        opts.update({'ras2': [soccpi]})
        opts.update({'ras3': [nmo - soccpi - doccpi - frozen_uocc]})
        opts.update({'frozen_uocc': [frozen_uocc]})
        opts.update({'ras4': [0]})
    elif(conf_space == "sdt"):
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
    psi4.set_options(opts)
    e_cas, wfn_cas = psi4.energy('detci', ref_wfn=wfn_rohf, return_wfn=True, 
                            molecule=mol)

    # printing useful info
    print("\n Root\tEnergy")
    print("-----------------------------------")
    n = 0
    while(psi4.core.has_variable("CI ROOT %i TOTAL ENERGY" % n)):
        n_str = "CI ROOT %i TOTAL ENERGY" % n
        e_n = psi4.core.variable(n_str)
        print(" %4i\t%6.12f" %(n, e_n))
        n = n + 1
    print("-----------------------------------\n")

    psi4.core.print_variables() # printing Psi4 variables

    # obtain eigenvectors if needed
    # partly based on Daniel Smith's answer on Psi4 forums
    if(write_ci_vects):
        wfn_cas_2 = psi4.core.CIWavefunction(wfn_rohf)
        C = np.zeros((wfn_cas_2.ndet(), n_roots))
        for i in range(n_roots):
            dvec = wfn_cas_2.new_civector(i+1, 53, True, True)
            dvec.set_nvec(i+1)
            dvec.init_io_files(True)
            dvec.read(i,0)
            C[:, i] = np.array(dvec)
        np.savetxt('ci_vect.txt', C)

    optstash.restore()

    print("Psi4FockCI complete. Have a good day!")

    return wfn_cas

