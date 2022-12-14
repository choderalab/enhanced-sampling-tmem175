import openmm, yaml
import openmm.app.topology as topology


def create_rmsd_restraint(positions,
                          atom_indicies,
                          spring_constant,
                          rmsd_max):
    """
    The goal is to make a flat-bottom, 1 side harmonic potential

    :param positions:
    :param atom_indicies:
    :param spring_constant:
    :param rmsd_max:
        Max RMSD value
    :return:
    """
    rmsd_cv = openmm.RMSDForce(positions, atom_indicies)
    energy_expression = f"(spring_constant/2)*max(0, RMSD-RMSDmax)^2"
    restraint_force = openmm.CustomCVForce(energy_expression)
    restraint_force.addCollectiveVariable('RMSD', rmsd_cv)
    restraint_force.addGlobalParameter('RMSDmax', rmsd_max)
    restraint_force.addGlobalParameter("spring_constant", spring_constant)
    return restraint_force


def build_rmsd_restraint_from_yaml(file_path,
                                   positions,
                                   topology):
    with open(file_path) as f:
        restraint_dict = yaml.safe_load(f)

    restraint_idx = get_openmm_idx(topology, restraint_dict["selection"])
    print(restraint_dict["spring_constant"], restraint_dict["rmsd_max"])
    rmsd_restraint_force = create_rmsd_restraint(positions=positions,
                                                 atom_indicies=restraint_idx,
                                                 spring_constant=restraint_dict["spring_constant"],
                                                 rmsd_max=restraint_dict["rmsd_max"]
                                                 )
    rmsd_restraint_force.setForceGroup(restraint_dict["force_group"])
    return rmsd_restraint_force


def create_harmonic_pulling_force(restraints_dict, spring_constant):
    """
    Expects a dictionary that looks like {atom_idx: (x0, y0, z0), ...}
    :return:
    """
    force = openmm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", spring_constant)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    for atom_idx, (x, y, z) in restraints_dict.items():
        force.addParticle(atom_idx, [x, y, z])
    return force


def get_openmm_idx(topology: topology.Topology, selection, res_list=False):
    """
    Filter based on selection and then further filter based on a passed in residue list.

    :param topology:
    :param selection:
    :param res_list:
    :return:
    """
    print(f"Using protein selection: {selection} and res_list: {res_list}")
    if selection == "protein":
        idx = [atom.index for atom in topology.atoms() if
               atom.residue.chain.index == 0 or atom.residue.chain.index == 1]
    elif selection == "protein_heavy":
        idx = [atom.index for atom in topology.atoms() if
               (atom.residue.chain.index == 0 or atom.residue.chain.index == 1) and atom.element.symbol != 'H']
    elif selection == "protein_ca":
        idx = [atom.index for atom in topology.atoms() if
               (atom.residue.chain.index == 0 or atom.residue.chain.index == 1) and atom.name == 'CA']
    else:
        raise NotImplementedError

    if res_list:
        atom_idx_list = [atom.index for atom in topology.atoms() if
                         atom.residue.id in res_list]
        idx = [index for index in idx if index in atom_idx_list]

    return idx
