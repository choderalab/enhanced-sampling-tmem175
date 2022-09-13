import openmm, openmmtools
import mdtraj

# def create_rmsd_restraint(positions, atom_indicies):
#     rmsd_cv = openmm.RMSDForce(positions, atom_indicies)
#     energy_expression = 'step(dRMSD) * (K_RMSD/2) * dRMSD^2; dRMSD = (RMSD-RMSD0);'
#     energy_expression += 'K_RMSD = %f;' % spring_constant.value_in_unit_system(md_unit_system)
#     energy_expression += 'RMSD0 = %f;' % restraint_distance.value_in_unit_system(md_unit_system)
#     restraint_force = CustomCVForce(energy_expression)
#     restraint_force.addCollectiveVariable('RMSD', rmsd_cv)
#     return restraint_force


def get_openmm_idx(topology: openmm.app.topology.Topology, selection, res_list=False):
    """
    Filter based on selection and then further filter based on a passed in residue list.

    :param topology:
    :param selection:
    :param res_list:
    :return:
    """
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
