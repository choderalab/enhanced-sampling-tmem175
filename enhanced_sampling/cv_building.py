import openmm, openmmtools

def create_rmsd_restraint(positions, atom_indicies):
    rmsd_cv = openmm.RMSDForce(positions, atom_indicies)
    energy_expression = 'step(dRMSD) * (K_RMSD/2) * dRMSD^2; dRMSD = (RMSD-RMSD0);'
    energy_expression += 'K_RMSD = %f;' % spring_constant.value_in_unit_system(md_unit_system)
    energy_expression += 'RMSD0 = %f;' % restraint_distance.value_in_unit_system(md_unit_system)
    restraint_force = CustomCVForce(energy_expression)
    restraint_force.addCollectiveVariable('RMSD', rmsd_cv)
    return restraint_force