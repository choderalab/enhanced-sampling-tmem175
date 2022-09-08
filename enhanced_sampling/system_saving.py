import openmm, bz2, os
from openmm.app import PDBxFile

def write_simulation_files(sim, output_dir):
    ## get state, system, integrator
    state = sim.context.getState(
        getPositions=True,
        getVelocities=True,
        getEnergy=True,
        getForces=True,
        enforcePeriodicBox=True,
        getParameters=True
    )
    system = sim.context.getSystem()
    integrator = sim.context.getIntegrator()

    output_state_path = os.path.join(output_dir, 'state.xml.bz2')
    output_pdb_path = os.path.join(output_dir, 'final_frame.cif')
    output_system_file = os.path.join(output_dir, 'system.xml.bz2')
    output_integrator_file = os.path.join(output_dir, 'integrator.xml.bz2')


    # Save and serialize the final state
    print("Saving state")
    with bz2.open(output_state_path, "wt") as outfile:
        xml = openmm.XmlSerializer.serialize(state)
        outfile.write(xml)

    # Save the final state as a PDBx File
    print("Saving cif")
    with open(output_pdb_path, "wt") as outfile:
        PDBxFile.writeFile(
            sim.topology,
            state.getPositions(),
            file=outfile,
            keepIds=True
        )

    # Save and serialize system
    print("Saving system")
    sim.system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
    with bz2.open(output_system_file, "wt") as outfile:
        xml = openmm.XmlSerializer.serialize(system)
        outfile.write(xml)

    print("Saving integrator")
    # Save and serialize integrator
    with bz2.open(output_integrator_file, "wt") as outfile:
        xml = openmm.XmlSerializer.serialize(integrator)
        outfile.write(xml)