## IMPORTS
from argparse import ArgumentParser
import os, sys
import openmm
import openmm.app
import logging

repo_path = f"{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}"
sys.path.append(repo_path)
from enhanced_sampling import utils, system_building as sb, system_saving as ss


## ARGUMENT PARSING
def get_args():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_dir", type=str,
                        help="Input directory", )
    parser.add_argument("-r", "--reference_dir", type=str,
                        help="Directory containing reference pose")
    parser.add_argument("-o", "--output_dir", type=str,
                        help="Output directory, should be outside of this repo")
    parser.add_argument("-c", "--charmm_param_dir", type=str, default="../forcefield_data/charmm-params",
                        help="Directory containing charmm parameter files")
    parser.add_argument("-p", "--params_file", type=str, default="../sim_params/defaults.yaml",
                        help="YAML file containing simulation running values")
    args = parser.parse_args()
    return args


## MAIN SCRIPT
def main():
    args = get_args()
    output_dir = utils.prep_output_dir(args.output_dir)

    logging.basicConfig(level=logging.DEBUG, format='%(message)s')
    logger = logging.getLogger()
    logger.addHandler(logging.FileHandler(os.path.join(output_dir, "log.txt"), 'a'))
    # sys.stderr.write = logger.error
    sys.stdout.write = logger.info

    utils.print_args(args)


    input_dict = sb.load_input_dir(args.input_dir, args.charmm_param_dir)
    print(input_dict.keys())
    psf = input_dict["psf"]
    print(psf.topology)

    params = sb.SimParams(args.params_file)
    print(params)

    system = psf.createSystem(input_dict["params"],
                              nonbondedMethod=params.nonbonded_method,
                              constraints=params.constraints,
                              removeCMMotion=False,
                              hydrogenMass=params.hydrogen_mass)

    integrator = openmm.LangevinMiddleIntegrator(params.temperature,
                                                 params.friction,
                                                 params.time_step)

    barostat = openmm.MonteCarloMembraneBarostat(params.pressure,
                                                 params.surface_tension,
                                                 params.temperature,
                                                 openmm.MonteCarloMembraneBarostat.XYIsotropic,
                                                 openmm.MonteCarloMembraneBarostat.ZFree
                                                 )

    # for some reason the __init__ won't accept it as an argument, but this works
    # the default is 25 timesteps, i've set it for 50
    barostat.setFrequency(50)

    system.addForce(barostat)

    vbond_force = sb.build_virtual_bond(psf, params)

    system.addForce(vbond_force)

    platform = sb.get_platform_from_params(params)

    ref_dict = sb.load_input_dir(args.reference_dir)
    ref_psf = ref_dict['psf']
    ref_positions = ref_dict['positions']

    assert len(ref_positions) == len(input_dict["positions"])

    idx = [atom.index for atom in ref_psf.topology.atoms() if
           atom.residue.chain.index == 0 or atom.residue.chain.index == 1]

    rmsd_force = openmm.RMSDForce(ref_positions, idx)

    rmsd_bias = openmm.app.metadynamics.BiasVariable(force=rmsd_force,
                                                     minValue=0.01,
                                                     maxValue=0.4,
                                                     biasWidth=0.02,
                                                     periodic=False)

    meta = openmm.app.Metadynamics(system=system,
                                   variables=[rmsd_bias],
                                   temperature=params.temperature,
                                   biasFactor=5,
                                   height=1,
                                   frequency=1,
                                   saveFrequency=10,
                                   biasDir=output_dir
                                   )

    sim = openmm.app.Simulation(psf.topology,
                                system=system,
                                integrator=integrator,
                                platform=platform)
    sim.context.setState(input_dict["state"])

    print(meta.getCollectiveVariables(sim))

    print(
        "  initial : %8.3f kcal/mol"
        % (
                sim.context.getState(getEnergy=True).getPotentialEnergy()
                / openmm.unit.kilocalories_per_mole
        )
    )

    sim.reporters.append(
        openmm.app.StateDataReporter(
            os.path.join(output_dir, "simulation_log.txt"),
            reportInterval=params.report_freq,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            temperature=True,
            speed=True,
            progress=True,
            remainingTime=True,
            totalSteps=params.nsteps,
            separator="\t",
        )
    )

    print("Running simulation")
    meta.step(sim, params.nsteps)
    print(meta.getCollectiveVariables(sim))

    print(f"Writing simulation files to {output_dir}")
    ss.write_simulation_files(sim, output_dir)

    print("Script cleanup")
    utils.save_env()
    utils.write_to_log(args,
                       os.path.basename(__file__))


## RUN COMMAND
if __name__ == "__main__":
    main()
