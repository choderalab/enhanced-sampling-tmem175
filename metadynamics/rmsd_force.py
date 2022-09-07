## IMPORTS
from argparse import ArgumentParser
import os, sys
import openmm as mm
import openmm.app

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

    utils.print_args(args)

    utils.prep_output_dir(args.output_dir)
    input_dict = sb.load_input_dir(args.input_dir, args.charmm_param_dir)
    print(input_dict)
    psf = input_dict["psf"]
    print(psf.topology)

    cif = input_dict['cif']

    param_dict_with_units = sb.load_simulation_params(args.params_file)
    print(param_dict_with_units)

    system = psf.createSystem(input_dict["params"],
                              nonbondedMethod=param_dict_with_units['nonbonded_method'],
                              constraints=param_dict_with_units['constraints'],
                              removeCMMotion=False,
                              hydrogenMass=param_dict_with_units['hydrogen_mass'])

    integrator = mm.LangevinMiddleIntegrator(param_dict_with_units['temperature'],
                                             param_dict_with_units['friction'],
                                             param_dict_with_units['time_step'])

    barostat = mm.MonteCarloMembraneBarostat(param_dict_with_units['pressure'],
                                             param_dict_with_units['surface_tension'],
                                             param_dict_with_units['temperature'],
                                             mm.MonteCarloMembraneBarostat.XYIsotropic,
                                             mm.MonteCarloMembraneBarostat.ZFree
                                             )
    barostat.setFrequency(50)  ## for some reason the __init__ won't accept it as an argument, but this works
    ## the default is 25 timesteps, i've set it for 50
    system.addForce(barostat)

    vbond_force = sb.build_virtual_bond(psf, param_dict_with_units)

    system.addForce(vbond_force)

    platform = sb.get_platform_from_params(param_dict_with_units)

    # ref_dict = sb.load_input_dir(args.reference_dir)
    # ref_psf = ref_dict['psf']
    # ref_cif = ref_dict['cif']
    #
    # assert len(ref_cif.positions) == len(cif.positions)
    #
    # idx = [atom.index for atom in ref_psf.topology.atoms() if
    #        atom.residue.chain.index == 0 or atom.residue.chain.index == 1]
    #
    # rmsd_force = openmm.RMSDForce(ref_cif.positions, idx)
    #
    # rmsd_bias = openmm.app.metadynamics.BiasVariable(force=rmsd_force,
    #                                                  minValue=0.01,
    #                                                  maxValue=0.4,
    #                                                  biasWidth=0.02,
    #                                                  periodic=False)
    #
    # meta = openmm.app.Metadynamics(system=system,
    #                                variables=[rmsd_bias],
    #                                temperature=param_dict_with_units['temperature'],
    #                                biasFactor=2,
    #                                height=1,
    #                                frequency=1,
    #                                saveFrequency=1,
    #                                biasDir=args.output_dir
    #                                )

    sim = openmm.app.Simulation(psf.topology,
                                system=system,
                                integrator=integrator,
                                platform=platform)

    sim.context.setPositions(cif.positions)
    sim.context.setVelocitiesToTemperature(param_dict_with_units['temperature'])
    # print(meta.getCollectiveVariables(sim))
    # print("Minimizing energy")
    # sim.minimizeEnergy()
    # meta.step(sim, 10)
    print("Running simulation")
    sim.step(10)
    # print(meta.getCollectiveVariables(sim))

    print(f"Writing simulation files to {args.output_dir}")
    ss.write_simulation_files(sim, args.output_dir)

    print("Script cleanup")
    utils.save_env()
    utils.write_to_log(args,
                       os.path.basename(__file__))


## RUN COMMAND
if __name__ == "__main__":
    main()
