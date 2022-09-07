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

    platform = mm.Platform.getPlatformByName(param_dict_with_units['platform'])

    platform.setPropertyDefaultValue('Precision', param_dict_with_units['precision'])
    platform.setPropertyDefaultValue('DeterministicForces', 'true')

    sim = openmm.app.Simulation(psf.topology,
                                system=system,
                                integrator=integrator,
                                platform=platform)

    sim.context.setPositions(input_dict['cif'].positions)
    sim.context.setVelocitiesToTemperature(param_dict_with_units['temperature'])

    ss.write_simulation_files(sim, args.output_dir)


    utils.save_env()
    utils.write_to_log(args,
                       os.path.basename(__file__))




## RUN COMMAND
if __name__ == "__main__":
    main()