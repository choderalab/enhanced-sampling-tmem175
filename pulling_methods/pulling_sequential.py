## IMPORTS
from argparse import ArgumentParser
import os, sys
import openmm, mdtraj
import openmm.app
import logging
import shutil

repo_path = f"{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}"
sys.path.append(repo_path)
from enhanced_sampling import (
    utils,
    system_building as sb,
    system_saving as ss,
    reporters,
    cv_building as cv,
    schema,
)


## ARGUMENT PARSING
def get_args():
    parser = ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_dir",
        type=str,
        help="Input directory",
    )
    parser.add_argument(
        "-e",
        "--ebdims_file",
        type=str,
        help="Path to file with ebdims trajectory",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        help="Output directory, should be outside of this repo",
    )
    parser.add_argument(
        "-c",
        "--charmm_param_dir",
        type=str,
        default="../forcefield_data/charmm-params",
        help="Directory containing charmm parameter files",
    )
    parser.add_argument(
        "-p",
        "--params_file",
        type=str,
        default="../sim_params/defaults.yaml",
        help="YAML file containing simulation running values",
    )
    parser.add_argument(
        "-y",
        "--pulling_params",
        type=str,
        default="./cas.yaml",
        help="Path to yaml file containing residue list and min_max values",
    )
    parser.add_argument(
        "--restart_from_frame",
        type=int,
        default=None,
        help="Frame being used to restart the pulling simulations.",
    )
    parser.add_argument(
        "-r",
        "--restart_dir",
        type=str,
        help="Restart input directory. Assumed to contain the state file for restarting",
    )
    args = parser.parse_args()
    return args


## MAIN SCRIPT
def main(args, input_dict, ref_positions, output_dir):
    print(input_dict.keys())
    psf = input_dict["psf"]
    print(psf.topology)
    # positions = input_dict["positions"]

    params = schema.SimParams(args.params_file)
    print(params)
    pulling_params = schema.PullingParams(args.pulling_params)
    print(pulling_params)
    system = psf.createSystem(
        input_dict["params"],
        nonbondedMethod=params.nonbonded_method,
        constraints=params.constraints,
        removeCMMotion=False,
        hydrogenMass=params.hydrogen_mass,
    )

    integrator = openmm.LangevinMiddleIntegrator(
        params.temperature, params.friction, params.time_step
    )

    barostat = openmm.MonteCarloMembraneBarostat(
        params.pressure,
        params.surface_tension,
        params.temperature,
        openmm.MonteCarloMembraneBarostat.XYIsotropic,
        openmm.MonteCarloMembraneBarostat.ZFree,
    )

    # for some reason the __init__ won't accept it as an argument, but this works
    # the default is 25 timesteps, i've set it for 50
    barostat.setFrequency(50)

    system.addForce(barostat)

    vbond_force = sb.build_virtual_bond(psf, params)

    system.addForce(vbond_force)

    platform = sb.get_platform_from_params(params)

    ## Implement Enhanced Sampling here! ###############################################################################
    ####################################################################################################################
    idx_list = cv.get_openmm_idx(
        psf.topology, selection=pulling_params.selection
    )

    assert len(ref_positions) == len(idx_list)

    restraint_positions = {
        idx_list[i]: ref_positions[i] for i in range(len(idx_list))
    }

    force = cv.create_harmonic_pulling_force(
        restraint_positions, pulling_params.spring_constant
    )
    force.setForceGroup(pulling_params.force_group)

    print(
        f"Adding {force.getName()} to system with {force.getNumParticles()} particles"
    )
    force_idx = system.addForce(force)

    ###################################################################################################################

    for force in system.getForces():
        if force.getForceGroup() > 6:
            print(force.getName(), force.getForceGroup())

    sim = openmm.app.Simulation(
        psf.topology, system=system, integrator=integrator, platform=platform
    )

    ## Set state (positions, box vectors, velocities)
    if input_dict["state"]:
        print(f"Setting state from state file...")
        sim.context.setState(input_dict["state"])
    elif input_dict["positions"]:
        print(f"Setting positions and velocities...")
        sim.context.setPositions(input_dict["positions"])
        sim.context.setVelocitiesToTemperature(params.temperature)
        print(
            "  before minimization : %8.3f kcal/mol"
            % (
                sim.context.getState(getEnergy=True).getPotentialEnergy()
                / openmm.unit.kilocalories_per_mole
            )
        )
        sim.minimizeEnergy()

    ## reset time and step count to 0
    sim.context.setTime(0)
    sim.context.setStepCount(0)

    print(
        "  before simulation : %8.3f kcal/mol"
        % (
            sim.context.getState(getEnergy=True).getPotentialEnergy()
            / openmm.unit.kilocalories_per_mole
        )
    )

    sim.reporters.append(
        openmm.app.StateDataReporter(
            os.path.join(output_dir, "trajectory.log"),
            reportInterval=params.report_freq,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            temperature=True,
            speed=True,
            progress=True,
            remainingTime=True,
            totalSteps=params.n_steps,
            separator="\t",
        )
    )

    sim.reporters.append(
        openmm.app.CheckpointReporter(
            file=os.path.join(output_dir, "trajectory.chk"),
            reportInterval=params.chk_freq,
        )
    )

    # Write out the trajectory
    sim.reporters.append(
        mdtraj.reporters.XTCReporter(
            file=os.path.join(output_dir, "trajectory.xtc"),
            reportInterval=params.traj_freq,
        )
    )

    sim.reporters.append(
        reporters.CustomForceReporter(
            file=os.path.join(output_dir, "forces.txt"),
            reportInterval=params.report_freq,
            force_group=pulling_params.force_group,
            atom_idx=idx_list,
        )
    )

    sim.reporters.append(
        reporters.CustomEnergyReporter(
            file=os.path.join(output_dir, "energy.txt"),
            reportInterval=params.report_freq,
            force_group=pulling_params.force_group,
        )
    )

    print("Running simulation")
    sim.step(params.n_steps)

    print(f"Writing simulation files to {output_dir}")
    ss.write_simulation_files(sim, output_dir)

    print("Script cleanup")
    utils.save_env()
    utils.write_to_log(args, os.path.basename(__file__))

    return output_dir


## RUN COMMAND
if __name__ == "__main__":
    args = get_args()

    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir)
    os.mkdir(args.output_dir)

    logging.basicConfig(level=logging.DEBUG, format="%(message)s")
    logger = logging.getLogger()
    logger.addHandler(
        logging.FileHandler(os.path.join(args.output_dir, "log.txt"), "a")
    )
    sys.stdout.write = logger.info

    print(f"Writing to {args.output_dir}")

    utils.print_args(args)

    t = mdtraj.load(args.ebdims_file)
    position_list = t.xyz

    input_dict = sb.load_input_dir(
        args.input_dir, args.charmm_param_dir, position_source="state"
    )

    ## options for if we're restarting from a frame
    if args.restart_from_frame:
        frames = range(len(position_list))[args.restart_from_frame + 1 :]
        new_input_dir = args.restart_dir
        print(f"Restarting from {new_input_dir} using frames {frames}")
    else:
        frames = range(len(position_list))
        new_input_dir = None

    for i in frames:
        print(f"Using positions from frame {i} from {args.ebdims_file}")
        ref_positions = position_list[i]

        output_dir = os.path.join(args.output_dir, f"frame{i:02d}")
        os.mkdir(output_dir)
        print(f"Writing to {output_dir}")

        if new_input_dir:
            new_input_dict = sb.load_input_dir(
                new_input_dir,
                load_psf=False,
                position_source="state",
            )
            input_dict["state"] = new_input_dict["state"]

        new_input_dir = main(
            args=args,
            input_dict=input_dict,
            ref_positions=ref_positions,
            output_dir=output_dir,
        )
        print(f"Trajectory info written to {new_input_dir}")
