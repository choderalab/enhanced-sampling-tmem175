## IMPORTS
from argparse import ArgumentParser
import os, sys

repo_path = f"{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}"
sys.path.append(repo_path)
from enhanced_sampling import utils, system_building as sb

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

    units_dict = sb.load_simulation_params(args.params_file)
    print(units_dict)

    utils.save_env()
    utils.write_to_log(args,
                       os.path.basename(__file__))


## RUN COMMAND
if __name__ == "__main__":
    main()