## IMPORTS
from argparse import ArgumentParser
import os, sys

## Add repo path
repo_path = f"{os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}"
sys.path.append(repo_path)
from enhanced_sampling import utils

## ARGUMENT PARSING
def get_args():
    parser = ArgumentParser()
    parser.add_argument("-y", "--param_yaml", type=str,
                        help="YAML file containing params for climber.", )

    args = parser.parse_args()
    return args

## OTHER FUNCTIONS



## MAIN SCRIPT
def main():
    args = get_args()
    utils.print_args(args)

    # utils.save_env()
    # utils.write_to_log(args,
    #                    os.path.basename(__file__))


## RUN COMMAND
if __name__ == "__main__":
    main()