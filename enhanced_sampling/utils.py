import argparse
import os, datetime, yaml


def args_to_dict(args):
    if type(args) == dict:
        arg_dict = args
    elif isinstance(args, argparse.Namespace):
        arg_dict = vars(args)
    else:
        raise NotImplementedError
    return arg_dict


def print_args(args):
    for k, v in args_to_dict(args).items():
        print(f"{k}: {v}")


def determine_log_file():
    log_path = os.path.join("../log", f"{datetime.date.today().isoformat()}_log.yaml")
    return log_path

def prep_output_dir(output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    assert os.path.exists(output_dir)


def load_yaml(file_path):
    file_path = os.path.join(file_path)
    with open(file_path) as file:
        yaml_dict = yaml.safe_load(file)
    if not yaml_dict:
        yaml_dict = {}
    return yaml_dict


def write_yaml(yaml_dict, file_path):
    with open(file_path, "w") as file:
        yaml.safe_dump(yaml_dict, file)


def run_lsf_script(lsf_script):
    command = f"bsub < {lsf_script}"
    print(f"Running '{command}'")
    os.system(command)


def save_env():
    os.system(f"conda env export --from-history > ../devtools/{datetime.date.today().isoformat()}_environment.yaml")


def write_to_log(args, python_script_name):
    print("Writing log")
    log_path = determine_log_file()
    if not os.path.exists(log_path):
        yaml_dict = {}
    else:
        yaml_dict = load_yaml(log_path)
    time = datetime.datetime.now().isoformat()
    yaml_dict[time] = {"args": args_to_dict(args),
                       "script": python_script_name,
                       "notes": '',
                       "success": ''}
    write_yaml(yaml_dict, log_path)

def convert_mdtraj_resid_to_resn(resid: int, n_residues=447):
    ## the h5 files have res 0 as ACE29 and res 1 as ILE30
    ## and residue 450 is ILE30 of chain B
    ## so the numbering is slightly different depending on the chain
    if resid <= n_residues:
        chain = 'A'
        res = resid + 29
    else:
        chain = 'B'
        res = resid - n_residues + 27
#     return f"{res}{chain}"
    return res

def get_res_array(np_array):
    res_array = []
    for pair in np_array:
        res_list = [convert_mdtraj_resid_to_resn(res) for res in pair]
        res_array.append(res_list)
    return res_array
