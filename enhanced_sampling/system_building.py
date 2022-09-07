import openmm
from openmm.app import CharmmParameterSet, CharmmPsfFile, PDBFile, PDBxFile
import json, os, yaml
from simtk import unit

def load_input_dir(input_dir, charmm_param_dir):
    assert os.path.exists(input_dir)

    ## load psf
    psf_path = os.path.join(input_dir, "step5_input.psf")
    psf = CharmmPsfFile(psf_path)

    param_paths = [os.path.join(charmm_param_dir, path) for path in os.listdir(charmm_param_dir)]
    params = CharmmParameterSet(*param_paths)

    cif_path = os.path.join(input_dir, "equilibrated.cif")
    cif = PDBxFile(cif_path)

    dat_file_path = os.path.join(input_dir, "sysinfo.dat")

    x, y, z = load_xyz_from_datfile(dat_file_path)

    psf.setBox(x, y, z)

    return {"psf": psf, "params": params, "cif": cif}

def load_xyz_from_datfile(dat_file_path):
    with open(dat_file_path) as file:
        data = json.load(file)
        x, y, z = map(float, data['dimensions'][:3]) * unit.angstroms
    print(f'Setting box size x: {x}, y: {y}, z:{z}')
    return (x,y,z)

def load_simulation_params(param_file):
    with open(param_file) as f:
        param_dict = yaml.safe_load(f)

    params_dict_with_units = {}
    params_dict_with_units['hydrogen_mass'] = param_dict['hydrogen_mass'] * unit.amu
    params_dict_with_units['temperature'] = param_dict['temperature'] * unit.kelvin
    params_dict_with_units['friction'] = param_dict['friction'] / unit.picoseconds
    params_dict_with_units['time_step'] = param_dict['time_step'] * unit.picoseconds
    params_dict_with_units['pressure'] = param_dict['pressure'] * unit.bar
    params_dict_with_units['surface_tension'] = param_dict['surface_tension']

    params_dict_with_units['nonbonded_method'] = getattr(openmm.app, param_dict['nonbonded_method'])
    params_dict_with_units['constraints'] = getattr(openmm.app, param_dict['constraints'])

    return params_dict_with_units
