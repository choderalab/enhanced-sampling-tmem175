import openmm, json, os, yaml, mdtraj, bz2
from openmm.app import CharmmParameterSet, CharmmPsfFile, PDBFile, PDBxFile
from simtk import unit


def load_input_dir(input_dir,
                   charmm_param_dir=None,
                   load_psf=True,
                   from_state=True
                   ):
    assert os.path.exists(input_dir)

    ## load psf
    if load_psf:
        psf_path = os.path.join(input_dir, "step5_input.psf")
        print(f"Loading psf from {psf_path}")
        psf = CharmmPsfFile(psf_path)
    else:
        psf = None

    if charmm_param_dir:
        print(f"Loading Charmm params from {charmm_param_dir}")
        param_paths = [os.path.join(charmm_param_dir, path) for path in os.listdir(charmm_param_dir)]
        params = CharmmParameterSet(*param_paths)
    else:
        params = None

    if from_state:
        input_state_path = os.path.join(input_dir, "state.xml.bz2")
        with bz2.open(input_state_path, 'rb') as infile:
            state = openmm.XmlSerializer.deserialize(infile.read().decode())

        if load_psf:
            x, y, z = state.getPeriodicBoxVectors()
            psf.setBox(x[0], y[1], z[2])

        positions = state.getPositions()
    else:
        state = None
        positions = None

    return {"psf": psf, "params": params, "state": state, "positions": positions}


def load_xyz_from_datfile(dat_file_path):
    with open(dat_file_path) as file:
        data = json.load(file)
        x, y, z = map(float, data['dimensions'][:3]) * unit.angstroms
    print(f'Setting box size x: {x}, y: {y}, z:{z}')
    return (x, y, z)


class SimParams():
    def __init__(self, param_file):
        print(f"Loading sim_params from {param_file}")
        with open(param_file) as f:
            param_dict = yaml.safe_load(f)

        ## Parameters with units
        self.hydrogen_mass = param_dict['hydrogen_mass'] * unit.amu
        self.temperature = param_dict['temperature'] * unit.kelvin
        self.friction = param_dict['friction'] / unit.picoseconds
        self.time_step = param_dict['time_step'] * unit.picoseconds
        self.pressure = param_dict['pressure'] * unit.bar

        ## Parameters without units
        self.surface_tension = param_dict['surface_tension']
        self.n_steps = param_dict["n_steps"]
        self.report_freq = param_dict["report_freq"]
        self.chk_freq = param_dict["chk_freq"]
        self.traj_freq = param_dict["traj_freq"]

        ## Not required parameters
        self.virtual_bond = param_dict.get("virtual_bond")
        self.platform = param_dict.get("platform")
        self.precision = param_dict.get("precision")

        ## Parameters that are actually objects
        self.nonbonded_method = getattr(openmm.app, param_dict['nonbonded_method'])
        self.constraints = getattr(openmm.app, param_dict['constraints'])

        self.steps_dict = {
            "Total Time": self.n_steps,
            "Report Freq": self.report_freq,
            "Checkpoint Freq": self.chk_freq,
            "Trajectory Freq": self.traj_freq
        }

        self.time_dict = {key: value * self.time_step
                          for key, value in self.steps_dict.items()}
        self.n_frames = int(self.n_steps / self.traj_freq)

    def __str__(self):
        repr_list = [f"{key:20}: {value}" for key, value, in self.__dict__.items()]
        repr_string = "\n".join(repr_list)
        return repr_string + '\n' + self.get_time_scales()

    def get_time_scales(self):
        repr_list = [f"{key:20}:\t{value.in_units_of(unit.nanoseconds).format('%04.2f'):20}{value.in_units_of(unit.picoseconds).format('%04.2f'):20}{value.in_units_of(unit.femtoseconds).format('%04.2f'):20}"
                     for key, value, in self.time_dict.items()]

        repr_list.append(f"{'Number of Frames':20}:\t{self.n_frames}")
        repr_str = "\n".join(repr_list)
        return repr_str

class MetaParams(object):
    """
    Parameter object for metadynamics
    """
    def __init__(self, param_file):
        print(f"Loading meta_params from {param_file}")

        with open(param_file) as f:
            param_dict = yaml.safe_load(f)
            self.res_list = param_dict.get('res_list')
            self.rmsd_sel = param_dict.get('selection')
            self.min_value = param_dict.get('min_value')
            self.max_value = param_dict.get('max_value')
            self.bias_width = param_dict.get('bias_width')
            self.spring_constant = param_dict.get('spring_constant')
            self.rmsd_max = param_dict.get('rmsd_max')

    def __str__(self):
        repr_list = [f"{key:20}: {value}" for key, value, in self.__dict__.items()]
        repr_string = "\n".join(repr_list)
        return repr_string + '\n'

def build_virtual_bond(psf, params):
    vbond_selections = params.virtual_bond
    print(vbond_selections)
    topology = mdtraj.Topology.from_openmm(psf.topology)
    print(topology)
    atom1 = topology.select(vbond_selections[0])[0]
    atom2 = topology.select(vbond_selections[1])[0]

    force = openmm.CustomBondForce("0*r")
    force.addBond(int(atom1), int(atom2))

    return force


def get_platform_from_params(params):
    platform_name = params.platform
    platform = openmm.Platform.getPlatformByName(platform_name)

    if not platform_name == 'CPU':
        print(f"Setting precision to {params.precision}")
        platform.setPropertyDefaultValue('Precision', params.precision)
        print(f"Setting DeterministicForces to 'true`")
        platform.setPropertyDefaultValue('DeterministicForces', 'true')

    return platform
