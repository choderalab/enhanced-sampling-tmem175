import yaml, openmm.unit as unit, openmm

class SimParams(object):
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
        repr_list = [
            f"{key:20}:\t{value.in_units_of(unit.nanoseconds).format('%04.2f'):20}{value.in_units_of(unit.picoseconds).format('%04.2f'):20}{value.in_units_of(unit.femtoseconds).format('%04.2f'):20}"
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
            self.selection = param_dict.get('selection')
            self.min_value = param_dict.get('min_value')
            self.max_value = param_dict.get('max_value')
            self.bias_width = param_dict.get('bias_width')

    def __str__(self):
        repr_list = [f"{key:20}: {value}" for key, value, in self.__dict__.items()]
        repr_string = "\n".join(repr_list)
        return repr_string + '\n'

class PullingParams(object):
    """
    Parameter object for pulling simulations.
    """