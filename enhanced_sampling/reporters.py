import openmm.unit as unit
import os
import numpy

class CustomReporter(object):
    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        """
        The return value should be a six element tuple, whose elements are as follows:

        The number of time steps until the next report.
        We calculate this as (report interval)-(current step)%(report interval).
        For example, if we want a report every 100 steps and the simulation is currently on step 530,
        we will return 100-(530%100) = 70.

        Whether the next report will need particle positions.

        Whether the next report will need particle velocities.

        Whether the next report will need forces.

        Whether the next report will need energies.

        Whether the positions should be wrapped to the periodic box.
        If None, it will automatically decide whether to wrap positions based on whether
        the System uses periodic boundary conditions.
        """
        pass

    def get_formated_str(self, out_list):
        formated_list = [f"{item:30}" for item in out_list]
        out_str = "\t".join(formated_list) + "\n"
        return out_str

    def write_header(self):
        self._out.write(self.get_formated_str(self.header_list))

class MetadynamicsReporter(CustomReporter):
    def __init__(self, collective_variable_file, reportInterval, meta):
        self._out = open(collective_variable_file, 'w')
        self._reportInterval = reportInterval
        self._meta = meta

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, None)

    def report(self, simulation, state):
        collective_variables = self._meta.getCollectiveVariables(simulation)
        cv_list = [str(cv) for cv in collective_variables]
        cv_str = ", ".join(cv_list)
        self._collective_variable_file.write(f"{cv_str}\n")


class CustomCVForceReporter(CustomReporter):
    """
    From <http://docs.openmm.org/latest/userguide/application/
    04_advanced_sim_examples.html#extracting-and-reporting-forces-and-other-data>
    """

    def __init__(self, file, reportInterval, force_group, force_idx):
        self._out = open(file, 'w')
        self._reportInterval = reportInterval
        self._force_group = force_group
        self._force_idx = force_idx
        self.header_list = ["Potential Energy", "Kinetic Energy", "CustomCV", "Forces"]
        self.write_header()

        ## Since I don't think I will every use bitmaps for this, enforce this to be a set
        if type(self._force_group) == int:
            self._force_group = {self._force_group}
        elif type(self._force_group) == list:
            self._force_group = set(self._force_group)

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, None)

    def report(self, simulation, state):
        state = simulation.context.getState(getForces=True,
                                            getEnergy=True,
                                            groups=self._force_group)
        system = simulation.context.getSystem()
        force = system.getForce(self._force_idx)

        out_list = [state.getPotentialEnergy().format('%.2f'),
                    state.getKineticEnergy().format('%.2f'),
                    str(force.getCollectiveVariableValues(simulation.context)),
                    str(state.getForces().value_in_unit(unit.kilojoules / unit.mole / unit.nanometer)[0])
                    ]
        self._out.write(self.get_formated_str(out_list))
        self._out.flush()


class CustomEnergyReporter(CustomReporter):
    """
    From <http://docs.openmm.org/latest/userguide/application/
    04_advanced_sim_examples.html#extracting-and-reporting-forces-and-other-data>
    """

    def __init__(self, file, reportInterval, force_group):
        self._out = open(file, 'w')
        self._reportInterval = reportInterval
        self._force_group = force_group
        self.header_list = [
            "Potential Energy (kJ/mol)",
            "Kinetic Energy (kJ/mol)",
        ]
        self.write_header()

        ## Since I don't think I will every use bitmaps for this, enforce this to be a set
        if type(self._force_group) == int:
            self._force_group = {self._force_group}
        elif type(self._force_group) == list:
            self._force_group = set(self._force_group)

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, None)

    def report(self, simulation, state):
        state = simulation.context.getState(getForces=True,
                                            getEnergy=True,
                                            groups=self._force_group)
        out_list = [
            state.getPotentialEnergy().format('%.2f'),
            state.getKineticEnergy().format('%.2f'),
        ]
        self._out.write(self.get_formated_str(out_list))
        self._out.flush()

class CustomForceReporter(CustomReporter):
    """
    From <http://docs.openmm.org/latest/userguide/application/
    04_advanced_sim_examples.html#extracting-and-reporting-forces-and-other-data>
    """

    def __init__(self, file, reportInterval, force_group, force_idx):
        self._out = open(file, 'w')
        self._reportInterval = reportInterval
        self._force_group = force_group
        self._force_idx = force_idx
        self.write_header()

        ## Since I don't think I will every use bitmaps for this, enforce this to be a set
        if type(self._force_group) == int:
            self._force_group = {self._force_group}
        elif type(self._force_group) == list:
            self._force_group = set(self._force_group)

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, None)

    def report(self, simulation, state):
        state = simulation.context.getState(getForces=True,
                                            getEnergy=True,
                                            groups=self._force_group)
        self._out.write(state.getForces())
        self._out.flush()

    def get_formated_str(self, out_list):
        formated_list = [f"{item:30}" for item in out_list]
        out_str = "\t".join(formated_list) + "\n"
        return out_str

def save_free_energies(output_dir, meta):
    print("Writing final free energies")
    with open(os.path.join(output_dir, "free_energies.npy"), "wb") as f:
        numpy.save(f, meta.getFreeEnergy())
