import openmm
import os
import numpy


class MetadynamicsReporter():
    def __init__(self, collective_variable_file, reportInterval, meta):
        self._collective_variable_file = open(collective_variable_file, 'w')
        self._reportInterval = reportInterval
        self._meta = meta

    def __del__(self):
        self._collective_variable_file.close()

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

        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False, None)

    def report(self, simulation, state):
        collective_variables = self._meta.getCollectiveVariables(simulation)
        cv_list = [str(cv) for cv in collective_variables]
        cv_str = ", ".join(cv_list)
        self._collective_variable_file.write(f"{cv_str}\n")

class ForceReporter(object):
    """
    From <http://docs.openmm.org/latest/userguide/application/
    04_advanced_sim_examples.html#extracting-and-reporting-forces-and-other-data>
    """
    def __init__(self, file, reportInterval):
        self._out = open(file, 'w')
        self._reportInterval = reportInterval

    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, False, True, False, None)

    def report(self, simulation, state):
        forces = state.getForces().value_in_unit(kilojoules/mole/nanometer)
        for f in forces:
            self._out.write('%g %g %g\n' % (f[0], f[1], f[2]))

def save_free_energies(output_dir, meta):
    print("Writing final free energies")
    with open(os.path.join(output_dir, "free_energies.npy"), "wb") as f:
        numpy.save(f, meta.getFreeEnergy())
