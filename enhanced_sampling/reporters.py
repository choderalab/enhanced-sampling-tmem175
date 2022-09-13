import openmm
import os
import numpy


class MetadynamicsReporter():
    def __init__(self, free_energy_file, collective_variable_file, reportInterval, meta):
        self._free_energy_file = open(free_energy_file, 'w')
        self._collective_variable_file = open(collective_variable_file, 'w')
        self._reportInterval = reportInterval
        self._meta = meta

    def __del__(self):
        self._free_energy_file.close()
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
        # free_energy = self._meta.getFreeEnergy()
        # self._free_energy_file.write(free_energy)

        collective_variables = self._meta.getCollectiveVariables(simulation)
        cv_list = [str(cv) for cv in collective_variables]
        cv_str = ", ".join(cv_list)
        # numpy.savez(self._collective_variable_file, collective_variables)
        self._collective_variable_file.write(f"{cv_str}\n")

def save_free_energies(output_dir, meta):
    print("Writing final free energies")
    with open(os.path.join(output_dir, "free_energies.npy"), "wb") as f:
        numpy.save(f, meta.getFreeEnergy())
