#! /bin/env python
#@+leo-ver=4-thin
#@+node:gcross.20100311125034.2323:@thin run.py
#@@first
#@<< Imports >>
#@+node:gcross.20100311125034.2324:<< Imports >>
import gc

import sys
sys.path.append("lib")

from vpi import *
import vpif

import itertools

from numpy import log

from scipy.misc import derivative
import __builtin__
from itertools import imap, combinations
#@-node:gcross.20100311125034.2324:<< Imports >>
#@nl
#@<< Parse command line >>
#@+node:gcross.20100311125034.2349:<< Parse command line >>
try:
    output_root_directory = sys.argv[1]
except IndexError:
    print "You must supply as an argument the directory in which to place the results."
    print "Example:  " + sys.argv[0] + " results"
    sys.exit()
#@-node:gcross.20100311125034.2349:<< Parse command line >>
#@nl

#@+others
#@-others

#@<< System Configuration >>
#@+node:gcross.20100311125034.2335:<< System Configuration >>
configuration = {
    # System parameters
    "number_of_particles": 10,
    "number_of_slices": 54,
    "lambda_": 0.5, # --> hbar/2m
    "initial_particle_distribution_size": 1,
    # Run parameters
    "total_number_of_observations": 10000,
    "number_of_prethermalization_steps": 1000,
    # Move parameters 
    # First move type is bridge, second is rigid, third is swap (unused).
    "dM": 20, # --> length of the path to move during a bridge move
    "move_type_probabilities": [0.9,0.1,0],
    "move_type_differentials": [0.1,1,0],
    "low_swap_dimension": 1,
    "high_swap_dimension": 3,
}
#@-node:gcross.20100311125034.2335:<< System Configuration >>
#@nl
#@<< Histogram Configuration >>
#@+node:gcross.20100311125034.2336:<< Histogram Configuration >>
_1d_density_histogram_left  = -2
_1d_density_histogram_right = +2
_1d_density_histogram_bin_count = 51

radial_densities_histogram_maximum_radius = 2.5
radial_densities_histogram_bin_count = 51

angular_densities_histogram_bin_count = 51

particle_separation_histogram_maximum_length = 2
particle_separation_histogram_bin_count = 51
#@nonl
#@-node:gcross.20100311125034.2336:<< Histogram Configuration >>
#@nl
#@<< System properties message >>
#@+node:gcross.20100311125034.2337:<< System properties message >>
system_properties_message = """\
Now simulating a system with
    *) {number_of_dimensions} dimension(s);
    *) and inside a harmonic oscillator trap with frequency {harmonic_oscillator_trap_frequency} in all dimensions"""
#@-node:gcross.20100311125034.2337:<< System properties message >>
#@nl

for number_of_dimensions in [1,2,3]:
    for harmonic_oscillator_trap_frequency in [0.5,1,2]:
        #@        << Run simulation for given parameters >>
        #@+node:gcross.20100311125034.2338:<< Run simulation for given parameters >>
        #@<< Add current parameters to the configuration >>
        #@+node:gcross.20100311125034.2343:<< Add current parameters to the configuration >>
        configuration["harmonic_oscillator_frequencies"] = array([harmonic_oscillator_trap_frequency]*number_of_dimensions)
        configuration["number_of_dimensions"] = number_of_dimensions
        #@-node:gcross.20100311125034.2343:<< Add current parameters to the configuration >>
        #@nl

        my_directory = "{output_root_directory}/{number_of_dimensions}D/trap-frequency={harmonic_oscillator_trap_frequency}".format(**vars()) 
        if (my_rank == 0):
            print
            print system_properties_message.format(**vars())

        system = System(**configuration)
        #@<< Initialize physics >>
        #@+node:gcross.20100311125034.2339:<< Initialize physics >>
        for physics in [
            HarmonicOscillator,
            FourthOrderGreensFunction,
            ]: system.add_physics(physics)
        #@-node:gcross.20100311125034.2339:<< Initialize physics >>
        #@nl
        #@<< Initialize observables >>
        #@+node:gcross.20100311125034.2340:<< Initialize observables >>
        for slice_name, slice_number in [("leftmost-path-slice",0),("center-path-slice",system.center_slice_number),("rightmost-slice",system.number_of_slices-1)]:
            density_slice_subdirectory = "{my_directory}/{slice_name}".format(**vars())
            for observable in [
                    PositionDensity1DHistogram(
                        slice_number,
                        [_1d_density_histogram_left]*number_of_dimensions,
                        [_1d_density_histogram_right]*number_of_dimensions,
                        _1d_density_histogram_bin_count,
                        [density_slice_subdirectory + "/1d-densities/" + label for label in ["x","y","z","w"][:system.number_of_dimensions]]
                    ),
                    RadialDensityHistogram(
                        slice_number,
                        radial_densities_histogram_maximum_radius,
                        radial_densities_histogram_bin_count,
                        density_slice_subdirectory + "/radial-density"
                    ),
                    ParticleSeparationHistogram(
                        slice_number,
                        particle_separation_histogram_maximum_length,
                        particle_separation_histogram_bin_count,
                        density_slice_subdirectory + "/particle-separation"
                    ),
                ]: system.add_observable(observable)

        center_slice = system.number_of_slices // 2
        averages_directory = "{my_directory}/../average-values-at-center-slice".format(**vars())
        for observable in [
            AverageRadiusEstimate(
                center_slice,
                averages_directory + "/distance-to-origin",
                harmonic_oscillator_trap_frequency
            ),
            AverageParticleSeparationEstimate(
                center_slice,
                averages_directory + "/particle-separation",
                harmonic_oscillator_trap_frequency
            ),
            TotalEnergyEstimate(
                averages_directory + "/total-energy",
                harmonic_oscillator_trap_frequency
            ),
            ]: system.add_observable(observable)
        #@-node:gcross.20100311125034.2340:<< Initialize observables >>
        #@nl
        system.run()
        system.total_and_write_observables()
        del system.observables
        del system
        gc.collect()
        #@-node:gcross.20100311125034.2338:<< Run simulation for given parameters >>
        #@nl
#@-node:gcross.20100311125034.2323:@thin run.py
#@-leo
