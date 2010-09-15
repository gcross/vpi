#! /bin/env python
#@+leo-ver=4-thin
#@+node:gcross.20100311125034.2384:@thin run.py
#@@first
#@<< Imports >>
#@+node:gcross.20100311125034.2385:<< Imports >>
import gc

import sys
sys.path.append("lib")

from vpi import *
import vpi.fortran as vpif

import itertools

from numpy import log

from scipy.misc import derivative
import __builtin__
from itertools import imap, combinations
#@nonl
#@-node:gcross.20100311125034.2385:<< Imports >>
#@nl
#@<< Parse command line >>
#@+node:gcross.20100311125034.2386:<< Parse command line >>
try:
    output_root_directory = sys.argv[1]
except IndexError:
    print "You must supply as an argument the directory in which to place the results."
    print "Example:  " + sys.argv[0] + " results"
    sys.exit()
#@-node:gcross.20100311125034.2386:<< Parse command line >>
#@nl

#@+others
#@-others

#@<< System Configuration >>
#@+node:gcross.20100311125034.2387:<< System Configuration >>
configuration = {
    # System parameters
    "number_of_dimensions": 3,
    "number_of_particles": 10,
    "number_of_slices": 102,
    "lambda_": 0.5, # --> hbar/2m
    "initial_particle_distribution_size": 4,
    # Physics parameters
    "harmonic_oscillator_frequencies": [1,1,1],
    # Run parameters
    "total_number_of_observations": 10000,
    "number_of_prethermalization_steps": 1000,
    # Move parameters 
    # First move type is bridge, second is rigid, third is swap (unused).
    "dM": 20, # --> length of the path to move during a bridge move
    "move_type_probabilities": [0.9,0.1,0],
    "move_type_differentials": [0.05,0.3,0],
    "low_swap_dimension": 1,
    "high_swap_dimension": 3,
}
#@-node:gcross.20100311125034.2387:<< System Configuration >>
#@nl
#@<< Histogram Configuration >>
#@+node:gcross.20100311125034.2388:<< Histogram Configuration >>
_1d_density_histogram_left  = -2
_1d_density_histogram_right = +2
_1d_density_histogram_bin_count = 51

radial_densities_histogram_maximum_radius = 2.5
radial_densities_histogram_bin_count = 51

angular_densities_histogram_bin_count = 51

particle_separation_histogram_maximum_length = 5
particle_separation_histogram_bin_count = 51
#@nonl
#@-node:gcross.20100311125034.2388:<< Histogram Configuration >>
#@nl
#@<< System properties message >>
#@+node:gcross.20100311125034.2389:<< System properties message >>
system_properties_message = """\
Now simulating a system with a hard sphere radius of {hard_sphere_radius}"""
#@-node:gcross.20100311125034.2389:<< System properties message >>
#@nl

for hard_sphere_radius in [0,0.1,0.2,0.4,0.8]:
    #@    << Run simulation for given parameters >>
    #@+node:gcross.20100311125034.2390:<< Run simulation for given parameters >>
    #@<< Add current parameters to the configuration >>
    #@+node:gcross.20100311125034.2391:<< Add current parameters to the configuration >>
    configuration["hard_sphere_radius"] = hard_sphere_radius
    #@-node:gcross.20100311125034.2391:<< Add current parameters to the configuration >>
    #@nl

    my_directory = "{output_root_directory}/radius={hard_sphere_radius}".format(**vars()) 
    if (my_rank == 0):
        print
        print system_properties_message.format(**vars())

    system = System(**configuration)
    #@<< Initialize physics >>
    #@+node:gcross.20100311125034.2392:<< Initialize physics >>
    for physics in [
        HarmonicOscillator,
        HardSphereInteraction,
        FourthOrderGreensFunction,
        ]: system.add_physics(physics)
    #@-node:gcross.20100311125034.2392:<< Initialize physics >>
    #@nl
    #@<< Initialize observables >>
    #@+node:gcross.20100311125034.2393:<< Initialize observables >>
    for slice_name, slice_number in [("leftmost-path-slice",0),("center-path-slice",system.center_slice_number),("rightmost-slice",system.number_of_slices-1)]:
        density_slice_subdirectory = "{my_directory}/{slice_name}".format(**vars())
        for observable in [
                PositionDensity1DHistogram(
                    slice_number,
                    [_1d_density_histogram_left]*configuration["number_of_dimensions"],
                    [_1d_density_histogram_right]*configuration["number_of_dimensions"],
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
    averages_directory = "{my_directory}/../averages-at-center-slice".format(**vars())
    for observable in [
        AverageRadiusEstimate(
            center_slice,
            averages_directory + "/distance-to-origin",
            hard_sphere_radius
        ),
        AverageParticleSeparationEstimate(
            center_slice,
            averages_directory + "/particle-separation",
            hard_sphere_radius
        ),
        TotalEnergyEstimate(
            averages_directory + "/total-energy",
            hard_sphere_radius
        ),
        ]: system.add_observable(observable)
    #@-node:gcross.20100311125034.2393:<< Initialize observables >>
    #@nl
    system.run()
    system.total_and_write_observables()
    del system.observables
    del system
    gc.collect()
    #@-node:gcross.20100311125034.2390:<< Run simulation for given parameters >>
    #@nl
#@-node:gcross.20100311125034.2384:@thin run.py
#@-leo
