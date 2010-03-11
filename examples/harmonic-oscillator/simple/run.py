#! /bin/env python
#@+leo-ver=4-thin
#@+node:gcross.20100311125034.2363:@thin run.py
#@@first
#@<< Imports >>
#@+node:gcross.20100311125034.2364:<< Imports >>
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
#@-node:gcross.20100311125034.2364:<< Imports >>
#@nl
#@<< Parse command line >>
#@+node:gcross.20100311125034.2365:<< Parse command line >>
try:
    output_root_directory = sys.argv[1]
except IndexError:
    print "You must supply as an argument the directory in which to place the results."
    print "Example:  " + sys.argv[0] + " results"
    sys.exit()
#@-node:gcross.20100311125034.2365:<< Parse command line >>
#@nl

#@+others
#@-others

#@<< System Configuration >>
#@+node:gcross.20100311125034.2366:<< System Configuration >>
configuration = {
    # System parameters
    "number_of_dimensions": 3,
    "number_of_particles": 10,
    "number_of_slices": 54,
    "lambda_": 0.5, # --> hbar/2m
    "initial_particle_distribution_size": 1,
    # Physics parameters
    "harmonic_oscillator_frequencies": [1,1,1],
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
#@-node:gcross.20100311125034.2366:<< System Configuration >>
#@nl
#@<< Histogram Configuration >>
#@+node:gcross.20100311125034.2367:<< Histogram Configuration >>
_1d_density_histogram_left  = -2
_1d_density_histogram_right = +2
_1d_density_histogram_bin_count = 51

radial_densities_histogram_maximum_radius = 2.5
radial_densities_histogram_bin_count = 51

angular_densities_histogram_bin_count = 51

particle_separation_histogram_maximum_length = 2
particle_separation_histogram_bin_count = 51
#@nonl
#@-node:gcross.20100311125034.2367:<< Histogram Configuration >>
#@nl
#@<< Run simulation >>
#@+node:gcross.20100311125034.2369:<< Run simulation >>
system = System(**configuration)
#@<< Initialize physics >>
#@+node:gcross.20100311125034.2371:<< Initialize physics >>
for physics in [
    HarmonicOscillator,
    FourthOrderGreensFunction,
    ]: system.add_physics(physics)
#@-node:gcross.20100311125034.2371:<< Initialize physics >>
#@nl
#@<< Initialize observables >>
#@+node:gcross.20100311125034.2372:<< Initialize observables >>
for slice_name, slice_number in [("leftmost-path-slice",0),("center-path-slice",system.center_slice_number),("rightmost-slice",system.number_of_slices-1)]:
    density_slice_subdirectory = "{output_root_directory}/{slice_name}".format(**vars())
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
averages_file = "{output_root_directory}/average-value-observables".format(**vars())
for observable in [
    AverageRadiusEstimate(
        center_slice,
        averages_file,
        "distance to origin:"
    ),
    AverageParticleSeparationEstimate(
        center_slice,
        averages_file,
        "particle separation:"
    ),
    TotalEnergyEstimate(
        averages_file,
        "total energy:"
    ),
    ]: system.add_observable(observable)
#@-node:gcross.20100311125034.2372:<< Initialize observables >>
#@nl
system.run()
system.total_and_write_observables()
del system.observables
del system
gc.collect()
#@-node:gcross.20100311125034.2369:<< Run simulation >>
#@nl
#@-node:gcross.20100311125034.2363:@thin run.py
#@-leo
