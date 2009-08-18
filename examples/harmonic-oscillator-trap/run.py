#@+leo-ver=4-thin
#@+node:gcross.20090818081913.1341:@thin run.py
import sys
sys.path.append("../../lib")

#@<< Imports >>
#@+node:gcross.20090818081913.1360:<< Imports >>
import vpi

from numpy import *
from numpy.random import rand
#@-node:gcross.20090818081913.1360:<< Imports >>
#@nl

#@<< Configuration >>
#@+node:gcross.20090818081913.1344:<< Configuration >>
# System parameters
n_particles = 2
n_slices = 54
lambda_ = 0.5

# Observable parameters
number_of_observations = 10000

# Move parameters
dM = 22
move_type_probabilities = [0.9,0.1,0]
move_type_differentials = [0.1,0.6,0]
low_swap_dimension = 1
high_swap_dimension= 3

# Potential parameters
harmonic_oscillator_coefficients = (1,1,1)

# Miscellaneous
#@-node:gcross.20090818081913.1344:<< Configuration >>
#@nl

#@+others
#@+node:gcross.20090818081913.1343:Functions
#@+node:gcross.20090818081913.1342:Potential function
def compute_potential(x,xij2):
    x_sq = x**2
    U = dot(x_sq,harmonic_oscillator_coefficients)/2.0
    gradU2 = sum(U,axis=1)
    return U, gradU2, False
#@-node:gcross.20090818081913.1342:Potential function
#@+node:gcross.20090818081913.1345:Trial function
def trial_function(x,xij2):
    return -sum(dot(x**2,harmonic_oscillator_coefficients))/2
#@-node:gcross.20090818081913.1345:Trial function
#@+node:gcross.20090818081913.1353:Trial derivatives
def trial_derivatives(x,xij2):
    gradient_of_log_trial_fn = -x*harmonic_oscillator_coefficients
    laplacian_of_log_trial_fn = -sum(harmonic_oscillator_coefficients)*x.shape[0]
    return gradient_of_log_trial_fn, laplacian_of_log_trial_fn
#@-node:gcross.20090818081913.1353:Trial derivatives
#@+node:gcross.20090818081913.1359:Observables
def compute_energy_estimates():
    estimates = []
    for i in [0,-1]:
        gradient_of_log_trial_fn, laplacian_of_log_trial_fn = trial_derivatives(x[i],xij2[i])
        estimates.append(
            vpi.observables.compute_local_energy_estimate(
                U[0],
                gradient_of_log_trial_fn, laplacian_of_log_trial_fn,
                lambda_,
            )
        )
    return estimates
#@-node:gcross.20090818081913.1359:Observables
#@-node:gcross.20090818081913.1343:Functions
#@-others

#@<< Initialization >>
#@+node:gcross.20090818081913.1349:<< Initialization >>
n_dimensions = 3

assert (n_slices % 2 == 0 and n_slices % 4 == 2)
c_slice = n_slices / 2

x = zeros((n_slices,n_particles,n_dimensions),dtype=double,order='Fortran')
x[...] = rand(1,n_particles,n_dimensions)

xij2 = zeros((n_slices,n_particles,n_particles),dtype=double,order='Fortran')
vpi.xij.update_xij(xij2,x)

U = zeros((n_slices,n_particles),dtype=double,order='Fortran')
gradU2 = zeros((n_slices),dtype=double,order='Fortran')

slice_move_attempted_counts = zeros((n_slices,),'i')
slice_move_accepted_counts = zeros((n_slices,),'i')

move_type_attempted_counts = zeros((3,),'i')
move_type_accepted_counts = zeros((3,),'i')

U_weights, gU2_weights = vpi.gfn.initialize_4th_order_weights(n_slices)

n_trials = n_particles*n_slices/dM
use_4th_order_green_function = True
#@-node:gcross.20090818081913.1349:<< Initialization >>
#@nl

#@<< Main iteration >>
#@+node:gcross.20090818081913.1346:<< Main iteration >>
vpi.thermalize.thermalize_path(
    x,xij2,
    U,gradU2,
    n_trials*1000,
    move_type_probabilities,move_type_differentials,
    dM,
    lambda_,
    low_swap_dimension,high_swap_dimension,
    slice_move_attempted_counts,move_type_attempted_counts,
    slice_move_accepted_counts,move_type_accepted_counts,
    compute_potential,trial_function,
    U_weights,gU2_weights,
    use_4th_order_green_function,
)

energy_estimates = zeros((2,),dtype=double)

for _ in xrange(number_of_observations):
    vpi.thermalize.thermalize_path(
        x,xij2,
        U,gradU2,
        n_trials,
        move_type_probabilities,move_type_differentials,
        dM,
        lambda_,
        low_swap_dimension,high_swap_dimension,
        slice_move_attempted_counts,move_type_attempted_counts,
        slice_move_accepted_counts,move_type_accepted_counts,
        compute_potential,trial_function,
        U_weights,gU2_weights,
        use_4th_order_green_function,
    )
    energy_estimates += compute_energy_estimates()

from mpi4py import MPI

comm = MPI.COMM_WORLD

total_energy_estimates = zeros((2,),dtype=double)

comm.Reduce((energy_estimates,MPI.DOUBLE),(total_energy_estimates,MPI.DOUBLE))

if comm.Get_rank() == 0:
    print "%.0f%%" % (float(move_type_accepted_counts[0])/float(move_type_attempted_counts[0])*100)
    print total_energy_estimates/(comm.Get_size()*number_of_observations)
#@-node:gcross.20090818081913.1346:<< Main iteration >>
#@nl
#@nonl
#@-node:gcross.20090818081913.1341:@thin run.py
#@-leo
