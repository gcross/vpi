#@+leo-ver=4-thin
#@+node:gcross.20100914200151.1700:@thin setup.py
#@@language python
#@@tabwidth -4
#@+others
#@+node:gcross.20100914200151.1701:setup declarations
sources = [
    "src/utils/constants.f95",
    "src/utils/erf.f",
    "src/utils/erfn.f95",
    "src/utils/timers.f95",
    "src/utils/rand_utils.f95",
    "src/utils/xij.f95",
    "src/utils/numeric_differentiation.f95",
    "src/measurement/observables.f95",
    "src/measurement/histograms.f95",
    "src/physics/gfn.f95",
    "src/physics/angular_momentum.f95",
    "src/physics/hard_sphere_interaction.f95",
    "src/physics/harmonic_oscillator.f95",
    "src/physics/lennard_jones_interaction.f95",
    "src/path-integral/lattice.f95",
    "src/path-integral/sample.f95",
    "src/path-integral/thermalize.f95",
    ]

#@-node:gcross.20100914200151.1701:setup declarations
#@+node:gcross.20100914200151.1702:configuration
def configuration():
    from numpy.distutils.misc_util import Configuration
    config = Configuration()
    config.add_extension(
        "vpi.fortran",
        sources=sources,
        libraries = ["blas"],
        )
    return config

#@-node:gcross.20100914200151.1702:configuration
#@-others
if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(
        name="vpi",
        version="1.0",
        description="Variational Path Integral Monte Carlo Simulator",
        author="Jonathan de Bois and Gregory Crosswhite",
        author_email="Gregory Crosswhite <gcross@phys.washington.edu>",
        package_dir = {"":"src"},
        packages = ["vpi"],
        requires = ["mpi4py"],
        **configuration().todict()
    )
#@-node:gcross.20100914200151.1700:@thin setup.py
#@-leo
