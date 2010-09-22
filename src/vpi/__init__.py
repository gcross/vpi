#@+leo-ver=4-thin
#@+node:gcross.20090902085220.2158:@thin vpi/__init__.py
#@@language python
#@@tabwidth -4

#@<< Imports >>
#@+node:gcross.20090902085220.2159:<< Imports >>
import vpi.fortran as vpif

from abc import ABCMeta, abstractmethod

from numpy import *
from numpy.random import rand

import os
import os.path

import sys

import itertools
from itertools import izip

import __builtin__
#@-node:gcross.20090902085220.2159:<< Imports >>
#@nl

#@<< MPI Initialization >>
#@+node:gcross.20090902085220.2160:<< MPI Initialization >>
from mpi4py import MPI

comm = MPI.COMM_WORLD

number_of_processors = comm.Get_size()
my_rank = comm.Get_rank()
#@-node:gcross.20090902085220.2160:<< MPI Initialization >>
#@nl

#@+others
#@+node:gcross.20090902085220.2286:Functions
#@+node:gcross.20090902085220.2287:ensure_path_to_file_exists
def ensure_path_to_file_exists(path):
    directory, _ = os.path.split(path)
    if not os.path.exists(directory):
        os.makedirs(directory)
#@-node:gcross.20090902085220.2287:ensure_path_to_file_exists
#@-node:gcross.20090902085220.2286:Functions
#@+node:gcross.20090902085220.2385:Classes
#@+node:gcross.20090902085220.2386:Exception classes
#@+node:gcross.20090902085220.2387:MoveRejected
class MoveRejected(Exception):
        pass
#@-node:gcross.20090902085220.2387:MoveRejected
#@+node:gcross.20090908085435.1634:FailureToContructValidSystemException
class FailureToContructValidSystemException(Exception):
    pass
#@-node:gcross.20090908085435.1634:FailureToContructValidSystemException
#@-node:gcross.20090902085220.2386:Exception classes
#@+node:gcross.20090902085220.2161:Observable classes
#@+others
#@+node:gcross.20090902085220.2162:Base classes
#@+node:gcross.20090902085220.2163:class Observable
class Observable(object):
    #@    @+others
    #@+node:gcross.20090902085220.2164:total_and_write
    def total_and_write(self):
        totals = self.compute_total()
        if my_rank == 0:
            self.write_out_totals(totals)
    #@-node:gcross.20090902085220.2164:total_and_write
    #@-others
#@-node:gcross.20090902085220.2163:class Observable
#@+node:gcross.20090902085220.2165:class AverageValuesEstimate
class AverageValuesEstimate(Observable):
    #@    @+others
    #@+node:gcross.20090902085220.2166:(fields)
    total_number_of_observations = 0
    #@-node:gcross.20090902085220.2166:(fields)
    #@+node:gcross.20090902085220.2167:__init__
    def __init__(self,*shape):
        self.estimate = zeros(shape,dtype=double)
        self.estimate_squared = zeros(shape,dtype=double)
    #@-node:gcross.20090902085220.2167:__init__
    #@+node:gcross.20090902085220.2168:write_out_totals
    def write_out_totals(self,totals_and_errors):
        means, means_of_squares = totals_and_errors
        variance = means_of_squares - means*means
        uncertainty_in_means = sqrt(variance / (self.total_number_of_observations * number_of_processors))
        self.write_out(means,uncertainty_in_means)
    #@-node:gcross.20090902085220.2168:write_out_totals
    #@+node:gcross.20090902085220.2169:compute_total
    def compute_total(self):
        total_estimate_and_square = zeros(2*prod(self.estimate.shape),dtype='d',order='Fortran')
        comm.Reduce((array([self.estimate,self.estimate_squared]).ravel(),MPI.DOUBLE),(total_estimate_and_square,MPI.DOUBLE))
        total_estimate_and_square /= (self.total_number_of_observations * number_of_processors)
        return total_estimate_and_square.reshape((2,)+self.estimate.shape)
    #@-node:gcross.20090902085220.2169:compute_total
    #@+node:gcross.20090902085220.2170:add
    def add(self,estimate):
        self.total_number_of_observations += 1
        self.estimate += estimate
        self.estimate_squared += estimate**2
    #@-node:gcross.20090902085220.2170:add
    #@+node:gcross.20090911091023.2447:reset
    def reset(self):
        self.estimate *= 0
        self.estimate_squared *= 0
    #@-node:gcross.20090911091023.2447:reset
    #@-others
#@-node:gcross.20090902085220.2165:class AverageValuesEstimate
#@+node:gcross.20090902085220.2171:class SingleAverageValueEstimate
class SingleAverageValueEstimate(AverageValuesEstimate):
    #@    @+others
    #@+node:gcross.20090902085220.2172:__init__
    def __init__(self):
        AverageValuesEstimate.__init__(self,1)
    #@-node:gcross.20090902085220.2172:__init__
    #@-others
#@-node:gcross.20090902085220.2171:class SingleAverageValueEstimate
#@+node:gcross.20090902085220.2173:class EstimatesAppendedToFile
class EstimatesAppendedToFile(Observable):
    #@    @+others
    #@+node:gcross.20090902085220.2174:__init__
    def __init__(self,filename,label):
        self.filename = filename
        self.label = label
    #@-node:gcross.20090902085220.2174:__init__
    #@+node:gcross.20090902085220.2175:write_out
    def write_out(self,total_estimate,total_estimate_error):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"a") as f:
            print >> f, self.label,
            for value, error in izip(total_estimate,total_estimate_error):
                print >> f, value, error,
            print >> f
    #@-node:gcross.20090902085220.2175:write_out
    #@-others
#@-node:gcross.20090902085220.2173:class EstimatesAppendedToFile
#@+node:gcross.20090902085220.2176:class SingleAverageValueEstimateAppendedToFile
class SingleAverageValueEstimateAppendedToFile(SingleAverageValueEstimate,EstimatesAppendedToFile):
    #@    @+others
    #@+node:gcross.20090902085220.2177:__init__
    def __init__(self,filename,label):
        SingleAverageValueEstimate.__init__(self)
        EstimatesAppendedToFile.__init__(self,filename,label)
    #@-node:gcross.20090902085220.2177:__init__
    #@-others
#@-node:gcross.20090902085220.2176:class SingleAverageValueEstimateAppendedToFile
#@+node:gcross.20090902085220.2178:class SingleAverageValueAtSliceEstimateAppendedToFile
class SingleAverageValueAtSliceEstimateAppendedToFile(SingleAverageValueEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090902085220.2179:__init__
    def __init__(self,slice_number,label,filename):
        SingleAverageValueEstimateAppendedToFile.__init__(self,label,filename)
        self.slice_number = slice_number
    #@-node:gcross.20090902085220.2179:__init__
    #@-others
#@-node:gcross.20090902085220.2178:class SingleAverageValueAtSliceEstimateAppendedToFile
#@-node:gcross.20090902085220.2162:Base classes
#@+node:gcross.20090902085220.2180:Histograms
#@+node:gcross.20090902085220.2181:class Histogram
class Histogram(Observable):
    #@    @+others
    #@+node:gcross.20090902085220.2182:compute_total
    def compute_total(self):
        total_histogram = zeros(self.histogram.shape,dtype='i',order='Fortran')
        comm.Reduce((self.histogram,MPI.INT),(total_histogram,MPI.INT))
        return total_histogram
    #@-node:gcross.20090902085220.2182:compute_total
    #@+node:gcross.20090902085220.2183:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        total_counts = float(sum(histogram))
        bin_width = float(self.right-self.left)/self.number_of_bins
        current = float(self.left)+bin_width/2
        with open(self.filename,"w") as f:
            for count in histogram:
                print >> f, "{0} {1}".format(current,count/total_counts)
                current += bin_width
    #@-node:gcross.20090902085220.2183:write_out_totals
    #@+node:gcross.20090911091023.2449:reset
    def reset(self):
        self.histogram *= 0
    #@-node:gcross.20090911091023.2449:reset
    #@-others
#@-node:gcross.20090902085220.2181:class Histogram
#@+node:gcross.20090902085220.2184:class PositionDensity1DHistogram
class PositionDensity1DHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090902085220.2185:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filenames):
        assert len(left) == len(right)
        self.left = array(left,dtype=double,order='Fortran')
        self.right = array(right,dtype=double,order='Fortran')
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,len(left)),dtype='i',order='Fortran')
        self.filenames = filenames
    #@-node:gcross.20090902085220.2185:__init__
    #@+node:gcross.20090902085220.2186:update
    def update(self):
        vpif.histograms.accumulate_1d_densities(
            self.system.x[...,self.slice_number],
            self.left,self.right,
            self.histogram
        )
    #@nonl
    #@-node:gcross.20090902085220.2186:update
    #@+node:gcross.20090902085220.2187:write_out_totals
    def write_out_totals(self,histograms):
        for filename, histogram, left, right in izip(self.filenames,histograms.transpose(),self.left,self.right):
            ensure_path_to_file_exists(filename)
            with open(filename,"w") as f:
                total_counts = float(sum(histogram))
                bin_width = float(right-left)/self.number_of_bins
                current = float(left)+bin_width/2
                for count in histogram:
                    print >> f, "{0} {1}".format(current,count/total_counts)
                    current += bin_width
    #@-node:gcross.20090902085220.2187:write_out_totals
    #@-others
#@-node:gcross.20090902085220.2184:class PositionDensity1DHistogram
#@+node:gcross.20090902085220.2188:class RadialDensityHistogram
class RadialDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090902085220.2189:(fields)
    left = 0
    #@-node:gcross.20090902085220.2189:(fields)
    #@+node:gcross.20090902085220.2190:__init__
    def __init__(self,slice_number,maximum_radius,number_of_bins,filename):
        self.right = self.maximum_radius = maximum_radius
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090902085220.2190:__init__
    #@+node:gcross.20090902085220.2191:update
    def update(self):
        vpif.histograms.accumulate_radial_densities(
            self.system.x[...,self.slice_number],
            self.maximum_radius,
            self.histogram
        )
    #@-node:gcross.20090902085220.2191:update
    #@+node:gcross.20090903090230.2077:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        total_counts = float(sum(histogram))
        bin_width = float(self.right-self.left)/self.number_of_bins
        current = float(self.left)+bin_width/2
        normalization_exponent = self.system.number_of_dimensions
        with open(self.filename,"w") as f:
            for count in histogram:
                normalization = (current+bin_width/2)**normalization_exponent \
                              - (current-bin_width/2)**normalization_exponent
                print >> f, "{0} {1}".format(current,count/(total_counts*normalization))
                current += bin_width
    #@-node:gcross.20090903090230.2077:write_out_totals
    #@-others
#@-node:gcross.20090902085220.2188:class RadialDensityHistogram
#@+node:gcross.20090902085220.2192:class PlaneRadialDensityHistogram
class PlaneRadialDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090902085220.2193:(fields)
    left = 0
    #@-node:gcross.20090902085220.2193:(fields)
    #@+node:gcross.20090902085220.2194:__init__
    def __init__(self,slice_number,maximum_radius,number_of_bins,filename):
        self.right = self.maximum_radius = maximum_radius
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090902085220.2194:__init__
    #@+node:gcross.20090902085220.2195:update
    def update(self):
        vpif.histograms.accumulate_plane_radial_densities(
            self.system.x[...,self.slice_number],
            self.maximum_radius,
            self.system.rotation_plane_axis_1,
            self.system.rotation_plane_axis_2,
            self.histogram
        )
    #@nonl
    #@-node:gcross.20090902085220.2195:update
    #@+node:gcross.20090903090230.2079:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        total_counts = float(sum(histogram))
        bin_width = float(self.right-self.left)/self.number_of_bins
        current = float(self.left)+bin_width/2
        normalization_exponent = 2
        with open(self.filename,"w") as f:
            for count in histogram:
                normalization = (current+bin_width/2)**normalization_exponent \
                              - (current-bin_width/2)**normalization_exponent
                print >> f, "{0} {1}".format(current,count/(total_counts*normalization))
                current += bin_width
    #@-node:gcross.20090903090230.2079:write_out_totals
    #@-others
#@-node:gcross.20090902085220.2192:class PlaneRadialDensityHistogram
#@+node:gcross.20090902085220.2196:class RecipricalRadiusSquaredDensityHistogram
class RecipricalRadiusSquaredDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090902085220.2197:(fields)
    left = 0
    #@-node:gcross.20090902085220.2197:(fields)
    #@+node:gcross.20090902085220.2198:__init__
    def __init__(self,slice_number,maximum_value,number_of_bins,filename):
        self.right = self.maximum_value = maximum_value
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090902085220.2198:__init__
    #@+node:gcross.20090902085220.2199:update
    def update(self):
        vpif.histograms.accumulate_reciprical_radius_squared_densities(
            self.system.x[...,self.slice_number],
            self.maximum_value,
            self.histogram
        )
    #@nonl
    #@-node:gcross.20090902085220.2199:update
    #@-others
#@-node:gcross.20090902085220.2196:class RecipricalRadiusSquaredDensityHistogram
#@+node:gcross.20090902085220.2200:class RecipricalPlaneRadiusSquaredDensityHistogram
class RecipricalPlaneRadiusSquaredDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090902085220.2201:(fields)
    left = 0
    #@-node:gcross.20090902085220.2201:(fields)
    #@+node:gcross.20090902085220.2202:__init__
    def __init__(self,slice_number,maximum_value,number_of_bins,filename):
        self.right = self.maximum_value = maximum_value
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090902085220.2202:__init__
    #@+node:gcross.20090902085220.2203:update
    def update(self):
        vpif.histograms.accumulate_recip_plane_r_sq_densities(
            self.system.x[...,self.slice_number],
            self.maximum_value,
            self.system.rotation_plane_axis_1,
            self.system.rotation_plane_axis_2,
            self.histogram
        )
    #@nonl
    #@-node:gcross.20090902085220.2203:update
    #@+node:gcross.20090902085220.2204:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"w") as f:
            total_counts = float(sum(histogram))
            bin_width = float(self.maximum_value)/self.number_of_bins
            current = bin_width/2
            for count in self.histogram:
                print >> f, "{0} {1}".format(current,count/total_counts)
                current += bin_width
    #@-node:gcross.20090902085220.2204:write_out_totals
    #@-others
#@-node:gcross.20090902085220.2200:class RecipricalPlaneRadiusSquaredDensityHistogram
#@+node:gcross.20090902085220.2205:class AngularSeparationDensityHistogram
class AngularSeparationDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090902085220.2206:(fields)
    left = 0
    right = 2*pi
    #@-node:gcross.20090902085220.2206:(fields)
    #@+node:gcross.20090902085220.2207:__init__
    def __init__(self,slice_number,number_of_bins,filename):
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090902085220.2207:__init__
    #@+node:gcross.20090902085220.2208:update
    def update(self):
        system = self.system
        x = system.x[...,self.slice_number]
        angles = arctan2(x[system.rotation_plane_axis_2],x[system.rotation_plane_axis_1])
        vpif.histograms.accumulate_angular_separation_densities(
            angles,
            self.histogram
        )
    #@nonl
    #@-node:gcross.20090902085220.2208:update
    #@+node:gcross.20090902085220.2209:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"w") as f:
            total_counts = float(sum(histogram))
            bin_width = float(2*pi)/self.number_of_bins
            current = bin_width/2
            for count in self.histogram:
                print >> f, "{0} {1}".format(current,count/total_counts)
                current += bin_width
    #@-node:gcross.20090902085220.2209:write_out_totals
    #@-others
#@-node:gcross.20090902085220.2205:class AngularSeparationDensityHistogram
#@+node:gcross.20090902085220.2210:class NeighborAngularSeparationDensityHistogram
class NeighborAngularSeparationDensityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090902085220.2211:(fields)
    left = 0
    right = 2*pi
    #@-node:gcross.20090902085220.2211:(fields)
    #@+node:gcross.20090902085220.2212:__init__
    def __init__(self,slice_number,number_of_bins,filename):
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090902085220.2212:__init__
    #@+node:gcross.20090902085220.2213:update
    def update(self):
        system = self.system
        x = system.x[...,self.slice_number]
        angles = arctan2(x[system.rotation_plane_axis_2],x[system.rotation_plane_axis_1])
        vpif.histograms.accumulate_neighbor_angular_separation_densities(
            angles,
            self.histogram
        )
    #@nonl
    #@-node:gcross.20090902085220.2213:update
    #@-others
#@-node:gcross.20090902085220.2210:class NeighborAngularSeparationDensityHistogram
#@+node:gcross.20090902085220.2214:class AngularVelocityHistogram
class AngularVelocityHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090902085220.2215:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filename):
        self.left = left
        self.right = right
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090902085220.2215:__init__
    #@+node:gcross.20090902085220.2216:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpif.angular_momentum.compute_angular_derivatives(
            system.x[...,self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        vpif.histograms.accumulate(first_derivatives,self.left,self.right,self.histogram)
    #@nonl
    #@-node:gcross.20090902085220.2216:update
    #@-others
#@-node:gcross.20090902085220.2214:class AngularVelocityHistogram
#@+node:gcross.20090902085220.2217:class AngularVelocitySquaredHistogram
class AngularVelocitySquaredHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090902085220.2218:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filename):
        self.left = left
        self.right = right
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090902085220.2218:__init__
    #@+node:gcross.20090902085220.2219:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpif.angular_momentum.compute_angular_derivatives(
            system.x[...,self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        first_derivatives **= 2
        vpif.histograms.accumulate(first_derivatives,self.left,self.right,self.histogram)
    #@nonl
    #@-node:gcross.20090902085220.2219:update
    #@-others
#@-node:gcross.20090902085220.2217:class AngularVelocitySquaredHistogram
#@+node:gcross.20090902085220.2220:class RotationQuadraticTermHistogram
class RotationQuadraticTermHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090902085220.2221:__init__
    def __init__(self,slice_number,left,right,number_of_bins,filename):
        self.left = left
        self.right = right
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090902085220.2221:__init__
    #@+node:gcross.20090902085220.2222:update
    def update(self):
        system = self.system
        x = system.x[...,self.slice_number]
        first_derivatives, _ = vpif.angular_momentum.compute_angular_derivatives(
            x,
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        term = (first_derivatives ** 2) / (x[system.rotation_plane_axis_2-1]**2+x[system.rotation_plane_axis_1-1]**2)
        vpif.histograms.accumulate(term,self.left,self.right,self.histogram)
    #@nonl
    #@-node:gcross.20090902085220.2222:update
    #@-others
#@-node:gcross.20090902085220.2220:class RotationQuadraticTermHistogram
#@+node:gcross.20090902085220.2223:class AngularVelocityAndRadiusHistogram
class AngularVelocityAndRadiusHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090902085220.2224:__init__
    def __init__(self,slice_number,maximum_angular_velocity,maximum_radius,number_of_bins,filename):
        self.maximum_angular_velocity = maximum_angular_velocity
        self.maximum_radius = maximum_radius
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,)*2,dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090902085220.2224:__init__
    #@+node:gcross.20090902085220.2225:update
    def update(self):
        system = self.system
        x = system.x[self.slice_number]
        first_derivatives, _ = vpif.angular_momentum.compute_angular_derivatives(
            x,
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        radii = sqrt(x[system.rotation_plane_axis_2-1]**2+x[system.rotation_plane_axis_1-1]**2)

        i_values = floor(first_derivatives/self.maximum_angular_velocity*self.number_of_bins)
        j_values = floor(radii/self.maximum_radius*self.number_of_bins)
        for (i,j) in izip(i_values,j_values):
            if (i >= 0) and (i < self.number_of_bins) and (j >= 0) and (j < self.number_of_bins):
                self.histogram[i,j] += 1
    #@nonl
    #@-node:gcross.20090902085220.2225:update
    #@+node:gcross.20090902085220.2226:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        total_counts = float(sum(histogram))
        with open(self.filename,"w") as f:
            for i in xrange(self.number_of_bins):
                for j in xrange(self.number_of_bins):
                    angular_velocity = (i+0.5)/self.number_of_bins*self.maximum_angular_velocity
                    radius = (j+0.5)/self.number_of_bins*self.maximum_radius
                    print >> f, "{0} {1} {2}".format(angular_velocity,radius,histogram[i,j]/total_counts)
                print >> f, ""
    #@-node:gcross.20090902085220.2226:write_out_totals
    #@-others
#@-node:gcross.20090902085220.2223:class AngularVelocityAndRadiusHistogram
#@+node:gcross.20090902085220.2227:class ParticleSeparationHistogram
class ParticleSeparationHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20090902085220.2228:(fields)
    left = 0
    #@-node:gcross.20090902085220.2228:(fields)
    #@+node:gcross.20090902085220.2229:__init__
    def __init__(self,slice_number,maximum_radius,number_of_bins,filename):
        self.right = self.maximum_value = maximum_radius
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20090902085220.2229:__init__
    #@+node:gcross.20090902085220.2230:update
    def update(self):
        vpif.histograms.accumulate_particle_separation_densities(
            self.system.xij2[...,self.slice_number],
            self.maximum_value,
            self.histogram
        )
    #@nonl
    #@-node:gcross.20090902085220.2230:update
    #@-others
#@-node:gcross.20090902085220.2227:class ParticleSeparationHistogram
#@+node:gcross.20100226131523.1663:class PositionDensity2DHistogram
class PositionDensity2DHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20100226131523.1664:__init__
    def __init__(self,slice_number,left_x,left_y,right_x,right_y,number_of_bins,filename):
        self.left_x = left_x
        self.left_y = left_y
        self.right_x = right_x
        self.right_y = right_y
        self.number_of_bins = number_of_bins
        self.slice_number = slice_number
        self.histogram = zeros((number_of_bins,number_of_bins),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20100226131523.1664:__init__
    #@+node:gcross.20100226131523.1665:update
    def update(self):
        vpif.histograms.accumulate_2d_density(
            self.system.x[...,self.slice_number],
            self.left_x,self.left_y,
            self.right_x,self.right_y,
            self.histogram
        )
    #@-node:gcross.20100226131523.1665:update
    #@+node:gcross.20100226131523.1674:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        total_counts = float(sum(histogram.ravel()))
        bin_width_x = float(self.right_x-self.left_x)/self.number_of_bins
        start_x = float(self.left_x)+bin_width_x/2
        bin_width_y = float(self.right_y-self.left_y)/self.number_of_bins
        start_y = float(self.left_y)+bin_width_y/2
        with open(self.filename,"w") as f:
            for i in xrange(self.number_of_bins):
                for j in xrange(self.number_of_bins):
                    print >> f, "{0} {1} {2}".format(start_x+i*bin_width_x,start_y+j*bin_width_y,histogram[i,j]/total_counts)
                print >> f
    #@-node:gcross.20100226131523.1674:write_out_totals
    #@-others
#@-node:gcross.20100226131523.1663:class PositionDensity2DHistogram
#@+node:gcross.20100226131523.1670:class PositionDensityMatrix2DHistogram
class PositionDensityMatrix2DHistogram(Histogram):
    #@    @+others
    #@+node:gcross.20100226131523.1671:__init__
    def __init__(self,number_of_slices,left_x,left_y,right_x,right_y,number_of_bins,filename):
        self.left_x = left_x
        self.left_y = left_y
        self.right_x = right_x
        self.right_y = right_y
        self.number_of_bins = number_of_bins
        self.center_slice_number = number_of_slices // 2
        self.histogram = zeros((number_of_bins,number_of_bins),dtype='i',order='Fortran')
        self.filename = filename
    #@-node:gcross.20100226131523.1671:__init__
    #@+node:gcross.20100226131523.1672:update
    def update(self):
        vpif.histograms.accumulate_2d_density_matrix(
            self.system.x[...,self.center_slice_number],
            self.system.x[...,self.center_slice_number+1],
            self.left_x,self.left_y,
            self.right_x,self.right_y,
            self.histogram
        )
    #@-node:gcross.20100226131523.1672:update
    #@+node:gcross.20100226131523.1676:write_out_totals
    def write_out_totals(self,histogram):
        ensure_path_to_file_exists(self.filename)
        with open(self.filename,"wb") as f:
            histogram.dump(f)
    #@-node:gcross.20100226131523.1676:write_out_totals
    #@-others
#@-node:gcross.20100226131523.1670:class PositionDensityMatrix2DHistogram
#@-node:gcross.20090902085220.2180:Histograms
#@+node:gcross.20090902085220.2231:Energy estimates
#@+node:gcross.20090902085220.2232:class TotalEnergyEstimate
class TotalEnergyEstimate(SingleAverageValueEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090902085220.2233:update
    def update(self):
        system = self.system
        for slice_number in [0,-1]:
            gradient_of_log_trial_fn, laplacian_of_log_trial_fn = system.compute_trial_derivatives(system.x[...,slice_number],system.xij2[...,slice_number])
            self.add(
                vpif.observables.compute_local_energy_estimate(
                    system.U[...,slice_number],
                    gradient_of_log_trial_fn, laplacian_of_log_trial_fn,
                    system.lambda_,
                )
            )
    #@-node:gcross.20090902085220.2233:update
    #@-others
#@-node:gcross.20090902085220.2232:class TotalEnergyEstimate
#@+node:gcross.20090902085220.2234:Slice estimates
#@+node:gcross.20090902085220.2235:class SliceEnergyEstimate
class SliceEnergyEstimate(SingleAverageValueEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090902085220.2236:__init__
    def __init__(self,slice_number,filename,label):
        SingleAverageValueEstimateAppendedToFile.__init__(self,filename,label)
        self.slice_number = slice_number
    #@-node:gcross.20090902085220.2236:__init__
    #@-others
#@-node:gcross.20090902085220.2235:class SliceEnergyEstimate
#@+node:gcross.20090902085220.2237:class EffectivePotentialSliceEnergyEstimate
class EffectivePotentialSliceEnergyEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090902085220.2238:update
    def update(self):
        self.add(self.system.compute_effective_potential(self.slice_number))
    #@-node:gcross.20090902085220.2238:update
    #@-others
#@-node:gcross.20090902085220.2237:class EffectivePotentialSliceEnergyEstimate
#@+node:gcross.20090902085220.2239:class PhysicalPotentialSliceEnergyEstimate
class PhysicalPotentialSliceEnergyEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090902085220.2240:update
    def update(self):
            self.add(self.system.compute_physical_potential(self.slice_number))
    #@-node:gcross.20090902085220.2240:update
    #@-others
#@-node:gcross.20090902085220.2239:class PhysicalPotentialSliceEnergyEstimate
#@+node:gcross.20090902085220.2241:class TotalPotentialSliceEnergyEstimate
class TotalPotentialSliceEnergyEstimate(SliceEnergyEstimate):
    #@    @+others
    #@+node:gcross.20090902085220.2242:update
    def update(self):
        self.add(self.system.compute_total_potential(self.slice_number))
    #@-node:gcross.20090902085220.2242:update
    #@-others
#@-node:gcross.20090902085220.2241:class TotalPotentialSliceEnergyEstimate
#@-node:gcross.20090902085220.2234:Slice estimates
#@+node:gcross.20090902085220.2243:Path estimates
#@+node:gcross.20090902085220.2244:class PathEnergyEstimates
class PathEnergyEstimates(AverageValuesEstimate):
    #@    @+others
    #@+node:gcross.20090902085220.2245:(fields)
    estimates = 0
    #@-node:gcross.20090902085220.2245:(fields)
    #@+node:gcross.20090902085220.2246:__init__
    def __init__(self,filename):
        self.filename = filename
    #@-node:gcross.20090902085220.2246:__init__
    #@+node:gcross.20090902085220.2247:write_out_totals
    def write_out_totals(self,total_estimates):
        ensure_path_to_file_exists(self.filename)
        center_slice_number = self.system.center_slice_number
        with open(self.filename,"w") as f:
            for slice_number, estimate in enumerate(total_estimates):
                print >> f, center_slice_number-abs(center_slice_number-slice_number), estimate, slice_number
    #@-node:gcross.20090902085220.2247:write_out_totals
    #@-others
#@-node:gcross.20090902085220.2244:class PathEnergyEstimates
#@+node:gcross.20090902085220.2252:class TotalPotentialPathEnergyEstimates
class TotalPotentialPathEnergyEstimates(PathEnergyEstimates):
    #@    @+others
    #@+node:gcross.20090902085220.2253:update
    def update(self):
        self.estimates += sum(self.system.U,axis=0)
    #@-node:gcross.20090902085220.2253:update
    #@-others
#@-node:gcross.20090902085220.2252:class TotalPotentialPathEnergyEstimates
#@-node:gcross.20090902085220.2243:Path estimates
#@-node:gcross.20090902085220.2231:Energy estimates
#@+node:gcross.20090902085220.2254:Position estimates
#@+node:gcross.20090902085220.2255:class AveragePositionEstimate
class AveragePositionEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    pass
#@-node:gcross.20090902085220.2255:class AveragePositionEstimate
#@+node:gcross.20090902085220.2256:class AverageAxialCoordinateEstimate
class AverageAxialCoordinateEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090902085220.2257:__init__
    def __init__(self,axis,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
        self.axis = axis
    #@-node:gcross.20090902085220.2257:__init__
    #@+node:gcross.20090902085220.2258:update
    def update(self):
        self.add(average(self.system.x[self.slice_number,:,self.axis]))
    #@-node:gcross.20090902085220.2258:update
    #@-others
#@-node:gcross.20090902085220.2256:class AverageAxialCoordinateEstimate
#@+node:gcross.20090902085220.2259:class AverageAxialDistanceEstimate
class AverageAxialDistanceEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090902085220.2260:__init__
    def __init__(self,axis,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
        self.axis = axis
    #@-node:gcross.20090902085220.2260:__init__
    #@+node:gcross.20090902085220.2261:update
    def update(self):
        self.add(average(abs(self.system.x[self.axis,:,self.slice_number])))
    #@-node:gcross.20090902085220.2261:update
    #@-others
#@-node:gcross.20090902085220.2259:class AverageAxialDistanceEstimate
#@+node:gcross.20090902085220.2262:class AverageRadiusEstimate
class AverageRadiusEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090902085220.2263:update
    def update(self):
        self.add(vpif.observables.compute_radius_average(self.system.x[...,self.slice_number]))
    #@nonl
    #@-node:gcross.20090902085220.2263:update
    #@-others
#@-node:gcross.20090902085220.2262:class AverageRadiusEstimate
#@+node:gcross.20090902085220.2264:class AveragePlaneRadiusEstimate
class AveragePlaneRadiusEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090902085220.2265:__init__
    def __init__(self,plane_axis_1,plane_axis_2,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
        assert plane_axis_1 >= 1
        assert plane_axis_2 >= 1
        assert not (plane_axis_1 == plane_axis_2)
        self.plane_axis_1 = plane_axis_1
        self.plane_axis_2 = plane_axis_2
    #@-node:gcross.20090902085220.2265:__init__
    #@+node:gcross.20090902085220.2266:update
    def update(self):
        self.add(
            vpif.observables.compute_plane_radius_average(
                self.system.x[...,self.slice_number],
                self.plane_axis_1,self.plane_axis_2
            )
        )
    #@-node:gcross.20090902085220.2266:update
    #@-others
#@-node:gcross.20090902085220.2264:class AveragePlaneRadiusEstimate
#@+node:gcross.20090902085220.2267:class AverageRecipricalPlaneRadiusSquaredEstimate
class AverageRecipricalPlaneRadiusSquaredEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090902085220.2268:__init__
    def __init__(self,plane_axis_1,plane_axis_2,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
        assert plane_axis_1 >= 1
        assert plane_axis_2 >= 1
        assert not (plane_axis_1 == plane_axis_2)
        self.plane_axis_1 = plane_axis_1
        self.plane_axis_2 = plane_axis_2
    #@-node:gcross.20090902085220.2268:__init__
    #@+node:gcross.20090902085220.2269:update
    def update(self):
        self.add(
            vpif.observables.compute_recip_plane_r_sq_average(
                self.system.x,
                self.plane_axis_1,self.plane_axis_2
            )
        )
    #@nonl
    #@-node:gcross.20090902085220.2269:update
    #@-others
#@-node:gcross.20090902085220.2267:class AverageRecipricalPlaneRadiusSquaredEstimate
#@+node:gcross.20090902085220.2270:class AverageParticleSeparationEstimate
class AverageParticleSeparationEstimate(AveragePositionEstimate):
    #@    @+others
    #@+node:gcross.20090902085220.2271:__init__
    def __init__(self,slice_number,filename,label):
        AveragePositionEstimate.__init__(self,slice_number,filename,label)
    #@-node:gcross.20090902085220.2271:__init__
    #@+node:gcross.20090902085220.2272:update
    def update(self):
        self.add(
            vpif.observables.compute_particle_separation_average(
                self.system.xij2[...,self.slice_number]
            )
        )
    #@nonl
    #@-node:gcross.20090902085220.2272:update
    #@-others
#@-node:gcross.20090902085220.2270:class AverageParticleSeparationEstimate
#@-node:gcross.20090902085220.2254:Position estimates
#@+node:gcross.20090902085220.2273:Rotation related estimates
#@+node:gcross.20090902085220.2274:class AverageAngularVelocityEstimate
class AverageAngularVelocityEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090902085220.2275:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpif.angular_momentum.compute_angular_derivatives(
            system.x[...,self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        self.add(average(first_derivatives))
    #@nonl
    #@-node:gcross.20090902085220.2275:update
    #@-others
#@-node:gcross.20090902085220.2274:class AverageAngularVelocityEstimate
#@+node:gcross.20090902085220.2276:class AverageAngularVelocitySquaredEstimate
class AverageAngularVelocitySquaredEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090902085220.2277:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpif.angular_momentum.compute_angular_derivatives(
            system.x[...,self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        first_derivatives **= 2
        self.add(average(first_derivatives))
    #@nonl
    #@-node:gcross.20090902085220.2277:update
    #@-others
#@-node:gcross.20090902085220.2276:class AverageAngularVelocitySquaredEstimate
#@+node:gcross.20090902085220.2278:class AverageRotationQuadraticTermEstimate
class AverageAngularVelocitySquaredEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090902085220.2279:update
    def update(self):
        system = self.system
        first_derivatives, _ = vpif.angular_momentum.compute_angular_derivatives(
            system.x[...,self.slice_number],
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        first_derivatives **= 2
        self.add(average(first_derivatives))
    #@nonl
    #@-node:gcross.20090902085220.2279:update
    #@-others
#@-node:gcross.20090902085220.2278:class AverageRotationQuadraticTermEstimate
#@+node:gcross.20090902085220.2280:class AverageAngularSeparationEstimate
class AverageAngularSeparationEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090902085220.2281:update
    def update(self):
        system = self.system
        x = system.x[...,self.slice_number]
        angles = arctan2(x[system.rotation_plane_axis_2],x[system.rotation_plane_axis_1])
        self.add(vpif.observables.compute_average_angular_separation(angles))
    #@nonl
    #@-node:gcross.20090902085220.2281:update
    #@-others
#@-node:gcross.20090902085220.2280:class AverageAngularSeparationEstimate
#@+node:gcross.20090902085220.2282:class AverageNeighborAngularSeparationEstimate
class AverageNeighborAngularSeparationEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090902085220.2283:update
    def update(self):
        system = self.system
        x = system.x[...,self.slice_number]
        angles = arctan2(x[system.rotation_plane_axis_2],x[system.rotation_plane_axis_1])
        self.add(vpif.observables.compute_avg_neighbor_angular_sep(angles))
    #@nonl
    #@-node:gcross.20090902085220.2283:update
    #@-others
#@-node:gcross.20090902085220.2282:class AverageNeighborAngularSeparationEstimate
#@+node:gcross.20090902085220.2284:class AverageRotationQuadraticTermEstimate
class AverageRotationQuadraticTermEstimate(SingleAverageValueAtSliceEstimateAppendedToFile):
    #@    @+others
    #@+node:gcross.20090902085220.2285:update
    def update(self):
        system = self.system
        x = system.x[...,self.slice_number]
        first_derivatives, _ = vpif.angular_momentum.compute_angular_derivatives(
            x,
            system.rotation_plane_axis_1, system.rotation_plane_axis_2,
            system.number_of_rotating_particles
        )
        term = (first_derivatives ** 2) / (x[system.rotation_plane_axis_2-1]**2+x[system.rotation_plane_axis_1-1]**2)
        self.add(average(term))
    #@nonl
    #@-node:gcross.20090902085220.2285:update
    #@-others
#@-node:gcross.20090902085220.2284:class AverageRotationQuadraticTermEstimate
#@-node:gcross.20090902085220.2273:Rotation related estimates
#@-others
#@-node:gcross.20090902085220.2161:Observable classes
#@+node:gcross.20090902085220.2355:Physics classes
#@+node:gcross.20100920135525.1710:Base classes
#@+node:gcross.20100920173658.1736:class MetaPhysics
class MetaPhysics(ABCMeta):
    #@    @+others
    #@+node:gcross.20100920232137.1735:(fields)
    known_physics_classes = {}
    #@-node:gcross.20100920232137.1735:(fields)
    #@+node:gcross.20100920232137.1736:__init__
    def __init__(self,name,bases,dict):
        ABCMeta.__init__(self,name,bases,dict)
        if name != "Physics":
            if Physics in bases:
                self.known_physics_classes[self] = [
                    method_name
                    for (method_name,method) in dict.iteritems()
                    if getattr(method,"__isabstractmethod__",False)
                ]
            for superclass, method_names in MetaPhysics.known_physics_classes.iteritems():
                if self == superclass: continue
                for method_name in method_names:
                    if method_name in dict:
                        if superclass in bases: continue
                        found = False
                        for base in bases:
                            if issubclass(base,superclass):
                                found = True
                                break
                        if found: continue
                        print "WARNING:  Class " + name + " has method '" + method_name + "', but since it is not a subclass of " + superclass.__name__ + " this method will be ignored by the simulator.  This was probably not what you intended, unless you did this on purpose."
                        break
    #@-node:gcross.20100920232137.1736:__init__
    #@-others
#@-node:gcross.20100920173658.1736:class MetaPhysics
#@+node:gcross.20100920135525.1716:class Physics
class Physics(object):
    __metaclass__ = MetaPhysics
    #@    @+others
    #@+node:gcross.20100920135525.1786:__init__
    def __init__(self,system):
        self.system = system
    #@-node:gcross.20100920135525.1786:__init__
    #@-others
#@-node:gcross.20100920135525.1716:class Physics
#@+node:gcross.20100920135525.1711:class PotentialFunction
class PotentialFunction(Physics):
    #@    @+others
    #@+node:gcross.20100920135525.1776:(fields)
    potential_function_registered = False
    #@nonl
    #@-node:gcross.20100920135525.1776:(fields)
    #@+node:gcross.20100920135525.1777:__init__
    def __init__(self,system,potentials_field_name):
        if not self.potential_function_registered:
            Physics.__init__(self,system)
            if isinstance(self,PhysicalPotentialFunction) and isinstance(self,EffectivePotentialFunction):
                raise TypeError("Error initiaizing " + self.__class__ + ": A potential function cannot be both an effective potential and a physical potential.  (FYI, these two superclasses are functionally the same, but they are distinguished so that the user can get an estimate of the physical portion of the potential.)")
            getattr(system,potentials_field_name).append(self)
            self.potential_function_registered = True
    #@-node:gcross.20100920135525.1777:__init__
    #@+node:gcross.20100920135525.1778:(abstract methods)
    @abstractmethod
    def accumulate_potential(self,x,xij2,U,gradU2): pass
    #@-node:gcross.20100920135525.1778:(abstract methods)
    #@-others
#@-node:gcross.20100920135525.1711:class PotentialFunction
#@+node:gcross.20100920135525.1738:class PhysicalPotentialFunction
class PhysicalPotentialFunction(PotentialFunction):
    def __init__(self,system):
        PotentialFunction.__init__(self,system,"physical_potentials")
#@-node:gcross.20100920135525.1738:class PhysicalPotentialFunction
#@+node:gcross.20100920135525.1740:class EffectivePotentialFunction
class EffectivePotentialFunction(PotentialFunction):
    def __init__(self,system):
        PotentialFunction.__init__(self,system,"effective_potentials")
#@-node:gcross.20100920135525.1740:class EffectivePotentialFunction
#@+node:gcross.20100920135525.1715:class TrialFunction
class TrialFunction(Physics):
    #@    @+others
    #@+node:gcross.20100920135525.1779:(fields)
    trial_function_registered = False
    #@-node:gcross.20100920135525.1779:(fields)
    #@+node:gcross.20100920135525.1780:__init__
    def __init__(self,system):
        if not self.trial_function_registered:
            Physics.__init__(self,system)
            system.trial_functions.append(self)
            self.trial_function_registered = True
    #@nonl
    #@-node:gcross.20100920135525.1780:__init__
    #@+node:gcross.20100920135525.1781:(abstract methods)
    @abstractmethod
    def compute_trial_weight(self,x,xij2): pass

    @abstractmethod
    def accumulate_trial_derivatives(self,x,xij2,gradient_of_log_trial_fn,laplacian_of_log_trial_fn): pass
    #@-node:gcross.20100920135525.1781:(abstract methods)
    #@-others
#@-node:gcross.20100920135525.1715:class TrialFunction
#@+node:gcross.20100920135525.1717:class GreensFunction
class GreensFunction(Physics):
    #@    @+others
    #@+node:gcross.20100920135525.1782:(fields)
    greens_function_registered = False
    #@-node:gcross.20100920135525.1782:(fields)
    #@+node:gcross.20100920135525.1783:__init__
    def __init__(self,system):
        if not self.greens_function_registered:
            Physics.__init__(self,system)
            system.greens_functions.append(self)
            self.greens_function_registered = True
    #@nonl
    #@-node:gcross.20100920135525.1783:__init__
    #@+node:gcross.20100920135525.1784:(abstract methods)
    @abstractmethod
    def compute_greens_function(self,x,xij2,U,gradU2,lam,dt,slice_start,slice_end,particle_number): pass
    #@-node:gcross.20100920135525.1784:(abstract methods)
    #@-others
#@-node:gcross.20100920135525.1717:class GreensFunction
#@-node:gcross.20100920135525.1710:Base classes
#@+node:gcross.20100920135525.1718:Harmonic Oscillator
#@+node:gcross.20090902085220.2356:class HarmonicOscillatorPotentialFunction
class HarmonicOscillatorPotentialFunction(PhysicalPotentialFunction):
    #@    @+others
    #@+node:gcross.20090902085220.2359:__init__
    def __init__(self,system,harmonic_oscillator_frequencies):
        assert len(harmonic_oscillator_frequencies) == system.number_of_dimensions
        PhysicalPotentialFunction.__init__(self,system)
        self.potential_coefficients = array(harmonic_oscillator_frequencies)**2
    #@-node:gcross.20090902085220.2359:__init__
    #@+node:gcross.20090902085220.2360:accumulate_potential
    def accumulate_potential(self,x,xij2,U,gradU2):
        vpif.harmonic_oscillator.accumulate_potential(
            x,
            self.potential_coefficients,
            U,
            gradU2
        )
    #@-node:gcross.20090902085220.2360:accumulate_potential
    #@-others
#@-node:gcross.20090902085220.2356:class HarmonicOscillatorPotentialFunction
#@+node:gcross.20100920135525.1725:class HarmonicOscillatorTrialFunction
class HarmonicOscillatorTrialFunction(TrialFunction):
    #@    @+others
    #@+node:gcross.20100920135525.1727:__init__
    def __init__(self,system,harmonic_oscillator_frequencies):
        assert len(harmonic_oscillator_frequencies) == system.number_of_dimensions
        TrialFunction.__init__(self,system)
        self.trial_coefficients = array(harmonic_oscillator_frequencies)
    #@-node:gcross.20100920135525.1727:__init__
    #@+node:gcross.20100920135525.1729:compute_trial_weight
    def compute_trial_weight(self,x,xij2):
        return vpif.harmonic_oscillator.compute_trial_weight(
            x,
            self.trial_coefficients
        )
    #@-node:gcross.20100920135525.1729:compute_trial_weight
    #@+node:gcross.20100920135525.1730:accumulate_trial_derivatives
    def accumulate_trial_derivatives(self,
            x,xij2,
            gradient_of_log_trial_fn,laplacian_of_log_trial_fn
        ):
        vpif.harmonic_oscillator.accumulate_trial_derivatives(
            x,self.trial_coefficients,
            gradient_of_log_trial_fn,laplacian_of_log_trial_fn
        )
    #@-node:gcross.20100920135525.1730:accumulate_trial_derivatives
    #@-others
#@-node:gcross.20100920135525.1725:class HarmonicOscillatorTrialFunction
#@+node:gcross.20100920135525.1741:class HarmonicOscillator
class HarmonicOscillator(HarmonicOscillatorPotentialFunction,HarmonicOscillatorTrialFunction):
    #@    @+others
    #@+node:gcross.20100920135525.1743:__init__
    def __init__(self,system,harmonic_oscillator_frequencies):
        HarmonicOscillatorPotentialFunction.__init__(self,system,harmonic_oscillator_frequencies)
        HarmonicOscillatorTrialFunction.__init__(self,system,harmonic_oscillator_frequencies)
    #@-node:gcross.20100920135525.1743:__init__
    #@-others
#@-node:gcross.20100920135525.1741:class HarmonicOscillator
#@-node:gcross.20100920135525.1718:Harmonic Oscillator
#@+node:gcross.20100920135525.1744:Hard Sphere Interaction
#@+node:gcross.20090902085220.2368:class HardSphereInteractionPotentialFunction
class HardSphereInteractionPotentialFunction(PhysicalPotentialFunction,GreensFunction):
    #@    @+others
    #@+node:gcross.20090902085220.2369:__init__
    def __init__(self,system,hard_sphere_radius):
        PhysicalPotentialFunction.__init__(self,system)
        GreensFunction.__init__(self,system)
        self.hard_sphere_radius = hard_sphere_radius
        self.hard_sphere_radius_squared = hard_sphere_radius ** 2
    #@-node:gcross.20090902085220.2369:__init__
    #@+node:gcross.20090902085220.2370:accumulate_potential
    def accumulate_potential(self,x,xij2,U,gradU):
        if vpif.hard_sphere_interaction.has_collision(
            xij2, self.hard_sphere_radius_squared
        ): raise MoveRejected
    #@-node:gcross.20090902085220.2370:accumulate_potential
    #@+node:gcross.20090902085220.2384:compute_greens_function
    def compute_greens_function(self,
            x,xij2,
            U,gradU2,
            lam,dt,
            slice_start,slice_end,
            particle_number
        ): return \
            vpif.hard_sphere_interaction.compute_greens_function(
                xij2,dt,
                self.hard_sphere_radius,
                slice_start,slice_end,particle_number
            )
    #@-node:gcross.20090902085220.2384:compute_greens_function
    #@-others
#@-node:gcross.20090902085220.2368:class HardSphereInteractionPotentialFunction
#@+node:gcross.20100920135525.1764:class HardSphereInteractionTrialFunction
class HardSphereInteractionTrialFunction(TrialFunction):
    #@    @+others
    #@+node:gcross.20100920135525.1766:__init__
    def __init__(self,system,hard_sphere_radius):
        TrialFunction.__init__(self,system)
        self.hard_sphere_radius = hard_sphere_radius
        self.hard_sphere_radius_squared = hard_sphere_radius ** 2
    #@-node:gcross.20100920135525.1766:__init__
    #@+node:gcross.20100920135525.1770:accumulate_trial_derivatives
    def accumulate_trial_derivatives(self,x,xij2,gradient_of_log_trial_fn,laplacian_of_log_trial_fn):
        vpif.hard_sphere_interaction.accumulate_trial_derivatives(
            x,xij2,self.hard_sphere_radius,
            gradient_of_log_trial_fn,laplacian_of_log_trial_fn
        )
    #@-node:gcross.20100920135525.1770:accumulate_trial_derivatives
    #@+node:gcross.20100920135525.1772:compute_trial_weight
    def compute_trial_weight(self,x,xij2):
        weight, reject_flag = \
            vpif.hard_sphere_interaction.compute_trial_weight(
                xij2,
                self.hard_sphere_radius
            )
        if reject_flag: raise MoveRejected
        return weight
    #@-node:gcross.20100920135525.1772:compute_trial_weight
    #@-others
#@-node:gcross.20100920135525.1764:class HardSphereInteractionTrialFunction
#@+node:gcross.20100920135525.1751:class HardSphereInteraction
class HardSphereInteraction(HardSphereInteractionPotentialFunction,HardSphereInteractionTrialFunction):
  #@  @+others
  #@+node:gcross.20100920135525.1752:__init__
  def __init__(self,system,hard_sphere_radius):
      HardSphereInteractionPotentialFunction.__init__(self,system,hard_sphere_radius)
      HardSphereInteractionTrialFunction.__init__(self,system,hard_sphere_radius)
  #@-node:gcross.20100920135525.1752:__init__
  #@-others
#@-node:gcross.20100920135525.1751:class HardSphereInteraction
#@-node:gcross.20100920135525.1744:Hard Sphere Interaction
#@+node:gcross.20100920135525.1773:Greens Functions
#@+node:gcross.20090902085220.2373:class SecondOrderGreensFunction
class SecondOrderGreensFunction(GreensFunction):
  #@  @+others
  #@+node:gcross.20100106123346.1702:__init__
  def __init__(self,system):
      GreensFunction.__init__(self,system)
      self.weights = vpif.gfn.initialize_2nd_order_weights(system.number_of_slices)
  #@-node:gcross.20100106123346.1702:__init__
  #@+node:gcross.20090902085220.2382:compute_greens_function
  def compute_greens_function(
          self,
          x,xij2,
          U,gradU2,
          lambda_,dt,
          slice_start,slice_end,
          particle_number
      ):  return vpif.gfn.gfn2_sp(slice_start,slice_end,U,self.weights,dt)
  #@-node:gcross.20090902085220.2382:compute_greens_function
  #@-others
#@-node:gcross.20090902085220.2373:class SecondOrderGreensFunction
#@+node:gcross.20091220132355.1699:class FourthOrderGreensFunction
class FourthOrderGreensFunction(GreensFunction):
  #@  @+others
  #@+node:gcross.20091220132355.1703:__init__
  def __init__(self,system):
      GreensFunction.__init__(self,system)
      self.weights = vpif.gfn.initialize_4th_order_weights(system.number_of_slices)
  #@-node:gcross.20091220132355.1703:__init__
  #@+node:gcross.20091220132355.1701:compute_greens_function
  def compute_greens_function(
          self,
          x,xij2,
          U,gradU2,
          lambda_,dt,
          slice_start,slice_end,
          particle_number
      ):  return vpif.gfn.gfn4_sp(slice_start,slice_end,U,gradU2,self.weights[0],self.weights[1],lambda_,dt)
  #@-node:gcross.20091220132355.1701:compute_greens_function
  #@-others
#@-node:gcross.20091220132355.1699:class FourthOrderGreensFunction
#@-node:gcross.20100920135525.1773:Greens Functions
#@-node:gcross.20090902085220.2355:Physics classes
#@+node:gcross.20090902085220.2288:class System
class System(object):
    #@    @+others
    #@+node:gcross.20100920232137.1737:__slots__
    __slots__ = [
        "number_of_slices",
        "number_of_particles",
        "number_of_dimensions",
        "center_slice_number",
        "lambda_",
        "initial_particle_distribution_size",
        "total_number_of_observations",
        "number_of_prethermalization_steps",
        "dM",
        "move_type_probabilities",
        "move_type_differentials",
        "low_swap_dimension",
        "high_swap_dimension",
        "x",
        "xij2",
        "U",
        "gradU2",
        "slice_move_attempted_counts",
        "slice_move_accepted_counts",
        "move_type_attempted_counts",
        "move_type_accepted_counts",
        "number_of_observations",
        "total_number_of_observations",
        "number_of_thermalizations_per_observation",
        "observables",
        "physical_potentials",
        "effective_potentials",
        "trial_functions",
        "greens_functions",
        ]
    #@-node:gcross.20100920232137.1737:__slots__
    #@+node:gcross.20090902085220.2333:Initialization
    #@+node:gcross.20090902085220.2303:__init__
    def __init__(
            self,
            number_of_slices,number_of_particles,number_of_dimensions,
            lambda_,initial_particle_distribution_size,
            total_number_of_observations,number_of_prethermalization_steps,
            dM,
            move_type_probabilities,move_type_differentials,
            low_swap_dimension,high_swap_dimension,
        ):
        self.number_of_slices = number_of_slices
        self.number_of_particles = number_of_particles
        self.number_of_dimensions = number_of_dimensions
        self.lambda_ = lambda_
        self.initial_particle_distribution_size = initial_particle_distribution_size
        self.total_number_of_observations = total_number_of_observations
        self.number_of_prethermalization_steps = number_of_prethermalization_steps
        self.dM = dM
        self.move_type_probabilities = move_type_probabilities
        self.move_type_differentials = move_type_differentials
        self.low_swap_dimension = low_swap_dimension
        self.high_swap_dimension = high_swap_dimension

        vpif.rand_utils.init_seed(my_rank)

        self.initialize_lattice()

        self.U = zeros((number_of_particles,number_of_slices),dtype=double,order='Fortran')
        self.gradU2 = zeros((number_of_slices),dtype=double,order='Fortran')

        self.slice_move_attempted_counts = zeros((number_of_slices,),'i')
        self.slice_move_accepted_counts = zeros((number_of_slices,),'i')

        self.move_type_attempted_counts = zeros((3,),'i')
        self.move_type_accepted_counts = zeros((3,),'i')

        self.number_of_observations = self.total_number_of_observations // number_of_processors + 1
        self.total_number_of_observations = self.number_of_observations * number_of_processors

        self.number_of_thermalizations_per_observation = number_of_particles * number_of_slices // self.dM

        assert (number_of_slices % 2 == 0 and number_of_slices % 4 == 2)
        self.center_slice_number = number_of_slices // 2

        self.observables = []
        self.physical_potentials = []
        self.effective_potentials = []
        self.trial_functions = []
        self.greens_functions = []
    #@-node:gcross.20090902085220.2303:__init__
    #@+node:gcross.20090902085220.2332:(re)initialize_lattice
    def initialize_lattice(self):
        self.x = vpif.lattice.make_lattice(
            self.initial_particle_distribution_size,
            self.number_of_slices,self.number_of_particles,self.number_of_dimensions
        )
        self.xij2 = vpif.xij.compute_xij(self.x)

    reinitialize_lattice = initialize_lattice
    #@-node:gcross.20090902085220.2332:(re)initialize_lattice
    #@-node:gcross.20090902085220.2333:Initialization
    #@+node:gcross.20090902085220.2300:Observable management
    #@+node:gcross.20090902085220.2301:add_observable
    def add_observable(self,observable):
        self.observables.append(observable)
        observable.system = self
    #@-node:gcross.20090902085220.2301:add_observable
    #@+node:gcross.20090902085220.2302:total_and_write_observables
    def total_and_write_observables(self):
        for observable in self.observables:
            observable.total_and_write()
    #@-node:gcross.20090902085220.2302:total_and_write_observables
    #@+node:gcross.20090911091023.2450:reset_observables
    def reset_observables(self):
        for observable in self.observables:
            observable.reset()
    #@-node:gcross.20090911091023.2450:reset_observables
    #@-node:gcross.20090902085220.2300:Observable management
    #@+node:gcross.20090902085220.2334:Physics
    #@+node:gcross.20090902085220.2342:Potential
    #@+node:gcross.20090902085220.2335:compute_potential
    def compute_potential(self,x,xij2):
        U = zeros(x.shape[-2:],dtype=double,order='Fortran')
        gradU2 = zeros(x.shape[-1:],dtype=double,order='Fortran')
        try:
            for potential in itertools.chain(
                    self.physical_potentials,
                    self.effective_potentials
                ): potential.accumulate_potential(x,xij2,U,gradU2)
            return U, gradU2, False
        except MoveRejected:
            return U, gradU2, True
    #@-node:gcross.20090902085220.2335:compute_potential
    #@+node:gcross.20090902085220.2337:compute_effective_potential
    def compute_effective_potential(self,slice_number):
        return (
             self.compute_total_potential(slice_number)
           - self.compute_physical_potential(slice_number)
        )
    #@-node:gcross.20090902085220.2337:compute_effective_potential
    #@+node:gcross.20090902085220.2341:compute_physical_potential
    def compute_physical_potential(self,slice_number):
        x = self.x[slice_number:slice_number+1]
        xij2 = self.xij2[slice_number:slice_number+1]
        U = zeros((1,x.shape[1]),dtype=double,order='Fortran')
        gradU2 = zeros((1,),dtype=double,order='Fortran')
        for potential in self.physical_potentials:
            potential.accumulate_potential(x,xij2,U,gradU2)
        return sum(U)
    #@-node:gcross.20090902085220.2341:compute_physical_potential
    #@+node:gcross.20090902085220.2339:compute_total_potential
    def compute_total_potential(self,slice_number):
        return sum(self.U[slice_number])
    #@-node:gcross.20090902085220.2339:compute_total_potential
    #@-node:gcross.20090902085220.2342:Potential
    #@+node:gcross.20090902085220.2347:Trial
    #@+node:gcross.20090902085220.2348:compute_trial_weight
    def compute_trial_weight(self,x,xij2):
        try:
            return __builtin__.sum(trial.compute_trial_weight(x,xij2) for trial in self.trial_functions), False
        except MoveRejected:
            return 0.0, True
    #@-node:gcross.20090902085220.2348:compute_trial_weight
    #@+node:gcross.20090902085220.2349:compute_trial_derivatives
    def compute_trial_derivatives(self,x,xij2):
        gradient_of_log_trial_fn = zeros(x.shape,dtype=double,order='Fortran')
        laplacian_of_log_trial_fn = zeros((),dtype=double,order='Fortran')
        for trial in self.trial_functions:
            trial.accumulate_trial_derivatives(
                x,xij2,
                gradient_of_log_trial_fn,laplacian_of_log_trial_fn
            )
        return gradient_of_log_trial_fn, laplacian_of_log_trial_fn
    #@-node:gcross.20090902085220.2349:compute_trial_derivatives
    #@-node:gcross.20090902085220.2347:Trial
    #@+node:gcross.20090902085220.2351:compute_greens_function
    def compute_greens_function(self,
        x,xij2,
        U,gradU2,
        lam,dt,
        slice_start,slice_end,
        particle_number
        ):
        return __builtin__.sum(
            greens_function.compute_greens_function(
                x,xij2,
                U,gradU2,
                lam,dt,
                slice_start,slice_end,
                particle_number
            ) for greens_function in self.greens_functions
        )

    #@-node:gcross.20090902085220.2351:compute_greens_function
    #@-node:gcross.20090902085220.2334:Physics
    #@+node:gcross.20090902085220.2304:run
    def run(self):
        #@    << Pre-run checks >>
        #@+node:gcross.20090902085220.2352:<< Pre-run checks >>
        if len(self.physical_potentials) == 0:
                raise Exception("No physical potentials have been specified!")
        if len(self.greens_functions) == 0:
                raise Exception("No Green's functions have been specified!")
        if len(self.trial_functions) == 0:
                raise Exception("No trial functions have been specified!")
        #@-node:gcross.20090902085220.2352:<< Pre-run checks >>
        #@nl
        #@    << Stash properties into local variables >>
        #@+node:gcross.20090902085220.2305:<< Stash properties into local variables >>
        x = self.x
        xij2 = self.xij2
        U = self.U
        gradU2 = self.gradU2
        number_of_prethermalization_steps = self.number_of_prethermalization_steps
        number_of_thermalizations_per_observation = self.number_of_thermalizations_per_observation
        move_type_probabilities = self.move_type_probabilities
        move_type_differentials = self.move_type_differentials
        dM = self.dM
        lambda_ = self.lambda_
        low_swap_dimension = self.low_swap_dimension
        high_swap_dimension = self.high_swap_dimension
        slice_move_attempted_counts = self.slice_move_attempted_counts
        move_type_attempted_counts = self.move_type_attempted_counts
        slice_move_accepted_counts = self.slice_move_accepted_counts
        move_type_accepted_counts = self.move_type_accepted_counts
        compute_potential = self.compute_potential
        compute_trial_weight = self.compute_trial_weight
        compute_greens_function = self.compute_greens_function
        observables = self.observables
        #@nonl
        #@-node:gcross.20090902085220.2305:<< Stash properties into local variables >>
        #@nl
        #@    << Ensure system is in valid state >>
        #@+node:gcross.20100920232137.1738:<< Ensure system is in valid state >>
        number_of_initialization_attempts = 1
        while compute_potential(x,xij2)[-1] and number_of_initialization_attempts < 1000:
            self.reinitialize_lattice()
            x = self.x
            xij2 = self.xij2
        if number_of_initialization_attempts == 1000:
            raise FailureToContructValidSystemException("Failed in 1000 attempts to construct a system which was not rejected by the potential function.")
        #@-node:gcross.20100920232137.1738:<< Ensure system is in valid state >>
        #@nl
        #@    << Prethermalize the system >>
        #@+node:gcross.20090902085220.2306:<< Prethermalize the system >>
        vpif.thermalize.thermalize_path(
            x,xij2,
            U,gradU2,
            number_of_prethermalization_steps,
            move_type_probabilities,move_type_differentials,
            dM,
            lambda_,
            low_swap_dimension,high_swap_dimension,
            slice_move_attempted_counts,move_type_attempted_counts,
            slice_move_accepted_counts,move_type_accepted_counts,
            compute_potential,compute_trial_weight,compute_greens_function,
        )
        #@-node:gcross.20090902085220.2306:<< Prethermalize the system >>
        #@nl
        #@    << Main iteration >>
        #@+node:gcross.20090902085220.2307:<< Main iteration >>
        decile = self.number_of_observations // 10
        for number_completed in xrange(self.number_of_observations):
            vpif.thermalize.thermalize_path(
                x,xij2,
                U,gradU2,
                number_of_thermalizations_per_observation,
                move_type_probabilities,move_type_differentials,
                dM,
                lambda_,
                low_swap_dimension,high_swap_dimension,
                slice_move_attempted_counts,move_type_attempted_counts,
                slice_move_accepted_counts,move_type_accepted_counts,
                compute_potential,compute_trial_weight,compute_greens_function,
            )

            for observable in observables:
                observable.update()
            if decile > 0 and (number_completed % decile == 0) and (my_rank == 0):
                print "{0:.0%} complete;  local bridge move acceptance rate = {1:.0%}, local rigid move acceptance rate = {2:.0%}".format(
                    float(number_completed)/self.number_of_observations,
                    float(move_type_accepted_counts[0])/move_type_attempted_counts[0],
                    float(move_type_accepted_counts[1])/move_type_attempted_counts[1],
                )
        #@-node:gcross.20090902085220.2307:<< Main iteration >>
        #@nl
    #@-node:gcross.20090902085220.2304:run
    #@-others
#@-node:gcross.20090902085220.2288:class System
#@-node:gcross.20090902085220.2385:Classes
#@-others
#@-node:gcross.20090902085220.2158:@thin vpi/__init__.py
#@-leo
