function [ net ] = getdefaultnet(  )
%GETDEFAULTNET Return a default network specification as an example
%   Returns a network definition that can be passed to spikingnet.m to run
%   a spiking net simulation. The network structure and parameters are
%   specified here.

net = struct();

net.group_sizes = [3, 2];
net.N = sum(net.group_sizes);
net.rand_seed = 1;

net.delays = [];
net.variance = [];
net.w = [];

net.rest = -65;
net.reset = -70;

net.sim_time_sec = 10;
net.inp = [];
net.ts = [];

net.voltages_to_save = [1 : net.N];



end

