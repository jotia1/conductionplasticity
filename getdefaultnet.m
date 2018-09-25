function [ net ] = getdefaultnet(  )
%GETDEFAULTNET Return a default network specification as an example
%   Returns a network definition that can be passed to spikingnet.m to run
%   a spiking net simulation. The network structure and parameters are
%   specified here.

net = struct();

net.group_sizes = [3, 2];
net.N = sum(net.group_sizes);
net.rand_seed = 1;

net.delays = zeros(net.N);
net.delays(1:3, 4:5) = 1;
net.variance = zeros(net.N);
net.variance(1:3, 4:5) = 0.1;
net.w = zeros(net.N);
net.w(1:3, 4:5) = 1;

net.fgi = 130;
net.neuron_tau = 20;
net.delay_max = 15;

net.v_rest = -65;
net.v_reset = -70;
net.v_thres = -55;
net.w_max = 10;
net.taupre = 20;
net.taupost = 20;
net.Apost = 0.1;
net.Apre = -0.1;

net.sim_time_sec = 1;
net.inp = [1, 2, 3, 1, 2 ,3];
net.ts = [1, 5, 8, 501, 505, 508];

net.voltages_to_save = [1 : net.N];



end

