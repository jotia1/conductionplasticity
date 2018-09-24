function [ outp ] = spikingnet( net )
%SPIKINGNET Computational engine of a spiking network
%   Given a network description run the simulation and compute the final
%   state variables of the network

outp = struct();
outp.timing_info.init_time = tic;
validatenet(net);
rng(net.rand_seed); 

ms_per_sec = 1000;

delays = sparse(net.delays);
variance = sparse(net.variance);
w = sparse(net.w);

N = net.N;
v = ones(N, 1) * net.rest;
vt = zeros(numel(net.voltages_to_save), net.sim_time_sec * ms_per_sec);

dApre = []; % TODO
dApost = []; % TODO

for sec = 1 :net.sim_time_sec
    
    % Trim data into seconds to speed searching later
    idxs = net.supplied_ts > (sec - 1) * 1000 & net.supplied_ts <= (sec * 1000);
    inp_trimmed = net.supplied_input(idxs);
    ts_trimmed = net.supplied_ts(idxs);
    
    for ms = 1 : ms_per_sec
        
        time = (sec - 1) * ms_per_sec + ms;
        vt(:, time) = v(net.voltages_to_save);
        
        %% Update membrane voltages  
        v = v + (net.v_rest + Iapp - v) / net.neuron_tau;
        
        %% Deal with neurons that just spiked
        
        
        
    end
    
end



end

