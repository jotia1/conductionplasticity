function [ outp ] = spikingnet( net )
%SPIKINGNET Computational engine of a spiking network
%   Given a network description run the simulation and compute the final
%   state variables of the network
%   In connectivity matrices, rows are the pre-synaptic (from) neurons and
%   columns are the post-synaptic (to) neurons.

outp = struct();
outp.timing_info.init_time = tic;
validatenet(net);
rng(net.rand_seed); 

ms_per_sec = 1000;

delays = sparse(net.delays);
variance = sparse(net.variance);
w = sparse(net.w);

N = net.N;
v = ones(N, 1) * net.v_rest;
vt = zeros(numel(net.voltages_to_save), net.sim_time_sec * ms_per_sec);
last_spike_time = zeros(net.N, 1) * -Inf;


conns = w ~= 0;
dApre = sparse(zeros(size(w)));
dApost = sparse(zeros(size(w)));
STDPdecaypre = exp(-1/net.taupre);
STDPdecaypost = exp(-1/net.taupost);
active_spikes = cell(net.delay_max, 1);  % To track when spikes arrive
active_idx = 1;

for sec = 1 :net.sim_time_sec
    
    spike_time_trace = [];
    
    % Trim data into seconds to speed searching later
    idxs = net.ts > (sec - 1) * 1000 & net.ts <= (sec * 1000);
    inp_trimmed = net.inp(idxs);
    ts_trimmed = net.ts(idxs);
    
    for ms = 1 : ms_per_sec
        
        time = (sec - 1) * ms_per_sec + ms;
        vt(:, time) = v(net.voltages_to_save);
        
        %% Calculate input
        Iapp = zeros(size(v));
        t0 = time - last_spike_time;  % TODO will need to reshape this (vector) to match delays (matrix)
        t0_negu = t0 - delays;
        p = net.fgi ./ sqrt(2 * pi * variance);
        g = p .* exp(- (t0_negu .^ 2) ./ (2 * variance));
        gaussian_values = w .* g;
        gaussian_values(isnan(gaussian_values)) = 0;
        Iapp = sum(gaussian_values, 1)';
        
        %% Update membrane voltages  
        v = v + (net.v_rest + Iapp - v) / net.neuron_tau;
        vt(:, time) = v(net.voltages_to_save);
        
        %% Deal with neurons that just spiked
        fired_naturally = find(v >= net.v_thres);
        fired_pixels = inp_trimmed(ts_trimmed == time);
        fired = unique([fired_naturally; fired_pixels';]); % unique in case of SDVL supervision later
        spike_time_trace = [spike_time_trace; time*ones(length(fired),1), fired];
        last_spike_time(fired) = time; 
        
        v(fired_naturally) = net.v_reset;
        % Any pre-synaptics weights should be increased
        w(:, fired) = w(:, fired) + (dApost(:, fired) .* conns(:, fired));
        % any post synaptic weights decreased
        w(fired, :) = w(fired, :) + (dApre(fired, :) .* conns(fired, :));
        dApost(fired, :) = dApost(fired, :) + net.Apost;
        dApre(:, fired) = dApre(:, fired) + net.Apre;

        for spike = 1 : length(fired_naturally)
            neuron_idx = fired(spike);
            
            
        end
        dApre = dApre * STDPdecaypre;
        dApost = dApost * STDPdecaypost;
        active_idx = mod(active_idx, net.delay_max) + 1;
    
        % Synaptic bounding - limit w to [0, w_max]
        w = max(0, min(net.w_max, w));        
        
    end

    
end

plot(vt(4:5, :)')
%set(gca, 'Color', 'k')
legend({'N4', 'N5'});

end















