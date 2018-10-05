
net = getdefaultnet();

% Make net larger (10 input, 2 outputs)
net.group_sizes = [10, 2];
N_inp = net.group_sizes(1);
net.N = sum(net.group_sizes);
ofgi = net.fgi;
%net.fgi = net.fgi * 0.55;
net.a1 = 3;
net.a2 = 2;
net.b1 = 5;
net.b2 = 5;
net.nu = 0.025;
net.nv = 0.025;
net.delays = zeros(net.N);
net.variance = zeros(net.N);
net.w = zeros(net.N);
net.w(1:N_inp, N_inp:net.N) = 1;
net.sim_time_sec = 100;
net.neuron_to_plot = N_inp + 1;

net.voltages_to_save = [1 : net.N];
net.variance_to_save = [11, 12];
net.delays_to_save = [11, 12];

output_folder = getValidFilename('task1');

num_trials = 10;
Ts = 12;

range = 1:10000; %0.01:0.02:0.2;
for di = 1 : numel(range)
    d = range(di);
    s1 = randi([1 Ts], 1, 10);
    s2 = randi([1 Ts], 1, 10);  
    
    % normalise
    s1 = s1 - min(s1(:)) + 1; s2 = s2 - min(s2(:)) + 1;
    diff = sum(abs(s1 - s2)) ./ (net.N * Ts);
    
    net.inp = repmat(1:N_inp, 1, net.sim_time_sec * 2);
    offsets = reshape(repmat(((0:(net.sim_time_sec * 2) -1 )' * 500)', N_inp, 1), 1, net.sim_time_sec * N_inp * 2);
    net.ts = repmat([s1, s2], 1, net.sim_time_sec) + offsets;
        
    fgi_range = 0.5 : 0.01 : 0.7;
    for fgi_multi = 1 : numel(fgi_range)
        fgi_mult = fgi_range(fgi_multi);
        
        net.fgi = ofgi * fgi_mult;
        
        hyperparams = struct();
        hyperparams.d = diff;
        hyperparams.fgi = net.fgi;
        hyperparams.description = 'Trial runs to get task 1 working and to test plotting etc.';
        
        net.delays(1:N_inp, N_inp:net.N) = 2 + (10-2) .* rand(N_inp, net.N - N_inp + 1);
        net.variance(1:N_inp, N_inp:net.N) = 1 + (3 - 1) .* rand(N_inp, net.N - N_inp + 1);
        
        out = spikingnet(net);
        
        save(sprintf('%s/%d-%d.mat', output_folder, di, fgi_multi), 'net', 'out', 'hyperparams');
        
    end
end


function [ name ] =  getValidFilename( ext )
%% Given a name (ext) return the name of a new empty folder.

% Get a valid folder to dump results to
count = 0;
exp_code = ext;
output_folder = sprintf('%s-%02d', exp_code, count);
while exist(output_folder, 'dir') == 7
    count = count + 1;
    output_folder = sprintf('%s-%02d', exp_code, count);
end

mkdir(output_folder);
name = output_folder;

end















